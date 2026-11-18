"""Build a release of the data format from a snapshot, by a recipe.

The format is defined in `docs/data-format.md` of the app repository `neherlab/pan-genome-visualization`. A release
is one folder of the data root, named by its id; `releases.json` of the data root lists it. The same recipe, snapshot,
and tool commit give the same bytes on every host that runs the Python of the data container: no clock or host value
enters an output, the gzip header is fixed, and every list has a defined order.
"""

import gzip
import hashlib
import json
import logging
import os
import re
import shutil
import sys
import tomllib
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from panx_data.archive import TAR_END, GzipWriter, check_zlib, gzip_member, json_bytes, tar_entry
from panx_data.catalog import CATALOG_CONFIG, Catalog, CatalogDataset, js_values, load_catalog
from panx_data.convert import (
    CONSENSUS,
    MetadataConfig,
    cluster_columns,
    cluster_entry,
    fasta_ids,
    gene_tree_json,
    gene_tree_newick,
    metadata_and_strains,
    metadata_tsv,
    newick_tips,
    non_root_nodes,
    older_fasta,
    strain_tree,
    tree_tips,
)
from panx_data.layout import BUCKET, Source, find_source
from panx_data.listings import load_listings
from panx_data.reader import DatasetFile, Job, OlderJob, job_files, overlay_modified, plan

LOG = logging.getLogger(__name__)
REPOSITORY = "https://github.com/neherlab/pan-genome-analysis"
DATA_ROOT = Path("out") / "data.pangenome.org"
RELEASES = "releases.json"
STAGING = ".staging"
RELEASE_ID = re.compile(r"\d{8}T\d{6}Z")
DATASET_ID = re.compile(r"[A-Za-z][A-Za-z0-9_]*")
CLUSTER_ID = re.compile(r"[A-Za-z0-9_]+")
PAGE_NAMES = frozenset({"tutorial", "pipeline", "pathogens", "index"})
IMMUTABLE = "public, max-age=31536000, immutable"
GENE_CLUSTER = "geneCluster"
PACK = "gene_clusters.tar.gz"
CORE_ALIGNMENTS = "core_gene_alignments.tar.gz"
ALL_ALIGNMENTS = "all_gene_alignments.tar.gz"


@dataclass(frozen=True)
class Recipe:
    """The definition of one release: `releases/<id>.toml` of data-tools, and the SHA-256 of its bytes."""

    id: str
    title: str
    snapshot: str
    datasets: tuple[str, ...]
    sha256: str


@dataclass(frozen=True)
class ToolVersion:
    """The commit of the build tool, and whether its folder had uncommitted changes."""

    commit: str
    uncommitted_changes: bool


@dataclass(frozen=True)
class StoredObject:
    """One object of a release folder, with the headers that the data host sends for it."""

    path: str
    size: int
    sha256: str
    # The field names of manifest.json
    contentType: str
    contentEncoding: str | None
    cacheControl: str
    contentDisposition: str | None


@dataclass(frozen=True)
class DatasetResult:
    """The counts and objects of one built dataset, and the gene-tree tips removed per cluster."""

    dataset: str
    strains: int
    clusters: int
    objects: tuple[StoredObject, ...]
    removed_tips: dict[str, int] = field(default_factory=dict)


def load_recipe(path: Path) -> Recipe:
    """Read and check a recipe file."""
    data = path.read_bytes()
    recipe = tomllib.loads(data.decode())
    if not RELEASE_ID.fullmatch(recipe["id"]) or path.stem != recipe["id"]:
        raise ValueError(f"'{path}': the id '{recipe['id']}' must be the file name and of the form YYYYMMDDThhmmssZ")
    datasets = tuple(recipe["datasets"])
    if len(set(datasets)) != len(datasets):
        raise ValueError(f"'{path}': a dataset is listed twice")
    for dataset in datasets:
        if not DATASET_ID.fullmatch(dataset) or dataset in PAGE_NAMES:
            raise ValueError(f"'{path}': '{dataset}' is no valid dataset id")
    return Recipe(recipe["id"], recipe["title"], recipe["snapshot"], datasets, hashlib.sha256(data).hexdigest())


def snapshot_checksum(snapshot: Path) -> str:
    """Return the SHA-256 of the snapshot's `src/MD5SUMS`, which names the content of every input."""
    return hashlib.sha256((snapshot / "src" / "MD5SUMS").read_bytes()).hexdigest()


def default_jobs() -> int:
    """Return the smaller of the CPU count minus 4 and the free memory in GB divided by 2, at least 1."""
    available = 0
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            available = int(line.split()[1]) // (1 << 20)
    return max(1, min((os.cpu_count() or 2) - 4, available // 2))


def release_build(
    recipe_path: Path,
    snapshot: Path,
    data_root: Path,
    jobs: int,
    tool: ToolVersion,
    python: str,
    zlib_version: str,
    only: frozenset[str] | None = None,
) -> Path:
    """Build the release of a recipe into the data root, and list it in the data root's `releases.json`.

    The datasets are built in a staging folder, which is renamed to the release id when every dataset is complete.
    """
    check_zlib()
    recipe = load_recipe(recipe_path)
    checksum = snapshot_checksum(snapshot)
    if checksum != recipe.snapshot:
        raise ValueError(f"the snapshot's src/MD5SUMS has the SHA-256 {checksum}, the recipe names {recipe.snapshot}")
    source = find_source(snapshot / "src")
    catalog = load_catalog(source.update_zip)
    selected = _selected(recipe, catalog, only)
    target = data_root / recipe.id
    if target.exists():
        raise FileExistsError(f"'{target}' exists; a release is written once")
    staging = data_root / STAGING
    if staging.exists():
        LOG.info("deleting the staging folder of an interrupted build, '%s'", staging)
        shutil.rmtree(staging)
    staging.mkdir(parents=True)
    results = _build_datasets(source, catalog, selected, staging, jobs)
    objects = [obj for name in selected for obj in results[name].objects]
    index = {
        "release": recipe.id,
        "categories": [
            {"id": category.id, "label": category.label}
            for category in catalog.categories
            if any(catalog.dataset(name).category == category.id for name in selected)
        ],
        "datasets": [_index_entry(catalog.dataset(name), results[name]) for name in selected],
    }
    objects.append(_write_object(staging, "index.json", json_bytes(index)))
    manifest = {
        "release": recipe.id,
        "build": {
            "repository": REPOSITORY,
            "commit": tool.commit,
            "uncommitted_changes": tool.uncommitted_changes,
            "python": python,
            "zlib": zlib_version,
            "snapshot": checksum,
            "recipe": recipe.sha256,
        },
        "objects": [obj.__dict__ for obj in sorted(objects, key=lambda obj: obj.path)],
    }
    _write_object(staging, "manifest.json", json_bytes(manifest))
    os.replace(staging, target)
    _list_release(data_root, recipe)
    removed = {name: len(result.removed_tips) for name, result in results.items() if result.removed_tips}
    for name, count in removed.items():
        LOG.info("%s: gene-tree tips without strain removed in %d clusters", name, count)
    return target


def build_dataset_release(
    job: Job | OlderJob, folder: Path, modified: Mapping[str, str | None], entry: CatalogDataset, older: Any | None
) -> DatasetResult:
    """Write the 8 files of one dataset into `folder`: `dataset.json`, `clusters.json`, `loci.json`, the pack, the two
    alignment archives, `strain_tree.nwk`, and `metadata.tsv`.

    `older` is the converted metadata configuration of an older dataset (`catalog/<dataset>.json`), None for a 2022
    dataset, whose `metaConfiguration.js` is read.
    """
    # The tree conversions recurse once per tree level, and gene trees of the older datasets are deeper than
    # Python's default limit
    sys.setrecursionlimit(100_000)
    files = {file.rel: file for file in job_files(job, modified)}
    dataset = job.dataset

    def read(name: str) -> bytes:
        if name not in files:
            raise ValueError(f"{dataset}: the pipeline output has no '{name}'")
        return files[name].content()

    if older is None:
        statements = js_values(read("metaConfiguration.js").decode())
        config = MetadataConfig(statements["meta_details"], statements["meta_display"])
        columns: list[dict[str, str]] = []
    else:
        config = MetadataConfig(older["meta_details"], older["meta_display"])
        columns = cluster_columns(older.get("new_columns_config", []))
    metadata, strains = metadata_and_strains(json.loads(read("strainMetainfo.json"))["data"], config)
    tree = strain_tree(json.loads(read("coreGenomeTree.json")))
    accessions = [strain["accession"] for strain in strains]
    _check_tips(dataset, "coreGenomeTree.json", tree_tips(tree), accessions)
    newick = read("strain_tree.nwk")
    _check_line_ends(dataset, "strain_tree.nwk", newick)
    _check_tips(dataset, "strain_tree.nwk", newick_tips(newick.decode()), accessions)
    rows = json.loads(read(f"{GENE_CLUSTER}.json"))
    context = _ClusterContext(dataset, files, frozenset(accessions), non_root_nodes(tree), older is not None)
    folder.mkdir(parents=True)
    clusters, loci, removed = _write_archives(folder, context, rows, [column["key"] for column in columns])
    document = {
        "id": dataset,
        "strain_tree": tree,
        "metadata": metadata,
        "strains": strains,
        "cluster_columns": columns,
        "accession_url": entry.accession_url,
    }
    objects = [
        _write_object(folder, "dataset.json", json_bytes(document), dataset),
        _write_object(folder, "clusters.json", json_bytes({"clusters": clusters}), dataset),
        _write_object(folder, "loci.json", json_bytes({"loci": loci}), dataset),
        _write_object(folder, "strain_tree.nwk", newick, dataset),
        _write_object(folder, "metadata.tsv", metadata_tsv(metadata, strains), dataset),
        *(_stored_file(folder, name, dataset) for name in (PACK, CORE_ALIGNMENTS, ALL_ALIGNMENTS)),
    ]
    return DatasetResult(dataset, len(strains), len(clusters), tuple(objects), removed)


def object_headers(path: str) -> tuple[str, str | None, str, str | None]:
    """Return Content-Type, Content-Encoding, Cache-Control, and Content-Disposition of a data path.

    The table "Object headers and bucket settings" of `docs/data-format.md` of the app repository.
    """
    folder, _, name = path.rpartition("/")
    dataset = folder.rpartition("/")[2]
    if path == RELEASES:
        return "application/json", "gzip", "public, max-age=300", None
    if name.endswith(".json"):
        return "application/json", "gzip", IMMUTABLE, None
    if name == "strain_tree.nwk":
        return "text/plain; charset=utf-8", "gzip", IMMUTABLE, f'attachment; filename="{dataset}_{name}"'
    if name == "metadata.tsv":
        return "text/tab-separated-values; charset=utf-8", "gzip", IMMUTABLE, f'attachment; filename="{dataset}_{name}"'
    if name.endswith(".tar.gz"):
        return "application/gzip", None, IMMUTABLE, f'attachment; filename="{dataset}_{name}"'
    raise ValueError(f"'{path}' is no path of the data format")


def cluster_entry_names(cluster_id: str, kind: str, has_tree: bool, has_newick: bool) -> list[str]:
    """Return the file names of a cluster's pack member, in member order."""
    names = [f"{cluster_id}_aa_aln.fa", f"{cluster_id}_aa_aln_reduced.fa"] if kind == "protein" else []
    names += [f"{cluster_id}_na_aln.fa", f"{cluster_id}_na_aln_reduced.fa"]
    names += [f"{cluster_id}_tree.json"] if has_tree else []
    names += [f"{cluster_id}.nwk"] if has_newick else []
    return [*names, f"{cluster_id}_patterns.json"]


@dataclass(frozen=True)
class _ClusterContext:
    dataset: str
    files: Mapping[str, DatasetFile]
    strains: frozenset[str]
    patterns_length: int
    older: bool


def _write_archives(
    folder: Path, context: _ClusterContext, rows: list[dict[str, Any]], extra_keys: list[str]
) -> tuple[list[dict[str, Any]], list[str], dict[str, int]]:
    """Write the pack and the two alignment archives in the order of the cluster table; return the clusters, the
    loci (the record identifiers of each cluster, joined by single spaces), and the removed gene-tree tips per
    cluster."""
    clusters: list[dict[str, Any]] = []
    loci: list[str] = []
    removed_tips: dict[str, int] = {}
    seen: set[str] = set()
    with (
        (folder / PACK).open("wb") as pack,
        (folder / CORE_ALIGNMENTS).open("wb") as core_file,
        (folder / ALL_ALIGNMENTS).open("wb") as all_file,
    ):
        core, every = GzipWriter(core_file), GzipWriter(all_file)
        for row in rows:
            cluster_id = row["msa"]
            if not CLUSTER_ID.fullmatch(cluster_id) or cluster_id in seen:
                raise ValueError(f"{context.dataset}: cluster id '{cluster_id}' is invalid or repeated")
            seen.add(cluster_id)
            entries, ids, removed = _cluster_entries(context, cluster_id)
            kind = "protein" if f"{cluster_id}_aa_aln.fa" in entries else "rna"
            start = pack.tell()
            pack.write(_member(context.dataset, entries))
            cluster = cluster_entry(row, kind, extra_keys)
            clusters.append(_with_range(cluster, start, pack.tell() - 1))
            loci.append(" ".join(ids))
            if removed:
                removed_tips[cluster_id] = removed
            if kind == "protein":
                single_copy = cluster["strains"] == len(context.strains) and not cluster["duplications"]
                for suffix in ("_aa_aln.fa", "_na_aln.fa"):
                    name = f"{cluster_id}{suffix}"
                    data, mtime = entries[name]
                    every.write(tar_entry(f"all_genes/{name}", data, mtime))
                    if single_copy:
                        core.write(tar_entry(f"core_genes/{name}", data, mtime))
        pack.write(gzip_member(TAR_END))
        for writer in (core, every):
            writer.write(TAR_END)
            writer.close()
    return clusters, loci, removed_tips


def _cluster_entries(context: _ClusterContext, cluster_id: str) -> tuple[dict[str, tuple[bytes, int]], list[str], int]:
    """Return the entries of a cluster's member by name, in member order, with their mtimes; the record identifiers of
    its full nucleotide alignment; and the number of removed gene-tree tips.

    A gene tree whose tips all name no strain is left out: the pipeline output of some older clusters holds the
    tree of another cluster."""

    def file(name: str) -> DatasetFile | None:
        return context.files.get(f"{GENE_CLUSTER}/{name}")

    def required(name: str) -> DatasetFile:
        found = file(name)
        if found is None:
            raise ValueError(f"{context.dataset}: cluster '{cluster_id}' has no '{name}'")
        return found

    protein = file(f"{cluster_id}_aa_aln.fa") is not None
    tree, newick = file(f"{cluster_id}_tree.json"), file(f"{cluster_id}.nwk")
    names = cluster_entry_names(cluster_id, "protein" if protein else "rna", tree is not None, newick is not None)
    entries: dict[str, tuple[bytes, int]] = {}
    removed = 0
    full_ids: list[str] = []
    for name in names:
        source = required(name)
        data = source.content()
        if name.endswith(".fa"):
            data, ids = _alignment(context, name, data)
            if name.endswith("_na_aln.fa"):
                full_ids = ids
            if name.endswith("_reduced.fa"):
                _check_reduced(context.dataset, name, ids, entries[name.replace("_reduced.fa", ".fa")][0])
        elif name.endswith((".nwk", "_tree.json")):
            tree_data, count = (_gene_tree_newick if name.endswith(".nwk") else _gene_tree_json)(context, name, data)
            removed = max(removed, count)
            if tree_data is None:
                LOG.warning("%s: '%s' has no tip of a strain and is left out", context.dataset, name)
                continue
            data = tree_data
        else:
            patterns = json.loads(data)["patterns"]
            if len(patterns) != context.patterns_length:
                raise ValueError(
                    f"{context.dataset}: '{name}' has {len(patterns)} patterns, the strain tree has "
                    f"{context.patterns_length} non-root nodes"
                )
        entries[name] = (data, source.mtime)
    return entries, full_ids, removed


def _alignment(context: _ClusterContext, name: str, data: bytes) -> tuple[bytes, list[str]]:
    _check_line_ends(context.dataset, name, data)
    text = data.decode()
    if context.older:
        text = older_fasta(text)
    ids = fasta_ids(text)
    for identifier in ids:
        if identifier != CONSENSUS and identifier.split("|", 1)[0] not in context.strains:
            raise ValueError(f"{context.dataset}: record '{identifier}' of '{name}' names no strain")
    return (text.encode() if context.older else data), ids


def _check_reduced(dataset: str, name: str, ids: list[str], full: bytes) -> None:
    records = ids[1:] if ids[:1] == [CONSENSUS] else ids
    if sorted(records) != sorted(fasta_ids(full.decode())):
        raise ValueError(f"{dataset}: '{name}' holds other records than its full alignment")


def _gene_tree_json(context: _ClusterContext, name: str, data: bytes) -> tuple[bytes | None, int]:
    original = json.loads(data)
    tree, removed = gene_tree_json(original, context.strains.__contains__)
    if tree is None:
        return None, removed
    return (data if tree == original else json.dumps(tree).encode()), removed


def _gene_tree_newick(context: _ClusterContext, name: str, data: bytes) -> tuple[bytes | None, int]:
    _check_line_ends(context.dataset, name, data)
    text, removed = gene_tree_newick(data.decode(), context.strains.__contains__)
    if text is None:
        return None, removed
    return (text.encode() if removed else data), removed


def _member(dataset: str, entries: Mapping[str, tuple[bytes, int]]) -> bytes:
    tar = b"".join(
        tar_entry(f"{dataset}/{GENE_CLUSTER}/{name}", data, mtime) for name, (data, mtime) in entries.items()
    )
    return gzip_member(tar)


def _with_range(cluster: dict[str, Any], first: int, last: int) -> dict[str, Any]:
    """Insert `range` before `extra`, in the field order of `clusters.json`."""
    extra = cluster.pop("extra")
    return {**cluster, "range": [first, last], "extra": extra}


def _check_tips(dataset: str, name: str, tips: list[str], accessions: list[str]) -> None:
    if len(set(accessions)) != len(accessions):
        raise ValueError(f"{dataset}: strainMetainfo.json names an accession twice")
    if len(set(tips)) != len(tips) or set(tips) != set(accessions):
        raise ValueError(f"{dataset}: the tips of '{name}' differ from the strain accessions")


def _check_line_ends(dataset: str, name: str, data: bytes) -> None:
    if b"\r" in data:
        raise ValueError(f"{dataset}: '{name}' has a carriage return; text files have LF line ends")


def _write_object(folder: Path, name: str, content: bytes, dataset: str | None = None) -> StoredObject:
    """Write a text object gzip-compressed and return its manifest entry."""
    stored = gzip_member(content)
    (folder / name).write_bytes(stored)
    return _object(f"{dataset}/{name}" if dataset else name, stored)


def _stored_file(folder: Path, name: str, dataset: str) -> StoredObject:
    digest = hashlib.sha256()
    with (folder / name).open("rb") as stored:
        for block in iter(lambda: stored.read(1 << 20), b""):
            digest.update(block)
    return _manifest_entry(f"{dataset}/{name}", (folder / name).stat().st_size, digest.hexdigest())


def _object(path: str, stored: bytes) -> StoredObject:
    return _manifest_entry(path, len(stored), hashlib.sha256(stored).hexdigest())


def _manifest_entry(path: str, size: int, sha256: str) -> StoredObject:
    return StoredObject(path, size, sha256, *object_headers(path))


def _index_entry(entry: CatalogDataset, result: DatasetResult) -> dict[str, Any]:
    return {
        "id": entry.id,
        "name": entry.name,
        "label": entry.label,
        "category": entry.category,
        "source": entry.source,
        "strains": result.strains,
        "clusters": result.clusters,
    }


def _selected(recipe: Recipe, catalog: Catalog, only: frozenset[str] | None) -> list[str]:
    """Return the datasets to build in catalog order: those of the recipe, or of `only`, which the recipe must list."""
    listed = set(recipe.datasets)
    unknown = listed - {dataset.id for dataset in catalog.datasets}
    if unknown:
        raise ValueError(f"the recipe lists datasets that the catalog lacks: {' '.join(sorted(unknown))}")
    outside = (only or frozenset()) - listed
    if outside:
        raise ValueError(f"not datasets of the recipe: {' '.join(sorted(outside))}")
    wanted = only or listed
    return [dataset.id for dataset in catalog.datasets if dataset.id in wanted]


def _build_datasets(
    source: Source, catalog: Catalog, selected: list[str], staging: Path, workers: int
) -> dict[str, DatasetResult]:
    listings = load_listings(source.listings, BUCKET)
    jobs = plan(source, listings, frozenset(selected))
    missing = set(selected) - {job.dataset for job in jobs}
    if missing:
        raise ValueError(f"the snapshot holds no files of: {' '.join(sorted(missing))}")
    results: dict[str, DatasetResult] = {}
    failures: list[str] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(
                build_dataset_release,
                job,
                staging / job.dataset,
                overlay_modified(job, listings),
                catalog.dataset(job.dataset),
                _older_config(job),
            ): job.dataset
            for job in jobs
        }
        for done, future in enumerate(as_completed(futures), start=1):
            name = futures[future]
            try:
                result = future.result()
            except Exception:  # report each failed dataset, finish the others, fail at the end
                LOG.exception("[%d/%d] %s failed", done, len(jobs), name)
                failures.append(name)
                continue
            results[name] = result
            LOG.info("[%d/%d] %s: %d strains, %d clusters", done, len(jobs), name, result.strains, result.clusters)
    if failures:
        raise RuntimeError(f"{len(failures)} datasets failed: {' '.join(sorted(failures))}")
    return results


def _older_config(job: Job | OlderJob) -> Any | None:
    if not isinstance(job, OlderJob):
        return None
    return json.loads((CATALOG_CONFIG.parent / f"{job.dataset}.json").read_text())


def _list_release(data_root: Path, recipe: Recipe) -> None:
    """Add the release to `releases.json`: create the file with the release as default, or add the entry and keep
    the other entries and the default."""
    path = data_root / RELEASES
    if path.exists():
        current = json.loads(gzip.decompress(path.read_bytes()))
        default, entries = current["default"], [entry for entry in current["releases"] if entry["id"] != recipe.id]
    else:
        default, entries = recipe.id, []
    entries.append({"id": recipe.id, "title": recipe.title})
    entries.sort(key=lambda entry: entry["id"], reverse=True)
    temporary = path.with_name(f".{RELEASES}.tmp")
    temporary.write_bytes(gzip_member(json_bytes({"default": default, "releases": entries})))
    os.replace(temporary, path)
