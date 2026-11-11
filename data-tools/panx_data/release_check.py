"""Check the releases of a data root against the pipeline output of the snapshot.

The format verifier of the app repository checks the rules of the format; this check compares a release with its
inputs: every row of each cluster table, the dataset counts, the gene trees of each cluster, and the metadata tables
that S3 rewrote. It also lists what the conversion left without a number or removed, for a reader to accept.
"""

import csv
import gzip
import io
import json
import logging
import sys
import tarfile
from collections.abc import Iterator, Mapping
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

from panx_data.convert import UNKNOWN, cluster_entry, gene_tree_newick
from panx_data.layout import BUCKET, find_source
from panx_data.listings import load_listings
from panx_data.reader import Job, OlderJob, Origin, job_files, overlay_modified, plan
from panx_data.release import GENE_CLUSTER, PACK, RELEASES, cluster_entry_names

LOG = logging.getLogger(__name__)


class Kind(Enum):
    """Kinds of findings; a finding of a failing kind fails the check."""

    CLUSTER_ROW = ("cluster-row", True)
    COUNT = ("count", True)
    PACK = ("pack", True)
    METADATA_DIFFERENCE = ("metadata-difference", True)
    NO_NUMBER = ("no-number", False)
    REMOVED_TIPS = ("removed-tips", False)
    REMOVED_TREE = ("removed-tree", False)

    @property
    def label(self) -> str:
        return self.value[0]

    @property
    def fails(self) -> bool:
        return self.value[1]


@dataclass(frozen=True)
class Finding:
    """One finding of the check, in one dataset of one release."""

    release: str
    dataset: str
    kind: Kind
    detail: str


def release_check(data_root: Path, snapshot: Path, workers: int) -> list[Finding]:
    """Check every release that `releases.json` of the data root lists; return the findings in dataset order."""
    releases = json.loads(gzip.decompress((data_root / RELEASES).read_bytes()))
    source = find_source(snapshot / "src")
    listings = load_listings(source.listings, BUCKET)
    findings: list[Finding] = []
    for release in releases["releases"]:
        folder = data_root / release["id"]
        index = json.loads(gzip.decompress((folder / "index.json").read_bytes()))
        names = [entry["id"] for entry in index["datasets"]]
        jobs = {job.dataset: job for job in plan(source, listings, frozenset(names))}
        results: dict[str, list[Finding]] = {}
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(
                    check_dataset,
                    release["id"],
                    jobs[entry["id"]],
                    folder,
                    entry,
                    overlay_modified(jobs[entry["id"]], listings),
                ): entry["id"]
                for entry in index["datasets"]
            }
            for done, future in enumerate(as_completed(futures), start=1):
                name = futures[future]
                results[name] = future.result()
                failing = sum(finding.kind.fails for finding in results[name])
                LOG.info("[%d/%d] %s/%s: %d failing findings", done, len(futures), release["id"], name, failing)
        findings += [finding for name in names for finding in results[name]]
    return findings


def check_dataset(
    release: str, job: Job | OlderJob, folder: Path, entry: Mapping[str, Any], modified: Mapping[str, str | None]
) -> list[Finding]:
    """Compare one dataset of a release with its pipeline output."""
    # The tree conversions recurse once per tree level, deeper than Python's default limit
    sys.setrecursionlimit(100_000)
    dataset_folder = folder / job.dataset
    document = _read_json(dataset_folder / "dataset.json")
    clusters = _read_json(dataset_folder / "clusters.json")["clusters"]
    strains = frozenset(strain["accession"] for strain in document["strains"])
    found: list[Finding] = []

    def report(kind: Kind, detail: str) -> None:
        found.append(Finding(release, job.dataset, kind, detail))

    rows: list[dict[str, Any]] = []
    cluster_files: set[str] = set()
    metainfo: bytes | None = None
    removed: dict[str, int] = {}
    dropped: set[str] = set()
    for file in job_files(job, modified):
        if file.rel == f"{GENE_CLUSTER}.json":
            rows = json.loads(file.content())
        elif file.rel == "metainfo.tsv" and file.origin is Origin.S3 and isinstance(job, Job):
            metainfo = file.content()
        elif file.rel.startswith(f"{GENE_CLUSTER}/"):
            name = file.rel.removeprefix(f"{GENE_CLUSTER}/")
            cluster_files.add(name)
            if name.endswith("_tree.json"):
                outside, tips = _tips_outside(json.loads(file.content()), strains)
                if outside:
                    removed[name.removesuffix("_tree.json")] = outside
                if outside == tips:
                    dropped.add(name)
            elif name.endswith(".nwk"):
                text, outside = gene_tree_newick(file.content().decode(), strains.__contains__)
                if text is None:
                    dropped.add(name)

    for key, count in (("strains", len(document["strains"])), ("clusters", len(clusters))):
        if entry[key] != count:
            report(Kind.COUNT, f"index.json gives {key} {entry[key]}, the dataset files {count}")
    extra_keys = [column["key"] for column in document["cluster_columns"]]
    found += [
        Finding(release, job.dataset, Kind.CLUSTER_ROW, detail)
        for detail in _row_differences(rows, clusters, cluster_files, extra_keys)
    ]
    found += [
        Finding(release, job.dataset, Kind.PACK, detail)
        for detail in _pack_differences(job.dataset, dataset_folder / PACK, clusters, cluster_files - dropped)
    ]
    if metainfo is not None:
        metadata = (dataset_folder / "metadata.tsv").read_bytes()
        found += [
            Finding(release, job.dataset, Kind.METADATA_DIFFERENCE, detail)
            for detail in metadata_differences(metainfo, gzip.decompress(metadata))
        ]
    for column in document["metadata"]:
        if column["type"] != "continuous":
            continue
        values = {strain["values"][column["key"]] for strain in document["strains"]} - {None}
        missing = sorted(values - set(column["numbers"]))
        if missing:
            report(Kind.NO_NUMBER, f"{column['key']}: {json.dumps(missing, ensure_ascii=False)}")
    # Files that no cluster of the table names are no input of the release
    ids = {cluster["id"] for cluster in clusters}
    for cluster_id, count in sorted(removed.items()):
        if cluster_id in ids:
            report(Kind.REMOVED_TIPS, f"{cluster_id}: {count} gene-tree tips without strain")
    for name in sorted(dropped):
        if name.removesuffix("_tree.json").removesuffix(".nwk") in ids:
            report(Kind.REMOVED_TREE, f"{name}: no tip of a strain, left out of the pack")
    return found


def metadata_differences(metainfo: bytes, metadata: bytes) -> Iterator[str]:
    """Compare the metadata table that S3 rewrote with `metadata.tsv`, by accession and column name; `unknown` in
    the S3 table equals an empty field."""
    theirs, ours = _table(metainfo), _table(metadata)
    their_columns = next(iter(theirs.values()), {}).keys() if theirs else []
    our_columns = next(iter(ours.values()), {}).keys() if ours else []
    for column in their_columns:
        if column not in our_columns:
            yield f"column '{column}' of the S3 metainfo.tsv is missing"
    for column in our_columns:
        if column not in their_columns:
            yield f"column '{column}' is not in the S3 metainfo.tsv"
    for accession in sorted(set(theirs) ^ set(ours)):
        yield f"accession '{accession}' is in one table only"
    for accession in sorted(set(theirs) & set(ours)):
        for column in their_columns:
            if column not in our_columns:
                continue
            their_value = theirs[accession][column]
            their_value = "" if their_value == UNKNOWN else their_value
            if their_value != ours[accession][column]:
                yield f"{accession} '{column}': S3 '{their_value}', release '{ours[accession][column]}'"


def _table(data: bytes) -> dict[str, dict[str, str]]:
    rows = list(csv.reader(io.StringIO(data.decode()), delimiter="\t"))
    header = rows[0]
    return {row[header.index("accession")]: dict(zip(header, row)) for row in rows[1:] if row}


def _row_differences(
    rows: list[dict[str, Any]], clusters: list[dict[str, Any]], files: set[str], extra_keys: list[str]
) -> Iterator[str]:
    if [row["msa"] for row in rows] != [cluster["id"] for cluster in clusters]:
        yield f"geneCluster.json has {len(rows)} rows, clusters.json {len(clusters)} clusters, or another order"
        return
    for row, cluster in zip(rows, clusters):
        kind = "protein" if f"{row['msa']}_aa_aln.fa" in files else "rna"
        expected = cluster_entry(row, kind, extra_keys)
        actual = {key: value for key, value in cluster.items() if key != "range"}
        if expected != actual:
            yield f"{row['msa']}: expected {json.dumps(expected)}, found {json.dumps(actual)}"


def _pack_differences(dataset: str, pack: Path, clusters: list[dict[str, Any]], files: set[str]) -> Iterator[str]:
    """Check that each member holds the entries that the pipeline output gives: `_tree.json` and `.nwk` exactly when
    the cluster has them."""
    with pack.open("rb") as stored:
        for cluster in clusters:
            first, last = cluster["range"]
            stored.seek(first)
            member = gzip.decompress(stored.read(last - first + 1))
            with tarfile.open(fileobj=io.BytesIO(member), mode="r:") as tar:
                names = tar.getnames()
            cluster_id = cluster["id"]
            expected = cluster_entry_names(
                cluster_id, cluster["kind"], f"{cluster_id}_tree.json" in files, f"{cluster_id}.nwk" in files
            )
            if names != [f"{dataset}/{GENE_CLUSTER}/{name}" for name in expected]:
                yield f"{cluster_id}: member holds {names}, the pipeline output gives {expected}"


def _tips_outside(tree: Mapping[str, Any], strains: frozenset[str]) -> tuple[int, int]:
    """Return the number of tips whose accession names no strain, and the number of tips."""
    accessions: list[str] = []
    stack = [tree]
    while stack:
        node = stack.pop()
        if node.get("children"):
            stack.extend(node["children"])
        else:
            accessions.append(node["accession"])
    return sum(accession not in strains for accession in accessions), len(accessions)


def _read_json(path: Path) -> Any:
    return json.loads(gzip.decompress(path.read_bytes()))
