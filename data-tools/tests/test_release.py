import gzip
import io
import json
import shutil
import tarfile
from pathlib import Path

import pytest

from panx_data.catalog import NCBI_NUCCORE, CatalogDataset
from panx_data.reader import Job, Upload
from panx_data.release import build_dataset_release, object_headers

ENTRY = CatalogDataset("Testus", "Testus", "Testus (2)", "refseq", "RefSeq", NCBI_NUCCORE)
META = (
    'var meta_details={"country":["CH","unknown"]},\n'
    'meta_display={"color_options":{"country":{"menuItem":"country","type":"discrete"}},'
    '"meta_display_order":["country"]};\nvar association_columns=[];\n'
)
ROWS = [
    {
        "geneId": 1,
        "geneLen": 6,
        "count": 2,
        "dupli": "no",
        "dup_detail": "None",
        "ann": "kinase",
        "msa": "p1_refined",
        "divers": 0.0,
        "event": "0",
        "allAnn": "kinase#2",
        "GName": "aroK",
        "allGName": "aroK#2",
        "locus": "A_1 B_1",
    },
    {
        "geneId": 2,
        "geneLen": 4,
        "count": 1,
        "dupli": "no",
        "dup_detail": "None",
        "ann": "tRNA",
        "msa": "RNA_r1",
        "divers": None,
        "event": "1",
        "allAnn": "tRNA#1",
        "GName": "None",
        "allGName": "None#1",
        "locus": "A_2",
    },
]
CORE_TREE = {"name": "N0", "children": [{"name": "A", "branch_length": 0.1}, {"name": "B", "branch_length": 0.2}]}
FILES = {
    "Testus/geneCluster.json": json.dumps(ROWS).encode(),
    "Testus/strainMetainfo.json": json.dumps(
        {
            "data": [
                {"accession": "A", "strain": "a", "country": "CH"},
                {"accession": "B", "strain": "b", "country": "unknown"},
            ]
        }
    ).encode(),
    "Testus/metaConfiguration.js": META.encode(),
    "Testus/coreGenomeTree.json": json.dumps(CORE_TREE).encode(),
    "Testus/strain_tree.nwk": b"(A:0.1,B:0.2)N0:0;\n",
    "Testus/geneCluster/p1_refined_aa_aln.fa.gz": gzip.compress(b">A|1\nMK\n>B|1\nMK\n"),
    "Testus/geneCluster/p1_refined_aa_aln_reduced.fa.gz": gzip.compress(b">consensus\nMK\n>A|1\n..\n>B|1\n..\n"),
    "Testus/geneCluster/p1_refined_na_aln.fa.gz": gzip.compress(b">A|1\nATGAAA\n>B|1\nATGAAA\n"),
    "Testus/geneCluster/p1_refined_na_aln_reduced.fa.gz": gzip.compress(
        b">consensus\nATGAAA\n>A|1\n......\n>B|1\n......\n"
    ),
    "Testus/geneCluster/p1_refined.nwk": b"(A|1:0.1,B|1:0.1):0.0;\n",
    "Testus/geneCluster/p1_refined_patterns.json": b'{"patterns": "33"}',
    "Testus/geneCluster/RNA_r1_na_aln.fa.gz": gzip.compress(b">A|2\nACGT\n"),
    "Testus/geneCluster/RNA_r1_na_aln_reduced.fa.gz": gzip.compress(b">consensus\nACGT\n>A|2\n....\n"),
    "Testus/geneCluster/RNA_r1_tree.json": b'{"name": "A|2", "accession": "A", "branch_length": 0.0}',
    "Testus/geneCluster/RNA_r1_patterns.json": b'{"patterns": "30"}',
    "Testus/geneCluster/orphan_na_aln.fa.gz": gzip.compress(b">A|9\nAC\n"),
}
NAMES = [
    "all_gene_alignments.tar.gz",
    "clusters.json",
    "core_gene_alignments.tar.gz",
    "dataset.json",
    "gene_clusters.tar.gz",
    "loci.json",
    "metadata.tsv",
    "strain_tree.nwk",
]


@pytest.mark.skipif(shutil.which("pigz") is None, reason="needs the pigz program")
def test_release_dataset_writes_the_files_of_format_1_with_the_same_bytes_twice(tmp_path: Path) -> None:
    tarball = _tarball(tmp_path / "Testus.tar.gz", FILES)
    job = Job("Testus", Upload.AUGUST_2022, tarball, None, None, tmp_path / "no-s3")

    result = build_dataset_release(job, tmp_path / "a", {}, ENTRY, None)
    build_dataset_release(job, tmp_path / "b", {}, ENTRY, None)

    folder = tmp_path / "a"
    assert sorted(path.name for path in folder.iterdir()) == NAMES
    assert all((folder / name).read_bytes() == (tmp_path / "b" / name).read_bytes() for name in NAMES)
    assert (result.strains, result.clusters) == (2, 2)
    clusters = _json(folder / "clusters.json")["clusters"]
    assert [(c["id"], c["kind"], c["diversity"]) for c in clusters] == [
        ("p1_refined", "protein", 0.0),
        ("RNA_r1", "rna", None),
    ]
    assert _json(folder / "loci.json") == {"format": 1, "loci": [["A|1", "B|1"], ["A|2"]]}
    pack = (folder / "gene_clusters.tar.gz").read_bytes()
    first, last = clusters[1]["range"]
    assert clusters[0]["range"] == [0, first - 1]
    with tarfile.open(fileobj=io.BytesIO(pack[first : last + 1]), mode="r:gz") as tar:
        assert [(info.name, info.mtime) for info in tar] == [
            ("Testus/geneCluster/RNA_r1_na_aln.fa", 1_657_000_000),
            ("Testus/geneCluster/RNA_r1_na_aln_reduced.fa", 1_657_000_000),
            ("Testus/geneCluster/RNA_r1_tree.json", 1_657_000_000),
            ("Testus/geneCluster/RNA_r1_patterns.json", 1_657_000_000),
        ]
    with tarfile.open(fileobj=io.BytesIO(pack), mode="r:gz") as tar:
        assert len(tar.getnames()) == 10
    with tarfile.open(folder / "core_gene_alignments.tar.gz", mode="r:gz") as tar:
        assert tar.getnames() == ["core_genes/p1_refined_aa_aln.fa", "core_genes/p1_refined_na_aln.fa"]
    assert gzip.decompress((folder / "metadata.tsv").read_bytes()) == b"accession\tstrain\tcountry\nA\ta\tCH\nB\tb\t\n"
    assert {obj.path for obj in result.objects} == {f"Testus/{name}" for name in NAMES}


@pytest.mark.skipif(shutil.which("pigz") is None, reason="needs the pigz program")
def test_release_dataset_fails_on_a_record_of_no_strain(tmp_path: Path) -> None:
    files = {**FILES, "Testus/geneCluster/RNA_r1_na_aln.fa.gz": gzip.compress(b">Z|2\nACGT\n")}
    job = Job("Testus", Upload.AUGUST_2022, _tarball(tmp_path / "Testus.tar.gz", files), None, None, tmp_path / "x")

    with pytest.raises(ValueError, match="names no strain"):
        build_dataset_release(job, tmp_path / "a", {}, ENTRY, None)


@pytest.mark.parametrize(
    ("path", "headers"),
    [
        ("releases.json", ("application/json", "gzip", "public, max-age=300", None)),
        ("index.json", ("application/json", "gzip", "public, max-age=31536000, immutable", None)),
        (
            "E_coli/metadata.tsv",
            (
                "text/tab-separated-values; charset=utf-8",
                "gzip",
                "public, max-age=31536000, immutable",
                'attachment; filename="E_coli_metadata.tsv"',
            ),
        ),
        (
            "E_coli/gene_clusters.tar.gz",
            (
                "application/gzip",
                None,
                "public, max-age=31536000, immutable",
                'attachment; filename="E_coli_gene_clusters.tar.gz"',
            ),
        ),
    ],
)
def test_object_headers_follow_the_table_of_the_format(path: str, headers: tuple) -> None:
    assert object_headers(path) == headers


def _json(path: Path) -> dict:
    return json.loads(gzip.decompress(path.read_bytes()))


def _tarball(path: Path, files: dict[str, bytes]) -> Path:
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 1_657_000_000
            tar.addfile(info, io.BytesIO(data))
    return path


@pytest.mark.skipif(shutil.which("pigz") is None, reason="needs the pigz program")
def test_release_dataset_leaves_out_a_gene_tree_without_strain_tips(tmp_path: Path) -> None:
    files = {**FILES, "Testus/geneCluster/RNA_r1_tree.json": b'{"name": "Z|2", "accession": "Z", "branch_length": 0.0}'}
    job = Job("Testus", Upload.AUGUST_2022, _tarball(tmp_path / "Testus.tar.gz", files), None, None, tmp_path / "x")

    result = build_dataset_release(job, tmp_path / "a", {}, ENTRY, None)

    first, last = _json(tmp_path / "a" / "clusters.json")["clusters"][1]["range"]
    member = (tmp_path / "a" / "gene_clusters.tar.gz").read_bytes()[first : last + 1]
    with tarfile.open(fileobj=io.BytesIO(member), mode="r:gz") as tar:
        assert "Testus/geneCluster/RNA_r1_tree.json" not in tar.getnames()
    assert result.removed_tips == {"RNA_r1": 1}
