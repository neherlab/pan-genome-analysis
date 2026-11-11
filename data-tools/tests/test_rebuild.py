import gzip
import hashlib
import io
import shutil
import tarfile
from pathlib import Path

import pytest

from panx_data.layout import Out
from panx_data.rebuild import Job, Upload, build_dataset, choose_gz

ORIGINAL = b"original-gz-bytes"
SAME_SIZE_COPY = b"new-gz-copy-bytes"
OTHER_SIZE_COPY = b"a longer new gz copy"


def test_rebuild_august_upload_keeps_the_original_of_the_same_size() -> None:
    assert choose_gz(SAME_SIZE_COPY, ORIGINAL, Upload.AUGUST_2022, None) is ORIGINAL


def test_rebuild_august_upload_takes_the_copy_of_another_size() -> None:
    assert choose_gz(OTHER_SIZE_COPY, ORIGINAL, Upload.AUGUST_2022, None) is OTHER_SIZE_COPY


def test_rebuild_november_upload_always_takes_the_copy() -> None:
    assert choose_gz(SAME_SIZE_COPY, ORIGINAL, Upload.NOVEMBER_2022, None) is SAME_SIZE_COPY


def test_rebuild_listed_etag_overrides_the_rule() -> None:
    etag = hashlib.md5(SAME_SIZE_COPY).hexdigest()

    assert choose_gz(SAME_SIZE_COPY, ORIGINAL, Upload.AUGUST_2022, etag) is SAME_SIZE_COPY


def test_rebuild_etag_of_neither_version_falls_back_to_the_rule() -> None:
    assert choose_gz(SAME_SIZE_COPY, ORIGINAL, Upload.AUGUST_2022, "0" * 32) is ORIGINAL


@pytest.mark.skipif(shutil.which("pigz") is None, reason="needs the pigz program")
def test_rebuild_dataset_gets_plain_files_gz_copies_zips_and_s3_objects(tmp_path: Path) -> None:
    alignment = b">strain_1\nMKV\n>strain_2\nMKI\n"
    tarball = _tarball(
        tmp_path / "Testus.tar.gz",
        {
            "Testus/geneCluster.json": b'[{"msa": "abc_refined"}]',
            "Testus/metainfo.tsv": b"old layout\n",
            "Testus/geneCluster/abc_refined_aa_aln.fa.gz": gzip.compress(alignment, mtime=0),
        },
    )
    zips = tmp_path / "zips"
    zips.mkdir()
    (zips / "all_gene_alignments.zip").write_bytes(b"zip-all")
    (zips / "core_gene_alignments.zip").write_bytes(b"zip-core")
    overlay = tmp_path / "overlay"
    overlay.mkdir()
    (overlay / "metainfo.tsv").write_bytes(b"new layout\n")
    out = Out(tmp_path / "out")
    for folder in (out.bucket / "dataset", out.staging, out.parts):
        folder.mkdir(parents=True)

    result = build_dataset(Job("Testus", Upload.AUGUST_2022, tarball, zips, None, overlay), out, {})

    folder = out.bucket / "dataset" / "Testus"
    files = sorted(path.relative_to(folder).as_posix() for path in folder.rglob("*") if path.is_file())
    assert files == [
        "all_gene_alignments.zip",
        "core_gene_alignments.zip",
        "geneCluster.json",
        "geneCluster.json.gz",
        "geneCluster/abc_refined_aa_aln.fa",
        "geneCluster/abc_refined_aa_aln.fa.gz",
        "metainfo.tsv",
        "metainfo.tsv.gz",
    ]
    assert (folder / "geneCluster/abc_refined_aa_aln.fa").read_bytes() == alignment
    assert gzip.decompress((folder / "geneCluster/abc_refined_aa_aln.fa.gz").read_bytes()) == alignment
    assert (folder / "metainfo.tsv").read_bytes() == b"new layout\n"
    manifest = dict(line.split("\t", 1) for line in (out.parts / "Testus.tsv").read_text().splitlines())
    assert manifest["dataset/Testus/metainfo.tsv"].endswith("\ts3")
    assert manifest["dataset/Testus/geneCluster.json"].endswith("\tscicore")
    assert (result.files, result.from_s3) == (8, 1)


def _tarball(path: Path, files: dict[str, bytes]) -> Path:
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 1_657_000_000
            tar.addfile(info, io.BytesIO(data))
    return path
