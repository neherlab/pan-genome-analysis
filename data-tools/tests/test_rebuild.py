import gzip
import hashlib
import io
import shutil
import tarfile
from pathlib import Path

import pytest

from panx_data.layout import Out
from panx_data.listings import S3Object
from panx_data.pigz_gzip import pigz_bytes
from panx_data.reader import Job, OlderJob, Origin, Upload, older_dataset_files
from panx_data.rebuild import build_dataset, build_older_dataset, choose_gz

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

    modified = {"dataset/Testus/metainfo.tsv": "2022-09-09T10:00:00+00:00"}

    result = build_dataset(Job("Testus", Upload.AUGUST_2022, tarball, zips, None, overlay), out, {}, modified)

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
    assert gzip.decompress((folder / "metainfo.tsv.gz").read_bytes()) == b"new layout\n"
    manifest = dict(line.split("\t", 1) for line in (out.parts / "Testus.tsv").read_text().splitlines())
    assert manifest["dataset/Testus/metainfo.tsv"].endswith("\ts3")
    assert manifest["dataset/Testus/geneCluster.json"].endswith("\tscicore")
    assert (result.files, result.from_s3) == (8, 1)


def test_rebuild_older_dataset_files_prefer_s3_objects_and_skip_working_folders(tmp_path: Path) -> None:
    old, new = b">a|1\nMKV\n", b">a|1\nMKI\n"
    tarball = _tarball(
        tmp_path / "Tus.tar.gz",
        {
            "Tus/vis/geneCluster.json": b"[]",
            "Tus/vis/geneCluster/a_aa_aln.fa.gz": gzip.compress(old, mtime=0),
            "Tus/vis/geneCluster/npy_bk/a.npy": b"working file",
            "Tus/input_GenBank/a.gbk": b"pipeline input",
        },
    )
    overlay = tmp_path / "s3" / "dataset" / "Testus"
    (overlay / "geneCluster").mkdir(parents=True)
    (overlay / "geneCluster" / "a_aa_aln.fa.gz").write_bytes(gzip.compress(new, mtime=0))
    (overlay / "newColumnConfig.js").write_bytes(b"var x = 1;")
    modified = {
        "dataset/Testus/geneCluster/a_aa_aln.fa.gz": "2022-08-25T15:44:03+00:00",
        "dataset/Testus/newColumnConfig.js": "2022-08-23T08:37:49+00:00",
    }

    files = list(older_dataset_files(OlderJob("Testus", tarball, "Tus", overlay), modified))

    by_name = {file.rel: file for file in files}
    assert sorted(by_name) == ["geneCluster.json", "geneCluster/a_aa_aln.fa", "newColumnConfig.js"]
    assert (by_name["geneCluster.json"].origin, by_name["geneCluster.json"].mtime) == (
        Origin.SCICORE_2018,
        1_657_000_000,
    )
    alignment = by_name["geneCluster/a_aa_aln.fa"]
    assert (alignment.origin, alignment.gzip, gzip.decompress(alignment.data)) == (Origin.S3, True, new)
    assert alignment.mtime == 1_661_442_243
    assert by_name["newColumnConfig.js"].gzip is False


def test_rebuild_older_dataset_object_without_listed_time_fails(tmp_path: Path) -> None:
    overlay = tmp_path / "Testus"
    overlay.mkdir()
    (overlay / "metaConfiguration.js").write_bytes(b"var meta_details = {};")

    with pytest.raises(ValueError, match="no LastModified"):
        list(older_dataset_files(OlderJob("Testus", None, "Testus", overlay), {}))


@pytest.mark.skipif(shutil.which("pigz") is None, reason="needs the pigz program")
def test_rebuild_older_dataset_writes_the_forms_that_the_bucket_holds(tmp_path: Path) -> None:
    table, alignment, extra = b'[{"msa": "a"}]', b">a|1\nMKV\n", b"(a:0.1);\n"
    stored = gzip.compress(alignment, mtime=0)
    tarball = _tarball(
        tmp_path / "Testus.tar.gz",
        {
            "Testus/vis/geneCluster.json": table,
            "Testus/vis/geneCluster/a_aa_aln.fa.gz": stored,
            "Testus/vis/geneCluster/a.nwk": extra,
        },
    )
    copy = pigz_bytes(table, Path("geneCluster.json"), 1_657_000_000)
    objects = {
        key: S3Object(key, hashlib.md5(data).hexdigest(), len(data))
        for key, data in (
            ("dataset/Testus/geneCluster.json", table),
            ("dataset/Testus/geneCluster.json.gz", copy),
            ("dataset/Testus/geneCluster/a_aa_aln.fa.gz", stored),
        )
    }
    out = Out(tmp_path / "out")
    for folder in (out.bucket / "dataset", out.staging, out.parts):
        folder.mkdir(parents=True)

    result = build_older_dataset(OlderJob("Testus", tarball, "Testus", tmp_path / "no-s3"), out, objects)

    folder = out.bucket / "dataset" / "Testus"
    files = {path.relative_to(folder).as_posix(): path.read_bytes() for path in folder.rglob("*") if path.is_file()}
    assert files == {
        "geneCluster.json": table,
        "geneCluster.json.gz": copy,
        "geneCluster/a_aa_aln.fa.gz": stored,
        "geneCluster/a.nwk": extra,
    }
    manifest = dict(line.split("\t", 1) for line in (out.parts / "Testus.tsv").read_text().splitlines())
    assert manifest["dataset/Testus/geneCluster/a.nwk"].endswith("\tscicore-2018")
    assert (result.files, result.from_s3) == (4, 0)


def _tarball(path: Path, files: dict[str, bytes]) -> Path:
    with tarfile.open(path, "w:gz") as tar:
        for name, data in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 1_657_000_000
            tar.addfile(info, io.BytesIO(data))
    return path
