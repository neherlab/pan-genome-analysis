import gzip
import hashlib
import json
from pathlib import Path

from panx_data.layout import Out, Source
from panx_data.listings import load_listings
from panx_data.verify import Deviation, classify, verify

BUCKET = "pangenome.org-data-master"


def test_verify_classify_names_s3_defects_and_real_differences() -> None:
    assert classify(100, 0) is Deviation.S3_EMPTY
    assert classify(100, 55) is Deviation.S3_SHORTER
    assert classify(100, 100) is Deviation.DIFFERS
    assert classify(100, None) is Deviation.DIFFERS


def test_verify_complete_scopes_decide_where_a_missing_object_counts(tmp_path: Path) -> None:
    listing = _listing_folder(tmp_path)
    _json(listing / "root.json", [])
    _json(listing / "top-level" / "Aus.json", [])
    (listing / "full").mkdir()
    (listing / "full" / "Bus.key-etag-modified.tsv").write_text("")

    listings = load_listings(tmp_path / "s3-listings", BUCKET)

    assert listings.covers("index.json")
    assert listings.covers("dataset/Aus/geneCluster.json")
    assert not listings.covers("dataset/Aus/geneCluster/x_tree.json")
    assert listings.covers("dataset/Bus/geneCluster/x_tree.json")
    assert not listings.covers("web/species-list-info.js")


def test_verify_gene_cluster_listing_completes_its_folder_and_keeps_modified_times(tmp_path: Path) -> None:
    listing = _listing_folder(tmp_path)
    (listing / "gene-cluster").mkdir()
    rows = 'dataset/Cus/geneCluster/x.nwk\t12\t"0123"\t2022-08-23T08:42:54+00:00\n'
    (listing / "gene-cluster" / "Cus.key-size-etag-modified.tsv").write_text(rows)

    listings = load_listings(tmp_path / "s3-listings", BUCKET)

    assert listings.covers("dataset/Cus/geneCluster/y_tree.json")
    assert not listings.covers("dataset/Cus/geneCluster.json")
    assert listings.objects["dataset/Cus/geneCluster/x.nwk"].modified == "2022-08-23T08:42:54+00:00"


def test_verify_sorts_each_difference_into_its_class(tmp_path: Path) -> None:
    good, broken = b"good bytes", b"correct bytes from the archive"
    listing = _listing_folder(tmp_path)
    _json(listing / "root.json", [])
    _json(
        listing / "full" / "Aus.json",
        [
            _object("dataset/Aus/good.json", good),
            _object("dataset/Aus/broken.json", b""),
            _object("dataset/Aus/all_gene_alignments.tar.gz", b"stale"),
        ],
    )
    _json(listing / "top-level" / "Old.json", [_object("dataset/Old/geneCluster.json", b"older dataset")])
    _json(
        listing / "full" / "Ous.json",
        [_object("dataset/Ous/a.fa.gz", b"gzip header 2022"), _object("dataset/Ous/b.fa.gz", b"gzip header 2022")],
    )
    out = Out(tmp_path / "out")
    out.root.mkdir()
    rows = [
        ("dataset/Aus/good.json", good, "scicore"),
        ("dataset/Aus/broken.json", broken, "scicore"),
        ("dataset/Aus/missing.json", good, "scicore"),
        ("dataset/Ous/a.fa.gz", b"gzip header 2017", "scicore-2018"),
        ("dataset/Ous/b.fa.gz", b"gzip header 2017", "s3"),
    ]
    with gzip.open(out.manifest, "wt") as manifest:
        manifest.write("path\tsize\tmd5\tsource\n")
        manifest.writelines(f"{p}\t{len(d)}\t{hashlib.md5(d).hexdigest()}\t{o}\n" for p, d, o in rows)
    source = Source(archives=tmp_path, s3_objects=tmp_path / "src-s3", listings=tmp_path / "s3-listings")
    (tmp_path / "src-s3").mkdir()

    counts = verify(source, out, load_listings(source.listings, BUCKET))

    assert counts == {
        Deviation.S3_EMPTY: 1,
        Deviation.MISSING_IN_S3: 1,
        Deviation.STALE_ARCHIVE: 1,
        Deviation.NOT_REBUILT: 1,
        Deviation.GZIP_HEADER: 1,
        Deviation.DIFFERS: 1,
    }


def _listing_folder(tmp_path: Path) -> Path:
    folder = tmp_path / "s3-listings" / "2026-11-11" / BUCKET
    folder.mkdir(parents=True)
    return folder


def _object(key: str, data: bytes) -> dict[str, object]:
    return {"Key": key, "Size": len(data), "ETag": f'"{hashlib.md5(data).hexdigest()}"'}


def _json(path: Path, contents: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"Contents": contents}))
