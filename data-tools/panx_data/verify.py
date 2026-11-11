"""Compare the rebuilt bucket and the copied S3 objects with the saved S3 listings.

Every difference goes to `out/deviations.tsv` with its class. The known classes describe the bucket, not the
rebuild: objects that the upload left out or broke, stale archives, and datasets that the archives lack. Any other
class means that the rebuild does not reproduce the bucket, and `verify` fails.
"""

import gzip
import hashlib
import logging
from collections import Counter
from collections.abc import Iterator
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from panx_data.layout import Out, Source
from panx_data.listings import Listings, S3Object

LOG = logging.getLogger(__name__)
MULTIPART_SIZES = (8 << 20, 5 << 20, 16 << 20, 1 << 20)
STALE_ARCHIVE_SUFFIX = "_gene_alignments.tar.gz"


def verify(source: Source, out: Out, listings: Listings) -> Counter["Deviation"]:
    """Write `out/deviations.tsv` and return the number of files in each class."""
    deviations = [*_check_out(out, listings), *_check_source(source, listings)]
    lines = [
        f"{d.path}\t{d.kind.value}\t{_text(d.s3_size)}\t{d.s3_etag or '-'}\t{_text(d.size)}\t{d.md5 or '-'}\n"
        for d in sorted(deviations, key=_path)
    ]
    out.deviations.write_text("path\tclass\ts3_size\ts3_etag\tlocal_size\tlocal_md5\n" + "".join(lines))
    return Counter(d.kind for d in deviations)


def classify(local_size: int, s3_size: int | None) -> "Deviation":
    """Name the kind of a file whose bytes differ from its S3 object."""
    if s3_size == 0:
        return Deviation.S3_EMPTY
    if s3_size is not None and s3_size < local_size:
        return Deviation.S3_SHORTER
    return Deviation.DIFFERS


def matches(obj: S3Object, md5: str, path: Path) -> bool:
    """Tell whether a file has the bytes of an S3 object, by MD5 or by the ETag of a multipart upload."""
    return _multipart_matches(path, obj.etag) if obj.is_multipart else md5 == obj.etag


class Deviation(Enum):
    """Kinds of differences between the rebuild and the S3 listings."""

    MISSING_IN_S3 = "missing-in-s3"
    S3_EMPTY = "s3-empty"
    S3_SHORTER = "s3-shorter"
    STALE_ARCHIVE = "only-in-s3:stale-archive"
    NOT_REBUILT = "only-in-s3:dataset-not-rebuilt"
    DIFFERS = "differs"
    ONLY_IN_S3 = "only-in-s3"
    SOURCE_DIFFERS = "src-differs-from-s3"

    @property
    def expected(self) -> bool:
        return self in EXPECTED


EXPECTED = frozenset(
    {
        Deviation.MISSING_IN_S3,
        Deviation.S3_EMPTY,
        Deviation.S3_SHORTER,
        Deviation.STALE_ARCHIVE,
        Deviation.NOT_REBUILT,
    }
)


@dataclass(frozen=True)
class Found:
    """One difference; local fields are None for an object that exists only in S3."""

    path: str
    kind: Deviation
    s3_size: int | None
    s3_etag: str | None
    size: int | None
    md5: str | None


def _check_out(out: Out, listings: Listings) -> Iterator[Found]:
    seen: set[str] = set()
    datasets: set[str] = set()
    for path, size, md5 in _manifest_rows(out.manifest):
        if path.startswith("dataset/"):
            datasets.add(path.split("/")[1])
        obj = listings.objects.get(path)
        if obj is not None:
            seen.add(path)
            if not matches(obj, md5, out.bucket / path):
                yield Found(path, classify(size, obj.size), obj.size, obj.etag, size, md5)
        elif listings.covers(path):
            yield Found(path, Deviation.MISSING_IN_S3, None, None, size, md5)
    for key, obj in listings.objects.items():
        if key not in seen:
            yield Found(key, _only_in_s3(key, datasets), obj.size, obj.etag, None, None)


def _check_source(source: Source, listings: Listings) -> Iterator[Found]:
    files = sorted(path for path in source.s3_objects.rglob("*") if path.is_file())
    for path in files:
        key = path.relative_to(source.s3_objects).as_posix()
        obj = listings.objects.get(key)
        if obj is None:
            LOG.warning("src object '%s' is in no listing", key)
            continue
        md5 = _md5_of(path)
        if not matches(obj, md5, path):
            yield Found(f"src:{key}", Deviation.SOURCE_DIFFERS, obj.size, obj.etag, path.stat().st_size, md5)
    LOG.info("checked %d objects of src/ against their listings", len(files))


def _only_in_s3(key: str, rebuilt: set[str]) -> Deviation:
    if key.endswith(STALE_ARCHIVE_SUFFIX):
        return Deviation.STALE_ARCHIVE
    if key.startswith("dataset/") and key.split("/")[1] not in rebuilt:
        return Deviation.NOT_REBUILT
    return Deviation.ONLY_IN_S3


def _manifest_rows(manifest: Path) -> Iterator[tuple[str, int, str]]:
    with gzip.open(manifest, "rt") as lines:
        next(lines)
        for line in lines:
            path, size, md5, _ = line.rstrip("\n").split("\t")
            yield path, int(size), md5


def _multipart_matches(path: Path, etag: str) -> bool:
    if not path.exists():
        return False
    digest, _, count = etag.partition("-")
    size = path.stat().st_size
    return any(_multipart_digest(path, part) == digest for part in MULTIPART_SIZES if -(-size // part) == int(count))


def _multipart_digest(path: Path, part: int) -> str:
    with path.open("rb") as source:
        parts = b"".join(hashlib.md5(chunk).digest() for chunk in iter(lambda: source.read(part), b""))
    return hashlib.md5(parts).hexdigest()


def _md5_of(path: Path) -> str:
    with path.open("rb") as source:
        return hashlib.file_digest(source, "md5").hexdigest()


def _path(found: Found) -> str:
    return found.path


def _text(value: int | None) -> str:
    return "-" if value is None else str(value)
