"""Rebuild the bucket's 2022 datasets from the sciCORE archives, in one pass over each dataset's tarball.

The upload scripts of `neherlab/pan-genome-visualization` filled the bucket:

- August 2022: every dataset from `13-07-2022/<dataset>.tar.gz`, the alignment zips from `zip_13-07-2022/<dataset>/`
- 2022-11-18: the datasets of `04-11-2022/` replaced their folders, with the zips and the `web/` files of
  `04-11-2022/zip_archives_update.zip`
- before each upload, `scripts/s3_compress.sh` decompressed every `.gz` file (`pigz -d`) and wrote a `.gz` copy of
  every file (`pigz -k`); the zips came later and have no `.gz` copy
- in August, `aws s3 sync` skipped a new `.gz` copy when the pipeline's original `.gz` file was already in the bucket
  with the same size, because the local file was older, so the bucket holds the original bytes there. A pass on
  2022-08-26 uploaded new copies of the same size for part of 17 datasets. Where a listing gives a `.gz` key's ETag,
  the ETag therefore decides between the two versions
- objects in `src/s3/` replace or add files: S3 changed them after the upload (the `metainfo.tsv` files of
  2022-09-09, a few `.gz` copies with another header time), or the archives lack them (the bucket root)

Each dataset is built in `out/.staging/` and moved into place when it is complete, so an interrupted run resumes.
"""

import gzip
import hashlib
import json
import logging
import os
import shutil
import subprocess
import tarfile
import time
import zipfile
from collections.abc import Iterator
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path

from panx_data.layout import (
    ALIGNMENT_ZIPS,
    JULY_RUN,
    JULY_ZIPS,
    NOVEMBER_RUN,
    UPDATE_ZIP_PREFIX,
    WEB_FILES,
    Out,
    Source,
)
from panx_data.listings import Listings
from panx_data.pigz_gzip import check_compressors, pigz_bytes

LOG = logging.getLogger(__name__)
ROOT_PART = "_root"
SUMMARY_COLUMNS = ("dataset", "upload", "files", "bytes", "originals_kept", "from_s3", "seconds")


def rebuild(
    source: Source, out: Out, listings: Listings, jobs: int, only: frozenset[str] | None = None
) -> list["Result"]:
    """Build every dataset of the bucket that the archives hold, then the bucket root, then the manifest."""
    if out.build_info.exists():
        raise FileExistsError(f"'{out.root}' holds a complete build; delete it to build again")
    check_compressors()
    planned = plan(source, listings, only)
    for folder in (out.bucket / "dataset", out.staging, out.parts):
        folder.mkdir(parents=True, exist_ok=True)
    pending = [job for job in planned if not _is_built(out, job.dataset)]
    LOG.info("%d of %d datasets already built, %d to build", len(planned) - len(pending), len(planned), len(pending))
    failures = _build_all(pending, out, listings.single_part_etags_by_dataset(), jobs)
    if failures:
        raise RuntimeError(f"{len(failures)} datasets failed: {' '.join(sorted(failures))}")
    build_root(source, out)
    names = [job.dataset for job in planned]
    _merge_manifest(out, names)
    results = [_read_result(out, name) for name in [*sorted(names), ROOT_PART]]
    _write_summary(out, results)
    shutil.rmtree(out.parts)
    shutil.rmtree(out.staging)
    return results


def plan(source: Source, listings: Listings, only: frozenset[str] | None) -> list["Job"]:
    """List the datasets to build: those in the bucket that have a tarball, the largest first."""
    july = _tarball_names(source.archives / JULY_RUN)
    november = _tarball_names(source.archives / NOVEMBER_RUN)
    unknown = (only or frozenset()) - listings.datasets
    if unknown:
        raise ValueError(f"not datasets of the bucket: {' '.join(sorted(unknown))}")
    names = (listings.datasets & july) & (only or listings.datasets)
    jobs = [
        Job(
            dataset=name,
            upload=(upload := Upload.NOVEMBER_2022 if name in november else Upload.AUGUST_2022),
            tarball=source.tarball(NOVEMBER_RUN if upload is Upload.NOVEMBER_2022 else JULY_RUN, name),
            zip_folder=None if upload is Upload.NOVEMBER_2022 else source.archives / JULY_ZIPS / name,
            update_zip=source.update_zip if upload is Upload.NOVEMBER_2022 else None,
            overlay=source.s3_objects / "dataset" / name,
        )
        for name in names
    ]
    return sorted(jobs, key=_tarball_size, reverse=True)


def build_dataset(job: "Job", out: Out, etags: dict[str, str]) -> "Result":
    """Write one dataset into `out/`: every file plain and as a `.gz` copy, the zips, then the objects from S3."""
    start = time.monotonic()
    work = out.staging / job.dataset
    shutil.rmtree(work, ignore_errors=True)
    rows: dict[str, Row] = {}
    writer = _Writer(work)
    kept = 0
    with _tar_members(job.tarball) as members:
        for member, data in members:
            top, _, rel = member.name.partition("/")
            if top != job.dataset or not rel:
                raise ValueError(f"'{job.tarball}': member '{member.name}' is outside folder '{job.dataset}/'")
            original = data if rel.endswith(".gz") else None
            plain = gzip.decompress(data) if original is not None else data
            key = f"dataset/{job.dataset}/{rel.removesuffix('.gz')}"
            if key in rows:
                raise ValueError(f"'{job.tarball}': two members give '{key}'")
            path = writer.write(key.split("/", 2)[2], plain, member.mtime)
            copy = pigz_bytes(plain, path, member.mtime)
            gz = choose_gz(copy, original, job.upload, etags.get(f"{key}.gz"))
            kept += gz is original
            writer.write(f"{key.split('/', 2)[2]}.gz", gz, member.mtime)
            rows |= {key: Row.of(plain, Origin.SCICORE), f"{key}.gz": Row.of(gz, Origin.SCICORE)}
    for name, data in _alignment_zips(job):
        writer.write(name, data, None)
        rows[f"dataset/{job.dataset}/{name}"] = Row.of(data, Origin.SCICORE)
    from_s3 = _apply_overlay(job.overlay, writer, f"dataset/{job.dataset}/", rows)
    target = out.bucket / "dataset" / job.dataset
    os.replace(work, target)
    result = Result(
        job.dataset,
        job.upload.value,
        len(rows),
        sum(row.size for row in rows.values()),
        kept,
        from_s3,
        round(time.monotonic() - start, 1),
    )
    _write_part(out, job.dataset, rows, result)
    return result


def build_root(source: Source, out: Out) -> "Result":
    """Write the files outside `dataset/`: the bucket root from `src/s3/` and `web/` from the update zip."""
    start = time.monotonic()
    writer = _Writer(out.bucket)
    with zipfile.ZipFile(source.update_zip) as update:
        web = {f"web/{name}": update.read(f"{UPDATE_ZIP_PREFIX}{name}") for name in WEB_FILES}
    rows = {key: Row.of(data, Origin.SCICORE) for key, data in web.items()}
    for key, data in web.items():
        writer.write(key, data, None)
    root_objects = sorted(
        path
        for path in source.s3_objects.rglob("*")
        if path.is_file() and path.relative_to(source.s3_objects).parts[0] != "dataset"
    )
    for path in root_objects:
        key = path.relative_to(source.s3_objects).as_posix()
        data = path.read_bytes()
        writer.write(key, data, int(path.stat().st_mtime))
        rows[key] = Row.of(data, Origin.S3)
    result = Result(
        ROOT_PART,
        "",
        len(rows),
        sum(row.size for row in rows.values()),
        0,
        len(root_objects),
        round(time.monotonic() - start, 1),
    )
    _write_part(out, ROOT_PART, rows, result)
    return result


def choose_gz(copy: bytes, original: bytes | None, upload: "Upload", etag: str | None) -> bytes:
    """Pick the `.gz` bytes that the bucket holds: the ETag decides, else the rule of the upload.

    Only the August upload can have kept the pipeline's `original`, and it kept it when both had the same size.
    """
    candidates = [copy] if original is None or upload is not Upload.AUGUST_2022 else [original, copy]
    by_etag = [data for data in candidates if etag is not None and hashlib.md5(data).hexdigest() == etag]
    if by_etag:
        return by_etag[0]
    return original if original is not None and upload is Upload.AUGUST_2022 and len(original) == len(copy) else copy


class Upload(Enum):
    """The upload that wrote a dataset into the bucket."""

    AUGUST_2022 = "2022-08"
    NOVEMBER_2022 = "2022-11-18"


class Origin(Enum):
    """Where the bytes of an output file come from."""

    SCICORE = "scicore"
    S3 = "s3"


@dataclass(frozen=True)
class Job:
    """One dataset to build, with its inputs."""

    dataset: str
    upload: Upload
    tarball: Path
    zip_folder: Path | None
    update_zip: Path | None
    overlay: Path


@dataclass(frozen=True)
class Result:
    """Counts of one built dataset, or of the bucket root."""

    dataset: str
    upload: str
    files: int
    bytes: int
    originals_kept: int
    from_s3: int
    seconds: float


@dataclass(frozen=True)
class Row:
    """One manifest row without its path."""

    size: int
    md5: str
    origin: Origin

    @staticmethod
    def of(data: bytes, origin: Origin) -> "Row":
        return Row(len(data), hashlib.md5(data).hexdigest(), origin)


def _tarball_names(folder: Path) -> frozenset[str]:
    return frozenset(path.name.removesuffix(".tar.gz") for path in folder.glob("*.tar.gz"))


def _tarball_size(job: Job) -> int:
    return job.tarball.stat().st_size


def _is_built(out: Out, dataset: str) -> bool:
    target = out.bucket / "dataset" / dataset
    if (out.parts / f"{dataset}.json").exists() and target.exists():
        return True
    if target.exists():  # moved into place, but its manifest part is missing: build it again
        shutil.rmtree(target)
    return False


def _build_all(jobs: list[Job], out: Out, etags: dict[str, dict[str, str]], workers: int) -> list[str]:
    failures: list[str] = []
    with ProcessPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(build_dataset, job, out, etags.get(job.dataset, {})): job for job in jobs}
        for done, future in enumerate(as_completed(futures), start=1):
            job = futures[future]
            try:
                result = future.result()
            except Exception:  # report each failed dataset, finish the others, fail at the end
                LOG.exception("[%d/%d] %s failed", done, len(jobs), job.dataset)
                failures.append(job.dataset)
                continue
            LOG.info(
                "[%d/%d] %s (%s): %d files, %.0f MB, %d original .gz kept, %d from S3, %.0f s",
                done,
                len(jobs),
                result.dataset,
                result.upload,
                result.files,
                result.bytes / 1e6,
                result.originals_kept,
                result.from_s3,
                result.seconds,
            )
    return failures


@contextmanager
def _tar_members(tarball: Path) -> Iterator[Iterator[tuple[tarfile.TarInfo, bytes]]]:
    """Stream the regular files of a gzipped tarball, decompressed by `pigz` in a second process."""
    with subprocess.Popen(["pigz", "-dc", "--", str(tarball)], stdout=subprocess.PIPE) as reader:
        assert reader.stdout is not None
        with tarfile.open(fileobj=reader.stdout, mode="r|") as tar:
            yield _regular_files(tar, tarball)
        reader.stdout.read()  # the padding after the end marker, so that pigz can finish writing
    if reader.returncode != 0:
        raise RuntimeError(f"'pigz -dc {tarball}' failed with exit code {reader.returncode}")


def _regular_files(tar: tarfile.TarFile, tarball: Path) -> Iterator[tuple[tarfile.TarInfo, bytes]]:
    for member in tar:
        if member.isdir():
            continue
        if not member.isfile():
            raise ValueError(f"'{tarball}': member '{member.name}' is not a regular file")
        stream = tar.extractfile(member)
        assert stream is not None
        yield member, stream.read()


def _alignment_zips(job: Job) -> Iterator[tuple[str, bytes]]:
    if job.update_zip is not None:
        with zipfile.ZipFile(job.update_zip) as update:
            yield from ((name, update.read(f"{UPDATE_ZIP_PREFIX}{job.dataset}/{name}")) for name in ALIGNMENT_ZIPS)
    else:
        assert job.zip_folder is not None
        yield from ((name, (job.zip_folder / name).read_bytes()) for name in ALIGNMENT_ZIPS)


def _apply_overlay(overlay: Path, writer: "_Writer", prefix: str, rows: dict[str, Row]) -> int:
    files = sorted(path for path in overlay.rglob("*") if path.is_file()) if overlay.exists() else []
    for path in files:
        rel = path.relative_to(overlay).as_posix()
        data = path.read_bytes()
        writer.write(rel, data, int(path.stat().st_mtime))
        rows[prefix + rel] = Row.of(data, Origin.S3)
    return len(files)


def _write_part(out: Out, name: str, rows: dict[str, Row], result: Result) -> None:
    lines = "".join(f"{key}\t{row.size}\t{row.md5}\t{row.origin.value}\n" for key, row in sorted(rows.items()))
    _write_atomic(out.parts / f"{name}.tsv", lines.encode())
    _write_atomic(out.parts / f"{name}.json", json.dumps(asdict(result)).encode())


def _read_result(out: Out, name: str) -> Result:
    return Result(**json.loads((out.parts / f"{name}.json").read_text()))


def _merge_manifest(out: Out, datasets: list[str]) -> None:
    """Join the parts into one gzipped manifest, sorted by path: `dataset/` sorts before the root files."""
    temporary = out.manifest.with_name(out.manifest.name + ".tmp")
    with temporary.open("wb") as target, subprocess.Popen(["pigz", "-c"], stdin=subprocess.PIPE, stdout=target) as gz:
        assert gz.stdin is not None
        gz.stdin.write(b"path\tsize\tmd5\tsource\n")
        for name in [*sorted(datasets), ROOT_PART]:
            gz.stdin.write((out.parts / f"{name}.tsv").read_bytes())
        gz.stdin.close()
    if gz.returncode != 0:
        raise RuntimeError(f"pigz failed with exit code {gz.returncode} while writing '{out.manifest}'")
    os.replace(temporary, out.manifest)


def _write_summary(out: Out, results: list[Result]) -> None:
    rows = ["\t".join(SUMMARY_COLUMNS), *("\t".join(str(getattr(r, c)) for c in SUMMARY_COLUMNS) for r in results)]
    _write_atomic(out.summary, ("\n".join(rows) + "\n").encode())


def _write_atomic(path: Path, data: bytes) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_bytes(data)
    os.replace(temporary, path)


class _Writer:
    """Write files below one folder, creating each subfolder once."""

    def __init__(self, root: Path) -> None:
        self._root = root
        self._folders: set[Path] = set()

    def write(self, rel: str, data: bytes, mtime: int | None) -> Path:
        path = self._root / rel
        if path.parent not in self._folders:
            path.parent.mkdir(parents=True, exist_ok=True)
            self._folders.add(path.parent)
        with path.open("wb") as target:
            target.write(data)
        if mtime is not None:
            os.utime(path, (mtime, mtime))
        return path
