"""Read the files of one dataset as the old data host serves them, from the archives and the objects of `src/s3/`.

A 2022 dataset comes from the tarball of its run (`plan`), and the plain objects of `src/s3/` replace or add a file:
the `metainfo.tsv` files that S3 rewrote on 2022-09-09. The `.gz` objects of `src/s3/` of a 2022 dataset differ from
the archive only in gzip bytes, so they are no input of the content; `rebuild` copies them.

An older dataset, which the 2022 run did not make, comes from the `vis/` folder of its 2018 archive on sciCORE and from
the objects of `src/s3/`, which replace a file of the archive: the copy chose S3 where its content is newer, or where
the archive lacks the file.

Both readers yield each file once, with its stored bytes, so a worker holds a dataset's files compressed. The mtime of
an object of `src/s3/` is the `LastModified` time of its key in the saved listings, never a file-system time, because
`src/MD5SUMS` covers the listings but not the file times.
"""

import gzip
import logging
import subprocess
import tarfile
from collections.abc import Callable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from enum import Enum
from pathlib import Path

from panx_data.layout import JULY_RUN, JULY_ZIPS, NOVEMBER_RUN, OLDER_FOLDER, Source, older_tarball_name
from panx_data.listings import Listings

LOG = logging.getLogger(__name__)


def plan(source: Source, listings: Listings, only: frozenset[str] | None) -> list["Job | OlderJob"]:
    """List the datasets to read, the largest first: those in the bucket that have a 2022 tarball, and the older
    datasets that have a 2018 archive or objects in `src/s3/`."""
    july = _tarball_names(source.archives / JULY_RUN)
    november = _tarball_names(source.archives / NOVEMBER_RUN)
    unknown = (only or frozenset()) - listings.datasets
    if unknown:
        raise ValueError(f"not datasets of the bucket: {' '.join(sorted(unknown))}")
    names = (listings.datasets & july) & (only or listings.datasets)
    jobs: list[Job | OlderJob] = [
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
    older = [
        OlderJob(name, source.older_tarball(name), older_tarball_name(name), source.s3_objects / "dataset" / name)
        for name in sorted((listings.datasets - july) & (only or listings.datasets))
    ]
    jobs += [job for job in older if job.tarball is not None or job.overlay.is_dir()]
    return sorted(jobs, key=_tarball_size, reverse=True)


def dataset_files(job: "Job", modified: Mapping[str, str | None]) -> Iterator["DatasetFile"]:
    """Yield each file of a 2022 dataset once, by its path below the dataset folder without `.gz`.

    The members of the tarball come first, in the tarball's order, with the tar member mtime; then the plain objects
    of `src/s3/`, each in place of the member of the same path. An object's mtime is the `LastModified` time of its key
    in `modified`.
    """
    replaced = _plain_objects(job.overlay)
    seen: set[str] = set()
    with tar_members(job.tarball) as members:
        for member, data in members:
            top, _, rel = member.name.partition("/")
            if top != job.dataset or not rel:
                raise ValueError(f"'{job.tarball}': member '{member.name}' is outside folder '{job.dataset}/'")
            name = rel.removesuffix(".gz")
            if name in seen:
                raise ValueError(f"'{job.tarball}': two members give '{name}'")
            seen.add(name)
            if name not in replaced:
                yield DatasetFile(name, data, rel.endswith(".gz"), member.mtime, Origin.SCICORE)
    for name, path in sorted(replaced.items()):
        key = f"dataset/{job.dataset}/{name}"
        yield DatasetFile(name, path.read_bytes(), False, listed_mtime(key, modified), Origin.S3)


def older_dataset_files(job: "OlderJob", modified: Mapping[str, str | None]) -> Iterator["DatasetFile"]:
    """Yield each file of an older dataset once, by its path below the dataset folder without `.gz`.

    A file comes from `src/s3/` when it has an object there, else from the `vis/` folder of the 2018 archive, with
    the tar member mtime. Of the objects of one file, the `.gz` object wins; when its content differs from the plain
    object, the file is reported. An object's mtime is the `LastModified` time of its key in `modified`.
    """
    overlay = _overlay_files(job.overlay)
    if job.tarball is not None:
        archived: dict[str, DatasetFile] = {}
        with tar_members(job.tarball, lambda name: _served_path(name, job.tar_folder) is not None) as members:
            for member, data in members:
                rel = _served_path(member.name, job.tar_folder)
                assert rel is not None
                name = rel.removesuffix(".gz")
                if name in overlay:
                    continue
                file = DatasetFile(name, data, rel.endswith(".gz"), member.mtime, Origin.SCICORE_2018)
                archived[name] = _prefer_gz(archived[name], file, job.tarball) if name in archived else file
        yield from archived.values()
    for name, path in sorted(overlay.items()):
        key = f"dataset/{job.dataset}/{path.relative_to(job.overlay).as_posix()}"
        yield DatasetFile(name, path.read_bytes(), path.name.endswith(".gz"), listed_mtime(key, modified), Origin.S3)


def job_files(job: "Job | OlderJob", modified: Mapping[str, str | None]) -> Iterator["DatasetFile"]:
    """Yield the files of any dataset with the reader of its kind."""
    return older_dataset_files(job, modified) if isinstance(job, OlderJob) else dataset_files(job, modified)


def overlay_modified(job: "Job | OlderJob", listings: Listings) -> dict[str, str | None]:
    """Return the `LastModified` time of each object of `src/s3/` of a dataset, by its key."""
    if not job.overlay.is_dir():
        return {}
    keys = (f"dataset/{job.dataset}/{path.relative_to(job.overlay).as_posix()}" for path in job.overlay.rglob("*"))
    return {key: obj.modified for key in keys if (obj := listings.objects.get(key)) is not None}


def listed_mtime(key: str, modified: Mapping[str, str | None]) -> int:
    """Return the `LastModified` time of an object of `src/s3/` as a Unix time."""
    listed = modified.get(key)
    if listed is None:
        raise ValueError(f"'{key}' of src/s3/ has no LastModified time in the saved listings")
    return int(datetime.fromisoformat(listed).timestamp())


class Upload(Enum):
    """The upload that wrote a dataset into the bucket."""

    AUGUST_2022 = "2022-08"
    NOVEMBER_2022 = "2022-11-18"


class Origin(Enum):
    """Where the bytes of a file come from."""

    SCICORE = "scicore"
    SCICORE_2018 = "scicore-2018"
    S3 = "s3"


@dataclass(frozen=True)
class Job:
    """One 2022 dataset, with its inputs."""

    dataset: str
    upload: Upload
    tarball: Path
    zip_folder: Path | None
    update_zip: Path | None
    overlay: Path


@dataclass(frozen=True)
class OlderJob:
    """One older dataset: its 2018 archive (None when it has none), the archive's folder name, the objects copied
    from S3."""

    dataset: str
    tarball: Path | None
    tar_folder: str
    overlay: Path


@dataclass(frozen=True)
class DatasetFile:
    """One file of a dataset: its path below the dataset folder without `.gz`, its stored bytes, whether they are
    gzip-compressed, its mtime, and its origin."""

    rel: str
    data: bytes
    gzip: bool
    mtime: int
    origin: Origin

    def content(self) -> bytes:
        """Return the plain content."""
        return gzip.decompress(self.data) if self.gzip else self.data


@contextmanager
def tar_members(
    tarball: Path, wanted: Callable[[str], bool] = lambda _: True
) -> Iterator[Iterator[tuple[tarfile.TarInfo, bytes]]]:
    """Stream the regular files of a gzipped tarball whose names are `wanted`, decompressed by `pigz` in a second
    process."""
    with subprocess.Popen(["pigz", "-dc", "--", str(tarball)], stdout=subprocess.PIPE) as reader:
        assert reader.stdout is not None
        with tarfile.open(fileobj=reader.stdout, mode="r|") as tar:
            yield _regular_files(tar, tarball, wanted)
        reader.stdout.read()  # the padding after the end marker, so that pigz can finish writing
    if reader.returncode != 0:
        raise RuntimeError(f"'pigz -dc {tarball}' failed with exit code {reader.returncode}")


def _regular_files(
    tar: tarfile.TarFile, tarball: Path, wanted: Callable[[str], bool]
) -> Iterator[tuple[tarfile.TarInfo, bytes]]:
    for member in tar:
        if member.isdir() or not wanted(member.name):
            continue
        if not member.isfile():
            raise ValueError(f"'{tarball}': member '{member.name}' is not a regular file")
        stream = tar.extractfile(member)
        assert stream is not None
        yield member, stream.read()


def _tarball_names(folder: Path) -> frozenset[str]:
    return frozenset(path.name.removesuffix(".tar.gz") for path in folder.glob("*.tar.gz"))


def _tarball_size(job: "Job | OlderJob") -> int:
    return 0 if job.tarball is None else job.tarball.stat().st_size


def _plain_objects(overlay: Path) -> dict[str, Path]:
    """Map each plain object of a 2022 dataset in `src/s3/` to its path below the dataset folder."""
    if not overlay.is_dir():
        return {}
    paths = (path for path in overlay.rglob("*") if path.is_file() and not path.name.endswith(".gz"))
    return {path.relative_to(overlay).as_posix(): path for path in paths}


def _served_path(name: str, tar_folder: str) -> str | None:
    """Return the path below the dataset folder of a 2018 archive member that the old layout serves, else None.

    The old layout has files at the top of `vis/` and in `vis/geneCluster/`; deeper folders hold working files.
    """
    prefix = f"{tar_folder}/{OLDER_FOLDER}/"
    if not name.startswith(prefix):
        return None
    rel = name.removeprefix(prefix)
    folder, _, rest = rel.rpartition("/")
    return rel if rest and folder in ("", "geneCluster") else None


def _overlay_files(overlay: Path) -> dict[str, Path]:
    """Map each file of an older dataset in `src/s3/` to the object to read: its `.gz` object when one exists."""
    files: dict[str, Path] = {}
    paths = sorted(path for path in overlay.rglob("*") if path.is_file()) if overlay.is_dir() else []
    for path in paths:
        name = path.relative_to(overlay).as_posix().removesuffix(".gz")
        if path.name.endswith(".gz"):
            plain = files.get(name)
            if plain is not None and gzip.decompress(path.read_bytes()) != plain.read_bytes():
                LOG.warning("'%s' and its .gz object differ in content; the .gz object is used", plain)
            files[name] = path
        else:
            files.setdefault(name, path)
    return files


def _prefer_gz(first: DatasetFile, second: DatasetFile, tarball: Path) -> DatasetFile:
    """Of a plain and a `.gz` member of one file, take the `.gz` member, and report a difference in content."""
    if first.gzip == second.gzip:
        raise ValueError(f"'{tarball}': two members give '{first.rel}'")
    plain, gz = (second, first) if first.gzip else (first, second)
    if gzip.decompress(gz.data) != plain.data:
        LOG.warning(
            "'%s': the plain and .gz members of '%s' differ in content; the .gz member is used", tarball, gz.rel
        )
    return gz
