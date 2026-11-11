"""Folders of a snapshot: `src/` holds the inputs at their source paths, `out/` the rebuilt bucket."""

from dataclasses import dataclass
from pathlib import Path

BUCKET = "pangenome.org-data-master"
JULY_RUN = "13-07-2022"
NOVEMBER_RUN = "04-11-2022"
JULY_ZIPS = "zip_13-07-2022"
UPDATE_ZIP = "zip_archives_update.zip"
UPDATE_ZIP_PREFIX = "zip_archives_update/"
ALIGNMENT_ZIPS = ("all_gene_alignments.zip", "core_gene_alignments.zip")
WEB_FILES = ("all_downloads_table.json", "species-list-info.js")


def find_source(src: Path) -> "Source":
    """Locate the inputs below `src`: the sciCORE folder that holds the run folders, the S3 objects, the listings."""
    runs = [path.parent for path in (src / "scicore").rglob(JULY_RUN) if path.is_dir()]
    if len(runs) != 1:
        raise FileNotFoundError(f"expected one folder '{JULY_RUN}/' below '{src / 'scicore'}', found {len(runs)}")
    return Source(archives=runs[0], s3_objects=src / "s3" / BUCKET, listings=src / "s3-listings")


@dataclass(frozen=True)
class Source:
    """Inputs: the sciCORE archives, objects copied from the bucket at their keys, and saved listings."""

    archives: Path
    s3_objects: Path
    listings: Path

    @property
    def update_zip(self) -> Path:
        return self.archives / NOVEMBER_RUN / UPDATE_ZIP

    def tarball(self, run: str, dataset: str) -> Path:
        return self.archives / run / f"{dataset}.tar.gz"


@dataclass(frozen=True)
class Out:
    """The rebuilt bucket and the files that describe it."""

    root: Path

    @property
    def bucket(self) -> Path:
        return self.root / BUCKET

    @property
    def manifest(self) -> Path:
        return self.root / "manifest.tsv.gz"

    @property
    def deviations(self) -> Path:
        return self.root / "deviations.tsv"

    @property
    def summary(self) -> Path:
        return self.root / "build-summary.tsv"

    @property
    def build_info(self) -> Path:
        return self.root / "BUILD_INFO.json"

    @property
    def staging(self) -> Path:
        return self.root / ".staging"

    @property
    def parts(self) -> Path:
        return self.root / ".manifest-parts"
