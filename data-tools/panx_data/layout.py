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
# The 2018 archive of the earlier panX server data: one `<name>.tar.gz` per older dataset, whose folder `<name>/vis/`
# holds the files in the old layout
OLDER_ARCHIVES = ("panX_archive", "leftover_data")
OLDER_FOLDER = "vis"
OLDER_TARBALL_NAMES = {"P_marinus_meta": "Pmnawg"}


def find_source(src: Path) -> "Source":
    """Locate the inputs below `src`: the sciCORE folder that holds the run folders, the S3 objects, the listings."""
    runs = [path.parent for path in (src / "scicore").rglob(JULY_RUN) if path.is_dir()]
    if len(runs) != 1:
        raise FileNotFoundError(f"expected one folder '{JULY_RUN}/' below '{src / 'scicore'}', found {len(runs)}")
    parent, name = OLDER_ARCHIVES
    older = [path for path in (src / "scicore").rglob(name) if path.is_dir() and path.parent.name == parent]
    if len(older) > 1:
        raise FileNotFoundError(f"expected at most one folder '{parent}/{name}/' below '{src / 'scicore'}'")
    return Source(
        archives=runs[0],
        s3_objects=src / "s3" / BUCKET,
        listings=src / "s3-listings",
        older_archives=older[0] if older else None,
    )


@dataclass(frozen=True)
class Source:
    """Inputs: the sciCORE archives, objects copied from the bucket at their keys, and saved listings."""

    archives: Path
    s3_objects: Path
    listings: Path
    older_archives: Path | None = None

    @property
    def update_zip(self) -> Path:
        return self.archives / NOVEMBER_RUN / UPDATE_ZIP

    def tarball(self, run: str, dataset: str) -> Path:
        return self.archives / run / f"{dataset}.tar.gz"

    def older_tarball(self, dataset: str) -> Path | None:
        """Return the 2018 archive of an older dataset, or None when it has none."""
        if self.older_archives is None:
            return None
        path = self.older_archives / f"{older_tarball_name(dataset)}.tar.gz"
        return path if path.is_file() else None


def older_tarball_name(dataset: str) -> str:
    """Name of the 2018 archive of an older dataset, which can differ from the dataset id."""
    return OLDER_TARBALL_NAMES.get(dataset, dataset)


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
