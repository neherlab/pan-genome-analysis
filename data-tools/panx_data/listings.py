"""Saved S3 listings of a bucket: its objects with ETags and sizes, and the places where a listing is complete.

Files below `src/s3-listings/<date>/<bucket>/`. The place of a file names the scope of its listing:

- `root.json`: `list-objects-v2 --delimiter /` at the bucket root, complete for keys without `/`
- `web.json`: `--prefix web/`, complete for `web/`
- `dataset.json`: `--prefix dataset/ --delimiter /`, the names of the datasets
- `top-level/<dataset>.json`: `--prefix dataset/<dataset>/ --delimiter /`, complete for the top level of the dataset
- `full/<dataset>.json`, `full/<dataset>.<columns>.tsv[.gz]`: every key of the dataset
- `samples/<dataset>.json`: one page of keys, complete for nothing
- `head-object.key-size-etag-encoding-modified.tsv`: sizes from `head-object`; failed requests are skipped

The `.json` files are `aws s3api list-objects-v2` output. The `.tsv` files are its `--output text` form, with the
columns that the file name gives, for example `key-etag-modified`. A listing of a later date replaces an object of
an earlier date.
"""

import gzip
import json
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass, replace
from pathlib import Path

HEAD_OBJECT = "head-object.key-size-etag-encoding-modified.tsv"


def load_listings(root: Path, bucket: str) -> "Listings":
    """Read every saved listing of `bucket` below `root`, the oldest date first."""
    folders = sorted(path / bucket for path in root.iterdir() if (path / bucket).is_dir())
    if not folders:
        raise FileNotFoundError(f"no listings of '{bucket}' below '{root}'")
    objects: dict[str, S3Object] = {}
    datasets: set[str] = set()
    top_level: set[str] = set()
    full: set[str] = set()
    for folder in folders:
        listed = [
            *_json_objects(folder / "root.json"),
            *_json_objects(folder / "web.json"),
            *(obj for path in sorted((folder / "top-level").glob("*.json")) for obj in _json_objects(path)),
            *(obj for path in sorted((folder / "full").glob("*")) for obj in _listing_objects(path)),
            *(obj for path in sorted((folder / "samples").glob("*.json")) for obj in _json_objects(path)),
        ]
        objects |= {obj.key: obj for obj in listed}
        objects |= {obj.key: _with_size(objects.get(obj.key), obj) for obj in _tsv_objects(folder / HEAD_OBJECT)}
        datasets |= {prefix.split("/")[1] for prefix in _json_prefixes(folder / "dataset.json")}
        top_level |= {path.stem for path in (folder / "top-level").glob("*.json")}
        full |= {path.name.split(".")[0] for path in (folder / "full").glob("*")}
    return Listings(
        objects=objects,
        datasets=frozenset(datasets),
        top_level_complete=frozenset(top_level),
        full_complete=frozenset(full),
        root_complete=any((folder / "root.json").exists() for folder in folders),
        web_complete=any((folder / "web.json").exists() for folder in folders),
    )


@dataclass(frozen=True)
class Listings:
    """Objects of the bucket, and the scopes where a missing key means a missing object."""

    objects: dict[str, "S3Object"]
    datasets: frozenset[str]
    top_level_complete: frozenset[str]
    full_complete: frozenset[str]
    root_complete: bool
    web_complete: bool

    def covers(self, key: str) -> bool:
        """Tell whether a complete listing includes the place of `key`."""
        head, _, rest = key.partition("/")
        if not rest:
            return self.root_complete
        if head == "web":
            return self.web_complete
        if head != "dataset":
            return False
        dataset, _, inner = rest.partition("/")
        return dataset in self.full_complete or ("/" not in inner and dataset in self.top_level_complete)

    def single_part_etags_by_dataset(self) -> dict[str, dict[str, str]]:
        """Group the ETags that are plain MD5 sums (single-part uploads) by dataset."""
        grouped: defaultdict[str, dict[str, str]] = defaultdict(dict)
        for key, obj in self.objects.items():
            if key.startswith("dataset/") and not obj.is_multipart:
                grouped[key.split("/")[1]][key] = obj.etag
        return dict(grouped)


@dataclass(frozen=True)
class S3Object:
    """One object of a listing; `size` is None when the listing has no sizes."""

    key: str
    etag: str
    size: int | None

    @property
    def is_multipart(self) -> bool:
        return "-" in self.etag


def _json_objects(path: Path) -> list[S3Object]:
    if not path.exists():
        return []
    contents = json.loads(path.read_text()).get("Contents") or []
    return [S3Object(item["Key"], item["ETag"].strip('"'), item["Size"]) for item in contents]


def _json_prefixes(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [item["Prefix"] for item in json.loads(path.read_text()).get("CommonPrefixes") or []]


def _listing_objects(path: Path) -> list[S3Object]:
    return _json_objects(path) if path.suffix == ".json" else list(_tsv_objects(path))


def _tsv_objects(path: Path) -> Iterator[S3Object]:
    """Read a text listing whose name gives its columns, for example `x.key-size-etag-modified.tsv.gz`."""
    if not path.exists():
        return
    columns = path.name.removesuffix(".gz").removesuffix(".tsv").split(".")[-1].split("-")
    with gzip.open(path, "rt") if path.suffix == ".gz" else path.open() as lines:
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(columns):
                continue  # a failed head-object request is logged as one error column
            row = dict(zip(columns, fields))
            yield S3Object(row["key"], row["etag"].strip('"'), int(row["size"]) if "size" in row else None)


def _with_size(known: S3Object | None, sized: S3Object) -> S3Object:
    return sized if known is None else replace(known, size=sized.size)
