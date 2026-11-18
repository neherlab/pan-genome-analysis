"""Publish a release of the data format to the data host, an S3 bucket folder, with the AWS CLI.

The steps follow "Publication of a release" of `docs/data-format.md` of the app repository:

- `upload` writes every object of the release folder with the headers of its manifest entry, then `manifest.json`.
  Each write is conditional on a missing key, so no object is replaced; a key that exists counts as done when its
  SHA-256 checksum equals the manifest's, so an interrupted upload resumes, and stops the upload otherwise
- `list` adds the release to `releases.json` and keeps the default; `promote` makes it the default, which is also a
  rollback to an earlier release; `unlist` takes a release that is not the default out of the list. On a data host
  without `releases.json`, `list` and `promote` create it with the release as its only entry and default. Every
  write of `releases.json` is conditional on the ETag that was read, so two publishes at the same time cannot drop
  a release from the list

With `dry_run`, every step prints its requests, their count, and their cost, and sends nothing.
"""

import base64
import gzip
import json
import logging
import shlex
import subprocess
import tempfile
import tomllib
from collections import Counter
from collections.abc import Callable, Sequence
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

from panx_data.archive import ZLIB_VERSION, gzip_member, json_bytes
from panx_data.release import RELEASE_ID, RELEASES, object_headers

LOG = logging.getLogger(__name__)
MANIFEST = "manifest.json"
STEPS = ("upload", "list", "promote", "unlist")
MISE_CONFIG = Path(__file__).resolve().parent.parent / ".config" / "mise.toml"
# S3 Standard in eu-central-1: PUT $5.40 and GET or HEAD $0.43 per 1M requests, storage $0.0245 per GB-month
PRICES = {"put-object": 5.40e-6, "get-object": 0.43e-6, "head-object": 0.43e-6}
STORAGE_PRICE = 0.0245
ETAG_OF_THE_READ = "<ETag of the get-object>"


@dataclass(frozen=True)
class Release:
    """A release folder of a local data root: its `releases.json` entry, its manifest, and the manifest's bytes."""

    id: str
    entry: dict[str, Any]
    manifest: dict[str, Any]
    manifest_sha256: str
    manifest_size: int


@dataclass(frozen=True)
class Request:
    """One AWS CLI call of `aws s3api`."""

    operation: str
    key: str
    options: tuple[str, ...] = ()
    body: str | None = None
    size: int = 0
    sha256: str = ""

    def argv(self, bucket: str) -> list[str]:
        body = ["--body", self.body] if self.operation == "put-object" and self.body else []
        out = [self.body] if self.operation == "get-object" and self.body else []
        return ["aws", "s3api", self.operation, "--bucket", bucket, "--key", self.key, *self.options, *body, *out]


def pinned_python() -> str:
    """Return the Python version of `.config/mise.toml`, the runtime of every release build."""
    tools = tomllib.loads(MISE_CONFIG.read_text())["tools"]
    return str(tools["http:python"]["version"]).partition("+")[0]


def load_release(data_root: Path, release_id: str) -> Release:
    """Read a release folder of a local data root and its entry in the data root's `releases.json`."""
    if not RELEASE_ID.fullmatch(release_id):
        raise ValueError(f"'{release_id}' is no release id of the form YYYYMMDDThhmmssZ")
    stored = (data_root / release_id / MANIFEST).read_bytes()
    listed = json.loads(gzip.decompress((data_root / RELEASES).read_bytes()))
    entries = [entry for entry in listed["releases"] if entry["id"] == release_id]
    if not entries:
        raise ValueError(f"'{data_root / RELEASES}' does not list the release {release_id}")
    return Release(release_id, entries[0], json.loads(gzip.decompress(stored)), sha256(stored).hexdigest(), len(stored))


def local_sizes(data_root: Path, release: Release) -> dict[str, int]:
    """Return the size of every manifest object that exists in the release folder."""
    folder = data_root / release.id
    return {
        obj["path"]: (folder / obj["path"]).stat().st_size
        for obj in release.manifest["objects"]
        if (folder / obj["path"]).is_file()
    }


def refusals(release: Release, sizes: dict[str, int], python: str, zlib_version: str) -> list[str]:
    """Return why the release must not be published: another release, another tool or runtime, or local files
    that the manifest does not describe."""
    found = []
    build = release.manifest["build"]
    if release.manifest["release"] != release.id:
        found.append(f"the manifest of the folder {release.id} names the release {release.manifest['release']}")
    if build["uncommitted_changes"]:
        found.append(f"the release was built by the tool commit {build['commit']} with uncommitted changes")
    if build["python"] != python:
        found.append(f"the release was built with Python {build['python']}, the data container has Python {python}")
    if build["zlib"] != zlib_version:
        found.append(f"the release was built with zlib {build['zlib']}, the data container has zlib {zlib_version}")
    for obj in release.manifest["objects"]:
        size = sizes.get(obj["path"])
        if size is None:
            found.append(f"'{obj['path']}' of the manifest is missing")
        elif size != obj["size"]:
            found.append(f"'{obj['path']}' has {size} bytes, the manifest names {obj['size']}")
    return found


def upload_requests(release: Release, prefix: str, folder: Path | None = None) -> list[Request]:
    """Return the writes of every object of the release, then of its manifest, each with its headers."""
    root = Path(release.id) if folder is None else folder
    requests = [
        _put(
            _key(prefix, f"{release.id}/{obj['path']}"),
            obj,
            str(root / obj["path"]),
        )
        for obj in release.manifest["objects"]
    ]
    content_type, encoding, cache, disposition = object_headers(MANIFEST)
    manifest = {
        "size": release.manifest_size,
        "sha256": release.manifest_sha256,
        "contentType": content_type,
        "contentEncoding": encoding,
        "cacheControl": cache,
        "contentDisposition": disposition,
    }
    requests.append(_put(_key(prefix, f"{release.id}/{MANIFEST}"), manifest, str(root / MANIFEST)))
    return requests


def releases_requests(prefix: str, release_id: str, step: str) -> list[Request]:
    """Return the requests of a change of `releases.json`, the read and the conditional write.

    `list` and `promote` first check that the upload of the release ended, by its `manifest.json`.
    """
    check = [] if step == "unlist" else [Request("head-object", _key(prefix, f"{release_id}/{MANIFEST}"))]
    read = Request("get-object", _key(prefix, RELEASES), body=RELEASES)
    return [*check, read, releases_put(prefix, RELEASES, 0, ETAG_OF_THE_READ)]


def releases_put(prefix: str, body: str, size: int, etag: str | None) -> Request:
    """Return the write of `releases.json`: replace the version with the ETag that was read, or create the file."""
    content_type, encoding, cache, _ = object_headers(RELEASES)
    condition = ("--if-match", etag) if etag else ("--if-none-match", "*")
    options = ("--content-type", content_type, "--content-encoding", str(encoding), "--cache-control", cache)
    return Request("put-object", _key(prefix, RELEASES), (*options, *condition), body, size)


def updated_releases(current: dict[str, Any] | None, step: str, entry: dict[str, Any]) -> dict[str, Any] | None:
    """Return `releases.json` after a step, or None when the step changes nothing."""
    if current is None:
        if step == "unlist":
            raise ValueError(f"the data host has no {RELEASES}, so no release is listed")
        return {"default": entry["id"], "releases": [entry]}
    listed = any(other["id"] == entry["id"] for other in current["releases"])
    if step == "unlist":
        if entry["id"] == current["default"]:
            raise ValueError(f"{entry['id']} is the default release; promote another release before the unlist")
        if not listed:
            return None
        releases = [other for other in current["releases"] if other["id"] != entry["id"]]
        return {**current, "releases": releases}
    default = entry["id"] if step == "promote" else current["default"]
    if listed and default == current["default"]:
        return None
    releases = current["releases"] if listed else [*current["releases"], entry]
    releases = sorted(releases, key=lambda other: other["id"], reverse=True)
    return {**current, "default": default, "releases": releases}


def head_matches(head: dict[str, Any], sha256_hex: str) -> bool:
    """Tell whether the `head-object --checksum-mode ENABLED` answer of a key names the manifest's SHA-256."""
    return head.get("ChecksumSHA256") == _base64(sha256_hex)


def summary(requests: Sequence[Request]) -> str:
    """Return the count and cost of the requests, and the size and monthly storage cost of the uploaded bytes."""
    counts = Counter(request.operation for request in requests)
    cost = sum(PRICES[operation] * count for operation, count in counts.items())
    kinds = ", ".join(f"{count} {operation}" for operation, count in sorted(counts.items()))
    text = f"{len(requests)} requests ({kinds}): {_dollars(cost)}"
    stored = sum(request.size for request in requests if request.operation == "put-object")
    if stored:
        text += f"; {stored / 1e9:.2f} GB stored: {_dollars(stored / 1e9 * STORAGE_PRICE)} per month"
    return text


def release_publish(
    step: str,
    data_root: Path,
    release_id: str,
    bucket: str,
    prefix: str,
    *,
    dry_run: bool,
    jobs: int,
    zlib_version: str = ZLIB_VERSION,
    aws: Callable[[list[str]], subprocess.CompletedProcess[str]] | None = None,
    print_line: Callable[[str], None] = print,
) -> None:
    """Run one publish step for a release of a local data root; with `dry_run`, print its requests instead."""
    release = load_release(data_root, release_id)
    if step != "unlist":
        found = refusals(release, local_sizes(data_root, release), pinned_python(), zlib_version)
        if found:
            raise ValueError("the release is refused:\n" + "\n".join(found))
    run = aws or _aws
    if step == "upload":
        requests = upload_requests(release, prefix, data_root / release.id)
        if dry_run:
            _print_requests(requests, bucket, print_line)
            print_line("plus one head-object per key that exists, when an upload resumes")
            return
        _upload(requests, bucket, jobs, run)
    elif dry_run:
        _print_requests(releases_requests(prefix, release.id, step), bucket, print_line)
        print_line(f"the put-object takes --if-none-match '*' when the data host has no {RELEASES}")
    else:
        _change_releases(step, release, bucket, prefix, run)


def _upload(
    requests: list[Request],
    bucket: str,
    jobs: int,
    run: Callable[[list[str]], subprocess.CompletedProcess[str]],
) -> None:
    def send(request: Request) -> bool:
        result = run(request.argv(bucket))
        if result.returncode == 0:
            return True
        if "PreconditionFailed" not in result.stderr:
            raise RuntimeError(f"'{request.key}': {result.stderr.strip()}")
        head = run(
            ["aws", "s3api", "head-object", "--bucket", bucket, "--key", request.key, "--checksum-mode", "ENABLED"]
        )
        if head.returncode != 0:
            raise RuntimeError(f"'{request.key}': {head.stderr.strip()}")
        if not head_matches(json.loads(head.stdout), request.sha256):
            raise RuntimeError(
                f"'{request.key}' exists with other bytes than the manifest names; a path is written once"
            )
        return False

    *objects, manifest = requests
    with ThreadPoolExecutor(jobs) as pool:
        written = sum(pool.map(send, objects))
    written += send(manifest)
    LOG.info("uploaded %d objects, %d existed with the same bytes", written, len(requests) - written)


def _change_releases(
    step: str,
    release: Release,
    bucket: str,
    prefix: str,
    run: Callable[[list[str]], subprocess.CompletedProcess[str]],
) -> None:
    *check, read, _ = releases_requests(prefix, release.id, step)
    for request in check:
        result = run(request.argv(bucket))
        if result.returncode != 0:
            raise RuntimeError(
                f"'{request.key}' cannot be read, so the upload of {release.id} has not ended: {result.stderr.strip()}"
            )
    with tempfile.TemporaryDirectory() as scratch:
        target = Path(scratch) / RELEASES
        result = run(Request(read.operation, read.key, body=str(target)).argv(bucket))
        if result.returncode == 0:
            current, etag = json.loads(gzip.decompress(target.read_bytes())), json.loads(result.stdout)["ETag"]
        elif "NoSuchKey" in result.stderr:
            current, etag = None, None
        else:
            raise RuntimeError(f"'{read.key}': {result.stderr.strip()}")
        updated = updated_releases(current, step, release.entry)
        if updated is None:
            LOG.info("%s of %s changes nothing in %s", step, release.id, RELEASES)
            return
        stored = gzip_member(json_bytes(updated))
        target.write_bytes(stored)
        result = run(releases_put(prefix, str(target), len(stored), etag).argv(bucket))
    if result.returncode != 0:
        raise RuntimeError(f"'{read.key}' was not written, it may have changed since the read: {result.stderr.strip()}")
    LOG.info(
        "%s: default %s, releases %s", RELEASES, updated["default"], [entry["id"] for entry in updated["releases"]]
    )


def _put(key: str, obj: dict[str, Any], body: str) -> Request:
    options = ["--content-type", obj["contentType"], "--cache-control", obj["cacheControl"]]
    if obj["contentEncoding"]:
        options += ["--content-encoding", obj["contentEncoding"]]
    if obj["contentDisposition"]:
        options += ["--content-disposition", obj["contentDisposition"]]
    options += ["--checksum-algorithm", "SHA256", "--checksum-sha256", _base64(obj["sha256"]), "--if-none-match", "*"]
    return Request("put-object", key, tuple(options), body, obj["size"], obj["sha256"])


def _key(prefix: str, path: str) -> str:
    folder = prefix.strip("/")
    return f"{folder}/{path}" if folder else path


def _base64(sha256_hex: str) -> str:
    return base64.b64encode(bytes.fromhex(sha256_hex)).decode()


def _dollars(amount: float) -> str:
    return "$" + (f"{float(f'{amount:.3g}'):.7f}".rstrip("0").rstrip(".") or "0")


def _print_requests(requests: Sequence[Request], bucket: str, print_line: Callable[[str], None]) -> None:
    for request in requests:
        print_line(shlex.join(request.argv(bucket)))
    print_line(summary(requests))


def _aws(argv: list[str]) -> subprocess.CompletedProcess[str]:
    command = [*argv[:3], "--output", "json", "--no-cli-pager", *argv[3:]]
    return subprocess.run(command, check=False, capture_output=True, text=True)
