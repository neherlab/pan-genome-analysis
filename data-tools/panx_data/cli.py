"""Command line: `rebuild` and `verify` work on a snapshot folder with `src/`; `fetch-s3` copies S3 objects into it;
`release-build` and `release-check` write and check a release of the data format; `release-publish` publishes it."""

import argparse
import json
import logging
import os
import platform
import subprocess
import sys
import time
import zlib
from datetime import UTC, datetime
from pathlib import Path

from panx_data.fetch_s3 import fetch_s3
from panx_data.layout import BUCKET, Out, find_source
from panx_data.listings import load_listings
from panx_data.pigz_gzip import pigz_version
from panx_data.publication import STEPS, release_publish
from panx_data.rebuild import rebuild
from panx_data.release import DATA_ROOT, ToolVersion, default_jobs, release_build
from panx_data.release_check import release_check
from panx_data.verify import verify

LOG = logging.getLogger("panx_data")


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s", stream=sys.stderr)
    return args.run(args)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="panx-data", description=__doc__)
    commands = parser.add_subparsers(required=True)

    build = commands.add_parser("rebuild", help="build out/ from src/ and write out/BUILD_INFO.json")
    build.add_argument("snapshot", type=Path, help="folder that holds src/")
    build.add_argument("--out", type=Path, help="write here instead of <snapshot>/out")
    build.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 4), help="datasets built in parallel")
    build.add_argument("--dataset", action="append", dest="datasets", help="build only this dataset (repeatable)")
    build.set_defaults(run=_rebuild)

    check = commands.add_parser(
        "verify", help="compare out/ and src/s3/ with the S3 listings; exit 1 on unknown differences"
    )
    check.add_argument("snapshot", type=Path, help="folder that holds src/")
    check.add_argument("--out", type=Path, help="read here instead of <snapshot>/out")
    check.set_defaults(run=_verify)

    fetch = commands.add_parser("fetch-s3", help="copy S3 objects to a folder, each at its key")
    fetch.add_argument("--bucket", default=BUCKET)
    fetch.add_argument("--keys", type=Path, required=True, help="file with one key per line")
    fetch.add_argument("--dest", type=Path, required=True, help="for example <snapshot>/src/s3/<bucket>")
    fetch.set_defaults(run=_fetch)

    release = commands.add_parser(
        "release-build", help="build the release of a recipe into a data root and list it in releases.json"
    )
    release.add_argument("recipe", type=Path, help="releases/<release id>.toml")
    release.add_argument("snapshot", type=Path, help="folder that holds src/")
    release.add_argument("--out", type=Path, help=f"data root (default <snapshot>/{DATA_ROOT})")
    release.add_argument("--jobs", type=int, default=default_jobs(), help="datasets built in parallel")
    release.add_argument("--dataset", action="append", dest="datasets", help="build only this dataset (repeatable)")
    release.set_defaults(run=_release_build)

    release_test = commands.add_parser(
        "release-check",
        help="compare the releases of a data root with the snapshot; print findings, exit 1 on a failing finding",
    )
    release_test.add_argument("data_root", type=Path, help="folder that holds releases.json")
    release_test.add_argument("snapshot", type=Path, help="folder that holds src/")
    release_test.add_argument("--jobs", type=int, default=default_jobs(), help="datasets checked in parallel")
    release_test.set_defaults(run=_release_check)

    publish = commands.add_parser(
        "release-publish",
        help="publish a release of a local data root to the data host with the AWS CLI, one step at a time",
        description="Steps in order: upload, list, promote. A rollback is a promote of an earlier release; unlist "
        "takes a release that is not the default out of releases.json. --dry-run prints the requests, their count, "
        "and their cost, and sends nothing.",
    )
    publish.add_argument("step", choices=STEPS)
    publish.add_argument("data_root", type=Path, help="folder that holds releases.json and the release folder")
    publish.add_argument("release", help="release id, such as 20261112T080000Z")
    publish.add_argument("--bucket", required=True, help="data bucket")
    publish.add_argument("--prefix", required=True, help="folder of the data host in the bucket")
    publish.add_argument("--dry-run", action="store_true", help="print the requests and their cost, send nothing")
    publish.add_argument("--jobs", type=int, default=8, help="uploads at the same time")
    publish.set_defaults(run=_release_publish)
    return parser


def _rebuild(args: argparse.Namespace) -> int:
    started = datetime.now(UTC)
    start = time.monotonic()
    source = find_source(args.snapshot / "src")
    out = Out(args.out or args.snapshot / "out")
    listings = load_listings(source.listings, BUCKET)
    results = rebuild(source, out, listings, args.jobs, frozenset(args.datasets) if args.datasets else None)
    info = {
        "tool": _tool_version(),
        "python": platform.python_version(),
        "zlib": zlib.ZLIB_RUNTIME_VERSION,
        "pigz": pigz_version(),
        "datasets": len(results) - 1,
        "files": sum(result.files for result in results),
        "bytes": sum(result.bytes for result in results),
        "jobs": args.jobs,
        "only": sorted(args.datasets or []),
        "started": started.isoformat(timespec="seconds"),
        "seconds": round(time.monotonic() - start),
    }
    out.build_info.write_text(json.dumps(info, indent=2) + "\n")
    LOG.info(
        "built %d datasets, %d files, %.1f GB in %d s",
        info["datasets"],
        info["files"],
        info["bytes"] / 1e9,
        info["seconds"],
    )
    return 0


def _verify(args: argparse.Namespace) -> int:
    source = find_source(args.snapshot / "src")
    out = Out(args.out or args.snapshot / "out")
    counts = verify(source, out, load_listings(source.listings, BUCKET))
    for kind, count in sorted(counts.items(), key=lambda item: item[0].value):
        LOG.info("%s %s: %d", "known" if kind.expected else "UNKNOWN", kind.value, count)
    unknown = sum(count for kind, count in counts.items() if not kind.expected)
    LOG.info("%s: %d differences, %d of them unknown", out.deviations, sum(counts.values()), unknown)
    return 1 if unknown else 0


def _release_build(args: argparse.Namespace) -> int:
    start = time.monotonic()
    version = _tool_version()
    target = release_build(
        args.recipe,
        args.snapshot,
        args.out or args.snapshot / DATA_ROOT,
        args.jobs,
        ToolVersion(str(version["commit"]), bool(version["uncommitted_changes"])),
        platform.python_version(),
        zlib.ZLIB_RUNTIME_VERSION,
        frozenset(args.datasets) if args.datasets else None,
    )
    LOG.info("wrote %s in %d s", target, round(time.monotonic() - start))
    return 0


def _release_check(args: argparse.Namespace) -> int:
    findings = release_check(args.data_root, args.snapshot, args.jobs)
    for finding in findings:
        print(f"{finding.release}/{finding.dataset}\t{finding.kind.label}\t{finding.detail}")
    failing = sum(finding.kind.fails for finding in findings)
    LOG.info("%d findings, %d of them failing", len(findings), failing)
    return 1 if failing else 0


def _release_publish(args: argparse.Namespace) -> int:
    release_publish(
        args.step, args.data_root, args.release, args.bucket, args.prefix, dry_run=args.dry_run, jobs=args.jobs
    )
    return 0


def _fetch(args: argparse.Namespace) -> int:
    keys = [line.strip() for line in args.keys.read_text().splitlines() if line.strip()]
    fetch_s3(args.bucket, keys, args.dest)
    return 0


def _tool_version() -> dict[str, str | bool]:
    """Name the commit of this tool, and whether its folder has uncommitted changes.

    Outside a git checkout the fields stay empty, so the git calls may fail.
    """
    folder = Path(__file__).resolve().parent.parent
    commit = subprocess.run(
        ["git", "-C", str(folder), "rev-parse", "HEAD"], check=False, capture_output=True, text=True
    )
    status = subprocess.run(
        ["git", "-C", str(folder), "status", "--porcelain", "--", "."], check=False, capture_output=True, text=True
    )
    remote = subprocess.run(
        ["git", "-C", str(folder), "remote", "get-url", "origin"], check=False, capture_output=True, text=True
    )
    return {
        "repository": remote.stdout.strip(),
        "commit": commit.stdout.strip(),
        "uncommitted_changes": bool(status.stdout.strip()),
    }
