"""Command line: `rebuild` and `verify` work on a snapshot folder with `src/`; `fetch-s3` copies S3 objects into it."""

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
from panx_data.rebuild import rebuild
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
