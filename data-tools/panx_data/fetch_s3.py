"""Copy single objects of a bucket to a folder, each at its key, with read-only `aws s3api get-object` requests.

The AWS CLI takes its credentials from the caller's environment. Each object costs one GET request.
"""

import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

LOG = logging.getLogger(__name__)


def fetch_s3(bucket: str, keys: list[str], dest: Path, workers: int = 8) -> None:
    """Download `keys` of `bucket` to `dest/<key>`; fail when any request fails."""
    with ThreadPoolExecutor(max_workers=workers) as pool:
        failed = [key for key, ok in zip(keys, pool.map(lambda key: _get(bucket, key, dest / key), keys)) if not ok]
    if failed:
        raise RuntimeError(f"{len(failed)} of {len(keys)} objects failed, first: {failed[0]}")
    LOG.info("fetched %d objects of '%s' to '%s'", len(keys), bucket, dest)


def _get(bucket: str, key: str, target: Path) -> bool:
    target.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            "aws",
            "s3api",
            "get-object",
            "--bucket",
            bucket,
            "--key",
            key,
            str(target),
            "--output",
            "text",
            "--query",
            "ETag",
        ],
        check=False,
        capture_output=True,
        text=True,  # failures are collected and reported together
    )
    if result.returncode != 0:
        LOG.error("get-object '%s': %s", key, result.stderr.strip())
    return result.returncode == 0
