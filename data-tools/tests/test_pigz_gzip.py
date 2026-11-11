import gzip
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from panx_data.pigz_gzip import PIGZ_BLOCK_SIZE, pigz_bytes

pytestmark = pytest.mark.skipif(shutil.which("pigz") is None, reason="needs the pigz program")


@pytest.mark.parametrize(
    "size", [0, 1, 4000, PIGZ_BLOCK_SIZE - 1, PIGZ_BLOCK_SIZE, PIGZ_BLOCK_SIZE + 1, 3 * PIGZ_BLOCK_SIZE]
)
def test_pigz_gzip_bytes_equal_pigz_keep_output(tmp_path: Path, size: int) -> None:
    data = (b">strain_1\nACGTTGCA-NNACGT\n" * (size // 26 + 1))[:size]
    path = tmp_path / "abc_refined_aa_aln.fa"
    path.write_bytes(data)
    os.utime(path, (1_657_000_000, 1_657_000_000))
    subprocess.run(["pigz", "-k", "--", str(path)], check=True)

    assert pigz_bytes(data, path, 1_657_000_000) == path.with_name(path.name + ".gz").read_bytes()


def test_pigz_gzip_copy_decompresses_to_the_input(tmp_path: Path) -> None:
    data = os.urandom(1000) * 300
    path = tmp_path / "geneCluster.json"
    path.write_bytes(data)

    assert gzip.decompress(pigz_bytes(data, path, 0)) == data
