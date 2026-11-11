"""Gzip bytes identical to those of `pigz` 2.8 with zlib 1.3 at the default level, the compressor of the 2022 upload.

pigz writes a gzip header with the file name and the modification time, and compresses an input of up to one block
(128 KiB) as a single raw deflate stream at level 6. A larger input is split into blocks that pigz compresses
independently, so this module hands such inputs to the `pigz` program itself.
"""

import struct
import subprocess
import zlib
from pathlib import Path

PIGZ_BLOCK_SIZE = 128 * 1024
REQUIRED_PIGZ = "pigz 2.8"
REQUIRED_ZLIB = "1.3"


def pigz_bytes(data: bytes, path: Path, mtime: int) -> bytes:
    """Return the bytes that `pigz -k` writes for `path`, whose content is `data` and modification time `mtime`."""
    if len(data) > PIGZ_BLOCK_SIZE:
        return subprocess.run(["pigz", "-c", "--", str(path)], capture_output=True, check=True).stdout
    header = b"\x1f\x8b\x08\x08" + struct.pack("<I", mtime) + b"\x00\x03" + path.name.encode() + b"\x00"
    compressor = zlib.compressobj(6, zlib.DEFLATED, -15, 8)
    body = compressor.compress(data) + compressor.flush()
    return header + body + struct.pack("<II", zlib.crc32(data), len(data) & 0xFFFFFFFF)


def check_compressors() -> None:
    """Fail unless pigz and zlib have the versions whose output the bucket holds."""
    pigz = pigz_version()
    if pigz != REQUIRED_PIGZ:
        raise RuntimeError(f"found '{pigz}', but byte-exact .gz copies need '{REQUIRED_PIGZ}'")
    if not zlib.ZLIB_RUNTIME_VERSION.startswith(REQUIRED_ZLIB):
        raise RuntimeError(
            f"found zlib {zlib.ZLIB_RUNTIME_VERSION}, but byte-exact .gz copies need zlib {REQUIRED_ZLIB}"
        )


def pigz_version() -> str:
    """Return the version line of the installed pigz, which prints it to stdout or stderr."""
    result = subprocess.run(["pigz", "--version"], capture_output=True, text=True, check=True)
    return (result.stdout + result.stderr).strip()
