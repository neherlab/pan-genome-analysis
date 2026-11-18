"""Byte encodings of the data format: gzip members with a fixed header, ustar entries, and compact JSON.

The build writes the gzip header and trailer itself, because Python's `gzip.compress` writes another OS byte in other
Python versions. Deflate bytes depend on the zlib build, so the build runs only with the zlib of `ZLIB_VERSION`, the
zlib of the pinned Python of the data container.
"""

import json
import struct
import tarfile
import zlib
from typing import Any, BinaryIO

# zlib.ZLIB_RUNTIME_VERSION of the Python of .config/mise.toml
ZLIB_VERSION = "1.3.2"
# No flags, mtime 0, XFL 2 (maximum compression), OS 255 (unknown)
GZIP_HEADER = bytes.fromhex("1f8b08000000000002ff")
TAR_BLOCK = 512
TAR_END = bytes(2 * TAR_BLOCK)


def check_zlib() -> None:
    """Fail unless the zlib of this Python gives the deflate bytes of the pinned runtime."""
    if zlib.ZLIB_RUNTIME_VERSION != ZLIB_VERSION:
        raise RuntimeError(
            f"found zlib {zlib.ZLIB_RUNTIME_VERSION}, but a release needs zlib {ZLIB_VERSION}, the zlib of the Python of "
            "the data container; run the build as `dev/docker/run just release-build ...`"
        )


def gzip_member(data: bytes) -> bytes:
    """Return one gzip member of `data`: the fixed header, raw deflate at level 9, the CRC32, and the size."""
    compressor = _compressor()
    body = compressor.compress(data) + compressor.flush()
    return GZIP_HEADER + body + struct.pack("<II", zlib.crc32(data), len(data) & 0xFFFFFFFF)


class GzipWriter:
    """Write one gzip member to a binary stream, from the data of several calls."""

    def __init__(self, target: BinaryIO) -> None:
        self._target = target
        self._compressor = _compressor()
        self._crc = 0
        self._size = 0
        target.write(GZIP_HEADER)

    def write(self, data: bytes) -> None:
        self._crc = zlib.crc32(data, self._crc)
        self._size += len(data)
        self._target.write(self._compressor.compress(data))

    def close(self) -> None:
        self._target.write(self._compressor.flush() + struct.pack("<II", self._crc, self._size & 0xFFFFFFFF))


def tar_entry(name: str, data: bytes, mtime: int) -> bytes:
    """Return the ustar header and the padded content of a regular file: mode 0644, owner 0:0 without names."""
    info = tarfile.TarInfo(name)
    info.size = len(data)
    info.mtime = mtime
    info.mode = 0o644
    info.uid = info.gid = 0
    info.uname = info.gname = ""
    info.type = tarfile.REGTYPE
    padding = -len(data) % TAR_BLOCK
    return info.tobuf(tarfile.USTAR_FORMAT, "utf-8", "strict") + data + bytes(padding)


def json_bytes(value: Any) -> bytes:
    """Encode compact UTF-8 JSON with finite numbers only."""
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode()


def _compressor() -> "zlib._Compress":
    return zlib.compressobj(9, zlib.DEFLATED, -15, 8, zlib.Z_DEFAULT_STRATEGY)
