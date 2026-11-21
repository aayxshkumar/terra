#!/usr/bin/env python3
"""Recompress terra's stored-deflate PNGs (standard library only).

terra writes PNGs with an uncompressed zlib stream so the engine needs no
zlib. This rewrites them with real compression and the Paeth-free 'up'
filter, which is plenty for the docs images.

  python3 tools/pngshrink.py in.png [out.png]
"""
import struct
import sys
import zlib


def chunks(data):
    pos = 8
    while pos < len(data):
        (n,) = struct.unpack(">I", data[pos:pos + 4])
        kind = data[pos + 4:pos + 8]
        yield kind, data[pos + 8:pos + 8 + n]
        pos += 12 + n


def chunk(kind, body):
    return struct.pack(">I", len(body)) + kind + body + struct.pack(">I", zlib.crc32(kind + body) & 0xFFFFFFFF)


def main():
    if len(sys.argv) not in (2, 3):
        print(__doc__)
        return 2
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) == 3 else src
    data = open(src, "rb").read()
    ihdr, idat = None, b""
    for kind, body in chunks(data):
        if kind == b"IHDR":
            ihdr = body
        elif kind == b"IDAT":
            idat += body
    w, h, depth, ctype = struct.unpack(">IIBB", ihdr[:10])
    if depth != 8 or ctype != 2:
        print("only 8-bit RGB is supported")
        return 1
    raw = zlib.decompress(idat)
    stride = w * 3
    rows, prev = [], bytes(stride)
    for y in range(h):
        line = raw[y * (stride + 1) + 1:(y + 1) * (stride + 1)]
        up = bytes((a - b) & 0xFF for a, b in zip(line, prev))
        rows.append(b"\x02" + up)
        prev = line
    out = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) + chunk(b"IDAT", zlib.compress(b"".join(rows), 9)) + chunk(b"IEND", b"")
    open(dst, "wb").write(out)
    print(f"{src}: {len(data)} -> {len(out)} bytes")
    return 0


if __name__ == "__main__":
    sys.exit(main())
