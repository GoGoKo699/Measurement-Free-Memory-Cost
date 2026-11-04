#!/usr/bin/env python3
"""Compare a fresh report build with its checked PDF, allowing PDF identifier drift."""
import argparse
import hashlib
import json
from pathlib import Path

import fitz


def digest(data):
    return hashlib.sha256(data).hexdigest()


def links(page):
    result = []
    for link in page.get_links():
        value = {}
        for key, item in link.items():
            if key in {"xref", "id"}:
                continue
            if isinstance(item, (fitz.Rect, fitz.Point)):
                item = tuple(item)
            value[key] = item
        result.append(value)
    return result


def compare(expected, fresh):
    left, right = fitz.open(expected), fitz.open(fresh)
    if len(left) != len(right):
        raise ValueError(f"Page count differs: {len(left)} versus {len(right)}")
    ignored = {"creationDate", "modDate"}
    metadata = lambda doc: {k: v for k, v in doc.metadata.items() if k not in ignored}
    if metadata(left) != metadata(right):
        raise ValueError("Report metadata differs beyond build timestamps")
    if left.metadata.get("author") or right.metadata.get("author"):
        raise ValueError("The explanatory report must have empty author metadata")
    pages = []
    for number, (a, b) in enumerate(zip(left, right), 1):
        pa = a.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), alpha=False)
        pb = b.get_pixmap(matrix=fitz.Matrix(1.6, 1.6), alpha=False)
        checks = {
            "page": number,
            "geometry_identical": tuple(a.rect) == tuple(b.rect),
            "text_identical": a.get_text() == b.get_text(),
            "decoded_content_streams_identical": a.read_contents() == b.read_contents(),
            "links_identical": links(a) == links(b),
            "rendered_pixels_identical": (pa.width, pa.height, pa.samples)
            == (pb.width, pb.height, pb.samples),
        }
        if not all(value for key, value in checks.items() if key != "page"):
            raise ValueError(f"Report page {number} differs: {checks}")
        pages.append(checks)
    return {
        "status": "PASS",
        "expected_sha256": digest(expected.read_bytes()),
        "fresh_sha256": digest(fresh.read_bytes()),
        "pdf_bytes_identical": expected.read_bytes() == fresh.read_bytes(),
        "pymupdf_version": fitz.VersionBind,
        "pages": pages,
        "scope": "Exact page content, text, geometry, links and raster comparison. PDF document identifiers and build timestamps may differ; this is not a scientific proof review.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("expected", type=Path)
    parser.add_argument("fresh", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = compare(args.expected, args.fresh)
    text = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
