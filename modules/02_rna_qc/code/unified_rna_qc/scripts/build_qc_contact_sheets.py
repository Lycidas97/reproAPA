#!/usr/bin/env python3
"""Build fixed-size contact sheets for visual review of per-library QC plots."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

from PIL import Image, ImageOps


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--index", type=Path, required=True)
    parser.add_argument("--per-sheet", type=int, default=16)
    args = parser.parse_args()
    paths = sorted(args.run_dir.glob("per_library/*/qc_diagnostic.png"))
    if not paths:
        raise ValueError("No per-library QC figures found")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tile = (600, 450)
    columns = 4
    rows_written = []
    for page, start in enumerate(range(0, len(paths), args.per_sheet), 1):
        batch = paths[start:start + args.per_sheet]
        rows = math.ceil(len(batch) / columns)
        canvas = Image.new("RGB", (columns * tile[0], rows * tile[1]), "white")
        for index, path in enumerate(batch):
            with Image.open(path) as source:
                image = ImageOps.contain(source.convert("RGB"), tile)
            x = (index % columns) * tile[0] + (tile[0] - image.width) // 2
            y = (index // columns) * tile[1] + (tile[1] - image.height) // 2
            canvas.paste(image, (x, y))
        output = args.output_dir / f"contact_sheet_{page:02d}.png"
        canvas.save(output, dpi=(150, 150))
        rows_written.append((page, output, len(batch), batch[0], batch[-1]))
    args.index.parent.mkdir(parents=True, exist_ok=True)
    with args.index.open("w", encoding="utf-8") as handle:
        handle.write("page\tpath\tn_libraries\tfirst_figure\tlast_figure\n")
        for page, output, count, first, last in rows_written:
            handle.write(f"{page}\t{output}\t{count}\t{first}\t{last}\n")


if __name__ == "__main__":
    main()
