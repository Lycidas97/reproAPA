"""Read matrix-matched barcode columns without importing analysis backends."""

from __future__ import annotations

import gzip
from pathlib import Path


def load_barcode_whitelist(path: Path) -> list[str]:
    """Load a one-barcode-per-line matrix column list without changing order."""
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        barcodes = [line.strip() for line in handle if line.strip()]
    if not barcodes:
        raise ValueError("Barcode whitelist is empty")
    if len(barcodes) != len(set(barcodes)):
        raise ValueError("Barcode whitelist contains duplicate entries")
    if any(len(barcode.split()) != 1 for barcode in barcodes):
        raise ValueError("Barcode whitelist must contain one barcode per line")
    return barcodes
