#!/usr/bin/env python3
"""Verify every immutable release member and cross-check release metadata artifacts."""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from common import sha256_file


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    release = root / "release"
    manifest_path = release / "release_manifest.tsv"
    expected_manifest_sha = (release / "release_manifest.sha256").read_text().split()[0]
    if sha256_file(manifest_path) != expected_manifest_sha:
        raise ValueError("Release manifest sidecar mismatch")
    manifest = pd.read_csv(manifest_path, sep="\t")
    rows = []
    for item in manifest.itertuples(index=False):
        path = root / item.path
        exists = path.is_file()
        size = path.stat().st_size if exists else -1
        digest = sha256_file(path) if exists else ""
        valid = exists and size == int(item.size_bytes) and digest == item.sha256
        rows.append({"path": item.path, "exists": exists, "expected_size_bytes": item.size_bytes, "observed_size_bytes": size, "expected_sha256": item.sha256, "observed_sha256": digest, "status": "PASS" if valid else "FAIL"})
    table = pd.DataFrame(rows)
    if table["status"].ne("PASS").any():
        raise ValueError("Release member verification failed")
    table_path = release / "release_verification.tsv"
    table.to_csv(table_path, sep="\t", index=False)
    (release / "release_verification.sha256").write_text(f"{sha256_file(table_path)}  {table_path.name}\n", encoding="utf-8")
    print(f"PASS: {len(table)} release members")


if __name__ == "__main__":
    main()
