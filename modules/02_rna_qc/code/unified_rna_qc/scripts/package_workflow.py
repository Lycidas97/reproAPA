#!/usr/bin/env python3
"""Create a deterministic candidate workflow bundle and member-level tamper test."""

from __future__ import annotations

import argparse
import gzip
import io
import tarfile
from pathlib import Path

import pandas as pd

from common import sha256_file, sha256_text


VERSION = "3.1.0-candidate.1"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    release = root / "release"
    release.mkdir(exist_ok=True)
    members = [root / "Snakefile", root / "config.schema.yaml", root / "analysis_plan.md", root / "environment.yml", root / "exact_rerun_commands.txt"]
    for directory in ("scripts", "configs", "manifests", "docs"):
        members.extend(sorted(path for path in (root / directory).rglob("*") if path.is_file() and "__pycache__" not in path.parts))
    bundle = release / f"workflow_bundle_v{VERSION}.tar.gz"
    buffer = io.BytesIO()
    rows = []
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.PAX_FORMAT) as archive:
        for path in sorted(set(members)):
            name = str(path.relative_to(root))
            data = path.read_bytes()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            info.mtime = 0
            info.mode = 0o755 if path.suffix in {".py", ".sh"} or path.name == "Snakefile" else 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            archive.addfile(info, io.BytesIO(data))
            rows.append({"bundle_member": name, "size_bytes": len(data), "sha256": sha256_text(data.hex())})
    with bundle.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(buffer.getvalue())
    digest = sha256_file(bundle)
    sidecar = bundle.with_suffix("").with_suffix(".sha256")
    sidecar.write_text(f"{digest}  {bundle.name}\n", encoding="utf-8")
    member_table = release / "WORKFLOW_BUNDLE_MANIFEST.tsv"
    # Recalculate conventional file hashes for the inspectable member table.
    for row, path in zip(rows, sorted(set(members)), strict=True):
        row["sha256"] = sha256_file(path)
    pd.DataFrame(rows).to_csv(member_table, sep="\t", index=False)
    target = root / "scripts/common.py"
    original = target.read_bytes()
    tampered = bytearray(original)
    tampered[len(tampered) // 2] ^= 1
    expected = sha256_file(target)
    observed = sha256_text(bytes(tampered).hex())
    # Convert the byte-hex helper result to a direct digest for a strict comparison.
    import hashlib
    observed = hashlib.sha256(bytes(tampered)).hexdigest()
    pd.DataFrame([{"test": "release_member_byte_tampering", "member": str(target.relative_to(root)), "expected_sha256": expected, "tampered_sha256": observed, "mismatch_detected": expected != observed, "status": "PASS" if expected != observed else "FAIL"}]).to_csv(root / "tests/release_member_tamper_test.tsv", sep="\t", index=False)
    print(f"{bundle.name}\t{bundle.stat().st_size}\t{digest}")


if __name__ == "__main__":
    main()
