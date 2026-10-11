#!/usr/bin/env python3
"""Package compact validation evidence for issue delivery without cell-level ledgers."""

from __future__ import annotations

import argparse
import gzip
import io
import tarfile
from pathlib import Path

from common import sha256_file


VERSION = "3.1.0-candidate.1"


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    members = []
    for relative in ["configs", "manifests", "docs", "tests"]:
        members.extend(path for path in (root / relative).rglob("*") if path.is_file() and not any(part.startswith(".") for part in path.relative_to(root).parts))
    members.extend([root / "analysis_plan.md", root / "exact_rerun_commands.txt", root / "release/WORKFLOW_BUNDLE_MANIFEST.tsv"])
    for run_name in ("gse178325_regression", "gse103221_small"):
        run = root / "runs" / run_name
        names = [
            "run_report.md", "run_summary.json", "tier_counts_and_hashes.tsv",
            "tier_counts_by_library.tsv", "sample_qc_decisions.tsv",
            "diagnostic_summary.tsv", "validation_summary.tsv",
            "independent_validation.tsv", "regression_membership_comparison.tsv",
            "input_inventory.tsv", "resource_usage.tsv", "aggregate_resource_usage.tsv",
            "figures/study_qc_summary.png", "figures/study_qc_summary.pdf",
        ]
        members.extend(run / name for name in names if (run / name).is_file())
        members.extend(sorted((run / "per_library").glob("*/qc_diagnostic.png")))
    bundle = root / "release" / f"validation_evidence_v{VERSION}.tar.gz"
    buffer = io.BytesIO()
    with tarfile.open(fileobj=buffer, mode="w", format=tarfile.PAX_FORMAT) as archive:
        for path in sorted(set(members)):
            data = path.read_bytes()
            info = tarfile.TarInfo(str(path.relative_to(root)))
            info.size = len(data)
            info.mtime = 0
            info.mode = 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            archive.addfile(info, io.BytesIO(data))
    with bundle.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(buffer.getvalue())
    digest = sha256_file(bundle)
    sidecar = bundle.with_suffix("").with_suffix(".sha256")
    sidecar.write_text(f"{digest}  {bundle.name}\n", encoding="utf-8")
    print(f"{bundle.name}\t{bundle.stat().st_size}\t{digest}")


if __name__ == "__main__":
    main()
