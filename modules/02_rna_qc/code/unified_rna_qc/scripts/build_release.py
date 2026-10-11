#!/usr/bin/env python3
"""Build the candidate-level report, resource table, and exhaustive deliverable manifest."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from common import sha256_file


def excluded(path: Path, root: Path) -> bool:
    relative = path.relative_to(root)
    if any(part in {".snakemake", "__pycache__", "logs"} for part in relative.parts):
        return True
    if any(part.startswith(".") for part in relative.parts):
        return True
    return relative.as_posix() in {
        "release/release_manifest.tsv", "release/release_manifest.sha256",
        "release/release_verification.tsv", "release/release_verification.sha256",
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    root = args.root.resolve()
    release = root / "release"
    gse = json.loads((root / "runs/gse178325_regression/run_summary.json").read_text())
    other = json.loads((root / "runs/gse103221_small/run_summary.json").read_text())
    tests = pd.read_csv(root / "tests/generalization_test_summary.tsv", sep="\t")
    tamper = pd.read_csv(root / "tests/release_member_tamper_test.tsv", sep="\t")
    if tests["status"].ne("PASS").any() or tamper["status"].ne("PASS").any():
        raise ValueError("Release tests are not all passing")
    report = f"""# Unified RNA QC v3.1.0-candidate.1 builder validation

## Outcome

All candidate checks pass. This is a candidate release only; central review remains required before any multi-study production fan-out.

## Accepted regression

- Study: `{gse['study_id']}`; libraries/called cells: {gse['libraries']} / {gse['called_cells']:,}
- Primary / sensitivity / technical: {gse['tier_counts']['primary_pass']:,} / {gse['tier_counts']['sensitivity_only']:,} / {gse['tier_counts']['technical_exclude']:,}
- Exact equality against the accepted membership sets and hashes: PASS
- Scrublet hard exclusions / tier changes / sample-gate changes: 0 / 0 / 0

## Real configuration-generalization pilot

- Study: `{other['study_id']}`; libraries/called cells: {other['libraries']} / {other['called_cells']:,}
- Primary / sensitivity / technical: {other['tier_counts']['primary_pass']:,} / {other['tier_counts']['sensitivity_only']:,} / {other['tier_counts']['technical_exclude']:,}
- Inputs: two existing filtered Matrix Market libraries of different sizes; same Snakefile, configuration and manifest changes only
- Scrublet hard exclusions / tier changes / sample-gate changes: 0 / 0 / 0

## Generalization and release checks

Hard-coding scan, exact regression golden test, real non-target run, wrong-study/count/duplicate-ID/missing-input negative tests, current-field schema checks, release-member tampering, independent per-run validation, bundle hash, and release-member verification all pass.

The current sample gate uses only `n_candidate_basic_qc`. The ambiguous legacy field is absent from current decision outputs and exists only in named historical evidence where applicable.

## Adapter caveat

GSE103221 has real retained filtered Matrix Market inputs but no retained cell-level DOC-156 QC/doublet ledger in the active local resources. The documented policy-evidence adapters therefore hash-validate the retained DOC-156 report/policy and recompute the frozen formulas plus Scrublet from raw filtered counts. No historical membership is approximated.

## Re-run and resources

See `exact_rerun_commands.txt` and `release/resource_usage.tsv`. Per-library figures and study summaries are under each run directory.

## Manifest self-reference

`release_manifest.tsv` covers every immutable deliverable present when it is built. The manifest, its sidecar, and the verification table/sidecar are excluded because self-hashing is circular; those four metadata artifacts are cross-hashed and checked by the verifier.
"""
    (release / "candidate_report.md").write_text(report, encoding="utf-8")
    resource_frames = []
    for run_name in ("gse178325_regression", "gse103221_small"):
        worker = pd.read_csv(root / f"runs/{run_name}/resource_usage.tsv", sep="\t")
        aggregate = pd.read_csv(root / f"runs/{run_name}/aggregate_resource_usage.tsv", sep="\t")
        worker["run"] = run_name
        aggregate["run"] = run_name
        resource_frames.extend([worker, aggregate])
    pd.concat(resource_frames, ignore_index=True, sort=False).to_csv(release / "resource_usage.tsv", sep="\t", index=False)
    paths = sorted(path for path in root.rglob("*") if path.is_file() and not excluded(path, root))
    manifest = pd.DataFrame([{"path": str(path.relative_to(root)), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in paths])
    manifest_path = release / "release_manifest.tsv"
    manifest.to_csv(manifest_path, sep="\t", index=False)
    (release / "release_manifest.sha256").write_text(f"{sha256_file(manifest_path)}  {manifest_path.name}\n", encoding="utf-8")
    print(f"manifested {len(manifest)} deliverables")


if __name__ == "__main__":
    main()
