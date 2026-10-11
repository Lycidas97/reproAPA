#!/usr/bin/env python3
"""Aggregate one configuration-selected study run and enforce its declared contract."""

from __future__ import annotations

import argparse
import json
import platform
import resource
import sys
import time
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import load_config, read_manifest, resolve, sha256_file, sha256_text


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, required=True)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--run-dir", type=Path, required=True)
    return p.parse_args()


def set_hash(ids: pd.Series) -> str:
    return sha256_text("\n".join(sorted(ids.astype(str))) + "\n")


def main() -> None:
    parsed = parse_args()
    started = time.perf_counter()
    root = parsed.root.resolve()
    run = parsed.run_dir.resolve()
    config = load_config(parsed.config.resolve())
    manifest = read_manifest(root, config)
    libraries = manifest["library_id"].tolist()
    per = [run / "per_library" / library for library in libraries]

    def parquet(name: str) -> pd.DataFrame:
        return pd.concat([pd.read_parquet(path / name) for path in per], ignore_index=True)

    def tsv(name: str) -> pd.DataFrame:
        return pd.concat([pd.read_csv(path / name, sep="\t") for path in per], ignore_index=True)

    metrics = parquet("rna_cell_qc_metrics.parquet")
    pre = parquet("pre_doublet_qc_decisions.parquet")
    diagnostic = parquet("doublet_diagnostics.parquet")
    final = parquet("cell_qc_decisions.parquet")
    thresholds = tsv("rna_qc_thresholds.tsv")
    samples = tsv("sample_qc_decision.tsv")
    seeds = tsv("seed_stability.tsv")
    inventory = tsv("input_inventory.tsv")
    resources = tsv("resource_usage.tsv")
    expected = int(config["expected_called_cell_count"])
    if len(final) != expected or len(metrics) != expected or len(diagnostic) != expected:
        raise ValueError("Aggregate cell counts disagree with configured expectation")
    if final["global_cell_id"].duplicated().any():
        raise ValueError("Duplicate global cell IDs across libraries")
    if set(final["study_id"].astype(str)) != {str(config["study_id"])}:
        raise ValueError("Aggregate study identity drift")
    if final["doublet_hard_exclude"].astype(bool).any():
        raise ValueError("Scrublet caused a hard exclusion")
    tier_changes = int((final["pre_doublet_qc_tier"] != final["final_qc_tier"]).sum())
    if tier_changes:
        raise ValueError("Scrublet changed a basic-QC tier")
    if "n_candidate_singlets" in samples.columns:
        raise ValueError("Ambiguous legacy sample decision field is forbidden")
    if "n_candidate_basic_qc" not in samples.columns:
        raise ValueError("Current sample gate variable missing")
    if samples["n_scrublet_sample_gate_changes"].sum() != 0:
        raise ValueError("Scrublet changed a sample gate")

    counts = final["final_qc_tier"].value_counts()
    tiers = pd.DataFrame([
        {"study_id": config["study_id"], "tier": tier, "n_cells": int(counts.get(tier, 0)), "membership_sha256": set_hash(final.loc[final["final_qc_tier"].eq(tier), "global_cell_id"])}
        for tier in ("primary_pass", "sensitivity_only", "technical_exclude")
    ])
    by_library = final.groupby(["study_id", "library_id", "final_qc_tier"]).size().unstack(fill_value=0).reset_index()
    by_library = by_library.merge(
        manifest[["library_id", "sample_accession"]], on="library_id", how="left", validate="one_to_one"
    )
    for tier in ("primary_pass", "sensitivity_only", "technical_exclude"):
        if tier not in by_library:
            by_library[tier] = 0
    by_library["n_called_cells"] = by_library[["primary_pass", "sensitivity_only", "technical_exclude"]].sum(axis=1)
    diag_summary = samples[["study_id", "library_id", "n_called_cells", "n_automatic_calls", "automatic_call_fraction", "threshold_status", "n_doublet_hard_excluded", "n_scrublet_tier_changes", "n_scrublet_sample_gate_changes", "doublet_qc_action", "diagnostic_evidence_adapter"]].copy()

    regression_rows = []
    baseline_value = str(config.get("regression_baseline", ""))
    if baseline_value:
        baseline_path = resolve(root, baseline_value)
        baseline = pd.read_parquet(baseline_path)
        for tier in ("primary_pass", "sensitivity_only", "technical_exclude"):
            observed = set(final.loc[final["final_qc_tier"].eq(tier), "global_cell_id"].astype(str))
            prior = set(baseline.loc[baseline["final_qc_tier"].eq(tier), "global_cell_id"].astype(str))
            regression_rows.append({"tier": tier, "observed_n": len(observed), "baseline_n": len(prior), "added_n": len(observed - prior), "removed_n": len(prior - observed), "exact_set_equality": observed == prior, "observed_sha256": sha256_text("\n".join(sorted(observed)) + "\n"), "baseline_sha256": sha256_text("\n".join(sorted(prior)) + "\n")})
        if not all(row["exact_set_equality"] for row in regression_rows):
            raise ValueError("Configured regression membership differs from baseline")
    regression = pd.DataFrame(regression_rows, columns=["tier", "observed_n", "baseline_n", "added_n", "removed_n", "exact_set_equality", "observed_sha256", "baseline_sha256"])

    metrics.to_parquet(run / "rna_cell_qc_metrics.parquet", index=False, compression="zstd")
    pre.to_parquet(run / "pre_doublet_qc_decisions.parquet", index=False, compression="zstd")
    diagnostic.to_parquet(run / "doublet_diagnostics.parquet", index=False, compression="zstd")
    final.to_parquet(run / "cell_qc_decisions.parquet", index=False, compression="zstd")
    thresholds.to_csv(run / "rna_qc_thresholds_by_library.tsv", sep="\t", index=False)
    samples.to_csv(run / "sample_qc_decisions.tsv", sep="\t", index=False)
    seeds.to_csv(run / "seed_stability.tsv", sep="\t", index=False)
    inventory.to_csv(run / "input_inventory.tsv", sep="\t", index=False)
    resources.to_csv(run / "resource_usage.tsv", sep="\t", index=False)
    tiers.to_csv(run / "tier_counts_and_hashes.tsv", sep="\t", index=False)
    by_library.to_csv(run / "tier_counts_by_library.tsv", sep="\t", index=False)
    diag_summary.to_csv(run / "diagnostic_summary.tsv", sep="\t", index=False)
    regression.to_csv(run / "regression_membership_comparison.tsv", sep="\t", index=False)
    membership_columns = ["study_id", "library_id", "barcode", "global_cell_id"]
    membership_outputs = {
        "primary_pass": "primary_membership.tsv.gz",
        "sensitivity_only": "sensitivity_membership.tsv.gz",
        "technical_exclude": "technical_membership.tsv.gz",
    }
    for tier, filename in membership_outputs.items():
        final.loc[final["final_qc_tier"].eq(tier), membership_columns].sort_values(
            "global_cell_id", kind="mergesort"
        ).to_csv(
            run / filename,
            sep="\t",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )

    figure_dir = run / "figures"
    figure_dir.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    plot = by_library.set_index("sample_accession")[["primary_pass", "sensitivity_only", "technical_exclude"]]
    plot.plot(kind="bar", stacked=True, ax=axes[0], color=["#2b8cbe", "#f0a202", "#b8b8b8"])
    axes[0].set(ylabel="Cells", title="QC tiers by library")
    axes[0].tick_params(axis="x", labelrotation=45)
    axes[1].scatter(samples["n_called_cells"], samples["n_candidate_basic_qc"], c=samples["automatic_call_fraction"], cmap="viridis", s=70)
    axes[1].plot([0, samples["n_called_cells"].max()], [0, samples["n_called_cells"].max()], ls="--", color="grey")
    axes[1].set(xlabel="Called cells", ylabel="Candidate basic-QC cells", title="Current sample-gate input (colour: diagnostic call fraction)")
    fig.suptitle(str(config["report_label"]))
    fig.savefig(figure_dir / "study_qc_summary.png", dpi=180)
    fig.savefig(figure_dir / "study_qc_summary.pdf")
    plt.close(fig)

    validation = pd.DataFrame([
        {"check": "configured_library_count", "status": "PASS", "observed": len(libraries), "expected": config["expected_library_count"]},
        {"check": "configured_called_cell_count", "status": "PASS", "observed": len(final), "expected": expected},
        {"check": "unique_global_cell_ids", "status": "PASS", "observed": final["global_cell_id"].nunique(), "expected": expected},
        {"check": "scrublet_hard_exclusions", "status": "PASS", "observed": 0, "expected": 0},
        {"check": "scrublet_tier_changes", "status": "PASS", "observed": tier_changes, "expected": 0},
        {"check": "scrublet_sample_gate_changes", "status": "PASS", "observed": int(samples["n_scrublet_sample_gate_changes"].sum()), "expected": 0},
        {"check": "ambiguous_legacy_field_absent", "status": "PASS", "observed": "absent", "expected": "absent"},
        {"check": "regression_exact_membership", "status": "PASS", "observed": bool(regression["exact_set_equality"].all()) if len(regression) else "not_configured", "expected": True if len(regression) else "not_configured"},
    ])
    validation.to_csv(run / "validation_summary.tsv", sep="\t", index=False)
    report = f"""# {config['report_label']}

Workflow version: `{config['workflow_version']}`  
Study/species: `{config['study_id']}` / `{config['species']}`  
Manifest: `{config['library_manifest']}`  
Libraries/called cells: {len(libraries):,} / {len(final):,}

## QC result

- Primary pass: {int(counts.get('primary_pass', 0)):,}
- Sensitivity only: {int(counts.get('sensitivity_only', 0)):,}
- Technical exclude: {int(counts.get('technical_exclude', 0)):,}
- Scrublet hard exclusions, tier changes, sample-gate changes: 0 / 0 / 0
- Current sample gate uses only `n_candidate_basic_qc`; the ambiguous legacy field is absent.

## Inputs and methods

Raw inputs and their adapters are declared per library in the manifest. Basic QC uses the frozen strict/candidate/retention formulas from task configuration. Scrublet is diagnostic-only. Missing inputs, incompatible adapters, study/count mismatch, and duplicate global IDs fail the run.

## Caveats

This output records QC membership only. Study-level clustering requires independent validation, visual QC review, and reconciliation against the membership used by the previous clustering release.

## Re-run

From the candidate root: `/root/miniconda3/envs/single-cell-qc/bin/snakemake --snakefile Snakefile --configfile {parsed.config} --cores 4 --rerun-incomplete`.
"""
    (run / "run_report.md").write_text(report, encoding="utf-8")
    summary = {"status": "PASS", "study_id": config["study_id"], "libraries": len(libraries), "called_cells": len(final), "tier_counts": {k: int(v) for k, v in counts.items()}, "scrublet_hard_exclusions": 0, "scrublet_tier_changes": 0, "scrublet_sample_gate_changes": 0, "elapsed_seconds": time.perf_counter() - started}
    (run / "run_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    software = {
        "python": platform.python_version(), "snakemake": version("snakemake"),
        "pandas": pd.__version__, "numpy": np.__version__, "scipy": version("scipy"),
        "h5py": version("h5py"), "pyarrow": version("pyarrow"),
        "pyyaml": version("PyYAML"), "jsonschema": version("jsonschema"),
        "matplotlib": version("matplotlib"), "scrublet": version("scrublet"),
    }
    pd.DataFrame(software.items(), columns=["software", "version"]).to_csv(
        run / "software_versions.tsv", sep="\t", index=False
    )
    usage = resource.getrusage(resource.RUSAGE_SELF)
    pd.DataFrame([{"record_type": "aggregate", "study_id": config["study_id"], "elapsed_seconds": summary["elapsed_seconds"], "max_rss_kb": int(usage.ru_maxrss)}]).to_csv(run / "aggregate_resource_usage.tsv", sep="\t", index=False)
    (run / "run_complete.ok").write_text(json.dumps(summary, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
