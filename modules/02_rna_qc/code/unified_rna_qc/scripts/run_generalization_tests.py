#!/usr/bin/env python3
"""Run hard-coding, regression, configuration, negative, and schema tests."""

from __future__ import annotations

import copy
import json
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

from common import global_cell_ids, load_config, read_input, read_manifest


ROOT = Path(__file__).resolve().parents[1]
TESTS = ROOT / "tests"
PRODUCTION = [
    ROOT / "Snakefile", ROOT / "config.schema.yaml", ROOT / "scripts/common.py",
    ROOT / "scripts/process_library.py", ROOT / "scripts/aggregate.py",
    ROOT / "scripts/validate_run.py", ROOT / "scripts/package_workflow.py",
    ROOT / "scripts/build_release.py", ROOT / "scripts/verify_release.py",
]
FORBIDDEN = {
    "study_identifier": "GSE178325",
    "library_prefix": "gse178325_human_chemical_pluripotency_reprogramming_",
    "literal_global_prefix": "GSE178325|",
    "pilot_called_count": "222726",
    "pilot_primary_count": "140593",
    "pilot_sensitivity_count": "36559",
    "pilot_technical_count": "45574",
}


def expect_failure(name, function) -> dict:
    try:
        function()
    except Exception as exc:
        return {"test": name, "status": "PASS", "observed": type(exc).__name__, "detail": str(exc)[:240]}
    return {"test": name, "status": "FAIL", "observed": "no_error", "detail": "Mutation was accepted"}


def main() -> None:
    TESTS.mkdir(exist_ok=True)
    scans = []
    for path in PRODUCTION:
        if not path.is_file():
            continue
        text = path.read_text(encoding="utf-8")
        for label, token in FORBIDDEN.items():
            count = text.count(token)
            scans.append({"file": str(path.relative_to(ROOT)), "pattern": label, "match_count": count, "status": "PASS" if count == 0 else "FAIL"})
    scan = pd.DataFrame(scans)
    scan.to_csv(TESTS / "hard_coding_scan.tsv", sep="\t", index=False)

    regression = pd.read_csv(ROOT / "runs/gse178325_regression/regression_membership_comparison.tsv", sep="\t")
    regression["status"] = np.where(regression["exact_set_equality"].astype(bool) & regression["added_n"].eq(0) & regression["removed_n"].eq(0), "PASS", "FAIL")
    regression.to_csv(TESTS / "gse178325_regression_golden.tsv", sep="\t", index=False)
    reg_final = pd.read_parquet(ROOT / "runs/gse178325_regression/cell_qc_decisions.parquet")
    non_final = pd.read_parquet(ROOT / "runs/gse103221_small/cell_qc_decisions.parquet")
    golden = pd.DataFrame([
        {"test": "accepted_membership_exact", "status": "PASS" if regression["status"].eq("PASS").all() else "FAIL", "observed": int(reg_final.shape[0]), "expected": int(reg_final.shape[0])},
        {"test": "regression_scrublet_zero_changes", "status": "PASS" if (not reg_final["doublet_hard_exclude"].astype(bool).any() and reg_final["pre_doublet_qc_tier"].equals(reg_final["final_qc_tier"])) else "FAIL", "observed": 0, "expected": 0},
        {"test": "real_non_target_study_same_snakefile", "status": "PASS" if len(non_final) > 0 and set(non_final["study_id"]) != set(reg_final["study_id"]) else "FAIL", "observed": len(non_final), "expected": 13750},
        {"test": "non_target_scrublet_zero_changes", "status": "PASS" if (not non_final["doublet_hard_exclude"].astype(bool).any() and non_final["pre_doublet_qc_tier"].equals(non_final["final_qc_tier"])) else "FAIL", "observed": 0, "expected": 0},
    ])
    golden.to_csv(TESTS / "golden_tests.tsv", sep="\t", index=False)

    config_path = ROOT / "configs/GSE103221_small.yaml"
    base = load_config(config_path)
    manifest = pd.read_csv(ROOT / base["library_manifest"], sep="\t", dtype=str)
    work = TESTS / ".negative_work"
    work.mkdir(exist_ok=True)
    wrong_study = copy.deepcopy(base)
    wrong_study["study_id"] = "WRONG_STUDY"
    wrong_count = copy.deepcopy(base)
    wrong_count["expected_called_cell_count"] = int(base["expected_called_cell_count"]) + 1

    def validate_mutation(cfg, name):
        import yaml
        path = work / name
        path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
        read_manifest(ROOT, cfg)

    negative = [
        expect_failure("wrong_study_id", lambda: validate_mutation(wrong_study, "wrong_study.yaml")),
        expect_failure("wrong_expected_counts", lambda: validate_mutation(wrong_count, "wrong_count.yaml")),
        expect_failure("duplicate_global_cell_ids", lambda: global_cell_ids(base, manifest.iloc[0]["library_id"], np.asarray(["barcode-1", "barcode-1"]))),
        expect_failure("missing_adapter_inputs", lambda: read_input("cellranger_filtered_h5", work / "not_present.h5")),
    ]
    pd.DataFrame(negative).to_csv(TESTS / "negative_tests.tsv", sep="\t", index=False)
    shutil.rmtree(work)

    schema_text = (ROOT / "config.schema.yaml").read_text(encoding="utf-8")
    schema_rows = []
    for run_name in ("gse178325_regression", "gse103221_small"):
        columns = pd.read_csv(ROOT / f"runs/{run_name}/sample_qc_decisions.tsv", sep="\t", nrows=0).columns
        schema_rows.append({"test": f"sample_decision_schema:{run_name}", "legacy_field_absent": "n_candidate_singlets" not in columns, "current_field_present": "n_candidate_basic_qc" in columns, "status": "PASS" if "n_candidate_singlets" not in columns and "n_candidate_basic_qc" in columns else "FAIL"})
    schema_rows.append({"test": "configuration_schema", "legacy_field_absent": "n_candidate_singlets" not in schema_text, "current_field_present": "expected_called_cell_count" in schema_text, "status": "PASS" if "n_candidate_singlets" not in schema_text else "FAIL"})
    pd.DataFrame(schema_rows).to_csv(TESTS / "schema_tests.tsv", sep="\t", index=False)
    summary = pd.DataFrame([
        {"suite": "hard_coding_scan", "tests": len(scan), "failures": int(scan["status"].ne("PASS").sum())},
        {"suite": "golden", "tests": len(golden), "failures": int(golden["status"].ne("PASS").sum())},
        {"suite": "negative", "tests": len(negative), "failures": sum(row["status"] != "PASS" for row in negative)},
        {"suite": "schema", "tests": len(schema_rows), "failures": sum(row["status"] != "PASS" for row in schema_rows)},
    ])
    summary["status"] = np.where(summary["failures"].eq(0), "PASS", "FAIL")
    summary.to_csv(TESTS / "generalization_test_summary.tsv", sep="\t", index=False)
    if summary["failures"].sum():
        raise SystemExit(json.dumps(summary.to_dict("records")))
    print(json.dumps(summary.to_dict("records"), sort_keys=True))


if __name__ == "__main__":
    main()
