"""Re-estimate Figure 5 P-M effects within common physical libraries.

This is a frozen-estimand sensitivity analysis.  It keeps the annotated PAS
models, candidate sites, state assignments, UMI thresholds, coverage threshold,
and downstream recurrence tests unchanged.  Only P-M state averaging changes:
P and M must pass the same gate in the same physical library before their
within-library usage difference contributes to the equal-library mean.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "analysis_code/release_2026-09-24/sources/F5/annotation-pas-global-v1"
RECURRENCE_SOURCE = ROOT / "analysis_code/release_2026-09-24/sources/F5/site-recurrence-v1"
OUT = Path(__file__).resolve().parent / "figure5-paired-pm"
ROUTES = {"mouse": ["CIP", "WOT"], "human": ["GSE118258", "GSE150311"]}


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def qualified_libraries(model: dict, state: str, cutoff: int) -> dict[str, dict]:
    return {
        row["library"]: row
        for row in model["profiles"].get(state, {}).get("libraries", [])
        if row["total"] >= cutoff and row["coverage"] >= 0.5
    }


def paired_state_stats(model: dict, cutoff: int) -> dict:
    """Return equal-library paired P-M usage differences and support counts."""
    m = qualified_libraries(model, "M", cutoff)
    p = qualified_libraries(model, "P", cutoff)
    common = sorted(set(m) & set(p))
    n_pas = len(model["pas"])
    if not common:
        missing = np.full(n_pas, np.nan)
        zeros = np.zeros(n_pas, dtype=float)
        return {
            "libraries": common,
            "delta_usage": missing,
            "m_counts": zeros,
            "p_counts": zeros,
            "m_qualified": len(m),
            "p_qualified": len(p),
        }
    differences = []
    m_counts = np.zeros(n_pas, dtype=float)
    p_counts = np.zeros(n_pas, dtype=float)
    for library in common:
        mu = np.asarray(m[library]["counts"], dtype=float)
        pu = np.asarray(p[library]["counts"], dtype=float)
        assert len(mu) == n_pas == len(pu)
        assert np.isclose(mu.sum(), m[library]["total"])
        assert np.isclose(pu.sum(), p[library]["total"])
        differences.append(pu / pu.sum() - mu / mu.sum())
        m_counts += mu
        p_counts += pu
    return {
        "libraries": common,
        "delta_usage": np.mean(differences, axis=0),
        "m_counts": m_counts,
        "p_counts": p_counts,
        "m_qualified": len(m),
        "p_qualified": len(p),
    }


def observed_both5(stats: dict, mask: np.ndarray) -> bool:
    if not stats["libraries"]:
        return False
    support = []
    for counts in [stats["m_counts"], stats["p_counts"]]:
        support.extend([float(counts @ mask), float(counts @ (1 - mask))])
    return min(support) >= 5


def main() -> None:
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "full-cohort-v2/data/mouse").mkdir(parents=True)
    (OUT / "full-cohort-v2/data/human").mkdir(parents=True)
    (OUT / "analysis").mkdir(parents=True)

    input_hashes: dict[str, str] = {}
    model_lookup: dict[tuple[str, str, str], dict] = {}
    paired_lookup: dict[tuple[str, str, str, int], dict] = {}
    coverage_rows = []
    for species, routes in ROUTES.items():
        for route in routes:
            path = SOURCE / f"full-cohort-v2/data/{species}/{route}/models.json"
            input_hashes[str(path.relative_to(ROOT))] = digest(path)
            models = json.loads(path.read_text())
            for gene_id, model in models.items():
                model_lookup[species, route, gene_id] = model
                if "MP" not in model.get("eligible_contrasts", []):
                    continue
                for cutoff in [20, 50]:
                    stats = paired_state_stats(model, cutoff)
                    paired_lookup[species, route, gene_id, cutoff] = stats
                    coverage_rows.append(
                        {
                            "species": species,
                            "route": route,
                            "gene_id": gene_id,
                            "gene_symbol": model["gene_symbol"],
                            "cutoff": cutoff,
                            "qualified_M_libraries": stats["m_qualified"],
                            "qualified_P_libraries": stats["p_qualified"],
                            "common_libraries": len(stats["libraries"]),
                            "common_library_ids": ";".join(stats["libraries"]),
                        }
                    )

    site_change_rows = []
    for species in ROUTES:
        source_path = SOURCE / f"full-cohort-v2/data/{species}/site_effects.tsv.gz"
        input_hashes[str(source_path.relative_to(ROOT))] = digest(source_path)
        sites = pd.read_csv(source_path, sep="\t", low_memory=False)
        original = sites.copy()
        mp_index = sites.index[sites.contrast.eq("MP")]
        for index in mp_index:
            row = sites.loc[index]
            mask_matrix = np.asarray(json.loads(row.mask_json), dtype=int)
            assert mask_matrix.ndim == 2 and mask_matrix.shape[1] == 2
            assert np.array_equal(mask_matrix[:, 0], mask_matrix[:, 1])
            mask = mask_matrix[:, 0]
            for cutoff in [20, 50]:
                stats = paired_lookup.get((species, row.route, row.gene_id, cutoff))
                effect = float(stats["delta_usage"] @ mask) if stats is not None else np.nan
                n_common = len(stats["libraries"]) if stats is not None else 0
                sites.at[index, f"eq{cutoff}_lo"] = effect
                sites.at[index, f"eq{cutoff}_hi"] = effect
                sites.at[index, f"eq{cutoff}_n_reference"] = n_common
                sites.at[index, f"eq{cutoff}_n_target"] = n_common
            site_change_rows.append(
                {
                    "species": species,
                    "route": row.route,
                    "gene_id": row.gene_id,
                    "gene_symbol": row.gene_symbol,
                    "category": row.category,
                    "start0": row.start0,
                    "end0": row.end0,
                    "marginal_effect20": original.at[index, "eq20_lo"],
                    "paired_effect20": sites.at[index, "eq20_lo"],
                    "marginal_effect50": original.at[index, "eq50_lo"],
                    "paired_effect50": sites.at[index, "eq50_lo"],
                }
            )
        # P-S is an invariant of this sensitivity analysis.
        unchanged = sites.contrast.ne("MP")
        effect_columns = [
            "eq20_lo", "eq20_hi", "eq50_lo", "eq50_hi",
            "eq20_n_reference", "eq20_n_target",
            "eq50_n_reference", "eq50_n_target",
        ]
        pd.testing.assert_frame_equal(
            sites.loc[unchanged, effect_columns],
            original.loc[unchanged, effect_columns],
            check_dtype=False,
        )
        sites.to_csv(
            OUT / f"full-cohort-v2/data/{species}/site_effects.tsv.gz",
            sep="\t",
            index=False,
        )

    region_path = SOURCE / "cis-regulatory-assessment-v1/analysis/all_coordinate_regions.tsv.gz"
    input_hashes[str(region_path.relative_to(ROOT))] = digest(region_path)
    regions = pd.read_csv(region_path, sep="\t")
    original_regions = regions.copy()
    for index in regions.index[regions.contrast.eq("MP")]:
        row = regions.loc[index]
        mask = np.asarray(json.loads(row["mask"]), dtype=int)
        for cutoff in [20, 50]:
            stats = paired_lookup.get((row.species, row.route, row.gene_id, cutoff))
            regions.at[index, f"effect{cutoff}"] = (
                float(stats["delta_usage"] @ mask) if stats is not None else np.nan
            )
        stats50 = paired_lookup.get((row.species, row.route, row.gene_id, 50))
        regions.at[index, "observed_both5"] = (
            observed_both5(stats50, mask) if stats50 is not None else False
        )
    pd.testing.assert_frame_equal(
        regions.loc[regions.contrast.ne("MP")],
        original_regions.loc[regions.contrast.ne("MP")],
        check_dtype=False,
    )
    paired_region_path = OUT / "all_coordinate_regions.paired_pm.tsv.gz"
    regions.to_csv(paired_region_path, sep="\t", index=False)
    original_support_regions = regions.copy()
    original_support_regions["observed_both5"] = original_regions["observed_both5"]
    original_support_region_path = OUT / "all_coordinate_regions.paired_pm_original_support.tsv.gz"
    original_support_regions.to_csv(original_support_region_path, sep="\t", index=False)
    pd.DataFrame(coverage_rows).to_csv(OUT / "paired_library_coverage.tsv", sep="\t", index=False)
    pd.DataFrame(site_change_rows).to_csv(OUT / "paired_vs_marginal_site_effects.tsv.gz", sep="\t", index=False)

    recurrence = load_module("paired_recurrence", RECURRENCE_SOURCE / "analyze.py")
    recurrence.P = OUT
    recurrence.OUT = OUT / "analysis"
    recurrence.REG = paired_region_path
    recurrence.main()
    consistency = load_module("paired_consistency", RECURRENCE_SOURCE / "consistency.py")
    consistency.ROOT = OUT
    consistency.OUT = OUT / "analysis"
    consistency.SOURCE = OUT / "analysis/exact_site_cross_study.tsv.gz"
    consistency.main()

    # Isolate whether attrition is caused by the paired estimator itself or by
    # also requiring the included/excluded support gate in the same libraries.
    recurrence_support = load_module("paired_recurrence_original_support", RECURRENCE_SOURCE / "analyze.py")
    recurrence_support.P = OUT
    recurrence_support.OUT = OUT / "analysis-original-support"
    recurrence_support.REG = original_support_region_path
    recurrence_support.main()
    consistency_support = load_module("paired_consistency_original_support", RECURRENCE_SOURCE / "consistency.py")
    consistency_support.ROOT = OUT
    consistency_support.OUT = OUT / "analysis-original-support"
    consistency_support.SOURCE = OUT / "analysis-original-support/exact_site_cross_study.tsv.gz"
    consistency_support.main()

    current = pd.read_csv(
        SOURCE / "site-recurrence-v1/analysis/region_consistency_statistics.tsv", sep="\t"
    )
    paired = pd.read_csv(OUT / "analysis/region_consistency_statistics.tsv", sep="\t")
    paired_original_support = pd.read_csv(
        OUT / "analysis-original-support/region_consistency_statistics.tsv", sep="\t"
    )
    primary = lambda d: d[d.primary_test].copy()
    comparison = primary(current).merge(
        primary(paired),
        on=["species", "contrast", "category", "unit", "primary_test"],
        suffixes=("_marginal", "_paired"),
        validate="one_to_one",
    )
    comparison.to_csv(OUT / "primary_region_statistics_comparison.tsv", sep="\t", index=False)
    three_way = primary(current).merge(
        primary(paired_original_support),
        on=["species", "contrast", "category", "unit", "primary_test"],
        suffixes=("_marginal", "_paired_original_support"),
        validate="one_to_one",
    ).merge(
        primary(paired),
        on=["species", "contrast", "category", "unit", "primary_test"],
        validate="one_to_one",
    )
    three_way.to_csv(OUT / "primary_region_statistics_three_way.tsv", sep="\t", index=False)

    coverage = pd.DataFrame(coverage_rows)
    audit = {
        "status": "PASS",
        "estimand": "mean within-physical-library P-minus-M inclusion difference",
        "fixed_elements": [
            "annotated PAS models",
            "candidate recognition intervals",
            "state assignments",
            "20 and 50 family-UMI thresholds",
            "family coverage >= 0.5",
            "strict 10 percentage-point change rule",
            "downstream recurrence and bootstrap definitions",
        ],
        "mp_gene_route_cutoff_rows": len(coverage),
        "rows_without_common_library": int(coverage.common_libraries.eq(0).sum()),
        "site_effect_rows": len(site_change_rows),
        "primary_result_rows": len(comparison),
        "input_sha256": input_hashes,
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(comparison.to_string(index=False))
    print(json.dumps({k: v for k, v in audit.items() if k != "input_sha256"}, indent=2))


if __name__ == "__main__":
    main()
