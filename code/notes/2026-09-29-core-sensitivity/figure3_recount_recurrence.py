"""Rerun Figure 3 P-S recurrence on the frozen endpoint recount.

The primary sensitivity views hold either the gene set or the contributing
metacells/libraries constant across all six counting rules.  The original
Figure 3 covariates and reference-WPP capacity matching are retained.  P-M is
out of scope because the endpoint recount froze P and documented-source S.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
INPUT = HERE / "fixed-catalog-results/route_wpp.tsv.gz"
F3 = Path(
    "/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot/outputs/figure3-baseline-wpp-null-2026-09-23"
)
OLD_SCRIPT = Path(
    "/Users/lycidas/work/科研/projects/2026-09-02_mouse-time-aware-integration-pilot/analysis/figure3-study-recurrence-null-2026-09-20/analyze.py"
)
OUT = HERE / "figure3-recount-recurrence"
MODES = [
    "original",
    "atlas_dedup",
    "transferred_observed",
    "transferred",
    "primary_clipped",
    "all_clipped",
]
VIEWS = ["fixed_genes", "common_contributors"]
THRESHOLDS = [0.02, 0.05]
N_PERM = 10_000
SEED = 20260929
EXCLUDED_STUDIES: set[str] = set()


def load_old():
    spec = importlib.util.spec_from_file_location("original_recurrence", OLD_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    module.N_PERM = N_PERM
    return module


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stratify(frame: pd.DataFrame, old) -> pd.DataFrame:
    result = old.add_strata(frame, include_pau=False)
    result["original_stratum"] = result.stratum
    threshold = float(result.threshold.iloc[0])
    result["can_negative"] = result.reference_WPP_min.ge(threshold - 1e-9)
    result["can_positive"] = result.reference_WPP_max.le(1 - threshold + 1e-9)
    assert result.loc[result.event.eq(-1), "can_negative"].all()
    assert result.loc[result.event.eq(1), "can_positive"].all()
    result["baseline_bin"] = result.groupby("unit").reference_WPP.transform(
        lambda x: old.qbin(x, q=5)
    )
    result["stratum"] = result[
        ["original_stratum", "baseline_bin", "can_negative", "can_positive"]
    ].astype(str).agg("|".join, axis=1)
    return result


def build_events(route: pd.DataFrame, covariates: pd.DataFrame, threshold: float) -> pd.DataFrame:
    rows = []
    for keys, frame in route.groupby(["species", "study_id", "gene_id"]):
        delta = frame.delta.to_numpy(float)
        event = 1 if np.all(delta >= threshold) else -1 if np.all(delta <= -threshold) else 0
        rows.append(
            dict(zip(["species", "unit", "gene_id"], keys))
            | {
                "event": event,
                "APA_delta_mean": float(delta.mean()),
                "n_routes": int(frame.route.nunique()),
                "reference_WPP": float(frame.reference_WPP.median()),
                "reference_WPP_min": float(frame.reference_WPP.min()),
                "reference_WPP_max": float(frame.reference_WPP.max()),
            }
        )
    events = pd.DataFrame(rows)
    events["contrast"] = "P_S"
    events["threshold"] = threshold
    events["level"] = "study"
    events["eligible_units"] = events.groupby(["species", "gene_id"])["unit"].transform("nunique")
    before = len(events)
    events = events.merge(
        covariates,
        on=["species", "contrast", "unit", "gene_id", "threshold"],
        how="left",
        validate="one_to_one",
    )
    assert len(events) == before
    assert events[["gene_symbol", "APA_coverage", "PAS_number"]].notna().all().all()
    return events


def recurrent_members(events: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for keys, frame in events.groupby(["species", "gene_id", "gene_symbol"]):
        positive = int((frame.event == 1).sum())
        negative = int((frame.event == -1).sum())
        support = max(positive, negative)
        if support >= 2:
            rows.append(
                dict(zip(["species", "gene_id", "gene_symbol"], keys))
                | {
                    "direction": "positive" if positive >= negative else "negative",
                    "same_direction_support": support,
                }
            )
    return pd.DataFrame(
        rows,
        columns=[
            "species",
            "gene_id",
            "gene_symbol",
            "direction",
            "same_direction_support",
        ],
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    assert INPUT.is_file(), f"Missing endpoint recount export: {INPUT}"
    old = load_old()
    route = pd.read_csv(INPUT, sep="\t")
    assert set(route["mode"]) == set(MODES)
    assert set(VIEWS) <= set(route.view)
    assert np.allclose(route.target_WPP - route.reference_WPP, route.delta, atol=1e-10)

    original_events = pd.read_csv(F3 / "event_baseline_matrix.tsv.gz", sep="\t")
    original_events = original_events[
        original_events.contrast.eq("P_S")
        & original_events.level.eq("study")
    ].copy()
    covariates = original_events[
        [
            "species",
            "contrast",
            "unit",
            "gene_id",
            "gene_symbol",
            "threshold",
            "mean_GEX",
            "APA_coverage",
            "PAS_number",
        ]
    ].drop_duplicates()

    summaries = []
    diagnostics = []
    member_parts = []
    event_parts = []
    transitions = []
    null_archive = {}
    original_overlap_checks = []

    # Recover the frozen original estimator before applying either sensitivity
    # view.  The fixed-gene and common-contributor views intentionally restrict
    # route contributors and therefore need not reproduce a multi-route study's
    # frozen study-level effect even under the original count mode.
    reapplied_route = route[
        route.view.eq("reapplied") & route["mode"].eq("original")
    ]
    for threshold in THRESHOLDS:
        recovered = build_events(reapplied_route, covariates, threshold)
        frozen = original_events[original_events.threshold.eq(threshold)]
        overlap = recovered.merge(
            frozen[
                [
                    "species",
                    "unit",
                    "gene_id",
                    "event",
                    "APA_delta_mean",
                    "reference_WPP",
                    "reference_WPP_min",
                    "reference_WPP_max",
                ]
            ],
            on=["species", "unit", "gene_id"],
            how="left",
            suffixes=("_recount", "_frozen"),
            validate="one_to_one",
        )
        assert len(overlap) == len(recovered) == len(frozen)
        assert overlap.event_recount.eq(overlap.event_frozen).all()
        assert np.allclose(
            overlap.APA_delta_mean_recount,
            overlap.APA_delta_mean_frozen,
            atol=1e-8,
        )
        assert np.allclose(
            overlap.reference_WPP_recount,
            overlap.reference_WPP_frozen,
            atol=1e-8,
        )
        for unit, unit_overlap in overlap.groupby("unit"):
            original_overlap_checks.append(
                {
                    "view": "reapplied",
                    "threshold": threshold,
                    "unit": unit,
                    "incomplete_route_excluded": False,
                    "rows": len(unit_overlap),
                    "event_mismatches": 0,
                    "effect_max_error": float(
                        (
                            unit_overlap.APA_delta_mean_recount
                            - unit_overlap.APA_delta_mean_frozen
                        )
                        .abs()
                        .max()
                    ),
                    "reference_WPP_max_error": float(
                        (
                            unit_overlap.reference_WPP_recount
                            - unit_overlap.reference_WPP_frozen
                        )
                        .abs()
                        .max()
                    ),
                }
            )

    for view in VIEWS:
        view_route = route[route.view.eq(view)]
        for threshold in THRESHOLDS:
            event_by_mode = {}
            member_by_mode = {}
            for mode in MODES:
                events = build_events(
                    view_route[view_route["mode"].eq(mode)], covariates, threshold
                )
                event_by_mode[mode] = events
                members = recurrent_members(events)
                member_by_mode[mode] = members
                event_parts.append(events.assign(view=view, mode=mode))
                member_parts.append(
                    members.assign(view=view, mode=mode, threshold=threshold)
                )
            baseline = event_by_mode["original"]
            frozen = original_events[original_events.threshold.eq(threshold)]
            overlap = baseline.merge(
                frozen[
                    [
                        "species",
                        "unit",
                        "gene_id",
                        "event",
                        "APA_delta_mean",
                        "reference_WPP",
                        "reference_WPP_min",
                        "reference_WPP_max",
                    ]
                ],
                on=["species", "unit", "gene_id"],
                how="left",
                suffixes=("_recount", "_frozen"),
                validate="one_to_one",
            )
            assert len(overlap) == len(baseline)
            for unit, unit_overlap in overlap.groupby("unit"):
                original_overlap_checks.append(
                    {
                        "view": view,
                        "threshold": threshold,
                        "unit": unit,
                        "incomplete_route_excluded": unit in EXCLUDED_STUDIES,
                        "rows": len(unit_overlap),
                        "event_mismatches": int(
                            unit_overlap.event_recount.ne(
                                unit_overlap.event_frozen
                            ).sum()
                        ),
                        "effect_max_error": float(
                            (
                                unit_overlap.APA_delta_mean_recount
                                - unit_overlap.APA_delta_mean_frozen
                            )
                            .abs()
                            .max()
                        ),
                        "reference_WPP_max_error": float(
                            (
                                unit_overlap.reference_WPP_recount
                                - unit_overlap.reference_WPP_frozen
                            )
                            .abs()
                            .max()
                        ),
                    }
                )
            baseline_members = member_by_mode["original"].set_index(["species", "gene_id"])
            frozen_members = recurrent_members(frozen).set_index(["species", "gene_id"])
            for mode in MODES:
                current = event_by_mode[mode]
                joined = baseline.merge(
                    current,
                    on=["species", "contrast", "unit", "gene_id", "threshold"],
                    suffixes=("_original", "_current"),
                    validate="one_to_one",
                )
                original_nonzero = joined.event_original.ne(0)
                transitions.append(
                    {
                        "view": view,
                        "mode": mode,
                        "threshold": threshold,
                        "study_gene_rows": len(joined),
                        "original_event_rows": int(original_nonzero.sum()),
                        "same_event_rows": int(
                            (original_nonzero & joined.event_current.eq(joined.event_original)).sum()
                        ),
                        "reversed_event_rows": int(
                            (original_nonzero & joined.event_current.eq(-joined.event_original)).sum()
                        ),
                        "became_subthreshold_rows": int(
                            (original_nonzero & joined.event_current.eq(0)).sum()
                        ),
                    }
                )
                current_members = member_by_mode[mode].set_index(["species", "gene_id"])
                for species in sorted(current.species.unique()):
                    original_set = set(baseline_members.loc[species].index) if species in baseline_members.index.get_level_values(0) else set()
                    current_set = set(current_members.loc[species].index) if species in current_members.index.get_level_values(0) else set()
                    frozen_set = set(frozen_members.loc[species].index) if species in frozen_members.index.get_level_values(0) else set()
                    summaries.append(
                        {
                            "view": view,
                            "mode": mode,
                            "species": species,
                            "threshold": threshold,
                            "observed_recurrent": len(current_set),
                            "published_full_recurrent": len(frozen_set),
                            "original_recurrent_same_view": len(original_set),
                            "original_view_retained_from_published": len(
                                frozen_set & original_set
                            ),
                            "retained_recurrent": len(original_set & current_set),
                            "lost_recurrent": len(original_set - current_set),
                            "gained_recurrent": len(current_set - original_set),
                        }
                    )

            for mode in MODES:
                events = event_by_mode[mode]
                for species, frame in events.groupby("species"):
                    stratified = stratify(frame, old)
                    key = f"{view}__{mode}__{species}__{threshold:.2f}"
                    seed = SEED + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
                    stats, _, diag, null = old.permute_analysis(
                        stratified, np.random.default_rng(seed)
                    )
                    row = stats[stats.metric.eq("N_ge2")].iloc[0].to_dict()
                    row.update(
                        {
                            "view": view,
                            "mode": mode,
                            "species": species,
                            "contrast": "P_S",
                            "threshold": threshold,
                            "seed": seed,
                            "permutations": N_PERM,
                        }
                    )
                    diagnostics.append(row)
                    null_archive[key] = null["N_ge2"]

    null_results = pd.DataFrame(diagnostics)
    null_results["q_within_view_mode_four_tests"] = np.nan
    for _, indexes in null_results.groupby(["view", "mode"]).groups.items():
        null_results.loc[indexes, "q_within_view_mode_four_tests"] = old.bh(
            null_results.loc[indexes, "empirical_p"]
        )
    summary = pd.DataFrame(summaries).merge(
        null_results[
            [
                "view",
                "mode",
                "species",
                "threshold",
                "null_mean",
                "null_median",
                "null_p025",
                "null_p975",
                "empirical_p",
                "q_within_view_mode_four_tests",
            ]
        ],
        on=["view", "mode", "species", "threshold"],
        validate="one_to_one",
    )
    summary.to_csv(OUT / "recurrence_robustness_summary.tsv", sep="\t", index=False)
    null_results.to_csv(OUT / "recurrence_null_results.tsv", sep="\t", index=False)
    pd.concat(event_parts, ignore_index=True).to_csv(
        OUT / "study_events.tsv.gz", sep="\t", index=False
    )
    pd.concat(member_parts, ignore_index=True).to_csv(
        OUT / "recurrent_members.tsv.gz", sep="\t", index=False
    )
    pd.DataFrame(transitions).to_csv(
        OUT / "study_event_transitions.tsv", sep="\t", index=False
    )
    pd.DataFrame(original_overlap_checks).to_csv(
        OUT / "original_overlap_checks.tsv", sep="\t", index=False
    )
    np.savez_compressed(OUT / "null_draws.npz", **null_archive)
    audit = {
        "status": "PROVISIONAL_EXCLUDED_STUDY_ROUTE" if EXCLUDED_STUDIES else "PASS",
        "scope": "P-S only; endpoint recount does not contain the complete P-M recurrence cohort",
        "excluded_studies_with_incomplete_routes": sorted(EXCLUDED_STUDIES),
        "views": VIEWS,
        "modes": MODES,
        "thresholds": THRESHOLDS,
        "permutations_per_test": N_PERM,
        "tests": len(null_results),
        "original_reapplied_event_recovery": "exact",
        "input_sha256": {
            str(INPUT): digest(INPUT),
            str(F3 / "event_baseline_matrix.tsv.gz"): digest(F3 / "event_baseline_matrix.tsv.gz"),
            str(OLD_SCRIPT): digest(OLD_SCRIPT),
        },
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--permutations", type=int, default=N_PERM)
    parser.add_argument("--excluded-study", action="append", default=[])
    args = parser.parse_args()
    INPUT = args.input
    OUT = args.out
    N_PERM = args.permutations
    EXCLUDED_STUDIES = set(args.excluded_study)
    main()
