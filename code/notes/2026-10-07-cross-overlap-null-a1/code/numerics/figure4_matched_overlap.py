"""Strengthen the original Figure 4A overlap null with WPP/capacity matching.

The original 18 P-S and 10 paired P-M ortholog sets are kept frozen.  This
script only changes the ortholog-label null by adding contrast-specific
measurement count, reference-WPP, and proximal-direction capacity to the
original measurement/PAS strata.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[2]
PROJECT = Path("/Users/lycidas/work/科研/projects/2026-09-13_reproapa-figure4-real")
PROGRAM = PROJECT / "data/figure4-program-v9"
ORIGINAL = PROJECT / "data/figure4-ortholog-null-v24"
OUT = Path(__file__).resolve().parent / "figure4-matched-overlap"
TAU = 0.02
B = 20_000
SEED = 20260929
CODES = ["W", "C", "H1", "H2", "H3", "H4"]
MOUSE = ["W", "C"]
HUMAN = ["H1", "H2", "H3", "H4"]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quantile_codes(values: pd.Series, bins: int = 5) -> pd.Series:
    result = pd.Series(-1, index=values.index, dtype=int)
    finite = values.notna()
    if finite.any():
        result.loc[finite] = pd.qcut(
            values.loc[finite], bins, labels=False, duplicates="drop"
        ).astype(int)
    return result


def range10(values: pd.Series) -> pd.Series:
    result = pd.Series(-1, index=values.index, dtype=int)
    finite = values.notna()
    result.loc[finite] = np.minimum(np.floor(values.loc[finite] * 10), 9).astype(int)
    return result


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    adjusted = np.empty(len(p), dtype=float)
    adjusted[order] = np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1].clip(0, 1)
    return adjusted


def analytic_stratum(n: int, nx: int, ny: int) -> tuple[float, float]:
    expected = nx * ny / n
    if n <= 1:
        return expected, 0.0
    variance = nx * ny * (n - nx) * (n - ny) / (n * n * (n - 1))
    return expected, variance


def matched_null(frame: pd.DataFrame, contrast: str, model: str, columns: list[str], draws=None) -> tuple[dict, pd.DataFrame]:
    x = frame[f"h_{contrast}"].to_numpy(bool)
    y = frame[f"m_{contrast}"].to_numpy(bool)
    labels = frame[columns].astype(str).agg("|".join, axis=1)
    groups = [np.asarray(index, dtype=int) for index in frame.groupby(labels, sort=True).indices.values()]
    rng = np.random.default_rng(SEED + sum(map(ord, contrast + model)))
    null = np.zeros(B, dtype=int)
    diagnostics = []
    expected = 0.0
    variance = 0.0
    fixed_overlap = 0
    for group in groups:
        n = len(group)
        nx = int(x[group].sum())
        ny = int(y[group].sum())
        mean, var = analytic_stratum(n, nx, ny)
        expected += mean
        variance += var
        draw = hypergeom.rvs(n, ny, nx, size=B, random_state=rng)
        null += draw
        is_fixed = var == 0
        if is_fixed:
            fixed_overlap += int(draw[0])
        diagnostics.append(
            {
                "contrast": contrast,
                "model": model,
                "stratum": labels.iloc[group[0]],
                "n": n,
                "human_events": nx,
                "mouse_events": ny,
                "analytic_expected": mean,
                "analytic_variance": var,
                "fixed_contribution": is_fixed,
            }
        )
    observed = int((x & y).sum())
    row = {
        "contrast": contrast,
        "model": model,
        "observed": observed,
        "analytic_null_mean": expected,
        "analytic_null_sd": float(np.sqrt(variance)),
        "null_mean": float(null.mean()),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "empirical_p": float((1 + np.count_nonzero(null >= observed)) / (B + 1)),
        "orthologs": len(frame),
        "strata": len(groups),
        "singleton_strata": sum(len(group) == 1 for group in groups),
        "fixed_overlap": fixed_overlap,
        "movable_null_sd": float(null.std()),
        "permutations": B,
        "seed": SEED + sum(map(ord, contrast + model)),
    }
    if draws is not None:
        draws[contrast + "|" + model] = null
    return row, pd.DataFrame(diagnostics)


def independent_null(
    frame: pd.DataFrame,
    contrast: str,
    model: str,
    human_columns: list[str],
    mouse_columns: list[str],
    draws=None,
) -> tuple[dict, pd.DataFrame]:
    """Shuffle each species' event labels within its own covariate strata."""
    human = frame[f"h_{contrast}"].to_numpy(bool)
    mouse = frame[f"m_{contrast}"].to_numpy(bool)
    human_labels = frame[human_columns].astype(str).agg("|".join, axis=1)
    mouse_labels = frame[mouse_columns].astype(str).agg("|".join, axis=1)
    human_groups = [
        np.asarray(index, dtype=int)
        for index in frame.groupby(human_labels, sort=True).indices.values()
    ]
    mouse_groups = [
        np.asarray(index, dtype=int)
        for index in frame.groupby(mouse_labels, sort=True).indices.values()
    ]
    human_probability = np.zeros(len(frame), dtype=float)
    mouse_probability = np.zeros(len(frame), dtype=float)
    diagnostics = []
    for side, labels, groups, events, probabilities in [
        ("human", human_labels, human_groups, human, human_probability),
        ("mouse", mouse_labels, mouse_groups, mouse, mouse_probability),
    ]:
        for group in groups:
            event_count = int(events[group].sum())
            probabilities[group] = event_count / len(group)
            diagnostics.append(
                {
                    "contrast": contrast,
                    "model": model,
                    "side": side,
                    "stratum": labels.iloc[group[0]],
                    "n": len(group),
                    "events": event_count,
                    "movable": 0 < event_count < len(group),
                }
            )
    analytic = float(np.sum(human_probability * mouse_probability))
    rng = np.random.default_rng(SEED + sum(map(ord, contrast + model)))
    null = np.empty(B, dtype=int)
    hp = np.empty(len(frame), dtype=bool)
    mp = np.empty(len(frame), dtype=bool)
    for draw in range(B):
        for group in human_groups:
            hp[group] = rng.permutation(human[group])
        for group in mouse_groups:
            mp[group] = rng.permutation(mouse[group])
        null[draw] = int(np.count_nonzero(hp & mp))
    observed = int(np.count_nonzero(human & mouse))
    diagnostic_frame = pd.DataFrame(diagnostics)
    row = {
        "contrast": contrast,
        "model": model,
        "observed": observed,
        "analytic_null_mean": analytic,
        "analytic_null_sd": np.nan,
        "null_mean": float(null.mean()),
        "null_median": float(np.median(null)),
        "null_q025": float(np.quantile(null, 0.025)),
        "null_q975": float(np.quantile(null, 0.975)),
        "empirical_p": float((1 + np.count_nonzero(null >= observed)) / (B + 1)),
        "orthologs": len(frame),
        "strata": len(human_groups) + len(mouse_groups),
        "singleton_strata": int(
            sum(len(group) == 1 for group in human_groups + mouse_groups)
        ),
        "fixed_overlap": np.nan,
        "movable_null_sd": float(null.std()),
        "human_movable_fraction": float(
            diagnostic_frame.loc[
                diagnostic_frame.side.eq("human") & diagnostic_frame.movable, "n"
            ].sum()
            / len(frame)
        ),
        "mouse_movable_fraction": float(
            diagnostic_frame.loc[
                diagnostic_frame.side.eq("mouse") & diagnostic_frame.movable, "n"
            ].sum()
            / len(frame)
        ),
        "permutations": B,
        "seed": SEED + sum(map(ord, contrast + model)),
    }
    if draws is not None:
        draws[contrast + "|" + model] = null
    return row, diagnostic_frame


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nodes_path = PROGRAM / "all_mapped_nodes.tsv.gz"
    ortholog_path = PROGRAM / "ortholog_universe.tsv"
    pm_path = PROGRAM / "paired_pm_study_effects.tsv.gz"
    original_feature_path = ORIGINAL / "ortholog_pair_features.tsv"
    original_result_path = ORIGINAL / "overlap_permutation_results.tsv"
    for path in [nodes_path, ortholog_path, pm_path, original_feature_path, original_result_path]:
        assert path.is_file(), path

    nodes = pd.read_csv(nodes_path, sep="\t")
    nodes = nodes[
        nodes.layer.eq("supported50")
        & nodes.measurement_pass
        & nodes.code.isin(CODES)
    ].copy()
    orthologs = pd.read_csv(ortholog_path, sep="\t")
    original_features = pd.read_csv(original_feature_path, sep="\t")
    original_results = pd.read_csv(original_result_path, sep="\t")

    def node_matrix(reference: str, target: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        pivot = nodes.pivot_table(
            index="gene_id", columns=["code", "node"], values="position_index", aggfunc="first"
        )
        effect = pd.DataFrame(index=pivot.index)
        ref = pd.DataFrame(index=pivot.index)
        for code in CODES:
            effect[code] = (
                pivot[(code, target)] - pivot[(code, reference)]
                if (code, target) in pivot and (code, reference) in pivot
                else np.nan
            )
            ref[code] = pivot[(code, reference)] if (code, reference) in pivot else np.nan
        return effect, ref

    ps_effect, ps_ref = node_matrix("S", "P")
    _, pm_ref = node_matrix("M", "P")
    pm = pd.read_csv(pm_path, sep="\t")
    pm = pm[pm.layer.eq("supported50") & pm.quality_pass].copy()
    pm = pm[
        ((pm.n_sample_units >= 2) & pm.code.isin(MOUSE))
        | ((pm.n_sample_units >= 1) & pm.code.isin(HUMAN))
    ]
    pm_effect = pm.pivot_table(index="gene_id", columns="code", values="mean_delta", aggfunc="first")
    pm_effect = pm_effect.reindex(columns=CODES)
    effects = {"PS": ps_effect, "PM": pm_effect}
    references = {"PS": ps_ref, "PM": pm_ref}

    frame = original_features.copy()
    assert len(frame) == 16_092
    assert int((frame.h_PS & frame.m_PS).sum()) == 18
    assert int((frame.h_PM & frame.m_PM).sum()) == 10

    for contrast in ["PS", "PM"]:
        effect = effects[contrast]
        reference = references[contrast]
        for prefix, ids, codes in [
            ("h", frame.human_gene_id, HUMAN),
            ("m", frame.mouse_gene_id, MOUSE),
        ]:
            values = []
            for gene_id in ids:
                if gene_id not in effect.index:
                    values.append((0, np.nan, np.nan, np.nan, 0))
                    continue
                e = effect.loc[gene_id, codes]
                r = reference.reindex(index=[gene_id], columns=codes).iloc[0]
                valid = e.notna() & r.notna()
                rv = r[valid]
                values.append(
                    (
                        int(valid.sum()),
                        float(rv.median()) if len(rv) else np.nan,
                        float(rv.min()) if len(rv) else np.nan,
                        float(rv.max()) if len(rv) else np.nan,
                        int((rv >= TAU).sum()),
                    )
                )
            values = pd.DataFrame(
                values,
                columns=[
                    f"{prefix}_{contrast}_effect_n",
                    f"{prefix}_{contrast}_ref_median",
                    f"{prefix}_{contrast}_ref_min",
                    f"{prefix}_{contrast}_ref_max",
                    f"{prefix}_{contrast}_prox_capacity",
                ],
            )
            for column in values:
                frame[column] = values[column].to_numpy()
            frame[f"{prefix}_{contrast}_ref_q5"] = quantile_codes(
                frame[f"{prefix}_{contrast}_ref_median"]
            )
            frame[f"{prefix}_{contrast}_pas_q3"] = quantile_codes(
                frame[f"{prefix}_pas"].replace(-1, np.nan), bins=3
            )
            frame[f"{prefix}_{contrast}_ref_min10"] = range10(
                frame[f"{prefix}_{contrast}_ref_min"]
            )
            frame[f"{prefix}_{contrast}_ref_max10"] = range10(
                frame[f"{prefix}_{contrast}_ref_max"]
            )

    base = ["h_n", "m_n", "h_pas", "m_pas"]
    results = []
    diagnostic_parts = []
    for contrast in ["PS", "PM"]:
        row, diagnostics = matched_null(frame, contrast, "published_joint_meas_pas", base)
        results.append(row)
        diagnostic_parts.append(diagnostics.assign(side="joint"))
        for model, human_columns, mouse_columns in [
            (
                "independent_q5_reference_capacity",
                [
                    f"h_{contrast}_effect_n",
                    f"h_{contrast}_pas_q3",
                    f"h_{contrast}_prox_capacity",
                    f"h_{contrast}_ref_q5",
                ],
                [
                    f"m_{contrast}_effect_n",
                    f"m_{contrast}_pas_q3",
                    f"m_{contrast}_prox_capacity",
                    f"m_{contrast}_ref_q5",
                ],
            ),
            (
                "independent_range10_reference_capacity",
                [
                    f"h_{contrast}_effect_n",
                    f"h_{contrast}_pas_q3",
                    f"h_{contrast}_prox_capacity",
                    f"h_{contrast}_ref_min10",
                    f"h_{contrast}_ref_max10",
                ],
                [
                    f"m_{contrast}_effect_n",
                    f"m_{contrast}_pas_q3",
                    f"m_{contrast}_prox_capacity",
                    f"m_{contrast}_ref_min10",
                    f"m_{contrast}_ref_max10",
                ],
            ),
        ]:
            row, diagnostics = independent_null(
                frame, contrast, model, human_columns, mouse_columns
            )
            results.append(row)
            diagnostic_parts.append(diagnostics)

    results = pd.DataFrame(results)
    results["q_within_model_two_contrasts"] = np.nan
    for model, indexes in results.groupby("model").groups.items():
        results.loc[indexes, "q_within_model_two_contrasts"] = bh(results.loc[indexes, "empirical_p"])

    # The rebuilt original stratification must recover the published null mean.
    published = original_results[original_results.null.eq("matched_meas_pas")]
    rebuilt = results[results.model.eq("published_joint_meas_pas")]
    check = rebuilt.merge(published, on="contrast", suffixes=("_rebuilt", "_published"))
    assert (check.observed_rebuilt == check.observed_published).all()
    assert np.max(np.abs(check.analytic_null_mean - check.expected)) < 0.05

    frame.to_csv(OUT / "ortholog_matched_features.tsv.gz", sep="\t", index=False)
    results.to_csv(OUT / "overlap_null_results.tsv", sep="\t", index=False)
    pd.concat(diagnostic_parts, ignore_index=True).to_csv(
        OUT / "stratum_diagnostics.tsv.gz", sep="\t", index=False
    )
    manifest = {
        "status": "PASS",
        "frozen_members": {"PS": 18, "PM": 10},
        "threshold": TAU,
        "permutations": B,
        "models": {
            "published_joint_meas_pas": "published joint measurement and supported-PAS strata",
            "independent_q5_reference_capacity": "independent species shuffles matched on contrast measurement count, PAS tercile, proximal capacity, and median reference-WPP quintile",
            "independent_range10_reference_capacity": "independent species shuffles matched on contrast measurement count, PAS tercile, proximal capacity, and min/max reference-WPP 0.1 bins",
        },
        "input_sha256": {
            str(path): digest(path)
            for path in [nodes_path, ortholog_path, pm_path, original_feature_path, original_result_path]
        },
    }
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
