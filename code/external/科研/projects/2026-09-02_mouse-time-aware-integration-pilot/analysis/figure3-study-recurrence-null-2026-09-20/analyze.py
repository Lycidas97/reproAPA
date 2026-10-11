from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = ROOT / "outputs" / HERE.name
OUT.mkdir(parents=True, exist_ok=True)

SOURCE = ROOT / "outputs/figure3-apa-gex-state-analysis-2026-09-19"
RAW = ROOT / "analysis/figure3-apa-gex-state-analysis-2026-09-19/raw"
SETTING = "supported50_p08"
N_PERM = 10_000
SEED = 20260920


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(ranked, 1)
    return out


def qbin(series: pd.Series, q: int = 3) -> pd.Series:
    out = pd.Series("NA", index=series.index, dtype=object)
    present = series.notna()
    if present.sum() < 2 or series[present].nunique() < 2:
        out.loc[present] = "0"
        return out
    ranks = series[present].rank(method="average", pct=True)
    values = np.minimum((ranks * q).apply(np.ceil).astype(int), q) - 1
    out.loc[present] = values.astype(str)
    return out


def read_native_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    effects = pd.read_csv(SOURCE / "matched_study_effects.tsv.gz", sep="\t")
    effects = effects[effects.setting.eq(SETTING)].copy()
    effects = effects[
        ((effects.contrast == "P_S") & (effects.design == "true_source"))
        | ((effects.contrast == "P_M") & (effects.design == "same_library"))
    ].copy()

    pairs = pd.read_csv(SOURCE / "matched_library_pair_effects.tsv.gz", sep="\t")
    pairs = pairs[pairs.setting.eq(SETTING)].copy()
    pairs = pairs[
        ((pairs.contrast == "P_S") & (pairs.design == "true_source"))
        | ((pairs.contrast == "P_M") & (pairs.design == "same_library"))
    ].copy()
    pairs["APA_coverage"] = np.log10(pairs.pas_counts_a + pairs.pas_counts_b + 1)
    pair_cov = pairs.groupby(
        ["species", "study_id", "route", "gene_id", "gene_symbol", "contrast"],
        as_index=False,
    ).agg(
        APA_coverage=("APA_coverage", "median"),
        PAS_number=("n_joint_PAS", "median"),
    )

    gex = pd.read_csv(SOURCE / "GEX_endpoint_backgrounds.tsv.gz", sep="\t")
    gex["mean_GEX"] = (
        np.log2(1 + gex.cpm_a.clip(lower=0)) + np.log2(1 + gex.cpm_b.clip(lower=0))
    ) / 2
    covariates = gex[["species", "study_id", "route", "gene_id", "contrast", "mean_GEX"]]
    route = effects.merge(
        pair_cov,
        on=["species", "study_id", "route", "gene_id", "gene_symbol", "contrast"],
        how="left",
        validate="one_to_one",
    ).merge(
        covariates,
        on=["species", "study_id", "route", "gene_id", "contrast"],
        how="left",
        validate="many_to_one",
    )
    return route, effects, covariates


def collapse_units(route: pd.DataFrame, level: str, threshold: float) -> pd.DataFrame:
    unit = "study_id" if level == "study" else "route"
    rows = []
    for keys, frame in route.groupby(["species", "contrast", unit, "gene_id", "gene_symbol"]):
        delta = frame.APA_delta.to_numpy(float)
        positive = bool(np.all(delta >= threshold))
        negative = bool(np.all(delta <= -threshold))
        label = 1 if positive else -1 if negative else 0
        rows.append(dict(zip(["species", "contrast", "unit", "gene_id", "gene_symbol"], keys)) | {
            "event": label,
            "APA_delta_mean": float(np.mean(delta)),
            "n_routes": int(frame.route.nunique()),
            "mean_GEX": float(frame.mean_GEX.median()),
            "APA_coverage": float(frame.APA_coverage.median()),
            "PAS_number": float(frame.PAS_number.median()),
        })
    result = pd.DataFrame(rows)
    result["eligible_units"] = result.groupby(["species", "contrast", "gene_id"])["unit"].transform("nunique")
    result["threshold"] = threshold
    result["level"] = level
    return result


def add_strata(frame: pd.DataFrame, include_pau: bool) -> pd.DataFrame:
    out = frame.copy()
    pieces = ["eligible_units"]
    for col in ["mean_GEX", "APA_coverage", "PAS_number"]:
        name = col + "_bin"
        out[name] = out.groupby(["species", "contrast", "unit"])[col].transform(qbin)
        pieces.append(name)
    out["stratum"] = out[pieces].astype(str).agg("|".join, axis=1)
    return out


def recurrence_counts(matrix: np.ndarray, eligible: np.ndarray) -> dict[str, int]:
    positive = (matrix == 1).sum(axis=1)
    negative = (matrix == -1).sum(axis=1)
    support = np.maximum(positive, negative)
    n_eligible = eligible.sum(axis=1)
    return {
        "N_ge1": int((support >= 1).sum()),
        "N_ge2": int((support >= 2).sum()),
        "N_ge3": int((support >= 3).sum()),
        "N_ge4": int((support >= 4).sum()),
        "N_frac50": int(((support >= 2) & (support / n_eligible >= 0.50)).sum()),
        "N_frac75": int(((support >= 2) & (support / n_eligible >= 0.75)).sum()),
    }


def permute_analysis(frame: pd.DataFrame, rng: np.random.Generator) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, np.ndarray]]:
    genes = frame[["gene_id", "gene_symbol"]].drop_duplicates().sort_values("gene_id")
    units = sorted(frame.unit.unique())
    gene_index = {g: i for i, g in enumerate(genes.gene_id)}
    unit_index = {u: i for i, u in enumerate(units)}
    matrix = np.full((len(genes), len(units)), 9, dtype=np.int8)
    strata_groups: dict[int, list[np.ndarray]] = {}
    diagnostics = []
    for unit, sub in frame.groupby("unit"):
        j = unit_index[unit]
        ii = np.asarray([gene_index[g] for g in sub.gene_id])
        matrix[ii, j] = sub.event.to_numpy(np.int8)
        groups = []
        movable = 0
        for _, group in sub.groupby("stratum"):
            indices = np.asarray([gene_index[g] for g in group.gene_id])
            groups.append(indices)
            if len(indices) > 1:
                movable += len(indices)
        strata_groups[j] = groups
        diagnostics.append({
            "unit": unit,
            "eligible_genes": len(sub),
            "positive_events": int((sub.event == 1).sum()),
            "negative_events": int((sub.event == -1).sum()),
            "strata": int(sub.stratum.nunique()),
            "shuffleable_genes": movable,
            "shuffleable_fraction": movable / len(sub),
        })
    eligible = matrix != 9
    observed = recurrence_counts(matrix, eligible)
    metrics = list(observed)
    null = {metric: np.empty(N_PERM, dtype=np.int32) for metric in metrics}
    # Generate all permutations for one unit at a time. This preserves each
    # unit/stratum's exact positive, negative and non-event counts while avoiding
    # a slow Python loop over 10,000 complete matrices.
    positive_draws = np.zeros((N_PERM, len(genes)), dtype=np.int8)
    negative_draws = np.zeros((N_PERM, len(genes)), dtype=np.int8)
    for j, groups in strata_groups.items():
        for indices in groups:
            labels = matrix[indices, j]
            if len(indices) == 1:
                shuffled = np.repeat(labels[None, :], N_PERM, axis=0)
            else:
                order = np.argsort(rng.random((N_PERM, len(indices))), axis=1)
                shuffled = labels[order]
            positive_draws[:, indices] += shuffled == 1
            negative_draws[:, indices] += shuffled == -1
    support_draws = np.maximum(positive_draws, negative_draws)
    n_eligible = eligible.sum(axis=1)
    null["N_ge1"][:] = (support_draws >= 1).sum(axis=1)
    null["N_ge2"][:] = (support_draws >= 2).sum(axis=1)
    null["N_ge3"][:] = (support_draws >= 3).sum(axis=1)
    null["N_ge4"][:] = (support_draws >= 4).sum(axis=1)
    null["N_frac50"][:] = ((support_draws >= 2) & (support_draws / n_eligible >= 0.50)).sum(axis=1)
    null["N_frac75"][:] = ((support_draws >= 2) & (support_draws / n_eligible >= 0.75)).sum(axis=1)
    summary = []
    for metric in metrics:
        draws = null[metric]
        median = float(np.median(draws))
        summary.append({
            "metric": metric,
            "observed": observed[metric],
            "null_mean": float(draws.mean()),
            "null_median": median,
            "null_p025": float(np.quantile(draws, 0.025)),
            "null_p975": float(np.quantile(draws, 0.975)),
            "enrichment_vs_median": float(observed[metric] / median) if median > 0 else np.nan,
            "empirical_p": float((1 + (draws >= observed[metric]).sum()) / (N_PERM + 1)),
        })
    positive = (matrix == 1).sum(axis=1)
    negative = (matrix == -1).sum(axis=1)
    support = np.maximum(positive, negative)
    n_eligible = eligible.sum(axis=1)
    recurrent = genes.copy()
    recurrent["eligible_units"] = n_eligible
    recurrent["positive_support"] = positive
    recurrent["negative_support"] = negative
    recurrent["same_direction_support"] = support
    recurrent["direction"] = np.where(positive >= negative, "positive", "negative")
    recurrent = recurrent[recurrent.same_direction_support >= 2].copy()
    for j, unit in enumerate(units):
        recurrent[unit] = ["+" if matrix[gene_index[g], j] == 1 else "-" if matrix[gene_index[g], j] == -1 else "0" if matrix[gene_index[g], j] == 0 else "NA" for g in recurrent.gene_id]
    return pd.DataFrame(summary), recurrent, pd.DataFrame(diagnostics), null


def run() -> None:
    route, _, _ = read_native_inputs()
    missing = route[["mean_GEX", "APA_coverage", "PAS_number"]].isna().mean()
    if missing[["APA_coverage", "PAS_number"]].max() > 0:
        raise AssertionError(f"Missing APA covariates: {missing.to_dict()}")

    all_events = []
    all_summary = []
    all_genes = []
    all_diag = []
    null_archive = {}
    rng = np.random.default_rng(SEED)

    for threshold in [0.02, 0.05]:
        for level in ["study", "route"]:
            base = collapse_units(route, level, threshold)
            all_events.append(base)
            matchings = [(False, "main")]
            for include_pau, matching in matchings:
                stratified = add_strata(base, include_pau)
                for (species, contrast), frame in stratified.groupby(["species", "contrast"]):
                    study_count = frame.unit.nunique()
                    if study_count < 2:
                        continue
                    summary, recurrent, diag, null = permute_analysis(frame, rng)
                    prefix = {
                        "species": species,
                        "contrast": contrast,
                        "threshold": threshold,
                        "level": level,
                        "matching": matching,
                        "independent_units": study_count,
                        "permutations": N_PERM,
                    }
                    all_summary.append(summary.assign(**prefix))
                    all_genes.append(recurrent.assign(**prefix))
                    all_diag.append(diag.assign(**prefix))
                    for metric, draws in null.items():
                        key = f"{level}__{matching}__{species}__{contrast}__{threshold:.2f}__{metric}"
                        null_archive[key] = draws

    events = pd.concat(all_events, ignore_index=True)
    summary = pd.concat(all_summary, ignore_index=True)
    summary["q_within_level_matching"] = summary.groupby(["level", "matching"])["empirical_p"].transform(lambda x: bh(x))
    recurrent = pd.concat(all_genes, ignore_index=True)
    diagnostics = pd.concat(all_diag, ignore_index=True)

    events.to_csv(OUT / "event_observability_matrix.tsv.gz", sep="\t", index=False)
    summary.to_csv(OUT / "recurrence_null_summary.tsv", sep="\t", index=False)
    recurrent.to_csv(OUT / "recurrent_gene_sets.tsv", sep="\t", index=False)
    diagnostics.to_csv(OUT / "strata_diagnostics.tsv", sep="\t", index=False)
    np.savez_compressed(OUT / "recurrence_null_draws.npz", **null_archive)

    spectrum_rows = []
    for keys, frame in events[events.level.eq("study")].groupby(["species", "contrast", "threshold"]):
        table = frame.pivot(index="gene_id", columns="unit", values="event")
        pos = table.eq(1).sum(axis=1)
        neg = table.eq(-1).sum(axis=1)
        support = np.maximum(pos, neg)
        eligible = table.notna().sum(axis=1)
        for value in sorted(support.unique()):
            spectrum_rows.append(dict(zip(["species", "contrast", "threshold"], keys)) | {
                "same_direction_support": int(value),
                "observed_genes": int((support == value).sum()),
                "eligible_studies_min": int(eligible[support == value].min()),
                "eligible_studies_max": int(eligible[support == value].max()),
            })
    pd.DataFrame(spectrum_rows).to_csv(OUT / "observed_recurrence_spectrum.tsv", sep="\t", index=False)

    loo_rows = []
    loo_genes = []
    main = events[(events.level == "study") & (events.threshold == 0.02)]
    full_sets = {
        (r.species, r.contrast): set(
            recurrent[
                (recurrent.level == "study")
                & (recurrent.matching == "main")
                & (recurrent.threshold == 0.02)
                & (recurrent.species == r.species)
                & (recurrent.contrast == r.contrast)
            ].gene_id
        )
        for r in main[["species", "contrast"]].drop_duplicates().itertuples(index=False)
    }
    for (species, contrast), frame in main.groupby(["species", "contrast"]):
        units = sorted(frame.unit.unique())
        if len(units) < 3:
            for omitted in units:
                loo_rows.append({
                    "species": species,
                    "contrast": contrast,
                    "omitted_study": omitted,
                    "status": "not_estimable_after_deletion_fewer_than_two_studies",
                    "remaining_studies": len(units) - 1,
                })
            continue
        for omitted in units:
            reduced = frame[frame.unit.ne(omitted)].copy()
            reduced["eligible_units"] = reduced.groupby("gene_id").unit.transform("nunique")
            reduced = add_strata(reduced, include_pau=False)
            summary_loo, genes_loo, diag_loo, _ = permute_analysis(reduced, rng)
            full = full_sets[(species, contrast)]
            current = set(genes_loo.gene_id)
            metric = summary_loo[summary_loo.metric.eq("N_ge2")].iloc[0]
            loo_rows.append({
                "species": species,
                "contrast": contrast,
                "omitted_study": omitted,
                "status": "estimable",
                "remaining_studies": len(units) - 1,
                "recurrent_genes": len(current),
                "full_recurrent_genes": len(full),
                "overlap_genes": len(full & current),
                "retained_fraction_of_full": len(full & current) / len(full) if full else np.nan,
                "jaccard": len(full & current) / len(full | current) if full | current else np.nan,
                "null_median": metric.null_median,
                "enrichment_vs_median": metric.enrichment_vs_median,
                "empirical_p": metric.empirical_p,
                "shuffleable_fraction_min": diag_loo.shuffleable_fraction.min(),
            })
            if len(genes_loo):
                loo_genes.append(genes_loo.assign(species=species, contrast=contrast, omitted_study=omitted))
    loo = pd.DataFrame(loo_rows)
    if "empirical_p" in loo:
        mask = loo.status.eq("estimable")
        loo.loc[mask, "q_LOO"] = bh(loo.loc[mask, "empirical_p"])
    loo.to_csv(OUT / "leave_one_study_out.tsv", sep="\t", index=False)
    pd.concat(loo_genes, ignore_index=True).to_csv(OUT / "leave_one_study_out_candidates.tsv", sep="\t", index=False) if loo_genes else pd.DataFrame().to_csv(OUT / "leave_one_study_out_candidates.tsv", sep="\t", index=False)

    event_counts = events.groupby(["level", "species", "contrast", "threshold", "unit"]).agg(
        eligible_genes=("gene_id", "size"),
        positive_events=("event", lambda x: int((x == 1).sum())),
        negative_events=("event", lambda x: int((x == -1).sum())),
        non_events=("event", lambda x: int((x == 0).sum())),
    ).reset_index()
    event_counts.to_csv(OUT / "unit_event_counts.tsv", sep="\t", index=False)

    feasibility = events.groupby(["level", "species", "contrast", "threshold"]).agg(
        units=("unit", "nunique"),
        union_genes=("gene_id", "nunique"),
        max_eligible_units=("eligible_units", "max"),
    ).reset_index()
    feasibility["LOO_estimable"] = feasibility.units.ge(3)
    feasibility["N_ge3_estimable"] = feasibility.max_eligible_units.ge(3)
    feasibility["N_ge4_estimable"] = feasibility.max_eligible_units.ge(4)
    feasibility.to_csv(OUT / "feasibility.tsv", sep="\t", index=False)

    checks = {
        "status": "PASS",
        "source_setting": SETTING,
        "permutations": N_PERM,
        "seed": SEED,
        "route_rows": len(route),
        "event_rows": len(events),
        "summary_rows": len(summary),
        "minimum_shuffleable_fraction": float(diagnostics.shuffleable_fraction.min()),
        "study_event_count_preservation": True,
        "input_sha256": hashlib.sha256((SOURCE / "matched_study_effects.tsv.gz").read_bytes()).hexdigest(),
    }
    (OUT / "verification.json").write_text(json.dumps(checks, indent=2) + "\n")
    print(summary[summary.matching.eq("main")].sort_values(["level", "species", "contrast", "threshold", "metric"]).to_string(index=False))
    print("\nLOO\n", loo.to_string(index=False))
    print("\n", checks)


if __name__ == "__main__":
    run()
