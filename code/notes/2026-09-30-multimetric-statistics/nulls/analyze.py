"""Metric-specific, study-level recurrence and documented-source sharing nulls."""
from pathlib import Path
import argparse
import hashlib
import json
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
THRESHOLDS = {"WPP": [.02, .05], "RWUI": [.02, .05], "MPRO": [.10, .20]}
HUMAN = ["GSE118258", "GSE178325"]
MOUSE = ["GSE103221", "GSE115943"]


def save(frame, name):
    frame.to_csv(HERE / name, sep="\t", index=False, na_rep="NA",
                 float_format="%.17g",
                 compression={"method": "gzip", "mtime": 0} if name.endswith(".gz") else None)


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def bh(values):
    p = np.asarray(values, float)
    assert np.isfinite(p).all() and ((p >= 0) & (p <= 1)).all()
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(q, 1)
    return result


def qbin(values, q=3):
    result = pd.Series("NA", index=values.index, dtype=object)
    present = values.notna()
    if values[present].nunique() < 2:
        result.loc[present] = "0"
    else:
        result.loc[present] = (np.minimum(np.ceil(values[present].rank(
            method="average", pct=True) * q).astype(int), q) - 1).astype(str)
    return result


def collapse(routes):
    """One study vote requires threshold agreement across all measurable series."""
    records = []
    keys = ["species", "contrast", "study_id", "gene_id"]
    for key, frame in routes.groupby(keys, sort=True):
        assert frame.route.nunique() == len(frame)
        base = dict(zip(["species", "contrast", "unit", "gene_id"], key))
        base["gene_symbol"] = frame.gene_symbol.iloc[0]
        base["n_routes"] = len(frame)
        for column in ["mean_GEX", "APA_coverage", "PAS_number"]:
            base[column] = frame[column].median()
        base["supported_catalogue_PAS_number"] = base["PAS_number"]
        for metric, thresholds in THRESHOLDS.items():
            effect_column = "MPRO" if metric == "MPRO" else "delta_" + metric
            effect = frame[effect_column].to_numpy(float)
            assert np.isfinite(effect).all()
            row = base | {"metric": metric, "effect_mean": effect.mean(),
                          "effect_min": effect.min(), "effect_max": effect.max()}
            if metric == "MPRO":
                row["baseline"] = frame.reference_concentration.median()
                row["baseline_min"] = frame.reference_concentration.min()
                row["baseline_max"] = frame.reference_concentration.max()
                # MPRO range(delta p) maximizes at a point mass on a minimum-p site.
                row["capacity_unsigned"] = (1 - frame.reference_min_usage + frame.reference_max_usage).min()
                row["capacity_negative"] = np.nan
                row["capacity_positive"] = np.nan
            else:
                reference = frame["reference_" + metric]
                row["baseline"] = reference.median()
                row["baseline_min"] = reference.min()
                row["baseline_max"] = reference.max()
                low = 0 if metric == "WPP" else frame.rwui_min_weight
                high = 1 if metric == "WPP" else frame.rwui_max_weight
                row["capacity_negative"] = (reference - low).min()
                row["capacity_positive"] = (high - reference).min()
                row["capacity_unsigned"] = np.nan
            for threshold in thresholds:
                event = int(np.all(effect >= threshold))
                if metric != "MPRO" and np.all(effect <= -threshold):
                    event = -1
                records.append(row | {"threshold": threshold, "event": event})
    result = pd.DataFrame(records)
    result["eligible_units"] = result.groupby(["species", "contrast", "metric", "threshold", "gene_id"]).unit.transform("nunique")
    return result


def stratify(frame, metric, threshold, pair=False):
    f = frame.reset_index(drop=True).copy()
    pieces = ["eligible_units"] + (["n_routes"] if pair else [])
    for column in ["mean_GEX", "APA_coverage", "PAS_number"]:
        f[column + "_bin"] = qbin(f[column])
        pieces.append(column + "_bin")
    if metric == "MPRO":
        f["concentration_bin"] = qbin(f.baseline)
        f["capacity_unsigned_bin"] = qbin(f.capacity_unsigned)
        f["can_event"] = f.capacity_unsigned.ge(threshold - 1e-9)
        assert f.loc[f.event.eq(1), "can_event"].all()
        pieces += ["concentration_bin", "capacity_unsigned_bin", "can_event"]
    else:
        f["baseline_bin"] = qbin(f.baseline, q=5)
        f["can_negative"] = f.capacity_negative.ge(threshold - 1e-9)
        f["can_positive"] = f.capacity_positive.ge(threshold - 1e-9)
        assert f.loc[f.event.eq(-1), "can_negative"].all()
        assert f.loc[f.event.eq(1), "can_positive"].all()
        pieces += ["baseline_bin", "can_negative", "can_positive"]
    f["stratum"] = f[pieces].astype(str).agg("|".join, axis=1)
    return f


def shuffle(frame, rng, n):
    """Uniform label permutation within each exact matching stratum."""
    labels = frame.event.to_numpy(np.int8)
    draws = np.broadcast_to(labels, (n, len(frame))).copy()
    probabilities = np.empty((len(frame), 3), float)
    diagnostics = []
    for stratum, group in frame.groupby("stratum", sort=True):
        idx = group.index.to_numpy()
        lab = labels[idx]
        frequencies = np.array([(lab == value).sum() for value in [-1, 0, 1]])
        probabilities[idx] = frequencies / len(idx)
        mixed = np.count_nonzero(frequencies) > 1
        if mixed:
            draws[:, idx] = rng.permuted(np.broadcast_to(lab, (n, len(idx))), axis=1)
        for value, count in zip([-1, 0, 1], frequencies):
            assert ((draws[:, idx] == value).sum(axis=1) == count).all()
        diagnostics.append(dict(stratum=stratum, genes=len(idx),
            negative=int(frequencies[0]), zero=int(frequencies[1]), positive=int(frequencies[2]),
            mixed_labels=mixed, singleton=len(idx) == 1,
            expected_labels_changed=float(len(idx) * (1 - np.sum((frequencies / len(idx)) ** 2))),
            baseline_min=float(group.baseline.min()), baseline_max=float(group.baseline.max())))
    if "can_event" in frame:
        assert not ((draws == 1) & ~frame.can_event.to_numpy()[None, :]).any()
    else:
        assert not ((draws == -1) & ~frame.can_negative.to_numpy()[None, :]).any()
        assert not ((draws == 1) & ~frame.can_positive.to_numpy()[None, :]).any()
    return draws, probabilities, pd.DataFrame(diagnostics)


def summarize(observed, null):
    expected = float(null.mean())
    return dict(observed=int(observed), null_mean=expected, null_median=float(np.median(null)),
                null_p025=float(np.quantile(null, .025)), null_p975=float(np.quantile(null, .975)),
                null_sd=float(null.std(ddof=1)), null_unique=len(np.unique(null)),
                enrichment=float(observed / expected) if expected else np.nan,
                empirical_p=float((1 + np.sum(null >= observed)) / (len(null) + 1)),
                inferentially_mobile=bool(null.std(ddof=1) > 0))


def recurrence(events, n, seed):
    summaries, members, strata, diagnostics, archive = [], [], [], [], {}
    for (metric, species, contrast, threshold), f in events.groupby(["metric", "species", "contrast", "threshold"]):
        if f.unit.nunique() < 2:
            raise ValueError(f"Cannot test recurrence: {(metric, species, contrast, threshold)}")
        context = dict(metric=metric, species=species, contrast=contrast, threshold=threshold)
        key = "recurrence__" + "__".join(map(str, context.values()))
        run_seed = seed + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(run_seed)
        genes = sorted(f.gene_id.unique())
        index = {gene: i for i, gene in enumerate(genes)}
        positives = np.zeros((n, len(genes)), np.int8)
        negatives = np.zeros_like(positives)
        observed_positive = np.zeros(len(genes), int)
        observed_negative = np.zeros(len(genes), int)
        for unit, sub in f.groupby("unit"):
            sf = stratify(sub, metric, threshold)
            draws, _, diag = shuffle(sf, rng, n)
            ii = np.array([index[gene] for gene in sf.gene_id])
            positives[:, ii] += draws == 1
            negatives[:, ii] += draws == -1
            observed_positive[ii] += sf.event.eq(1).to_numpy()
            observed_negative[ii] += sf.event.eq(-1).to_numpy()
            strata.append(sf.assign(analysis="recurrence"))
            diagnostics.append(diag.assign(analysis="recurrence", unit=unit, **context))
        # Original Figure 3 permits opposing study labels; retain this exact rule.
        support = np.maximum(observed_positive, observed_negative)
        null = (np.maximum(positives, negatives) >= 2).sum(axis=1)
        summaries.append(context | summarize((support >= 2).sum(), null) | {"seed": run_seed, "permutations": n})
        labels = f.pivot(index="gene_id", columns="unit", values="event").reindex(genes)
        eligible = labels.notna().sum(axis=1).to_numpy()
        for i, gene in enumerate(genes):
            if support[i] < 2:
                continue
            row = dict(gene_id=gene, gene_symbol=f.loc[f.gene_id.eq(gene), "gene_symbol"].iloc[0],
                       eligible_units=int(eligible[i]), positive_support=int(observed_positive[i]),
                       negative_support=int(observed_negative[i]), support=int(support[i]),
                       opposing_support=int(min(observed_positive[i], observed_negative[i])),
                       direction="unsigned" if metric == "MPRO" else "positive" if observed_positive[i] >= observed_negative[i] else "negative")
            members.append(row | context)
        archive[key] = null
        print(key, summaries[-1], flush=True)
    summary = pd.DataFrame(summaries)
    summary["q_metric_8"] = summary.groupby("metric").empirical_p.transform(bh)
    summary["q_all_24"] = bh(summary.empirical_p)
    save(summary, "recurrence_summary.tsv")
    save(pd.DataFrame(members), "recurrence_members.tsv.gz")
    save(pd.concat(strata, ignore_index=True), "recurrence_strata.tsv.gz")
    return diagnostics, archive


def cross_species(events, ortholog, n, seed):
    summaries, backgrounds, members, strata, diagnostics, archive = [], [], [], [], [], {}
    for metric, thresholds in THRESHOLDS.items():
        for h in HUMAN:
            for m in MOUSE:
                for threshold in thresholds:
                    context = dict(metric=metric, human_study=h, mouse_study=m, threshold=threshold)
                    selected = events[events.metric.eq(metric) & events.threshold.eq(threshold) & events.contrast.eq("P_S")]
                    hs = selected[selected.unit.eq(h)]
                    ms = selected[selected.unit.eq(m)]
                    f = ortholog.merge(hs.add_prefix("h_"), left_on="human_gene_id", right_on="h_gene_id", validate="one_to_one")
                    f = f.merge(ms.add_prefix("m_"), left_on="mouse_gene_id", right_on="m_gene_id", validate="one_to_one")
                    f = f.sort_values("human_gene_id").reset_index(drop=True)
                    assert len(f) > 0
                    key = "sharing__" + "__".join(map(str, context.values()))
                    run_seed = seed + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
                    rng = np.random.default_rng(run_seed)
                    draws, probabilities = [], []
                    for species, prefix in [("human", "h_"), ("mouse", "m_")]:
                        sub = f[[column for column in f if column.startswith(prefix)]].rename(columns=lambda x: x[len(prefix):])
                        sf = stratify(sub, metric, threshold, pair=True)
                        shuffled, prob, diag = shuffle(sf, rng, n)
                        draws.append(shuffled)
                        probabilities.append(prob)
                        strata.append(sf.assign(human_gene_id=f.human_gene_id, mouse_gene_id=f.mouse_gene_id, analysis="sharing", **context))
                        diagnostics.append(diag.assign(analysis="sharing", species=species, **context))
                    tests = ["shared_unsigned"] if metric == "MPRO" else ["same_direction", "shared_negative", "shared_positive"]
                    backgrounds.append(f.assign(**context))
                    shared = f.h_event.eq(f.m_event) & f.h_event.ne(0)
                    members.append(f[shared].assign(**context))
                    for test in tests:
                        values = [1] if test == "shared_unsigned" or test == "shared_positive" else [-1] if test == "shared_negative" else [-1, 1]
                        observed = sum(((f.h_event == value) & (f.m_event == value)).sum() for value in values)
                        null = sum(((draws[0] == value) & (draws[1] == value)).sum(axis=1) for value in values)
                        analytic = sum((probabilities[0][:, value + 1] * probabilities[1][:, value + 1]).sum() for value in values)
                        assert abs(null.mean() - analytic) <= 6 * null.std(ddof=1) / np.sqrt(n) + 1e-10
                        summaries.append(context | {"test": test, "common_orthologs": len(f), "analytic_null_mean": analytic,
                                                    "seed": run_seed, "permutations": n} | summarize(observed, null))
                        archive[key + "__" + test] = null
                    print(key, len(f), summaries[-len(tests)], flush=True)
    summary = pd.DataFrame(summaries)
    summary["q_metric_all"] = summary.groupby("metric").empirical_p.transform(bh)
    for metric, f in summary.groupby("metric"):
        primary = f[f.test.isin(["same_direction", "shared_unsigned"])]
        assert len(primary) == 8
        summary.loc[primary.index, "q_primary_8"] = bh(primary.empirical_p)
        directional = f[~f.test.isin(["same_direction", "shared_unsigned"])]
        if len(directional):
            assert len(directional) == 16
            summary.loc[directional.index, "q_directional_16"] = bh(directional.empirical_p)
    primary = summary[summary.test.isin(["same_direction", "shared_unsigned"])]
    summary.loc[primary.index, "q_joint_primary_24"] = bh(primary.empirical_p)
    summary["q_joint_all_56"] = bh(summary.empirical_p)
    save(summary, "sharing_summary.tsv")
    save(pd.concat(backgrounds, ignore_index=True), "sharing_backgrounds.tsv.gz")
    save(pd.concat(members, ignore_index=True), "sharing_members.tsv.gz")
    save(pd.concat(strata, ignore_index=True), "sharing_strata.tsv.gz")
    return diagnostics, archive


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--effects", type=Path, required=True)
    parser.add_argument("--ortholog", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260930)
    parser.add_argument("--legacy-events", type=Path)
    args = parser.parse_args()
    started = time.monotonic()
    routes = pd.read_csv(args.effects, sep="\t")
    required = ["species", "study_id", "route", "gene_id", "gene_symbol", "contrast", "delta_WPP", "reference_WPP",
                "delta_RWUI", "reference_RWUI", "rwui_min_weight", "rwui_max_weight", "MPRO", "reference_concentration",
                "reference_min_usage", "reference_max_usage", "mean_GEX", "APA_coverage", "PAS_number"]
    missing = set(required) - set(routes.columns)
    if missing:
        raise ValueError(f"Missing input columns: {sorted(missing)}")
    routes = routes[routes.contrast.isin(["P_S", "P_M"])].copy()
    assert not routes.duplicated(["species", "study_id", "route", "gene_id", "contrast"]).any()
    assert routes[["APA_coverage", "PAS_number"]].notna().all().all()
    assert routes.PAS_number.ge(2).all()
    assert routes.MPRO.between(-1e-9, 2 + 1e-9).all()
    assert routes.reference_concentration.between(0, 1 + 1e-9).all()
    assert routes.reference_min_usage.between(-1e-9, .5 + 1e-9).all()
    ortholog = pd.read_csv(args.ortholog, sep="\t")
    if "homology_type" in ortholog:
        assert ortholog.homology_type.eq("ortholog_one2one").all()
    ortholog = ortholog[["human_gene_id", "mouse_gene_id"]]
    assert not ortholog.human_gene_id.duplicated().any()
    assert not ortholog.mouse_gene_id.duplicated().any()
    events = collapse(routes)
    legacy_audit = None
    if args.legacy_events:
        old = pd.read_csv(args.legacy_events, sep="\t")
        if "level" in old:
            old = old[old.level.eq("study")]
        keys = ["species", "contrast", "unit", "gene_id", "threshold"]
        actual = events[events.metric.eq("WPP")]
        check = actual.merge(old, on=keys, how="outer", suffixes=("_new", "_old"), indicator=True,
                             validate="one_to_one")
        assert check._merge.eq("both").all(), check._merge.value_counts().to_dict()
        assert check.event_new.eq(check.event_old).all()
        assert check.n_routes_new.eq(check.n_routes_old).all()
        differences = {col: float((check[col + "_new"] - check[col + "_old"]).abs().max())
                       for col in ["mean_GEX", "APA_coverage", "PAS_number", "eligible_units"]}
        baseline_difference = float((check.baseline - check.reference_WPP).abs().max())
        baseline_min_difference = float((check.baseline_min - check.reference_WPP_min).abs().max())
        baseline_max_difference = float((check.baseline_max - check.reference_WPP_max).abs().max())
        cov_keys = ["species", "contrast", "unit", "gene_id"]
        cov_columns = ["mean_GEX", "APA_coverage", "PAS_number"]
        assert old.groupby(cov_keys)[cov_columns].nunique(dropna=False).le(1).all().all()
        covariates = old[cov_keys + cov_columns].drop_duplicates(cov_keys)
        save(covariates, "original_matching_covariates.tsv.gz")
        events = events.rename(columns={col: "reconstructed_" + col for col in cov_columns})
        events = events.merge(covariates, on=cov_keys, how="left", validate="many_to_one", indicator=True)
        assert events._merge.eq("both").all()
        events = events.drop(columns="_merge")
        save(events[cov_keys + ["metric", "threshold", "supported_catalogue_PAS_number"] + cov_columns +
                    ["reconstructed_" + col for col in cov_columns]], "covariate_adaptation.tsv.gz")
        legacy_audit = dict(rows=len(check), all_WPP_events_identical=True, all_series_counts_identical=True,
                            maximum_raw_covariate_differences=differences, path=str(args.legacy_events),
                            sha256=sha(args.legacy_events),
                            covariate_adaptation="all metrics use unchanged original study-level mean_GEX/APA_coverage/PAS_number; catalogue count retained separately",
                            maximum_WPP_reference_difference=baseline_difference,
                            maximum_WPP_reference_min_difference=baseline_min_difference,
                            maximum_WPP_reference_max_difference=baseline_max_difference,
                            WPP_null_scope="conditional reconstructed baseline comparator; not a rerun of joint six-rule consensus")
        (HERE / "legacy_event_validation.json").write_text(json.dumps(legacy_audit, indent=2) + "\n")
    save(routes, "route_effects_used.tsv.gz")
    save(events, "study_event_covariates.tsv.gz")
    diag1, archive1 = recurrence(events, args.permutations, args.seed)
    diag2, archive2 = cross_species(events, ortholog, args.permutations, args.seed)
    diagnostics = pd.concat(diag1 + diag2, ignore_index=True)
    save(diagnostics, "stratum_diagnostics.tsv.gz")
    diagnostics["mixed_genes"] = diagnostics.genes * diagnostics.mixed_labels
    diagnostics["singleton_genes"] = diagnostics.genes * diagnostics.singleton
    keys = ["analysis", "metric", "species", "contrast", "unit", "human_study", "mouse_study", "threshold"]
    mobility = diagnostics.groupby(keys, dropna=False).agg(genes=("genes", "sum"), strata=("stratum", "size"),
        mixed_genes=("mixed_genes", "sum"), singleton_genes=("singleton_genes", "sum"),
        expected_labels_changed=("expected_labels_changed", "sum")).reset_index()
    mobility["mixed_label_fraction"] = mobility.mixed_genes / mobility.genes
    mobility["expected_changed_fraction"] = mobility.expected_labels_changed / mobility.genes
    save(mobility, "null_mobility.tsv")
    np.savez_compressed(HERE / "null_draws.npz", **(archive1 | archive2))
    audit = dict(status="computed", inputs_sha256={str(path): sha(path) for path in [args.effects, args.ortholog]},
                 script_sha256=sha(__file__), rows=len(routes), event_rows=len(events),
                 permutations=args.permutations, seed=args.seed, thresholds=THRESHOLDS,
                 minutes=(time.monotonic() - started) / 60,
                 legacy_event_validation=legacy_audit,
                 recurrence_rule=">=2 matching signed study labels; opposing studies allowed as in original Figure 3; MPRO unsigned >=2",
                 inference_unit="study-level set label conditional on observed effects; not gene-level or biological-replicate inference",
                 metric_definitions="RWUI i/n; MPRO max_i(delta usage_i)-min_i(delta usage_i); concentration sum_i(reference usage_i squared)")
    (HERE / "audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
