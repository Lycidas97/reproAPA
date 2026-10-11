"""Fixed recurrence cohorts on the common-contributor P/M/S measurements."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables"
SETTING = "supported50_p08"
N_DRAW = 10000
N_BOOT = 10000
SOURCES = {}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, **kwargs):
    SOURCES[str(path.relative_to(ROOT))] = digest(path)
    return pd.read_csv(path, sep="\t", float_precision="round_trip", **kwargs)


def save(frame, name):
    compression = {"method": "gzip", "mtime": 0} if name.endswith(".gz") else None
    pd.DataFrame(frame).to_csv(HERE / name, sep="\t", index=False, na_rep="NA",
                              float_format="%.17g", compression=compression)


def seed(*args):
    return 20261010 + int(hashlib.sha256("|".join(map(str, args)).encode()).hexdigest()[:12], 16)


def bh(values):
    values = np.asarray(values, float)
    assert np.isfinite(values).all()
    order = np.argsort(values)
    out = np.empty(len(values))
    out[order] = np.minimum(1, np.minimum.accumulate(
        (values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1])
    return out


def load_weights():
    path = ROOT / "analysis_code/pas_rerun/polyapipe_effects.py"
    SOURCES[str(path.relative_to(ROOT))] = digest(path)
    tree = ast.parse(path.read_text())
    functions = [x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "weights"]
    assert len(functions) == 1
    namespace = dict(np=np, pd=pd)
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace["weights"]


def validate_vectors(frame, triplets):
    checks = []
    for metric in ["DWUI", "RWUI", "PDUI", "PPUI"]:
        error = np.max(np.abs(frame[metric + "_PS"] - frame[metric + "_PM"] - frame[metric + "_MS"]))
        assert error < 1e-12
        checks.append(dict(check=metric + "_additivity_all_settings", n=len(frame), max_abs_error=error))
    vector_error = 0.
    for row in frame.itertuples():
        states = {state: np.asarray(json.loads(getattr(row, state + "_usage"))) for state in "PMS"}
        assert all(np.isclose(x.sum(), 1, atol=1e-12, rtol=0) for x in states.values())
        for co, a, b in [("PS", "P", "S"), ("PM", "P", "M"), ("MS", "M", "S")]:
            vector_error = max(vector_error, np.max(np.abs(
                np.asarray(json.loads(getattr(row, "delta_usage_" + co))) - states[a] + states[b])))
            assert abs(getattr(row, "DWUI_" + co) - getattr(row, a + "_DWUI") + getattr(row, b + "_DWUI")) < 1e-12
    assert vector_error < 1e-12
    checks.append(dict(check="state_and_full_usage_vectors_all_settings", n=len(frame), max_abs_error=vector_error))
    weights = load_weights()
    keys = ["species", "study_id", "route", "setting", "gene_id"]
    primary = frame[frame.setting.eq(SETTING)].set_index(keys)
    triplets = triplets[triplets.setting.eq(SETTING)]
    assert set(triplets.set_index(keys).index) == set(primary.index)
    error, n_checked = 0., 0
    for key, part in triplets.groupby(keys, sort=False):
        row = primary.loc[key]
        w = weights(part)
        assert np.allclose(w, part.triplet_weight, rtol=0, atol=1e-12)
        assert len(part) == row.n_triplets
        assert part.library_key_a.nunique() == row.n_target_libraries
        assert part.library_key_b.nunique() == row.n_source_libraries
        assert part.time_a.nunique() == row.n_target_dates
        assert part.MS_joint_PAS.ge(2).all()
        assert part.PAS_ids.eq(part.M_PAS_ids).all()
        assert part.target_usage.eq(part.PM_P_usage).all()
        for state, column in [("P", "target_usage"), ("M", "M_usage"), ("S", "reference_usage")]:
            estimate = w @ np.stack(part[column].map(json.loads))
            error = max(error, np.max(np.abs(estimate - np.asarray(json.loads(row[state + "_usage"])))))
        n_checked += 1
    assert error < 1e-12
    checks.append(dict(check="primary_common_weights_and_contributors_replayed", n=n_checked, max_abs_error=error))
    return checks


def dual_reference(matched, members, threshold):
    rows = []
    for species, routes in [("human", {"GSE118258"}), ("mouse", {"CIP", "WOT"})]:
        fixed = members[members.species.eq(species) & members.contrast.eq("P_M") & members.threshold.eq(threshold)]
        for row in fixed.itertuples():
            q = matched[matched.species.eq(species) & matched.gene_id.eq(row.gene_id) & matched.route.isin(routes)]
            assert not q.route.duplicated().any()
            if set(q.route) != routes:
                continue
            d = 1 if row.direction == "distal" else -1
            if (q[["DWUI_PS", "DWUI_PM"]] * d).ge(threshold).all().all():
                rows.append(dict(species=species, threshold=threshold, gene_id=row.gene_id,
                                 gene_symbol=row.gene_symbol, direction=row.direction,
                                 minimum_aligned_effect=float((q[["DWUI_PS", "DWUI_PM"]] * d).min().min()),
                                 route_scope=";".join(sorted(routes))))
    return pd.DataFrame(rows)


def bins(values, count):
    edges = np.unique(np.quantile(values, np.linspace(0, 1, count + 1)))
    # Equal values are never split using gene IDs or row ordering.
    inner = edges[1:-1]
    return np.searchsorted(inner, values, side="left"), edges.tolist()


def matching(frame, strength_bins, include_native_strength=False):
    frame = frame.copy()
    edges = {}
    columns = ["native_sign", "availability"]
    continuous = ["abs_PS", "S_DWUI", "M_DWUI", "mean_log_GEX", "mean_log_APA",
                  "native_PAS_number", "minimum_target_libraries"]
    if include_native_strength:
        continuous.append("native_mean_abs_effect")
    for column in continuous:
        assert np.isfinite(frame[column]).all(), column
        count = strength_bins if column == "abs_PS" else 2
        name = column + "_bin"
        frame[name], edges[column] = bins(frame[column].to_numpy(), count)
        columns.append(name)
    cells = list(frame.groupby(columns, dropna=False, sort=True).indices.values())
    return frame, cells, edges


def label_null(values, labels, cells, n_draw, rng):
    """Joint labels preserve both direction counts in every matched stratum."""
    counts = {d: int((labels == d).sum()) for d in [-1, 1]}
    assert all(counts.values())
    expected = {d: sum((labels[ii] == d).sum() * values[ii].mean() for ii in cells) / counts[d]
                for d in [-1, 1]}
    result = []
    for start in range(0, n_draw, 256):
        n = min(256, n_draw - start)
        shuffled = np.broadcast_to(labels, (n, len(labels))).copy()
        for ii in cells:
            if len(np.unique(labels[ii])) > 1:
                shuffled[:, ii] = rng.permuted(np.broadcast_to(labels[ii], (n, len(ii))), axis=1)
        assert all(np.all((shuffled == d).sum(1) == counts[d]) for d in [-1, 1])
        result.append(np.column_stack([(shuffled == d) @ values / counts[d] for d in [-1, 1]]))
    return np.concatenate(result), expected


def inference(observed, draws):
    sd = float(np.std(draws, ddof=1))
    informative = np.ptp(draws) > 1e-12
    # Dense dot products and filtered means can differ at round-off for the same labels.
    tail_count = int(((draws >= observed) | np.isclose(draws, observed, atol=1e-12, rtol=0)).sum())
    return dict(observed=float(observed), null_mean=float(draws.mean()),
                null_p025=float(np.quantile(draws, .025)), null_p975=float(np.quantile(draws, .975)),
                null_sd=sd, excess=float(observed - draws.mean()), informative=bool(informative),
                p_upper=(1 + tail_count) / (len(draws) + 1) if informative else np.nan,
                tail_count=tail_count if informative else np.nan)


def describe(values, rng):
    values = np.asarray(values, float)
    boot = np.concatenate([values[rng.integers(0, len(values), size=(min(256, N_BOOT - start), len(values)))].mean(1)
                           for start in range(0, N_BOOT, 256)])
    quantiles = np.quantile(values, [.1, .25, .5, .75, .9])
    return dict(n=len(values), mean=float(values.mean()), q10=quantiles[0], q25=quantiles[1],
                median=quantiles[2], q75=quantiles[3], q90=quantiles[4],
                gene_bootstrap_mean_lo=float(np.quantile(boot, .025)),
                gene_bootstrap_mean_hi=float(np.quantile(boot, .975)))


def add_bh(frame):
    frame = frame.copy()
    frame["q_BH_per_threshold"] = np.nan
    frame["q_BH_joint_thresholds"] = np.nan
    frame["BH_family_size_per_threshold"] = 0
    frame["BH_family_size_joint_thresholds"] = 0
    for model, part in frame.groupby("null_model", sort=False):
        finite = part[part.p_upper.notna()]
        # Uninformative positions consume a planned family slot, without a p/q claim.
        joint = pd.Series(bh(part.p_upper.fillna(1).to_numpy()), index=part.index)
        frame.loc[finite.index, "q_BH_joint_thresholds"] = joint.loc[finite.index]
        frame.loc[part.index, "BH_family_size_joint_thresholds"] = len(part)
        for threshold, q in part.groupby("threshold", sort=False):
            adjusted = pd.Series(bh(q.p_upper.fillna(1).to_numpy()), index=q.index)
            keep = q[q.p_upper.notna()].index
            frame.loc[keep, "q_BH_per_threshold"] = adjusted.loc[keep]
            frame.loc[q.index, "BH_family_size_per_threshold"] = len(q)
    return frame


def main():
    for path in [HERE / "protocol.json", HERE / "analyze.py", ROOT / "analysis_code/pas_rerun/downstream_cohorts.py",
                 ROOT / "analysis_code/pas_rerun/downstream_candidate_screen.py"]:
        SOURCES[str(path.relative_to(ROOT))] = digest(path)
    matched_all = read(BASE / "matched_pms/matched_state_vectors.tsv.gz")
    members = read(ROOT / "notes/2026-10-10-directional-recurrence/observed_members.tsv.gz")
    covariates = read(BASE / "recurrence_controls/recurrence_frozen_covariates.tsv.gz")
    original = read(BASE / "recurrence/recurrence_members.tsv.gz")
    original = original[original.metric.eq("WPP") & original.variant.eq("primary") & original.series_matching]
    assert set(zip(members.species, members.contrast, members.threshold, members.gene_id, members.direction)) == set(zip(
        original.species, original.contrast, original.threshold, original.gene_id,
        original.direction.map({"negative": "proximal", "positive": "distal"})))
    assert not matched_all.duplicated(["species", "route", "setting", "gene_id"]).any()
    keys = ["species", "study_id", "route", "setting", "gene_id"]
    triplet_columns = keys + ["PAS_ids", "target_usage", "reference_usage", "M_usage", "M_PAS_ids", "PM_P_usage",
                              "MS_joint_PAS", "library_key_a", "library_key_b", "time_a", "time_b",
                              "sample_unit_a", "sample_unit_b", "triplet_weight"]
    triplets = read(BASE / "matched_pms/weighted_triplets.tsv.gz", usecols=triplet_columns)
    checks = validate_vectors(matched_all, triplets)
    del triplets
    matched = matched_all[matched_all.setting.eq(SETTING)].copy()
    matched["A"] = matched.DWUI_MS - matched.DWUI_PM
    dual = pd.concat([dual_reference(matched, members, t) for t in [.02, .05]], ignore_index=True)
    exported = read(BASE / "candidate_screen/dual_reference_candidates.tsv")
    assert set(zip(dual[dual.threshold.eq(.02)].species, dual[dual.threshold.eq(.02)].gene_id,
                   dual[dual.threshold.eq(.02)].direction)) == set(zip(exported.species, exported.gene_id,
                   exported.direction.map({"negative": "proximal", "positive": "distal"})))
    checks.append(dict(check="dual_reference_primary_exact_exported_roster", n=len(exported), max_abs_error=0))
    save(dual, "dual_reference_members.tsv")
    scopes = [("human", "GSE118258", {"GSE118258"}), ("mouse", "CIP", {"CIP"}),
              ("mouse", "WOT", {"WOT"}), ("mouse", "CIP_WOT_common", {"CIP", "WOT"})]
    coverages, effects, summary, tests, interactions, strata_rows, nulls, edge_rows = [], [], [], [], [], [], {}, []
    for species, scope, routes in scopes:
        q = matched[matched.species.eq(species) & matched.route.isin(routes)]
        complete = q.groupby("gene_id").route.nunique().eq(len(routes))
        q = q[q.gene_id.isin(complete.index[complete])]
        mean = q.groupby("gene_id").agg(gene_symbol=("gene_symbol", "first"),
            DWUI_PS=("DWUI_PS", "mean"), DWUI_PM=("DWUI_PM", "mean"), DWUI_MS=("DWUI_MS", "mean"),
            A=("A", "mean"), P_DWUI=("P_DWUI", "mean"), M_DWUI=("M_DWUI", "mean"), S_DWUI=("S_DWUI", "mean"),
            minimum_target_libraries=("n_target_libraries", "min"), minimum_target_dates=("n_target_dates", "min"),
            minimum_triplets=("n_triplets", "min"), PAS_number=("PAS_number", "mean")).reset_index()
        assert np.allclose(mean.DWUI_PS, mean.DWUI_PM + mean.DWUI_MS, rtol=0, atol=1e-12)
        for contrast in ["P_S", "P_M"]:
            for threshold in [.02, .05]:
                context = dict(species=species, scope=scope, selection_contrast=contrast, threshold=threshold)
                fixed = members[members.species.eq(species) & members.contrast.eq(contrast) & members.threshold.eq(threshold)]
                fixed = fixed[["gene_id", "direction"]]
                observed = mean.merge(fixed, on="gene_id", validate="one_to_one")
                observed["direction_sign"] = observed.direction.map({"proximal": -1, "distal": 1})
                observed["aligned_allocation"] = observed.direction_sign * observed.A
                observed["aligned_PS"] = observed.direction_sign * observed.DWUI_PS
                observed["aligned_PM"] = observed.direction_sign * observed.DWUI_PM
                observed["aligned_MS"] = observed.direction_sign * observed.DWUI_MS
                effects.append(observed.assign(**context))
                for direction, d in [("proximal", -1), ("distal", 1)]:
                    x = observed[observed.direction.eq(direction)]
                    coverage = context | dict(direction=direction, fixed_genes=int(fixed.direction.eq(direction).sum()),
                        matched_genes=len(x), matched_background_genes=len(mean),
                        PS_direction_retained=int((d * x.DWUI_PS).ge(threshold).sum()),
                        PM_direction_retained=int((d * x.DWUI_PM).ge(threshold).sum()),
                        MS_direction_retained=int((d * x.DWUI_MS).ge(threshold).sum()),
                        MS_dominates_aligned=int(x.aligned_allocation.gt(0).sum()),
                        PM_dominates_aligned=int(x.aligned_allocation.lt(0).sum()),
                        opposing_components=int((x.DWUI_MS * x.DWUI_PM).lt(0).sum()),
                        matched_same_direction_components=int((d*x.DWUI_MS).gt(0).mul((d*x.DWUI_PM).gt(0)).sum()))
                    coverages.append(coverage)
                    if len(x):
                        for column in ["DWUI_PS", "DWUI_PM", "DWUI_MS", "A", "aligned_PS", "aligned_PM", "aligned_MS", "aligned_allocation"]:
                            summary.append(context | dict(direction=direction, metric=column) |
                                           describe(x[column], np.random.default_rng(seed(*context.values(), direction, column, "bootstrap"))))
                if scope not in ["GSE118258", "CIP_WOT_common"]:
                    continue
                native = covariates[covariates.species.eq(species) & covariates.contrast.eq(contrast) & covariates.threshold.eq(threshold)].copy()
                native["log_GEX"] = np.log1p(native.mean_GEX)
                native["log_APA"] = np.log1p(native.APA_coverage)
                aggregate = native.groupby("gene_id").agg(eligible_studies=("unit", "nunique"),
                    native_mean_effect=("effect_mean", "mean"), native_mean_abs_effect=("effect_mean", lambda a: a.abs().mean()),
                    mean_log_GEX=("log_GEX", "mean"),
                    mean_log_APA=("log_APA", "mean"), native_PAS_number=("PAS_number", "mean"))
                signatures = native.groupby("gene_id").apply(
                    lambda a: ";".join(f"{r.unit}:{r.series_pattern}" for r in a.sort_values("unit").itertuples()), include_groups=False)
                aggregate["availability"] = signatures
                universe = mean.merge(aggregate[aggregate.eligible_studies.ge(2)], on="gene_id", validate="one_to_one")
                universe = universe.merge(fixed, on="gene_id", how="left", validate="one_to_one")
                universe["label"] = universe.direction.map({"proximal": -1, "distal": 1}).fillna(0).astype(int)
                universe["abs_PS"] = universe.DWUI_PS.abs()
                universe["native_sign"] = np.sign(universe.native_mean_effect).astype(int)
                assert set(universe[universe.label.ne(0)].gene_id) == set(observed.gene_id)
                values, labels = universe.A.to_numpy(), universe.label.to_numpy()
                for level, strength_bins, native_strength in [("median_strength", 2, False),
                        ("quintile_strength", 5, False), ("native_selection_strength", 2, True)]:
                    model, cells, edges = matching(universe, strength_bins, native_strength)
                    key = "|".join(map(str, [species, contrast, threshold, level]))
                    run_seed = seed(key, "labels")
                    draws, analytic = label_null(values, labels, cells, N_DRAW, np.random.default_rng(run_seed))
                    nulls[key] = draws
                    edge_rows.append(context | dict(null_model=level, edges=edges, matching_columns=[c for c in model if c.endswith("_bin")] + ["native_sign", "availability"]))
                    model["stratum"] = -1
                    for cell_id, ii in enumerate(cells):
                        model.loc[ii, "stratum"] = cell_id
                        strata_rows.append(context | dict(null_model=level, stratum=cell_id, genes=len(ii),
                            proximal=int((labels[ii] == -1).sum()), distal=int((labels[ii] == 1).sum()),
                            non_recurrent=int((labels[ii] == 0).sum()), mixed_labels=len(np.unique(labels[ii])) > 1))
                    save(model.assign(**context, null_model=level), f"background-{species}-{contrast}-{str(threshold).replace('.', 'p')}-{level}.tsv.gz")
                    for index, (direction, d) in enumerate([("proximal", -1), ("distal", 1)]):
                        estimate = float(values[labels == d].mean())
                        result = inference(estimate, draws[:, index])
                        if result["null_sd"] > 1e-14:
                            assert abs(result["null_mean"] - analytic[d]) < 6 * result["null_sd"] / np.sqrt(N_DRAW) + 1e-12
                        selected_movable = sum(int((labels[ii] == d).sum()) for ii in cells if len(np.unique(labels[ii])) > 1)
                        tests.append(context | dict(null_model=level, direction=direction, n=int((labels == d).sum()),
                            universe_genes=len(universe), strata=len(cells), singleton_strata=sum(len(ii) == 1 for ii in cells),
                            movable_selected_genes=selected_movable, movable_selected_fraction=selected_movable / (labels == d).sum(),
                            analytic_null_mean=analytic[d], aligned_observed=d*estimate, aligned_null_mean=d*analytic[d],
                            aligned_test_tail="upper" if d == 1 else "lower", permutations=N_DRAW, seed=run_seed) | result)
                    # Distal aligned minus proximal aligned = distal raw A + proximal raw A.
                    interaction_null = draws[:, 1] + draws[:, 0]
                    interaction_observed = float(values[labels == 1].mean() + values[labels == -1].mean())
                    interactions.append(context | dict(null_model=level, statistic="distal_aligned_minus_proximal_aligned",
                        universe_genes=len(universe), n_distal=int((labels == 1).sum()), n_proximal=int((labels == -1).sum()),
                        permutations=N_DRAW, seed=run_seed) | inference(interaction_observed, interaction_null))
    save(pd.concat(effects, ignore_index=True), "fixed_cohort_component_effects.tsv.gz")
    save(coverages, "coverage_and_component_counts.tsv")
    save(summary, "component_distributions.tsv")
    save(add_bh(pd.DataFrame(tests)), "matched_component_tests.tsv")
    save(add_bh(pd.DataFrame(interactions)), "matched_interaction_tests.tsv")
    save(strata_rows, "matched_null_stratum_margins.tsv.gz")
    save(checks, "verification.tsv")
    np.savez_compressed(HERE / "matched_component_null_draws.npz", **nulls)
    (HERE / "matching_edges.json").write_text(json.dumps(edge_rows, indent=2) + "\n")
    (HERE / "input_sha256.json").write_text(json.dumps(SOURCES, indent=2, sort_keys=True) + "\n")
    (HERE / "validation.json").write_text(json.dumps(dict(status="complete", vector_checks=len(checks),
        validation_rows=len(matched_all), primary_matched_rows=len(matched), null_contexts=len(nulls),
        permutations_per_context=N_DRAW, all_fixed_cohorts_preserved=True, manuscript_changed=False), indent=2) + "\n")
    print(json.dumps(dict(status="complete", primary_matched_rows=len(matched), null_contexts=len(nulls),
        dual_reference_primary=dual[dual.threshold.eq(.02)].groupby(["species", "direction"]).size().reset_index(name="n").to_dict("records"))))


if __name__ == "__main__":
    main()
