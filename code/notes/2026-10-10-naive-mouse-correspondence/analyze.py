"""Conditional ortholog overlap between fixed human hN and mouse P responses."""
import hashlib
import json
import math
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
TABLES = ROOT / "notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables"
SETTING = "supported50_p08"
B = 20000
INPUTS = {}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(relative, **kwargs):
    path = TABLES / relative
    INPUTS[str(path.relative_to(ROOT))] = digest(path)
    return pd.read_csv(path, sep="\t", float_precision="round_trip", **kwargs)


def save(frame, name):
    frame.to_csv(HERE / name, sep="\t", index=False, na_rep="NA",
                 float_format="%.17g",
                 compression={"method": "gzip", "mtime": 0} if name.endswith(".gz") else None)


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    out = np.empty(len(p))
    out[order] = np.minimum(1, np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return out


def seed(*parts):
    return 20261010 + int(hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()[:12], 16)


def bins(values, n):
    out = pd.Series(-1, index=values.index, dtype=int)
    valid = values.notna()
    if valid.any():
        out.loc[valid] = pd.qcut(values[valid], n, labels=False, duplicates="drop").fillna(0).astype(int)
    return out


def cells(frame, columns):
    return list(frame.groupby(columns, dropna=False, sort=True).indices.values())


def membership_probabilities(labels, strata):
    result = np.zeros(len(labels))
    for ii in strata:
        result[ii] = labels[ii].mean()
    return result


def binary_draw(labels, strata, n, random):
    draw = np.broadcast_to(labels, (n, len(labels))).copy()
    for ii in strata:
        count = int(labels[ii].sum())
        if 0 < count < len(ii):
            draw[:, ii] = random.permuted(np.broadcast_to(labels[ii], (n, len(ii))), axis=1)
        assert np.all(draw[:, ii].sum(1) == count)
    return draw


def simulate(human, mouse, human_strata, mouse_strata, n, random):
    draws = []
    for start in range(0, n, 250):
        size = min(250, n - start)
        h = binary_draw(human, human_strata, size, random)
        m = binary_draw(mouse, mouse_strata, size, random)
        draws.append((h & m).sum(1))
    return np.concatenate(draws)


def hypergeom_tail(population, human, mouse, observed):
    """One-sided hypergeometric tail, for the unmatched descriptive reference."""
    def logchoose(n, k):
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)
    low = max(0, human + mouse - population)
    high = min(human, mouse)
    if observed <= low:
        return 1.
    return min(1., math.fsum(math.exp(logchoose(human, k) +
                   logchoose(population - human, mouse - k) - logchoose(population, mouse))
                   for k in range(observed, high + 1)))


def human_frames(pairs, axes):
    pairs = pairs[pairs.version.eq("v2_primary") & pairs.setting.eq(SETTING)].copy()
    a = axes[axes.study_id.eq("GSE150311") & axes.setting.eq(SETTING)]
    assert a.atlas_pas_id.is_unique
    weight = a.set_index("atlas_pas_id").WPP
    for column in ["reference", "target", "axis_min", "axis_max", "supported_baseline"]:
        pairs[column] = np.nan
    for index, row in pairs.iterrows():
        ids = json.loads(row.PAS_ids)
        w = weight.loc[ids].to_numpy()
        reference = np.asarray(json.loads(row.usage_b)) @ w
        target = np.asarray(json.loads(row.usage_a)) @ w
        assert abs(target - reference - row.delta_DWUI) < 1e-12
        baseline = (reference - w.min()) / (w.max() - w.min())
        pairs.loc[index, ["reference", "target", "axis_min", "axis_max", "supported_baseline"]] = reference, target, w.min(), w.max(), baseline
    pairs["coverage"] = (pairs.pas_counts_a / pairs.n_mc_a + pairs.pas_counts_b / pairs.n_mc_b) / 2
    pairs["capacity"] = pairs.reference - pairs.axis_min
    pairs["min_metacells"] = pairs[["n_mc_a", "n_mc_b"]].min(axis=1)
    q = pairs[pairs.node_a.eq("N") & pairs.node_b.eq("P")].copy()
    q["a"] = q.library_a.str[-1]
    q["b"] = q.library_b.str[-1]
    pooled = q[q.a.isin(["6", "8"]) & q.b.eq("7")].copy()
    counts = pooled.groupby("gene_id").size()
    pooled = pooled[pooled.gene_id.isin(counts[counts.eq(2)].index)]
    assert not pooled.duplicated(["gene_id", "a"]).any()
    wide = pooled.pivot(index="gene_id", columns="a", values="delta_DWUI")
    common = pooled.groupby("gene_id").agg(
        h_gene_symbol=("gene_symbol", "first"), h_PAS=("PAS_number", "first"),
        h_baseline=("supported_baseline", "mean"), h_capacity=("capacity", "min"),
        h_coverage=("coverage", "mean"), h_min_metacells=("min_metacells", "min"),
        h_tandem=("tandem_utr_compatible", "first"))
    assert pooled.groupby("gene_id").PAS_number.nunique().eq(1).all()
    assert pooled.groupby("gene_id").reference.agg(lambda x: x.max() - x.min()).lt(1e-12).all()
    common["h_naive_delta"] = wide["6"]
    common["h_RSeT_delta"] = wide["8"]
    frames = {}
    for label, effect in [("pooled_two_hN_vs_primed", wide.mean(1)),
                          ("both_hN_conditions_proximal", wide.max(1)),
                          ("naive_condition_vs_primed", wide["6"]),
                          ("RSeT_condition_vs_primed", wide["8"])]:
        f = common.copy()
        f["h_delta"] = effect
        # Single-condition matching uses its actual measured coverage and metacells.
        if label in ["naive_condition_vs_primed", "RSeT_condition_vs_primed"]:
            k = "6" if label.startswith("naive_") else "8"
            f["h_coverage"] = pooled[pooled.a.eq(k)].set_index("gene_id").coverage
            f["h_min_metacells"] = pooled[pooled.a.eq(k)].set_index("gene_id").min_metacells
        frames[label] = f.reset_index().rename(columns={"gene_id": "human_gene_id"})
    same = q[q.a.eq("8") & q.b.eq("8")].copy()
    assert same.gene_id.is_unique
    frames["RSeT_same_library_N_vs_P"] = same.rename(columns={
        "gene_id": "human_gene_id", "gene_symbol": "h_gene_symbol", "PAS_number": "h_PAS",
        "supported_baseline": "h_baseline", "capacity": "h_capacity", "coverage": "h_coverage",
        "min_metacells": "h_min_metacells", "tandem_utr_compatible": "h_tandem",
        "delta_DWUI": "h_delta"})[["human_gene_id", "h_gene_symbol", "h_PAS", "h_baseline",
            "h_capacity", "h_coverage", "h_min_metacells", "h_tandem", "h_delta"]]
    return frames, pairs


def complete_human_axes(pairs, axes, eligibility):
    """Recover translation-invariant human weights from signed cuts and full span.

    Native absolute PAS origin is not recoverable for all naive-only genes.
    The test therefore matches supported-axis-relative baseline consistently.
    Every stored delta is subsequently verified against these weights.
    """
    q = pairs[pairs.version.eq("v2_primary") & pairs.setting.eq(SETTING)]
    spans = eligibility[eligibility.study_id.eq("GSE150311")].set_index("gene_id").native_axis_span_nt
    original = axes[axes.study_id.eq("GSE150311") & axes.setting.eq(SETTING)].set_index("atlas_pas_id").WPP
    out = []
    for gene, rows in q.groupby("gene_id"):
        lists = rows.PAS_ids.map(json.loads)
        ids = lists.iloc[0]
        assert all(value == ids for value in lists)
        signed = []
        for pas in ids:
            match = re.fullmatch(r"polyApipe:GSE150311:(chr[^_]+)_(\d+)_([fr])", pas)
            assert match, pas
            signed.append(int(match[2]) * (1 if match[3] == "f" else -1))
        signed = np.asarray(signed, float)
        weights = (signed - signed.min()) / spans.loc[gene]
        assert np.ptp(signed) <= spans.loc[gene] and spans.loc[gene] > 0
        if all(p in original.index for p in ids):
            old = original.loc[ids].to_numpy()
            assert np.max(np.abs((old - old.min()) - weights)) < 1e-12
        out.extend(dict(atlas_pas_id=pas, species="human", study_id="GSE150311",
                        setting=SETTING, WPP=float(weight)) for pas, weight in zip(ids, weights))
    return pd.concat([axes[~axes.study_id.eq("GSE150311")], pd.DataFrame(out)], ignore_index=True)


def mouse_frames(frozen, recurrence, vectors, axes, threshold):
    f = frozen[frozen.species.eq("mouse") & frozen.contrast.eq("P_M") &
               frozen.threshold.eq(threshold)].copy()
    assert set(f.unit) == {"GSE103221", "GSE115943"}
    measured = f.groupby("gene_id").unit.nunique()
    f = f[f.gene_id.isin(measured[measured.eq(2)].index)]
    assert not f.duplicated(["gene_id", "unit"]).any()
    votes = f.pivot(index="gene_id", columns="unit", values="event")
    events = votes.eq(-1).all(1)
    expected = recurrence[recurrence.species.eq("mouse") & recurrence.contrast.eq("P_M") &
                          recurrence.threshold.eq(threshold) & recurrence.direction.eq("negative")]
    assert set(events[events].index) == set(expected.gene_id)
    cov = f.groupby("gene_id").agg(
        m_gene_symbol=("gene_symbol", "first"), m_GEX=("mean_GEX", "median"),
        m_coverage=("APA_coverage", "median"), m_PAS=("PAS_number", "median"),
        m_baseline=("baseline", "median"), m_capacity=("capacity_negative", "min"))
    cov["m_selected"] = events
    out = {"P_M_proximal_recurrent": cov.reset_index().rename(columns={"gene_id": "mouse_gene_id"})}
    v = vectors[vectors.species.eq("mouse") & vectors.setting.eq(SETTING) &
                vectors.route.isin(["CIP", "WOT"])].copy()
    complete = v.groupby("gene_id").route.nunique()
    v = v[v.gene_id.isin(complete[complete.eq(2)].index) & v.gene_id.isin(cov.index)]
    assert not v.duplicated(["gene_id", "route"]).any()
    weights = axes[axes.species.eq("mouse") & axes.setting.eq(SETTING)].set_index("atlas_pas_id").WPP
    assert weights.index.is_unique
    v["axis_min"] = v.PAS_ids.map(lambda x: float(weights.loc[json.loads(x)].min()))
    v["dual_baseline"] = (v.M_DWUI + v.S_DWUI) / 2
    v["dual_capacity"] = v[["M_DWUI", "S_DWUI"]].min(axis=1) - v.axis_min
    passes = v.DWUI_PM.le(-threshold) & v.DWUI_PS.le(-threshold)
    v["dual_pass"] = passes
    dual = cov.loc[v.gene_id.unique()].copy()
    dual["m_baseline"] = v.groupby("gene_id").dual_baseline.median()
    dual["m_capacity"] = v.groupby("gene_id").dual_capacity.min()
    dual["m_selected"] &= v.groupby("gene_id").dual_pass.all()
    if threshold == .02:
        assert int(cov.m_selected.sum()) == 372
        assert int(dual.m_selected.sum()) == 207
    out["P_M_PS_same_contributor_dual_proximal"] = dual.reset_index().rename(columns={"gene_id": "mouse_gene_id"})
    return out


def make_strata(frame, side, threshold, model):
    frame = frame.copy()
    frame[side + "_capacity_ok"] = frame[side + "_capacity"].ge(threshold)
    columns = [side + "_capacity_ok"]
    for column, n in [("PAS", 3), ("coverage", 3), ("baseline", 3 if model == "coarse_assay" else 5)]:
        name = side + "_" + column + "_bin"
        frame[name] = bins(frame[side + "_" + column], n)
        columns.append(name)
    if side == "m" and model == "fine_assay_mouse_expression":
        frame["m_GEX_bin"] = bins(frame.m_GEX, 3)
        columns.append("m_GEX_bin")
    strata = cells(frame, columns)
    label = frame[side + "_selected"].to_numpy(bool)
    assert not (label & ~frame[side + "_capacity_ok"].to_numpy()).any()
    probabilities = membership_probabilities(label, strata)
    assert probabilities[~frame[side + "_capacity_ok"].to_numpy()].sum() == 0
    return frame, strata, probabilities


def main():
    INPUTS[str((HERE / "protocol.json").relative_to(ROOT))] = digest(HERE / "protocol.json")
    orth = read("cross_overlap/ortholog_matched_features.tsv.gz")
    assert len(orth) == 16092 and orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    assert orth.human_partner_count.eq(1).all() and orth.mouse_partner_count.eq(1).all()
    orth = orth[["human_gene_id", "mouse_gene_id"]]
    pairs = read("naive_state/library_pair_effects.tsv.gz")
    axes = read("state_vectors/PAS_axes.primary.tsv.gz", usecols=["atlas_pas_id", "species", "study_id", "setting", "WPP"])
    eligibility = read("tandem/tandem_utr_eligibility.tsv.gz")
    axes = complete_human_axes(pairs, axes, eligibility)
    frozen = read("recurrence_controls/recurrence_frozen_covariates.tsv.gz")
    recurrence = read("recurrence/recurrence_members.tsv.gz")
    recurrence = recurrence[recurrence.metric.eq("WPP") & recurrence.series_matching & recurrence.variant.eq("primary")]
    vectors = read("matched_pms/matched_state_vectors.tsv.gz")
    human, reconstructed = human_frames(pairs, axes)
    contrast = read("naive_state/contrast_gene_effects.tsv.gz")
    q = contrast[contrast.version.eq("v2_primary") & contrast.setting.eq(SETTING) &
                 contrast.contrast.eq("hN_P_primed_condition") & contrast.complete_contributors]
    h = human["pooled_two_hN_vs_primed"].set_index("human_gene_id").h_delta
    assert set(h.index) == set(q.gene_id)
    assert np.max(np.abs(h.loc[q.gene_id].to_numpy() - q.delta_DWUI.to_numpy())) < 1e-12
    rows, memberships, inventories, diagnostics, balances, draws = [], [], [], [], [], {}
    for threshold in [.02, .05]:
        mouse = mouse_frames(frozen, recurrence, vectors, axes, threshold)
        for human_label, hf in human.items():
            hf = hf.assign(h_selected=hf.h_delta.le(-threshold))
            for target, mf in mouse.items():
                f = orth.merge(hf, on="human_gene_id", validate="one_to_one").merge(
                    mf, on="mouse_gene_id", validate="one_to_one").sort_values("human_gene_id").reset_index(drop=True)
                context = dict(human_label=human_label, mouse_target=target, threshold=threshold)
                x, y = f.h_selected.to_numpy(bool), f.m_selected.to_numpy(bool)
                obs = int((x & y).sum())
                n, hn, mn = len(f), int(x.sum()), int(y.sum())
                unmatched = hn * mn / n
                b, c, d = hn - obs, mn - obs, n - hn - mn + obs
                odds = obs * d / (b * c) if b * c else np.inf
                inv = context | dict(human_measurable=len(hf), human_full_selected=int(hf.h_selected.sum()),
                    human_mapped_measurable=int(hf.human_gene_id.isin(orth.human_gene_id).sum()),
                    mouse_measurable=len(mf), mouse_full_selected=int(mf.m_selected.sum()),
                    mouse_mapped_measurable=int(mf.mouse_gene_id.isin(orth.mouse_gene_id).sum()),
                    paired_measurable_orthologs=n, human_selected_in_pair_background=hn,
                    mouse_selected_in_pair_background=mn, observed_overlap=obs,
                    unmatched_expected=unmatched, raw_odds_ratio=odds,
                    unmatched_hypergeometric_p=hypergeom_tail(n, hn, mn, obs))
                inventories.append(inv)
                memberships.append(f.assign(**context))
                for model in ["coarse_assay", "fine_assay_mouse_expression"]:
                    _, hs, hp = make_strata(f, "h", threshold, model)
                    _, ms, mp = make_strata(f, "m", threshold, model)
                    run_seed = seed("naive_mouse_overlap", human_label, target, threshold, model)
                    null = simulate(x, y, hs, ms, B, np.random.default_rng(run_seed))
                    analytic = float(hp @ mp)
                    null_sd = float(null.std(ddof=1))
                    assert abs(null.mean() - analytic) < 6 * null_sd / np.sqrt(B) + 1e-12
                    key = "|".join(map(str, (human_label, target, threshold, model)))
                    draws[key] = null
                    row = context | dict(model=model, n_orthologs=n, human_selected=hn, mouse_selected=mn,
                        observed=obs, analytic_matched_expectation=analytic, null_mean=float(null.mean()),
                        null_sd=null_sd, null_p025=float(np.quantile(null, .025)), null_p975=float(np.quantile(null, .975)),
                        fold_vs_matched_expectation=obs / analytic if analytic else np.nan,
                        excess_vs_matched=obs - analytic, p=(1 + (null >= obs).sum()) / (B + 1),
                        tail_count=int((null >= obs).sum()), permutations=B, seed=run_seed,
                        unmatched_expected=unmatched, raw_odds_ratio=odds,
                        unmatched_hypergeometric_p=inv["unmatched_hypergeometric_p"])
                    for side, labels, strata, probabilities in [("human", x, hs, hp), ("mouse", y, ms, mp)]:
                        row[side + "_strata"] = len(strata)
                        row[side + "_movable_fraction"] = float(((probabilities > 0) & (probabilities < 1)).mean())
                        row[side + "_expected_changed_fraction"] = float((2 * probabilities * (1 - probabilities)).mean())
                        row[side + "_selected_movable_fraction"] = float(((probabilities > 0) & (probabilities < 1))[labels].mean()) if labels.any() else 0.
                        row[side + "_fixed_selected"] = int((probabilities == 1).sum())
                        for j, ii in enumerate(strata):
                            diagnostics.append(context | dict(model=model, side=side, stratum=j,
                                genes=len(ii), selected=int(labels[ii].sum()), probability=float(probabilities[ii[0]])))
                        prefix = "h" if side == "human" else "m"
                        for column in ["baseline", "capacity", "PAS", "coverage"] + (["GEX"] if side == "mouse" else ["min_metacells"]):
                            values = f[prefix + "_" + column].to_numpy(float)
                            valid = np.isfinite(values)
                            selected = labels & valid
                            total = probabilities[valid].sum()
                            balances.append(context | dict(model=model, side=side, covariate=column,
                                observed_selected_mean=float(values[selected].mean()) if selected.any() else np.nan,
                                null_selected_mean=float(probabilities[valid] @ values[valid] / total) if total else np.nan,
                                observed_selected_missing=int((labels & ~valid).sum()),
                                expected_selected_missing=float(probabilities[~valid].sum())))
                    rows.append(row)
                    print(key, "N", n, "selected", hn, mn, "overlap", obs,
                          "expected", round(analytic, 4), "p", row["p"], flush=True)
    result = pd.DataFrame(rows)
    result["primary_family"] = result.threshold.eq(.02) & result.human_label.eq("pooled_two_hN_vs_primed")
    result["q_primary_or_sensitivity_within_model"] = result.groupby(
        ["model", "primary_family"]).p.transform(bh)
    result["q_all20_within_model"] = result.groupby("model").p.transform(bh)
    result["q_all40_exploratory"] = bh(result.p)
    assert result.primary_family.sum() == 4 and len(result) == 40
    save(result, "overlap_results.tsv")
    save(pd.DataFrame(inventories), "measurability_and_unmatched_reference.tsv")
    save(pd.concat(memberships, ignore_index=True), "ortholog_membership_and_covariates.tsv.gz")
    save(pd.DataFrame(diagnostics), "stratum_margins.tsv.gz")
    save(pd.DataFrame(balances), "covariate_balance.tsv")
    np.savez_compressed(HERE / "overlap_null_draws.npz", **draws)
    assert all(digest(ROOT / path) == value for path, value in INPUTS.items())
    outputs = ["overlap_results.tsv", "measurability_and_unmatched_reference.tsv",
               "ortholog_membership_and_covariates.tsv.gz", "stratum_margins.tsv.gz",
               "covariate_balance.tsv", "overlap_null_draws.npz"]
    manifest = dict(status="PASS", input_sha256=INPUTS, script_sha256=digest(Path(__file__)),
        output_sha256={name: digest(HERE / name) for name in outputs},
        python_executable=sys.executable, python=sys.version, numpy=np.__version__, pandas=pd.__version__,
        human_pair_rows_with_exact_DWUI_reconstruction=len(reconstructed),
        pooled_contrast_replay=True, mouse_primary_roster_replay=True, mouse_207_dual_reference_replay=True,
        ortholog_universe=16092, exact_permutation_margins_asserted=True,
        analytic_matched_expectation_MC_verified=True,
        limitation="Human N GEX and absolute full-native PAS origin unavailable: human baseline is relative to supported PAS, proximal capacity is attainable movement on frozen supported axis. Mouse dual-target expression/coverage are original recurrence contributor proxies; libraries fixed, not replicate inference.")
    (HERE / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
