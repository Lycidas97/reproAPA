"""Audit Figure 4C uncertainty with unique orthologs as resampling units.

This reads the frozen v26 matches. It never rematches or edits the source results.
The bootstrap is conditional on the selected matching graph; each reused control
receives one draw, shared by every set in which that gene occurs.
"""

from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/figure4-recurrence-validation-v26/matched_comparators.tsv.gz"
OLD = ROOT / "data/figure4-recurrence-validation-v26/primary_reciprocal_results.tsv"
OLD_SENSITIVITY = ROOT / "data/figure4-recurrence-validation-v26/threshold_layer_sensitivity.tsv"
OUT = ROOT / "analysis/figure4c-reused-controls-v53"
OUT.mkdir(exist_ok=True)
B_PRIMARY = 20000
B_SENSITIVITY = 5000
SEED = 20260923


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    order = np.argsort(p)
    ranks = np.arange(1, len(p) + 1)
    q_sorted = np.minimum.accumulate((p[order] * len(p) / ranks)[::-1])[::-1]
    q = np.empty_like(p)
    q[order] = np.minimum(q_sorted, 1.0)
    return q


def odds(case_rate, control_rate, n_sets):
    a = n_sets * case_rate
    b = n_sets * control_rate
    return ((a + 0.5) / (n_sets - a + 0.5)) / ((b + 0.5) / (n_sets - b + 0.5))


def analysis(q, seed, draws):
    rng = np.random.default_rng(seed)
    cols = ["human_gene_id", "role", "signed_validation", "validation_same_direction",
            "validation_strict_recurrence", "stratum", "disc_direction",
            "disc_eligible", "val_eligible"]
    u = q[cols].drop_duplicates("human_gene_id").reset_index(drop=True)
    assert q.groupby("human_gene_id")["role"].nunique().max() == 1
    for col in ["signed_validation", "validation_same_direction", "validation_strict_recurrence", "stratum"]:
        assert q.groupby("human_gene_id")[col].nunique().max() == 1
    assert np.allclose(q.groupby("match_id").apply(
        lambda x: x.loc[x.role.eq("control"), "match_weight"].sum(), include_groups=False), 1.0)
    n = q.match_id.nunique()
    case_ids = q.loc[q.role.eq("recurrent"), ["match_id", "human_gene_id"]].drop_duplicates()
    assert len(case_ids) == n and case_ids.human_gene_id.is_unique
    control_edges = q[q.role.eq("control")]
    control_weight = control_edges.groupby("human_gene_id").match_weight.sum()
    assert np.isclose(control_weight.sum(), n)
    case = u[u.role.eq("recurrent")].set_index("human_gene_id")
    ctrl = u[u.role.eq("control")].set_index("human_gene_id").loc[control_weight.index]
    weights = control_weight.to_numpy(dtype=float)

    yc = case.signed_validation.to_numpy(dtype=float)
    yt = ctrl.signed_validation.to_numpy(dtype=float)
    case_rate = case.validation_same_direction.to_numpy(dtype=float)
    ctrl_rate = ctrl.validation_same_direction.to_numpy(dtype=float)
    case_strict = case.validation_strict_recurrence.to_numpy(dtype=float)
    ctrl_strict = ctrl.validation_strict_recurrence.to_numpy(dtype=float)
    observed = yc.mean() - np.average(yt, weights=weights)

    # Sample each distinct case/control gene once per draw. The aggregate
    # control weight carries every reuse of that gene in the fixed match graph.
    draw_case = rng.integers(0, len(case), size=(draws, len(case)))
    draw_ctrl = rng.integers(0, len(ctrl), size=(draws, len(ctrl)))
    sampled_weights = weights[draw_ctrl]
    denominator = sampled_weights.sum(axis=1)
    case_cont_boot = yc[draw_case].mean(axis=1)
    ctrl_cont_boot = (yt[draw_ctrl] * sampled_weights).sum(axis=1) / denominator
    delta_boot = case_cont_boot - ctrl_cont_boot
    case_bin_boot = case_rate[draw_case].mean(axis=1)
    ctrl_bin_boot = (ctrl_rate[draw_ctrl] * sampled_weights).sum(axis=1) / denominator
    case_strict_boot = case_strict[draw_case].mean(axis=1)
    ctrl_strict_boot = (ctrl_strict[draw_ctrl] * sampled_weights).sum(axis=1) / denominator

    # Preserve the exact matching variables when resampling unique genes.
    # This avoids adding variation from shifts in the direction/study-count mix.
    exact_cols = ["disc_direction", "disc_eligible", "val_eligible"]
    for frame in (case, ctrl):
        frame["exact_stratum"] = frame[exact_cols].astype(str).agg("|".join, axis=1)
    strat_case_cont = np.zeros(draws)
    strat_ctrl_cont = np.zeros(draws)
    strat_case_bin = np.zeros(draws)
    strat_ctrl_bin = np.zeros(draws)
    strat_case_strict = np.zeros(draws)
    strat_ctrl_strict = np.zeros(draws)
    for stratum, cases in case.groupby("exact_stratum"):
        controls = ctrl[ctrl.exact_stratum.eq(stratum)]
        assert not controls.empty
        w = control_weight.loc[controls.index].to_numpy(dtype=float)
        factor = len(cases) / n
        assert np.isclose(w.sum(), len(cases))
        ic = rng.integers(0, len(cases), size=(draws, len(cases)))
        it = rng.integers(0, len(controls), size=(draws, len(controls)))
        wb = w[it]
        den = wb.sum(axis=1)
        strat_case_cont += factor * cases.signed_validation.to_numpy(dtype=float)[ic].mean(axis=1)
        strat_ctrl_cont += factor * (controls.signed_validation.to_numpy(dtype=float)[it] * wb).sum(axis=1) / den
        strat_case_bin += factor * cases.validation_same_direction.to_numpy(dtype=float)[ic].mean(axis=1)
        strat_ctrl_bin += factor * (controls.validation_same_direction.to_numpy(dtype=float)[it] * wb).sum(axis=1) / den
        strat_case_strict += factor * cases.validation_strict_recurrence.to_numpy(dtype=float)[ic].mean(axis=1)
        strat_ctrl_strict += factor * (controls.validation_strict_recurrence.to_numpy(dtype=float)[it] * wb).sum(axis=1) / den
    strat_delta_boot = strat_case_cont - strat_ctrl_cont

    # Under conditional exchangeability, keep the matching graph and its
    # repeated-control links intact and permute one outcome per unique gene.
    gene_weight = pd.Series(0.0, index=u.human_gene_id)
    gene_weight.loc[case.index] = 1.0 / n
    gene_weight.loc[ctrl.index] = -weights / n
    coefficients = u.human_gene_id.map(gene_weight).to_numpy(dtype=float)
    outcomes = u.signed_validation.to_numpy(dtype=float)
    assert np.isclose(coefficients @ outcomes, observed)
    nulls = {}
    for mode, key_cols in {
        "exact": ["disc_direction", "disc_eligible", "val_eligible"],
        "effect_bin": ["stratum"],
    }.items():
        groups = [np.asarray(list(indices), dtype=int) for indices in u.groupby(key_cols).indices.values()]
        null = np.empty(draws)
        for b in range(draws):
            permuted = outcomes.copy()
            for indices in groups:
                permuted[indices] = outcomes[rng.permutation(indices)]
            null[b] = coefficients @ permuted
        center = float(null.mean())
        p = (1 + np.count_nonzero(np.abs(null - center) >= abs(observed - center))) / (draws + 1)
        group_keys = u[key_cols].astype(str).agg("|".join, axis=1)
        role_counts = pd.crosstab(group_keys, u.role)
        mixed = role_counts[(role_counts.get("recurrent", 0) > 0) & (role_counts.get("control", 0) > 0)]
        nulls[mode] = {
            f"{mode}_null_mean": center,
            f"{mode}_null_lo": float(np.quantile(null, 0.025)),
            f"{mode}_null_hi": float(np.quantile(null, 0.975)),
            f"{mode}_permutation_p_two_sided": float(p),
            f"{mode}_strata": len(groups),
            f"{mode}_swappable_genes": int(sum(len(g) for g in groups if len(g) > 1)),
            f"{mode}_mixed_role_strata": len(mixed),
            f"{mode}_mixed_role_genes": int(mixed.sum().sum()),
        }
    reuse = control_edges.groupby("human_gene_id").match_id.nunique()
    result = {
        "n_sets": n,
        "n_unique_cases": len(case),
        "n_unique_controls": len(ctrl),
        "control_edge_count": int(reuse.sum()),
        "n_reused_controls": int((reuse > 1).sum()),
        "max_control_reuse": int(reuse.max()),
        "control_effective_n_by_weights": float(weights.sum() ** 2 / np.sum(weights ** 2)),
        "observed_continuous_difference": float(observed),
        "gene_bootstrap_continuous_lo": float(np.quantile(delta_boot, 0.025)),
        "gene_bootstrap_continuous_hi": float(np.quantile(delta_boot, 0.975)),
        "stratified_gene_bootstrap_continuous_lo": float(np.quantile(strat_delta_boot, 0.025)),
        "stratified_gene_bootstrap_continuous_hi": float(np.quantile(strat_delta_boot, 0.975)),
        "observed_binary_or": float(odds(case_rate.mean(), np.average(ctrl_rate, weights=weights), n)),
        "gene_bootstrap_binary_or_lo": float(np.quantile(odds(case_bin_boot, ctrl_bin_boot, n), 0.025)),
        "gene_bootstrap_binary_or_hi": float(np.quantile(odds(case_bin_boot, ctrl_bin_boot, n), 0.975)),
        "stratified_gene_bootstrap_binary_or_lo": float(np.quantile(odds(strat_case_bin, strat_ctrl_bin, n), 0.025)),
        "stratified_gene_bootstrap_binary_or_hi": float(np.quantile(odds(strat_case_bin, strat_ctrl_bin, n), 0.975)),
        "observed_strict_or": float(odds(case_strict.mean(), np.average(ctrl_strict, weights=weights), n)),
        "gene_bootstrap_strict_or_lo": float(np.quantile(odds(case_strict_boot, ctrl_strict_boot, n), 0.025)),
        "gene_bootstrap_strict_or_hi": float(np.quantile(odds(case_strict_boot, ctrl_strict_boot, n), 0.975)),
        "stratified_gene_bootstrap_strict_or_lo": float(np.quantile(odds(strat_case_strict, strat_ctrl_strict, n), 0.025)),
        "stratified_gene_bootstrap_strict_or_hi": float(np.quantile(odds(strat_case_strict, strat_ctrl_strict, n), 0.975)),
        **nulls["exact"], **nulls["effect_bin"],
    }
    return result


matched = pd.read_csv(SOURCE, sep="\t")
old = pd.read_csv(OLD, sep="\t")
selected = matched[matched.layer.eq("supported50") & matched.tau.eq(0.02)]
rows = []
for i, ((contrast, discovery), q) in enumerate(selected.groupby(["contrast", "discovery"])):
    rows.append({"layer": "supported50", "contrast": contrast, "discovery": discovery,
                 "tau": 0.02, **analysis(q, SEED + i, B_PRIMARY)})
result = pd.DataFrame(rows).merge(old, on=["layer", "contrast", "discovery", "tau"], validate="one_to_one")
assert np.allclose(result.observed_continuous_difference, result.continuous_difference)
assert np.allclose(result.observed_binary_or, result.binary_OR)
assert np.allclose(result.observed_strict_or, result.strict_OR)
for mode in ["exact", "effect_bin"]:
    result[f"{mode}_q_four"] = bh_adjust(result[f"{mode}_permutation_p_two_sided"])
result.to_csv(OUT / "primary_unique_gene_audit.tsv", sep="\t", index=False)

old_sensitivity = pd.read_csv(OLD_SENSITIVITY, sep="\t")
primary_keys = set(zip(result.layer, result.contrast, result.discovery, result.tau))
sensitivity_rows = []
for i, (key, q) in enumerate(matched.groupby(["layer", "contrast", "discovery", "tau"])):
    if key in primary_keys:
        continue
    if q.match_id.nunique() < 5:
        continue
    sensitivity_rows.append(dict(zip(["layer", "contrast", "discovery", "tau"], key)) |
                            analysis(q, SEED + 100 + i, B_SENSITIVITY))
all_audit = pd.concat([result, pd.DataFrame(sensitivity_rows)], ignore_index=True)
all_audit = all_audit.merge(old_sensitivity[["layer", "contrast", "discovery", "tau",
                                             "continuous_difference", "continuous_signflip_p"]],
                            on=["layer", "contrast", "discovery", "tau"],
                            how="inner", validate="one_to_one", suffixes=("", "_sensitivity"))
assert len(all_audit) == len(old_sensitivity) == 30
assert np.allclose(all_audit.observed_continuous_difference, all_audit.continuous_difference_sensitivity)
for mode in ["exact", "effect_bin"]:
    all_audit[f"{mode}_q_thirty"] = bh_adjust(all_audit[f"{mode}_permutation_p_two_sided"])
all_audit.to_csv(OUT / "thirty_setting_unique_gene_audit.tsv", sep="\t", index=False)
manifest = {
    "source_match_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    "source_primary_sha256": hashlib.sha256(OLD.read_bytes()).hexdigest(),
    "source_sensitivity_sha256": hashlib.sha256(OLD_SENSITIVITY.read_bytes()).hexdigest(),
    "seed": SEED, "primary_draws_per_setting": B_PRIMARY,
    "sensitivity_draws_per_setting": B_SENSITIVITY,
    "primary_universe": "supported50; absolute Delta WPP >= 0.02; four reciprocal analyses",
    "bootstrap": "separate unique-gene resamples in recurrent and control pools, fixed match graph; primary audit stratifies by exact discovery direction and both eligible-study counts; a reused control has one aggregate weight",
    "permutation": "shuffle each unique gene's aligned validation effect once within discovery-side strata; fixed graph; two-sided deviation from empirical null mean",
    "scope": "conditional on frozen matching; gene-level sampling uncertainty, not biological library/study replication",
}
(OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(result[["contrast", "discovery", "n_sets", "n_unique_controls", "n_reused_controls",
              "control_effective_n_by_weights", "continuous_lo", "continuous_hi",
              "stratified_gene_bootstrap_continuous_lo", "stratified_gene_bootstrap_continuous_hi",
              "continuous_signflip_p", "exact_permutation_p_two_sided", "exact_q_four",
              "effect_bin_permutation_p_two_sided", "effect_bin_q_four"]].to_string(index=False))
print("Sensitivity family:", len(all_audit), "settings; minimum exact q =",
      all_audit.exact_q_thirty.min(), "; minimum effect-bin q =",
      all_audit.effect_bin_q_thirty.min())
