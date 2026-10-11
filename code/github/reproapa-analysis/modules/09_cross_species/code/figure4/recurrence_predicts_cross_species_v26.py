"""Independent within-species discovery followed by reciprocal validation."""
from pathlib import Path
import hashlib, json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

P = Path(__file__).resolve().parents[1]
D = P / "data/figure4-recurrence-validation-v26"
A = P / "analysis/figure4-recurrence-validation-v26"
D.mkdir(exist_ok=True)
A.mkdir(exist_ok=True)
NODES_PATH = P / "data/figure4-program-v9/all_mapped_nodes.tsv.gz"
PM_PATH = P / "data/figure4-program-v9/paired_pm_study_effects.tsv.gz"
nodes = pd.read_csv(NODES_PATH, sep="\t")
pm = pd.read_csv(PM_PATH, sep="\t")
codes_mouse = ["W", "C"]
codes_human = ["H1", "H2", "H3", "H4"]
taus = [0.01, 0.02, 0.03, 0.05]
rng = np.random.default_rng(260920)
B = 500


def build_effects(layer):
    q = nodes[nodes.layer.eq(layer) & nodes.measurement_pass & nodes.code.isin(codes_mouse + codes_human)].copy()
    w = q.pivot_table(index=["human_gene_id", "code"], columns="node", values="position_index", aggfunc="first")
    ps = (w.P - w.S).rename("effect").dropna().reset_index()
    # Measurement covariates are minimum endpoint support for P-S.
    m = q[q.node.isin(["S", "P"])].groupby(["human_gene_id", "code"]).agg(
        pas=("supported_pas_in_route", "min"), support=("qualified_metacells", "min"),
        libraries=("n_libraries", "min"), span=("supported_span", "min"),
    ).reset_index()
    ps = ps.merge(m, on=["human_gene_id", "code"], validate="one_to_one")
    ps["contrast"] = "P-S"

    z = pm[pm.layer.eq(layer) & pm.quality_pass & pm.code.isin(codes_mouse + codes_human)].copy()
    z = z[((z.code.isin(codes_mouse)) & z.n_sample_units.ge(2)) | ((z.code.isin(codes_human)) & z.n_sample_units.ge(1))]
    z = z.groupby(["human_gene_id", "code"], as_index=False).agg(
        effect=("mean_delta", "first"), pas=("supported_sites", "first"),
        support=("n_sample_units", "first"), libraries=("n_libraries", "first"),
        span=("supported_span", "first"),
    )
    z["contrast"] = "P-M"
    return pd.concat([ps, z], ignore_index=True)


effects = pd.concat([build_effects(x).assign(layer=x) for x in ["supported50", "expanded_primary"]], ignore_index=True)
effects.to_csv(D / "study_effects_and_covariates.tsv.gz", sep="\t", index=False)


def species_summary(frame, codes, prefix, tau):
    x = frame[frame.code.isin(codes)].copy()
    e = x.pivot(index="human_gene_id", columns="code", values="effect").reindex(columns=codes)
    out = pd.DataFrame(index=e.index)
    out[f"{prefix}_eligible"] = e.notna().sum(axis=1)
    out[f"{prefix}_meta"] = e.median(axis=1, skipna=True)
    out[f"{prefix}_direction"] = np.sign(out[f"{prefix}_meta"]).replace(0, np.nan)
    signed = e.mul(out[f"{prefix}_direction"], axis=0)
    out[f"{prefix}_supporting"] = signed.ge(tau).sum(axis=1)
    out[f"{prefix}_opposing"] = signed.le(-tau).sum(axis=1)
    out[f"{prefix}_score"] = out[f"{prefix}_supporting"] / out[f"{prefix}_eligible"]
    out[f"{prefix}_recurrent"] = out[f"{prefix}_eligible"].ge(2) & out[f"{prefix}_supporting"].ge(2) & out[f"{prefix}_opposing"].eq(0)
    out[f"{prefix}_active"] = out[f"{prefix}_eligible"].ge(2) & e.abs().ge(tau).any(axis=1)
    feat = x.groupby("human_gene_id").agg(
        pas=("pas", "median"), support=("support", "median"),
        libraries=("libraries", "median"), span=("span", "median"),
    )
    feat = feat.add_prefix(prefix + "_")
    return out.join(feat)


def validation_table(layer, contrast, discovery, tau):
    frame = effects[effects.layer.eq(layer) & effects.contrast.eq(contrast)]
    if discovery == "human":
        ds = species_summary(frame, codes_human, "disc", tau)
        vs = species_summary(frame, codes_mouse, "val", tau)
    else:
        ds = species_summary(frame, codes_mouse, "disc", tau)
        vs = species_summary(frame, codes_human, "val", tau)
    x = ds.join(vs, how="inner")
    x = x[x.disc_active & x.val_eligible.ge(2) & x.disc_direction.notna()].copy()
    x["signed_validation"] = x.val_meta * x.disc_direction
    x["validation_same_direction"] = x.signed_validation.ge(tau)
    # Strict outcome asks whether the validation species itself meets the same
    # recurrence rule in the discovery direction, without redefining direction.
    v = frame[frame.code.isin(codes_mouse if discovery == "human" else codes_human)].pivot(
        index="human_gene_id", columns="code", values="effect"
    )
    signed = v.mul(x.disc_direction, axis=0)
    x["validation_strict_recurrence"] = signed.ge(tau).sum(axis=1).ge(2) & signed.le(-tau).sum(axis=1).eq(0)
    x["layer"] = layer
    x["contrast"] = contrast
    x["discovery"] = discovery
    x["tau"] = tau
    return x.reset_index()


def add_bins(x):
    x = x.copy()
    # Effect magnitude is the primary matching variable. Quantile bins plus
    # exact direction/eligible counts keep the comparator close without using
    # validation outcomes.
    x["effect_bin"] = pd.qcut(x.disc_meta.abs(), q=min(10, x.disc_meta.abs().nunique()), duplicates="drop", labels=False)
    x["stratum"] = x[["disc_direction", "disc_eligible", "val_eligible", "effect_bin"]].astype(str).agg("|".join, axis=1)
    return x


def match_sets(x):
    x = add_bins(x)
    recurrent = x[x.disc_recurrent]
    controls = x[x.disc_active & ~x.disc_recurrent]
    technical_covs = ["disc_pas", "disc_support", "disc_libraries", "disc_span"]
    log_covs = np.log1p(x[technical_covs].clip(lower=0))
    z_covs = (log_covs - log_covs.mean()) / log_covs.std().replace(0, 1)
    rows = []
    for ridx, r in recurrent.iterrows():
        pool = controls[
            controls.disc_direction.eq(r.disc_direction)
            & controls.disc_eligible.eq(r.disc_eligible)
            & controls.val_eligible.eq(r.val_eligible)
        ]
        target = abs(r.disc_meta)
        caliper = max(0.005, 0.25 * target)
        lower = pool[pool.disc_meta.abs().le(target) & (target - pool.disc_meta.abs()).le(caliper)].copy()
        upper = pool[pool.disc_meta.abs().ge(target) & (pool.disc_meta.abs() - target).le(caliper)].copy()
        if lower.empty or upper.empty:
            continue
        lower = lower.loc[(lower.disc_meta.abs() - target).abs().nsmallest(10).index]
        upper = upper.loc[(upper.disc_meta.abs() - target).abs().nsmallest(10).index]
        best = None
        for li, lo_candidate in lower.iterrows():
            for hi_i, hi_candidate in upper.iterrows():
                lo_abs, hi_abs = abs(lo_candidate.disc_meta), abs(hi_candidate.disc_meta)
                if np.isclose(lo_abs, hi_abs):
                    lo_weight, hi_weight = 1.0, 0.0
                else:
                    lo_weight = (hi_abs - target) / (hi_abs - lo_abs)
                    hi_weight = (target - lo_abs) / (hi_abs - lo_abs)
                weighted = z_covs.loc[li] * lo_weight + z_covs.loc[hi_i] * hi_weight
                distance = float(((weighted - z_covs.loc[ridx]) ** 2).sum())
                if best is None or distance < best[0]:
                    best = (distance, lo_candidate, hi_candidate, lo_weight, hi_weight)
        _, lo, hi, lo_weight, hi_weight = best
        rows.append({"match_id": int(ridx), "role": "recurrent", "match_weight": 1.0, **r.to_dict()})
        if lo.name == hi.name or np.isclose(abs(lo.disc_meta), abs(hi.disc_meta)):
            rows.append({"match_id": int(ridx), "role": "control", "match_weight": 1.0, **lo.to_dict()})
        else:
            rows.append({"match_id": int(ridx), "role": "control", "match_weight": lo_weight, **lo.to_dict()})
            rows.append({"match_id": int(ridx), "role": "control", "match_weight": hi_weight, **hi.to_dict()})
    return pd.DataFrame(rows)


def matched_stats(m):
    sets = []
    for mid, q in m.groupby("match_id"):
        case = q[q.role.eq("recurrent")].iloc[0]
        ctrl = q[q.role.eq("control")]
        weights = ctrl.match_weight.to_numpy()
        sets.append({
            "match_id": mid,
            "case_binary": float(case.validation_same_direction),
            "control_binary": np.average(ctrl.validation_same_direction, weights=weights),
            "case_strict": float(case.validation_strict_recurrence),
            "control_strict": np.average(ctrl.validation_strict_recurrence, weights=weights),
            "case_continuous": case.signed_validation,
            "control_continuous": np.average(ctrl.signed_validation, weights=weights),
        })
    s = pd.DataFrame(sets)
    def odds(a, b):
        return ((a + 0.5) / (len(s) - a + 0.5)) / ((b + 0.5) / (len(s) - b + 0.5))
    def boot(metric_case, metric_control, fn, n=2000):
        vals = []
        for _ in range(n):
            q = s.iloc[rng.integers(0, len(s), len(s))]
            vals.append(fn(q[metric_case], q[metric_control]))
        return np.quantile(vals, [0.025, 0.975])
    a = s.case_binary.sum(); b = s.control_binary.sum(); binary_or = odds(a, b)
    strict_a = s.case_strict.sum(); strict_b = s.control_strict.sum(); strict_or = odds(strict_a, strict_b)
    binary_ci = boot("case_binary", "control_binary", lambda x, y: odds(x.sum(), y.sum()))
    strict_ci = boot("case_strict", "control_strict", lambda x, y: odds(x.sum(), y.sum()))
    delta = (s.case_continuous - s.control_continuous).mean()
    continuous_ci = boot("case_continuous", "control_continuous", lambda x, y: (x - y).mean())
    signs = rng.choice([-1, 1], size=(B, len(s)))
    p_cont = (1 + (np.abs((signs * (s.case_continuous - s.control_continuous).to_numpy()).mean(axis=1)) >= abs(delta)).sum()) / (B + 1)
    return s, {
        "n_recurrent_matched": len(s), "mean_controls_per_set": len(m[m.role.eq("control")]) / len(s),
        "case_binary_rate": s.case_binary.mean(), "control_binary_rate": s.control_binary.mean(),
        "binary_OR": binary_or, "binary_OR_lo": binary_ci[0], "binary_OR_hi": binary_ci[1],
        "case_strict_rate": s.case_strict.mean(), "control_strict_rate": s.control_strict.mean(),
        "strict_OR": strict_or, "strict_OR_lo": strict_ci[0], "strict_OR_hi": strict_ci[1],
        "case_signed_validation": s.case_continuous.mean(), "control_signed_validation": s.control_continuous.mean(),
        "continuous_difference": delta, "continuous_lo": continuous_ci[0], "continuous_hi": continuous_ci[1],
        "continuous_signflip_p": p_cont,
    }


all_tables = []
match_outputs = []
result_rows = []
balance_rows = []
dose_rows = []
for layer in ["supported50", "expanded_primary"]:
    for contrast in ["P-S", "P-M"]:
        for discovery in ["human", "mouse"]:
            for tau in taus:
                x = validation_table(layer, contrast, discovery, tau)
                all_tables.append(x)
                m = match_sets(x)
                if m.empty or m.match_id.nunique() < 5:
                    continue
                sets, stats = matched_stats(m)
                result_rows.append({"layer": layer, "contrast": contrast, "discovery": discovery, "tau": tau, **stats})
                m["analysis"] = f"{layer}|{contrast}|{discovery}|{tau}"
                match_outputs.append(m)
                # Balance uses recurrent rows and controls weighted equally per set.
                for cov in ["disc_meta", "disc_pas", "disc_support", "disc_libraries", "disc_span"]:
                    case = m[m.role.eq("recurrent")].drop_duplicates("match_id")[cov].abs() if cov == "disc_meta" else m[m.role.eq("recurrent")].drop_duplicates("match_id")[cov]
                    ctrl = m[m.role.eq("control")].groupby("match_id").apply(
                        lambda q: np.average(q[cov].abs() if cov == "disc_meta" else q[cov], weights=q.match_weight), include_groups=False
                    )
                    pooled = np.sqrt((case.var() + ctrl.var()) / 2)
                    balance_rows.append({"layer": layer, "contrast": contrast, "discovery": discovery, "tau": tau, "covariate": cov, "standardized_difference": (case.mean() - ctrl.mean()) / pooled if pooled else 0})
                # Dose response is descriptive and uses all APA-active genes.
                x["score_group"] = pd.cut(x.disc_score, [-0.001, 0.49, 0.99, 1.001], labels=["low", "medium", "high"])
                for group, q in x.groupby("score_group", observed=True):
                    dose_rows.append({"layer": layer, "contrast": contrast, "discovery": discovery, "tau": tau, "score_group": group, "n": len(q), "validation_rate": q.validation_same_direction.mean(), "mean_signed_validation": q.signed_validation.mean()})

all_gene = pd.concat(all_tables, ignore_index=True)
all_gene.to_csv(D / "reciprocal_gene_universe.tsv.gz", sep="\t", index=False)
pd.concat(match_outputs, ignore_index=True).to_csv(D / "matched_comparators.tsv.gz", sep="\t", index=False)
results = pd.DataFrame(result_rows)
primary = results[results.layer.eq("supported50") & results.tau.eq(0.02)].copy()
primary["continuous_q_four"] = false_discovery_control(primary.continuous_signflip_p)
primary.to_csv(D / "primary_reciprocal_results.tsv", sep="\t", index=False)
results.to_csv(D / "threshold_layer_sensitivity.tsv", sep="\t", index=False)
pd.DataFrame(balance_rows).to_csv(D / "matching_balance.tsv", sep="\t", index=False)
pd.DataFrame(dose_rows).to_csv(D / "recurrence_score_dose_response.tsv", sep="\t", index=False)

# Dose-response test adjusted by the same discovery-only strata used for the
# orthology permutation. Validation outcomes never define the strata.
dose_tests = []
for (layer, contrast, discovery, tau_value), universe in all_gene.groupby(["layer", "contrast", "discovery", "tau"]):
    if layer != "supported50" or tau_value != 0.02:
        continue
    universe = add_bins(universe).reset_index(drop=True)
    groups = [np.array(idx, dtype=int) for idx in universe.groupby("stratum").groups.values()]
    score = universe.disc_score.to_numpy()
    for outcome in ["signed_validation", "validation_same_direction"]:
        values = universe[outcome].astype(float).to_numpy()
        observed = pd.Series(score).corr(pd.Series(values), method="spearman")
        null = np.empty(B)
        for b in range(B):
            shuffled = values.copy()
            for idx in groups:
                shuffled[idx] = values[rng.permutation(idx)]
            null[b] = pd.Series(score).corr(pd.Series(shuffled), method="spearman")
        dose_tests.append({
            "layer": layer, "contrast": contrast, "discovery": discovery, "tau": tau_value,
            "outcome": outcome, "n": len(universe), "spearman": observed,
            "null_lo": np.quantile(null, 0.025), "null_hi": np.quantile(null, 0.975),
            "p_greater": (1 + (null >= observed).sum()) / (B + 1),
        })
dose_test = pd.DataFrame(dose_tests)
dose_test["q_eight"] = false_discovery_control(dose_test.p_greater)
dose_test.to_csv(D / "recurrence_score_trend_tests.tsv", sep="\t", index=False)

# Leave-one-study sensitivity for the discovery or validation side. These are
# raw recurrent-set summaries; matched inference remains the primary analysis.
leave_rows = []
for contrast in ["P-S", "P-M"]:
    frame = effects[effects.layer.eq("supported50") & effects.contrast.eq(contrast)]
    # Human discovery after dropping one human study; validate in both mice.
    for omitted in codes_human:
        disc_codes = [c for c in codes_human if c != omitted]
        ds = species_summary(frame, disc_codes, "disc", 0.02)
        vs = species_summary(frame, codes_mouse, "val", 0.02)
        x = ds.join(vs, how="inner")
        x = x[x.disc_recurrent & x.val_eligible.ge(2) & x.disc_direction.notna()].copy()
        signed = x.val_meta * x.disc_direction
        leave_rows.append({"contrast": contrast, "direction": "human_to_mouse", "side": "discovery", "omitted": omitted, "n_recurrent_testable": len(x), "validation_rate": signed.ge(0.02).mean(), "mean_signed_validation": signed.mean()})
    # Full human discovery; validate separately in W and C.
    ds = species_summary(frame, codes_human, "disc", 0.02)
    for kept in codes_mouse:
        v = frame[frame.code.eq(kept)].set_index("human_gene_id").effect
        x = ds.join(v.rename("validation"), how="inner")
        x = x[x.disc_recurrent & x.disc_direction.notna() & x.validation.notna()].copy()
        signed = x.validation * x.disc_direction
        leave_rows.append({"contrast": contrast, "direction": "human_to_mouse", "side": "validation", "omitted": ",".join(c for c in codes_mouse if c != kept), "n_recurrent_testable": len(x), "validation_rate": signed.ge(0.02).mean(), "mean_signed_validation": signed.mean()})
    # Full mouse discovery; validate after dropping each human study.
    ds = species_summary(frame, codes_mouse, "disc", 0.02)
    for omitted in codes_human:
        val_codes = [c for c in codes_human if c != omitted]
        vs = species_summary(frame, val_codes, "val", 0.02)
        x = ds.join(vs, how="inner")
        x = x[x.disc_recurrent & x.val_eligible.ge(2) & x.disc_direction.notna()].copy()
        signed = x.val_meta * x.disc_direction
        leave_rows.append({"contrast": contrast, "direction": "mouse_to_human", "side": "validation", "omitted": omitted, "n_recurrent_testable": len(x), "validation_rate": signed.ge(0.02).mean(), "mean_signed_validation": signed.mean()})
    leave_rows.append({"contrast": contrast, "direction": "mouse_to_human", "side": "discovery", "omitted": "W or C", "n_recurrent_testable": np.nan, "validation_rate": np.nan, "mean_signed_validation": np.nan, "note": "not estimable: one remaining mouse study cannot define cross-study recurrence"})
pd.DataFrame(leave_rows).to_csv(D / "leave_one_study_sensitivity.tsv", sep="\t", index=False)

# Orthology x recurrence interaction: shuffle validation outcomes within the
# exact matching strata, keeping discovery recurrence and matching fixed.
interaction = []
for key, m in pd.concat(match_outputs, ignore_index=True).groupby("analysis"):
    layer, contrast, discovery, tau_s = key.split("|")
    if layer != "supported50" or float(tau_s) != 0.02:
        continue
    observed = []
    for _, q in m.groupby("match_id"):
        ctrl = q[q.role.eq("control")]
        observed.append(q[q.role.eq("recurrent")].signed_validation.iloc[0] - np.average(ctrl.signed_validation, weights=ctrl.match_weight))
    observed = np.mean(observed)
    universe = all_gene[(all_gene.layer.eq(layer)) & (all_gene.contrast.eq(contrast)) & (all_gene.discovery.eq(discovery)) & all_gene.tau.eq(float(tau_s))]
    universe = add_bins(universe).reset_index(drop=True)
    gid_to_index = dict(zip(universe.human_gene_id, universe.index))
    set_indices = []
    for _, q in m.groupby("match_id"):
        case_index = gid_to_index[q[q.role.eq("recurrent")].human_gene_id.iloc[0]]
        ctrl = q[q.role.eq("control")]
        control_indices = np.array([gid_to_index[g] for g in ctrl.human_gene_id], dtype=int)
        set_indices.append((case_index, control_indices, ctrl.match_weight.to_numpy()))
    groups = [np.array(idx, dtype=int) for idx in universe.groupby("stratum").groups.values()]
    base = universe.signed_validation.to_numpy()
    null = np.empty(B)
    for b in range(B):
        shuffled = base.copy()
        for idx in groups:
            shuffled[idx] = base[rng.permutation(idx)]
        null[b] = np.mean([shuffled[case] - np.average(shuffled[controls], weights=weights) for case, controls, weights in set_indices])
    interaction.append({"layer": layer, "contrast": contrast, "discovery": discovery, "tau": float(tau_s), "observed_recurrence_orthology_interaction": observed, "null_lo": np.quantile(null, 0.025), "null_hi": np.quantile(null, 0.975), "p_greater": (1 + (null >= observed).sum()) / (B + 1)})
interaction = pd.DataFrame(interaction)
interaction["q_four"] = false_discovery_control(interaction.p_greater)
interaction.to_csv(D / "orthology_recurrence_interaction.tsv", sep="\t", index=False)

manifest = {
    "discovery": "within one species only; >=2 same-direction strong studies, no strong opposite study",
    "validation": "other species, >=2 eligible studies; binary median effect aligned >= tau; strict >=2 aligned and no opposite",
    "matching": "common-support interpolation: exact direction and eligible-study counts; lower/upper controls weighted to exactly match discovery absolute effect; pair chosen to minimize PAS/support/library/span distance",
    "expression_matching": "not available on a harmonized all-ortholog background; omitted and disclosed",
    "permutations": B,
    "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in [NODES_PATH, PM_PATH]},
}
(A / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(primary.to_string(index=False))
print(pd.DataFrame(interaction).to_string(index=False))
