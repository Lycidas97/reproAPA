"""Outcome-blind disjoint matching and conditional set-level PAS validation."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment


HERE = Path(__file__).resolve().parent
OUT = HERE / "matched_results"
COVARIATES = ["mean_log2_last_exon_CPM", "log2_min_shared_count",
              "log2_external_last_exon_sites", "min_shared_fraction",
              "log2_shared_span_bp", "MEF_shared_axis_usage"]
SAMPLES = ["MEF_R1", "MEF_R2", "Untransfected_D3_R1", "Untransfected_D3_R2",
           "Untransfected_D6_R1", "Untransfected_D6_R2", "iPSC_R1", "iPSC_R2"]
TARGETS = {"iPSC": ["iPSC_R1", "iPSC_R2"], "D3": SAMPLES[2:4], "D6": SAMPLES[4:6]}


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


def save(frame, name):
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.12g")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build_measurements():
    genes = read(HERE / "results/all_gene_feasibility.tsv").set_index("gene_id")
    sites = read(HERE / "results/site_mapping.tsv")
    sites = sites[sites.usable].copy()
    raw = read(HERE / "inputs/GSE99922_readcounts-table.txt.gz", usecols=["COORDS", *SAMPLES])
    raw = raw.drop_duplicates("COORDS")
    assert raw.COORDS.is_unique
    sites = sites.merge(raw, on="COORDS", validate="one_to_one")
    sites["signed_cut"] = sites.mapped_cut * np.where(sites.ref_strand.eq("+"), 1, -1)
    sites["local_cut_bp"] = sites.signed_cut - sites.groupby("gene_id").signed_cut.transform("min")
    for sample in SAMPLES:
        numerator = (sites.local_cut_bp * sites[sample]).groupby(sites.gene_id).sum()
        denominator = sites.groupby("gene_id")[sample].sum()
        genes[f"local_mean_bp__{sample}"] = numerator / denominator.replace(0, np.nan)
    for target, samples in TARGETS.items():
        delta = genes[[f"local_mean_bp__{s}" for s in samples]].mean(axis=1) - genes[[f"local_mean_bp__{s}" for s in SAMPLES[:2]]].mean(axis=1)
        old = genes[f"delta_bp__{target}"]
        both = delta.notna() & old.notna()
        assert np.allclose(delta[both], old[both], rtol=1e-8, atol=1e-6)
        genes[f"precise_delta_bp__{target}"] = delta
    depths = read(HERE / "results/library_qc.tsv").set_index("sample").count_sum_all_unique_sites
    endpoint = SAMPLES[:2] + TARGETS["iPSC"]
    cpm = pd.DataFrame({s: genes[f"last_exon_count__{s}"] / depths[s] * 1e6 for s in endpoint})
    genes["mean_log2_last_exon_CPM"] = np.log2(cpm + .5).mean(axis=1)
    genes["log2_min_shared_count"] = np.log2(genes.min_shared_count__iPSC.clip(lower=1))
    genes["log2_external_last_exon_sites"] = np.log2(genes.n_external_last_exon_sites)
    genes["min_shared_fraction"] = genes.min_shared_fraction__iPSC
    genes["log2_shared_span_bp"] = np.log2(genes.shared_span_bp.clip(lower=1))
    genes["MEF_shared_axis_usage"] = genes[[f"local_mean_bp__{s}" for s in SAMPLES[:2]]].mean(axis=1) / genes.shared_span_bp.replace(0, np.nan)
    genes["pas_stratum"] = np.where(genes.n_shared_sites.eq(2), "2", np.where(genes.n_shared_sites.le(4), "3-4", "5+"))
    eligible = genes.measurable20__iPSC
    assert genes.loc[eligible, COVARIATES].notna().all().all()
    assert genes.loc[eligible, "MEF_shared_axis_usage"].between(0, 1).all()
    return genes


def match_cases(case_covariates, control_covariates, k, protocol, covariates=None):
    """This function accepts only IDs, matching covariates and PAS strata."""
    covariates = COVARIATES if covariates is None else covariates
    allowed_columns = covariates + ["pas_stratum"]
    cases = case_covariates[allowed_columns].sort_index()
    controls = control_covariates[allowed_columns].sort_index()
    pooled = pd.concat([cases[covariates], controls[covariates]])
    scales = pooled.std(ddof=0).replace(0, 1)
    a, b = cases[covariates].to_numpy(), controls[covariates].to_numpy()
    difference = a[:, None, :] - b[None, :, :]
    standardized = difference / scales.to_numpy()[None, None, :]
    limits = protocol["calipers"]
    valid = (cases.pas_stratum.to_numpy()[:, None] == controls.pas_stratum.to_numpy()[None, :])
    for index, name in enumerate(covariates):
        key = name + ("_sd" if index < 3 else "_absolute")
        current = standardized[:, :, index] if index < 3 else difference[:, :, index]
        valid &= np.abs(current) <= limits[key] + 1e-12
    cost = np.square(standardized).sum(axis=2)
    rows = np.repeat(np.arange(len(cases)), k)
    constrained = np.where(valid[rows], cost[rows], 1e6)
    augmented = np.concatenate([constrained, np.full((len(rows), len(rows)), 1e4)], axis=1)
    row_indices, assignment = linear_sum_assignment(augmented)
    assert np.array_equal(row_indices, np.arange(len(rows)))
    assignment = assignment.reshape(len(cases), k)
    complete = (assignment < len(controls)).all(axis=1)
    matched = []
    for i in np.flatnonzero(complete):
        assert valid[i, assignment[i]].all()
        matched.append([cases.index[i], *controls.index[assignment[i]].tolist()])
    used = [gene for group in matched for gene in group[1:]]
    assert len(used) == len(set(used)), "Controls must be disjoint across matched sets"
    balance = []
    if matched:
        selected_cases = cases.loc[[group[0] for group in matched], covariates]
        selected_controls = controls.loc[used, covariates]
        for name in covariates:
            balance.append({"covariate": name, "case_mean": selected_cases[name].mean(),
                "control_mean": selected_controls[name].mean(),
                "smd_before": (selected_cases[name].mean() - controls[name].mean()) / scales[name],
                "smd_after": (selected_cases[name].mean() - selected_controls[name].mean()) / scales[name]})
    status = pd.DataFrame({"gene_id": cases.index, "caliper_qualified_controls": valid.sum(axis=1),
                           "complete_disjoint_match": complete})
    return matched, pd.DataFrame(balance), status


def residuals(aligned):
    k = aligned.shape[-1] - 1
    signs = (aligned > 0).astype(float)
    direction = signs[..., 0] - signs[..., 1:].mean(axis=-1)
    physical = aligned[..., 0] - aligned[..., 1:].mean(axis=-1)
    return direction, physical


def evaluate(effects, directions, rng, permutations):
    aligned = effects * directions[:, None]
    excess, physical = residuals(aligned)
    observed = np.array([excess.mean(), np.median(physical)])
    n, members = aligned.shape
    selections = rng.integers(0, members, size=(permutations, n))
    chosen = aligned[np.arange(n)[None, :], selections]
    signs = (aligned > 0).astype(float)
    chosen_sign = signs[np.arange(n)[None, :], selections]
    permuted_excess = (members * chosen_sign - signs.sum(axis=1)[None, :]) / (members - 1)
    permuted_physical = (members * chosen - aligned.sum(axis=1)[None, :]) / (members - 1)
    null = np.column_stack([permuted_excess.mean(axis=1), np.median(permuted_physical, axis=1)])
    p = (1 + (null >= observed[None, :] - 1e-12).sum(axis=0)) / (permutations + 1)
    return observed, p, null


def bootstrap_intervals(genes, groups, directions, rng, draws):
    n = len(groups)
    effects = genes.loc[np.array(groups).ravel(), "precise_delta_bp__iPSC"].to_numpy().reshape(n, -1)
    aligned = effects * directions[:, None]
    excess, physical = residuals(aligned)
    sampled = rng.integers(0, n, size=(draws, n))
    conditional = np.column_stack([excess[sampled].mean(axis=1), np.median(physical[sampled], axis=1)])
    # A whole-library draw is shared across all genes and controls.
    target_weight = rng.binomial(2, .5, size=draws) / 2
    source_weight = rng.binomial(2, .5, size=draws) / 2
    values = {s: genes.loc[np.array(groups).ravel(), f"local_mean_bp__{s}"].to_numpy().reshape(n, -1) for s in SAMPLES[:2] + TARGETS["iPSC"]}
    effect_draws = (target_weight[:, None, None] * values["iPSC_R1"][None, :, :]
        + (1 - target_weight[:, None, None]) * values["iPSC_R2"][None, :, :]
        - source_weight[:, None, None] * values["MEF_R1"][None, :, :]
        - (1 - source_weight[:, None, None]) * values["MEF_R2"][None, :, :])
    joint_excess, joint_physical = residuals(effect_draws * directions[None, :, None])
    row = np.arange(draws)[:, None]
    joint = np.column_stack([joint_excess[row, sampled].mean(axis=1), np.median(joint_physical[row, sampled], axis=1)])
    return np.quantile(conditional, [.025, .975], axis=0), np.quantile(joint, [.025, .975], axis=0)


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(adjusted, 1)
    return result


def run():
    OUT.mkdir(exist_ok=True)
    protocol_path = HERE / "matched-protocol.json"
    protocol = json.loads(protocol_path.read_text())
    frozen = read(HERE / "inputs/frozen_candidates.tsv")
    assert frozen.groupby("set_name").size().to_dict() == {"primary_fixed130": 130, "secondary_PS488": 488}
    genes = build_measurements()
    save(genes.reset_index(), "measurement_covariates.tsv")
    excluded = set(frozen.gene_id)
    configs = [("primary", 20, 3), ("count50", 50, 3), ("count100", 100, 3), ("controls2", 20, 2)]
    summaries, matches, balances, statuses, crosspairs, timepoints = [], [], [], [], [], []
    archives = {}
    for config_index, (configuration, threshold, k) in enumerate(configs):
        eligible = genes[genes[f"measurable{threshold}__iPSC"]].copy()
        controls = eligible[~eligible.index.isin(excluded)]
        for set_index, (name, roster) in enumerate(frozen.groupby("set_name", sort=True)):
            cases = eligible[eligible.index.isin(roster.gene_id)]
            groups, balance, status = match_cases(cases[COVARIATES + ["pas_stratum"]], controls[COVARIATES + ["pas_stratum"]], k, protocol)
            balance["set_name"], balance["configuration"] = name, configuration
            balances.append(balance)
            complete_roster = roster.merge(status, on="gene_id", how="left", validate="one_to_one")
            complete_roster["configuration"] = configuration
            complete_roster["matching_status"] = np.where(complete_roster.complete_disjoint_match.fillna(False), "matched", np.where(complete_roster.gene_id.isin(cases.index), "insufficient_disjoint_caliper_controls", "not_measurable_under_frozen_pilot_rules"))
            statuses.append(complete_roster)
            # Persist outcome-blind matching status before calculating effect statistics.
            save(pd.concat(statuses), "candidate_matching_status.tsv")
            direction_lookup = roster.set_index("gene_id").direction_sign
            directions = direction_lookup.loc[[group[0] for group in groups]].to_numpy()
            for group in groups:
                for role, gene in enumerate(group):
                    matches.append({"configuration": configuration, "set_name": name, "candidate_id": group[0], "role": "candidate" if role == 0 else "control", "gene_id": gene, "candidate_direction": int(direction_lookup[group[0]])})
            retained = len(groups) / len(cases) if len(cases) else 0
            max_smd = balance.smd_after.abs().max() if len(balance) else np.inf
            gate = retained >= protocol["design_quality_gate"]["minimum_fraction_of_measurable_candidates_retained"] and max_smd <= protocol["design_quality_gate"]["maximum_absolute_standardized_mean_difference"]
            assert len(groups) >= 3, "Insufficient matched sets; no significance result is defined"
            flat = np.array(groups).ravel()
            effects = genes.loc[flat, "precise_delta_bp__iPSC"].to_numpy().reshape(len(groups), k + 1)
            rng = np.random.default_rng(protocol["seed"] + config_index * 100 + set_index)
            observed, p, null = evaluate(effects, directions, rng, protocol["permutations"])
            conditional_ci, joint_ci = bootstrap_intervals(genes, groups, directions, rng, protocol["permutations"])
            aligned = effects * directions[:, None]
            for i, statistic in enumerate(["excess_direction_concordance", "median_aligned_excess_bp"]):
                summaries.append({"configuration": configuration, "set_name": name, "statistic": statistic,
                    "frozen_n": len(roster), "measurable_n": len(cases), "matched_n": len(groups), "controls_n": len(groups) * k,
                    "fraction_measurable_retained": retained, "max_abs_smd": max_smd, "design_quality_gate_passed": bool(gate),
                    "estimate": observed[i], "p_one_sided": p[i],
                    "conditional_ci_low": conditional_ci[0, i], "conditional_ci_high": conditional_ci[1, i],
                    "gene_and_library_ci_low": joint_ci[0, i], "gene_and_library_ci_high": joint_ci[1, i],
                    "candidate_concordance": (aligned[:, 0] > 0).mean(), "matched_control_concordance": (aligned[:, 1:] > 0).mean(),
                    "null_low": np.quantile(null[:, i], .025), "null_high": np.quantile(null[:, i], .975)})
            archives[f"{configuration}__{name}"] = null
            archives[f"{configuration}__{name}__observed"] = observed
            if configuration == "primary":
                for source in SAMPLES[:2]:
                    for target in TARGETS["iPSC"]:
                        delta = (genes.loc[flat, f"local_mean_bp__{target}"] - genes.loc[flat, f"local_mean_bp__{source}"]).to_numpy().reshape(len(groups), k + 1)
                        exc, bp = residuals(delta * directions[:, None])
                        crosspairs.append({"set_name": name, "MEF_library": source, "iPSC_library": target,
                            "matched_n": len(groups), "excess_direction_concordance": exc.mean(), "median_aligned_excess_bp": np.median(bp)})
                common = [j for j, group in enumerate(groups) if genes.loc[group, "measurable20__D3"].all() and genes.loc[group, "measurable20__D6"].all()]
                for target in ["iPSC", "D3", "D6"]:
                    if common:
                        common_flat = np.array(groups)[common].ravel()
                        delta = genes.loc[common_flat, f"precise_delta_bp__{target}"].to_numpy().reshape(len(common), k + 1)
                        align = delta * directions[common, None]
                        exc, bp = residuals(align)
                        timepoints.append({"set_name": name, "target": target, "complete_matched_n": len(common),
                            "candidate_concordance": (align[:, 0] > 0).mean(), "matched_control_concordance": (align[:, 1:] > 0).mean(),
                            "excess_direction_concordance": exc.mean(), "median_aligned_excess_bp": np.median(bp)})
    result = pd.DataFrame(summaries)
    primary = result.configuration.eq("primary")
    assert primary.sum() == 4
    result.loc[primary, "q_primary_four_test_family"] = bh(result.loc[primary, "p_one_sided"])
    save(result, "matched_summary.tsv")
    save(pd.DataFrame(matches), "matched_members.tsv")
    save(pd.concat(balances), "matching_balance.tsv")
    save(pd.DataFrame(crosspairs), "single_library_sensitivity.tsv")
    save(pd.DataFrame(timepoints), "common_set_timepoints.tsv")
    np.savez_compressed(OUT / "permutation_draws.npz", **archives)
    inputs = [protocol_path, HERE / "inputs/frozen_candidates.tsv", HERE / "inputs/GSE99922_readcounts-table.txt.gz",
              HERE / "results/all_gene_feasibility.tsv", HERE / "results/site_mapping.tsv", HERE / "results/library_qc.tsv"]
    receipt = {"status": "completed", "protocol_sha256": digest(protocol_path), "script_sha256": digest(Path(__file__)),
        "input_sha256": {str(p.relative_to(HERE)): digest(p) for p in inputs},
        "primary_test_family_size": int(primary.sum()), "eligible_background_genes_count20": int((genes.measurable20__iPSC & ~genes.index.isin(excluded)).sum()),
        "inference_scope": "Conditional observational matched-set test, limited to measurable matched frozen members; baseline discovery members and directions unchanged; no per-gene FDR",
        "uncertainty_scope": "Descriptive bootstrap intervals; two endpoint libraries constrain biological generalization; crosspairs are sensitivity measurements, not independent replicates"}
    (OUT / "run.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(result.loc[primary].to_string(index=False), flush=True)


if __name__ == "__main__":
    run()
