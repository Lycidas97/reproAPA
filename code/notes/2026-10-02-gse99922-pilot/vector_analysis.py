"""Unsigned redistribution and signed support with whole-library diagnostics."""
from itertools import combinations, product
from pathlib import Path
import json

import numpy as np
import pandas as pd

from expanded_measurements import HERE, OUT, ENDPOINT, profile_metrics, save
from matched_analysis import digest, read, match_cases

METRICS = ["tvd", "mpro", "aligned_shared_axis", "aligned_bp", "direction_concordance"]
CONFIGS = [("primary", 20, 3), ("count50", 50, 3), ("count100", 100, 3), ("controls2", 20, 2)]
ASSIGNMENTS = list(combinations(range(4), 2))


def statistics(effects, groups, signs):
    values = effects.loc[np.asarray(groups).ravel()].to_numpy().reshape(len(groups), -1, 4)
    aligned = values[:, :, 2:] * signs[:, None, None]
    individual = np.concatenate([values[:, :, :2], aligned,
        (aligned[:, :, 1:2] > 0).astype(float)], axis=2)
    excess = individual[:, 0] - individual[:, 1:].mean(axis=1)
    result = np.r_[np.median(excess[:, :4], axis=0), excess[:, 4].mean()]
    return result, individual, excess


def run():
    protocol_path = HERE / "vector-protocol.json"
    protocol = json.loads(protocol_path.read_text())
    covariates = protocol["matching_covariates"]
    genes = read(OUT / "all_gene_measurements.tsv.gz").set_index("gene_id")
    vectors = read(OUT / "shared_pas_vectors.tsv.gz")
    eligible20 = genes.index[genes.measurable20]
    vectors = vectors[vectors.gene_id.isin(eligible20)]
    profiles = {g: (v[[f"usage__{s}" for s in ENDPOINT]].to_numpy().T,
                    v.local_cut_bp.to_numpy()) for g, v in vectors.groupby("gene_id")}
    effects = {}
    for source in ASSIGNMENTS:
        effects[source] = pd.DataFrame({g: profile_metrics(p, x, source) for g, (p, x) in profiles.items()}).T
        effects[source].columns = ["tvd", "mpro", "delta_shared_axis", "delta_bp"]
    # Coverage, covariates and thus all matching decisions are label-invariant.
    assert genes.loc[eligible20, covariates].notna().all().all()
    frozen = read(HERE / "inputs/frozen_candidates.tsv")
    excluded = set(frozen.gene_id)
    summaries, members, balances, statuses, nulls, sensitivities = [], [], [], [], [], []
    for configuration, threshold, k in CONFIGS:
        measured = genes[genes[f"measurable{threshold}"]]
        controls = measured[~measured.index.isin(excluded)]
        for name, roster in frozen.groupby("set_name"):
            cases = measured[measured.index.isin(roster.gene_id)]
            groups, balance, matching = match_cases(cases[covariates + ["pas_stratum"]],
                controls[covariates + ["pas_stratum"]], k, protocol, covariates=covariates)
            balance["configuration"], balance["set_name"] = configuration, name
            balances.append(balance)
            status = roster.merge(matching, on="gene_id", how="left", validate="one_to_one")
            status["configuration"] = configuration
            status["measurable"] = status.gene_id.isin(cases.index)
            status["match_status"] = np.where(~status.measurable, "not_measurable",
                np.where(status.complete_disjoint_match.eq(True), "matched", "no_complete_match"))
            statuses.append(status)
            if not groups:
                raise ValueError(f"No matched sets for {configuration}/{name}; do not silently omit a design")
            signs = roster.set_index("gene_id").loc[[group[0] for group in groups], "direction_sign"].to_numpy()
            observed, individual, excess = statistics(effects[(0, 1)], groups, signs)
            null = []
            for source in ASSIGNMENTS:
                stat, _, _ = statistics(effects[source], groups, signs)
                null.append(stat)
                for index, metric in enumerate(METRICS):
                    nulls.append({"configuration": configuration, "set_name": name, "statistic": metric,
                        "source_libraries": ";".join(ENDPOINT[i] for i in source), "value": stat[index]})
            null = np.asarray(null)
            tail = (null >= observed[None, :] - 1e-12).sum(axis=0) / len(ASSIGNMENTS)
            weight_values = []
            selected = {g for group in groups for g in group}
            for sw, tw in product([0., .5, 1.], repeat=2):
                records = {}
                for gene in selected:
                    p, x = profiles[gene]
                    source = sw * p[0] + (1 - sw) * p[1]
                    target = tw * p[2] + (1 - tw) * p[3]
                    difference = target - source
                    records[gene] = [np.abs(difference).sum() / 2, np.ptp(difference),
                                    (difference @ x) / np.ptp(x), difference @ x]
                values, _, _ = statistics(pd.DataFrame(records).T, groups, signs)
                weight_values.append(values)
                for index, metric in enumerate(METRICS):
                    sensitivities.append({"configuration": configuration, "set_name": name,
                        "statistic": metric, "MEF_R1_weight": sw, "iPSC_R1_weight": tw,
                        "value": values[index]})
            weights = np.asarray(weight_values)
            maximum_smd = balance.smd_after.abs().max()
            retention = len(groups) / len(cases)
            passed = retention >= .8 and maximum_smd <= .2
            for index, metric in enumerate(METRICS):
                summaries.append({"configuration": configuration, "set_name": name, "statistic": metric,
                    "frozen_n": len(roster), "measurable_n": len(cases), "matched_n": len(groups),
                    "controls_n": len(groups) * k, "retention": retention, "max_abs_smd": maximum_smd,
                    "design_gate_pass": bool(passed), "estimate": observed[index],
                    "candidate_median_or_rate": individual[:, 0, index].mean() if index == 4 else np.median(individual[:, 0, index]),
                    "control_median_or_rate": individual[:, 1:, index].mean() if index == 4 else np.median(individual[:, 1:, index]),
                    "whole_library_upper_tail_diagnostic": tail[index], "label_assignments": 6,
                    "whole_library_sensitivity_min": weights[:, index].min(),
                    "whole_library_sensitivity_max": weights[:, index].max()})
            for group_index, group in enumerate(groups):
                for role, gene in enumerate(group):
                    row = {"configuration": configuration, "set_name": name,
                        "case_gene": group[0], "gene_id": gene,
                        "role": "candidate" if role == 0 else "control", "direction_sign": signs[group_index]}
                    row.update({metric: individual[group_index, role, index] for index, metric in enumerate(METRICS)})
                    members.append(row)
    for frame, name in [(pd.DataFrame(summaries), "vector_summary.tsv"),
        (pd.DataFrame(members), "vector_matched_members.tsv"), (pd.concat(balances), "vector_matching_balance.tsv"),
        (pd.concat(statuses), "vector_candidate_status.tsv"), (pd.DataFrame(nulls), "whole_library_label_assignments.tsv"),
        (pd.DataFrame(sensitivities), "whole_library_sensitivity.tsv")]:
        save(frame, name)
    record = {"status": "completed", "script_sha256": digest(Path(__file__)),
        "matching_script_sha256": digest(HERE / "matched_analysis.py"), "protocol_sha256": digest(protocol_path),
        "input_sha256": {str(p.relative_to(HERE)): digest(p) for p in
            [OUT / "all_gene_measurements.tsv.gz", OUT / "shared_pas_vectors.tsv.gz", HERE / "inputs/frozen_candidates.tsv"]},
        "diagnostic_scope": "All six whole-library assignments preserve gene dependence; exchangeability is not established and n=2 prevents fine permutation resolution; no confirmatory p/q or gene-as-replicate interval",
        "old_analysis_preserved": True}
    (OUT / "vector_run.json").write_text(json.dumps(record, indent=2) + "\n")
    print(pd.DataFrame(summaries).query("configuration == 'primary'").to_string(index=False))


if __name__ == "__main__":
    run()
