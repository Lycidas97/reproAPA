"""Matched state-specific human--mouse APA correspondence and technical audit."""
from pathlib import Path
import hashlib, json

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, pearsonr, spearmanr

P = Path(__file__).resolve().parents[1]
UP = P.parent / "2026-09-02_mouse-time-aware-integration-pilot"
D = P / "data/figure4-state-recurrence-v25"
A = P / "analysis/figure4-state-recurrence-v25"
D.mkdir(exist_ok=True)
A.mkdir(exist_ok=True)

NODES_PATH = P / "data/figure4-program-v9/all_mapped_nodes.tsv.gz"
nodes = pd.read_csv(NODES_PATH, sep="\t")
codes = ["W", "C", "H1", "H2", "H3", "H4"]
study_names = dict(nodes.drop_duplicates("code")[["code", "study"]].itertuples(index=False, name=None))
tau = 0.02
B = 1000
rng = np.random.default_rng(250920)

valid = nodes[nodes.code.isin(codes) & nodes.measurement_pass].copy()
wide = valid.pivot(index=["layer", "human_gene_id", "code"], columns="node", values="position_index")
for node in ["S", "R", "E", "P"]:
    if node not in wide:
        wide[node] = np.nan
wide = wide[["S", "R", "E", "P"]]


def effect_frame(layer, mouse, human, state, matched_state=None):
    """Return true ortholog effects; optionally require state and P together."""
    m = wide.xs((layer, mouse), level=("layer", "code"))
    h = wide.xs((layer, human), level=("layer", "code"))
    needed = ["S", state] if matched_state is None else ["S", state, matched_state]
    m = m.dropna(subset=needed)
    h = h.dropna(subset=needed)
    ids = m.index.intersection(h.index)
    out = pd.DataFrame(index=ids)
    for x in sorted(set([state] + ([matched_state] if matched_state else []))):
        out[f"m_{x}"] = m.loc[ids, x] - m.loc[ids, "S"]
        out[f"h_{x}"] = h.loc[ids, x] - h.loc[ids, "S"]
    return out


def metrics(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    strong = (np.abs(x) >= tau) & (np.abs(y) >= tau)
    same = strong & (np.sign(x) == np.sign(y))
    opposite = strong & (np.sign(x) != np.sign(y))
    return {
        "n": len(x),
        "spearman": spearmanr(x, y).statistic,
        "pearson": pearsonr(x, y).statistic,
        "n_both_strong": int(strong.sum()),
        "direction_concordance": float(same.sum() / strong.sum()) if strong.any() else np.nan,
        "direction_score_all": float((same.sum() - opposite.sum()) / len(x)),
    }


state_pairs = {
    "R": [(m, h) for m in ["W", "C"] for h in ["H1", "H3"]],
    "E": [(m, h) for m in ["W", "C"] for h in ["H1", "H2"]],
    "P": [(m, h) for m in ["W", "C"] for h in ["H1", "H2", "H3", "H4"]],
}

# Unselected state-specific backgrounds.
state_rows = []
state_null_rows = []
for layer in ["supported50", "expanded_primary"]:
    for state, pairs in state_pairs.items():
        pair_data = []
        for mouse, human in pairs:
            q = effect_frame(layer, mouse, human, state)
            pair_data.append(q)
            state_rows.append({"layer": layer, "state": state, "pair": f"{mouse}-{human}", **metrics(q[f"m_{state}"], q[f"h_{state}"])})
        observed = {
            k: np.nanmean([metrics(q[f"m_{state}"], q[f"h_{state}"])[k] for q in pair_data])
            for k in ["spearman", "pearson", "direction_concordance", "direction_score_all"]
        }
        null = {k: np.empty(B) for k in observed}
        for b in range(B):
            vals = {k: [] for k in observed}
            for q in pair_data:
                perm = rng.permutation(len(q))
                met = metrics(q[f"m_{state}"].to_numpy(), q[f"h_{state}"].to_numpy()[perm])
                for k in observed:
                    vals[k].append(met[k])
            for k in observed:
                null[k][b] = np.nanmean(vals[k])
        for k in observed:
            state_null_rows.append({
                "layer": layer, "state": state, "metric": k, "observed": observed[k],
                "null_mean": float(null[k].mean()), "null_lo": float(np.quantile(null[k], 0.025)),
                "null_hi": float(np.quantile(null[k], 0.975)),
                "p_greater": float((1 + (null[k] >= observed[k]).sum()) / (B + 1)),
            })
pd.DataFrame(state_rows).to_csv(D / "state_specific_pair_metrics.tsv", sep="\t", index=False)
state_null = pd.DataFrame(state_null_rows)
state_null["q_within_layer_12"] = state_null.groupby("layer").p_greater.transform(
    lambda x: false_discovery_control(x.to_numpy())
)
state_null["q_all_24"] = false_discovery_control(state_null.p_greater)
state_null.to_csv(D / "state_specific_permutation.tsv", sep="\t", index=False)

# Matched P-vs-R and P-vs-E. The same human permutation is applied to both
# states so that study architecture, state distributions, and missingness stay fixed.
matched_rows = []
permutation_rows = []
observation_rows = []
for layer in ["supported50", "expanded_primary"]:
    for state in ["R", "E"]:
        pair_data = []
        for mouse, human in state_pairs[state]:
            q = effect_frame(layer, mouse, human, state, "P")
            pair_data.append((mouse, human, q))
            for x in [state, "P"]:
                met = metrics(q[f"m_{x}"], q[f"h_{x}"])
                matched_rows.append({"layer": layer, "comparison": f"P_vs_{state}", "state": x, "pair": f"{mouse}-{human}", **met})
                strong = (q[f"m_{x}"].abs() >= tau) & (q[f"h_{x}"].abs() >= tau)
                relation = np.where(~strong, "weak_or_mixed", np.where(np.sign(q[f"m_{x}"]) == np.sign(q[f"h_{x}"]), "concordant", "discordant"))
                for gid, rel in zip(q.index, relation):
                    observation_rows.append({
                        "layer": layer, "comparison": f"P_vs_{state}", "state": x,
                        "pair": f"{mouse}-{human}", "mouse_code": mouse, "human_code": human,
                        "human_gene_id": gid, "mouse_effect": q.loc[gid, f"m_{x}"],
                        "human_effect": q.loc[gid, f"h_{x}"], "relation": rel,
                    })

        observed = {}
        for metric in ["spearman", "pearson", "direction_concordance", "direction_score_all"]:
            observed[metric] = {
                x: np.nanmean([metrics(q[f"m_{x}"], q[f"h_{x}"])[metric] for _, _, q in pair_data])
                for x in [state, "P"]
            }
        null = {k: np.empty(B) for k in observed}
        for b in range(B):
            vals = {k: {state: [], "P": []} for k in observed}
            for _, _, q in pair_data:
                perm = rng.permutation(len(q))
                for x in [state, "P"]:
                    met = metrics(q[f"m_{x}"].to_numpy(), q[f"h_{x}"].to_numpy()[perm])
                    for k in observed:
                        vals[k][x].append(met[k])
            for k in observed:
                null[k][b] = np.nanmean(vals[k]["P"]) - np.nanmean(vals[k][state])
        for k in observed:
            delta = observed[k]["P"] - observed[k][state]
            permutation_rows.append({
                "layer": layer, "comparison": f"P_vs_{state}", "metric": k,
                "state_value": observed[k][state], "P_value": observed[k]["P"],
                "observed_difference": delta, "null_mean": float(null[k].mean()),
                "null_lo": float(np.quantile(null[k], 0.025)),
                "null_hi": float(np.quantile(null[k], 0.975)),
                "p_P_stronger": float((1 + (null[k] >= delta).sum()) / (B + 1)),
            })

pd.DataFrame(matched_rows).to_csv(D / "matched_pair_metrics.tsv", sep="\t", index=False)
matched_perm = pd.DataFrame(permutation_rows)
matched_perm["q_within_layer_8"] = matched_perm.groupby("layer").p_P_stronger.transform(
    lambda x: false_discovery_control(x.to_numpy())
)
matched_perm["q_all_16"] = false_discovery_control(matched_perm.p_P_stronger)
matched_perm.to_csv(D / "matched_state_permutation.tsv", sep="\t", index=False)
obs = pd.DataFrame(observation_rows)

# Technical covariates from node summaries and node-level library dispersion.
meta = valid.set_index(["layer", "human_gene_id", "code", "node"])[
    ["n_libraries", "qualified_metacells", "supported_pas_in_route", "supported_span"]
]
dispersion_parts = []
input_hashes = {str(NODES_PATH): hashlib.sha256(NODES_PATH.read_bytes()).hexdigest()}
for code in codes:
    study = study_names[code]
    for layer in ["supported50", "expanded_primary"]:
        path = UP / "data/apa-majority-shared-v1/profiles" / study / f"{layer}_library_positions.tsv.gz"
        input_hashes[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        z = pd.read_csv(path, sep="\t", usecols=["gene_id", "node", "position_index"])
        mapping = valid[valid.code.eq(code)][["gene_id", "human_gene_id"]].drop_duplicates().set_index("gene_id").human_gene_id
        z["human_gene_id"] = z.gene_id.map(mapping)
        z = z[z.node.isin(["S", "R", "E", "P"]) & z.human_gene_id.notna()]
        g = z.groupby(["human_gene_id", "node"]).position_index.agg(["std", "count"]).reset_index()
        g["code"] = code
        g["layer"] = layer
        dispersion_parts.append(g)
dispersion = pd.concat(dispersion_parts).set_index(["layer", "human_gene_id", "code", "node"])

tech_rows = []
for r in obs.itertuples(index=False):
    vals = []
    disps = []
    for code in [r.mouse_code, r.human_code]:
        for node in ["S", r.state]:
            key = (r.layer, r.human_gene_id, code, node)
            vals.append(meta.loc[key])
            disps.append(dispersion.loc[key, "std"] if key in dispersion.index else np.nan)
    frame = pd.DataFrame(vals)
    tech_rows.append({
        **r._asdict(),
        "min_libraries": frame.n_libraries.min(),
        "min_metacells": frame.qualified_metacells.min(),
        "min_pas": frame.supported_pas_in_route.min(),
        "min_span": frame.supported_span.min(),
        "max_library_sd": np.nanmax(disps) if np.isfinite(disps).any() else np.nan,
        "all_nodes_multi_library": bool((frame.n_libraries >= 2).all()),
    })
tech = pd.DataFrame(tech_rows)

# Cross-layer relation stability for the same matched observation.
keycols = ["comparison", "state", "pair", "human_gene_id"]
relwide = tech.pivot_table(index=keycols, columns="layer", values="relation", aggfunc="first")
relwide["relation_stable"] = relwide.get("supported50").eq(relwide.get("expanded_primary"))
tech = tech.merge(relwide[["relation_stable"]].reset_index(), on=keycols, how="left", validate="many_to_one")
tech.to_csv(D / "technical_observation_audit.tsv.gz", sep="\t", index=False)

# Association audit: within each study-pair/state stratum, permute concordant /
# discordant labels and compare mean log covariates. This is diagnostic because
# genes and studies recur across rows.
strong = tech[tech.relation.isin(["concordant", "discordant"]) & tech.layer.eq("supported50")].copy()
covariates = ["min_libraries", "min_metacells", "min_pas", "min_span", "max_library_sd"]
audit_rows = []
for cov in covariates:
    q = strong.dropna(subset=[cov]).copy()
    q["value"] = np.log1p(q[cov])
    observed = q.groupby("relation").value.mean()
    delta = observed.get("discordant", np.nan) - observed.get("concordant", np.nan)
    null = np.empty(B)
    groups = list(q.groupby(["comparison", "state", "pair"]).indices.values())
    labels = q.relation.to_numpy()
    values = q.value.to_numpy()
    for b in range(B):
        perm_labels = labels.copy()
        for g in groups:
            perm_labels[g] = labels[rng.permutation(g)]
        null[b] = values[perm_labels == "discordant"].mean() - values[perm_labels == "concordant"].mean()
    audit_rows.append({
        "covariate": cov, "n": len(q), "n_concordant": int((labels == "concordant").sum()),
        "n_discordant": int((labels == "discordant").sum()), "observed_log_difference": delta,
        "null_lo": np.quantile(null, 0.025), "null_hi": np.quantile(null, 0.975),
        "p_two_sided": (1 + (np.abs(null) >= abs(delta)).sum()) / (B + 1),
    })
audit = pd.DataFrame(audit_rows)
audit["q_five_covariates"] = false_discovery_control(audit.p_two_sided)
audit.to_csv(D / "technical_confound_permutation.tsv", sep="\t", index=False)

binary_rows = []
for cov in ["all_nodes_multi_library", "relation_stable"]:
    q = strong.dropna(subset=[cov]).copy()
    labels = q.relation.to_numpy()
    values = q[cov].astype(float).to_numpy()
    observed = values[labels == "discordant"].mean() - values[labels == "concordant"].mean()
    groups = list(q.groupby(["comparison", "state", "pair"]).indices.values())
    null = np.empty(B)
    for b in range(B):
        perm_labels = labels.copy()
        for g in groups:
            perm_labels[g] = labels[rng.permutation(g)]
        null[b] = values[perm_labels == "discordant"].mean() - values[perm_labels == "concordant"].mean()
    binary_rows.append({
        "covariate": cov, "n": len(q),
        "concordant_fraction": values[labels == "concordant"].mean(),
        "discordant_fraction": values[labels == "discordant"].mean(),
        "observed_fraction_difference": observed,
        "null_lo": np.quantile(null, 0.025), "null_hi": np.quantile(null, 0.975),
        "p_two_sided": (1 + (np.abs(null) >= abs(observed)).sum()) / (B + 1),
    })
binary = pd.DataFrame(binary_rows)
binary["q_two_covariates"] = false_discovery_control(binary.p_two_sided)
binary.to_csv(D / "technical_binary_audit.tsv", sep="\t", index=False)

summary = {
    "state_pairs": state_pairs,
    "threshold": tau,
    "permutations": B,
    "matched_rule": "same study pair and same genes measured at S, intermediate state, and P in both species",
    "permutation_rule": "permute human ortholog labels jointly for intermediate and P within each study pair",
    "technical_rule": "within pair/state label permutation; diagnostic, not biological-replicate inference",
    "input_sha256": input_hashes,
}
(A / "manifest.json").write_text(json.dumps(summary, indent=2) + "\n")
print(pd.DataFrame(permutation_rows).query("layer == 'supported50'").to_string(index=False))
print(audit.to_string(index=False))
print(binary.to_string(index=False))
