"""Joint conditional overlap and direct specificity of the shared proximal module."""
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD_PATH = ROOT / "notes/2026-10-10-naive-mouse-correspondence/analyze.py"
spec = importlib.util.spec_from_file_location("naive_base", OLD_PATH)
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
TABLES = old.TABLES
INPUTS = {}
B = 20000
LABELS = {
    "naive_direction": ["pooled_hN_P7"],
    "Liu_primed_source": ["pooled_hN_P7", "primed_P7_S5_pool"],
    "Liu_same_library_controls": ["pooled_hN_P7", "primed_P7_S7", "M7_S7", "E7_S7", "primed_P7_M7"],
    "human_recurrent_controls": ["pooled_hN_P7", "human_P_S_recurrent", "human_P_M_recurrent"],
}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(relative, **kwargs):
    path = TABLES / relative
    INPUTS[str(path.relative_to(ROOT))] = digest(path)
    return pd.read_csv(path, sep="\t", float_precision="round_trip", **kwargs)


def save(frame, name):
    frame.to_csv(HERE / name, sep="\t", index=False, na_rep="NA", float_format="%.17g",
                 compression={"method": "gzip", "mtime": 0} if name.endswith(".gz") else None)


def joint_draw(labels, strata, n, rng):
    """Permute rows, preserving the entire response signature, not its columns."""
    out = np.broadcast_to(labels, (n,) + labels.shape).copy()
    for ii in strata:
        if len(ii) < 2 or np.unique(labels[ii], axis=0).shape[0] == 1:
            continue
        perm = rng.permuted(np.broadcast_to(np.arange(len(ii)), (n, len(ii))), axis=1)
        out[:, ii] = labels[ii][perm]
        assert np.array_equal(out[:, ii].sum(1), np.broadcast_to(labels[ii].sum(0), (n, labels.shape[1])))
    return out


def probabilities(labels, strata):
    out = np.zeros(labels.shape, float)
    for ii in strata:
        out[ii] = labels[ii].mean(0)
    return out


def simulate(h, m, hs, ms, n, rng):
    out = []
    for start in range(0, n, 250):
        size = min(250, n-start)
        hd, md = joint_draw(h, hs, size, rng), joint_draw(m, ms, size, rng)
        out.append(np.einsum("bnh,bnm->bhm", hd.astype(np.int16), md.astype(np.int16)))
    return np.concatenate(out)


def strata(frame, side, names, threshold, model):
    z = pd.DataFrame(index=frame.index)
    if side == "h":
        for name in names:
            if name + "_observation_signature" in frame:
                z[name+"_observation_signature"] = frame[name+"_observation_signature"]
            for direction in ["proximal", "distal"]:
                z[name + "_" + direction + "_possible"] = frame[name + "_capacity_" + direction].ge(threshold)
    else:
        for direction in ["proximal", "distal"]:
            z[direction + "_possible"] = frame["m_capacity_" + direction].ge(threshold)
    for covariate, number in [("PAS", 3), ("coverage", 3), ("baseline", 5 if "fine" in model else 3)]:
        z[covariate] = old.bins(frame[side + "_" + covariate], number)
    if side == "m" and "fine" in model or side == "h" and "human_expression" in model:
        z["GEX"] = old.bins(frame[side + "_GEX"], 3)
    return old.cells(z, list(z.columns)), z


def direction_labels(delta, threshold):
    delta = np.asarray(delta, float)
    return np.column_stack([delta <= -threshold, delta >= threshold])


def human_base(pairs, axes, state_values, control_pairs):
    frames, reconstructed = old.human_frames(pairs, axes)
    naive = frames["pooled_two_hN_vs_primed"].copy()
    q = reconstructed[reconstructed.node_a.eq("N") & reconstructed.node_b.eq("P") &
                      reconstructed.library_a.str.endswith("6") & reconstructed.library_b.str.endswith("7")].set_index("gene_id")
    naive["capacity_proximal"] = naive.h_capacity
    naive["capacity_distal"] = naive.human_gene_id.map(q.axis_max-q.reference)
    naive = naive.rename(columns={"h_delta": "delta", "h_PAS": "PAS", "h_coverage": "coverage", "h_baseline": "baseline"})
    keep = ["human_gene_id", "h_gene_symbol", "delta", "PAS", "coverage", "baseline", "capacity_proximal", "capacity_distal"]
    out = {"pooled_hN_P7": naive[keep].copy()}
    p = reconstructed[reconstructed.node_a.eq("P") & reconstructed.node_b.eq("S") &
                      reconstructed.library_a.str.endswith("7") & reconstructed.library_b.str.endswith("5")].copy()
    assert p.gene_id.is_unique
    p["capacity_proximal"] = p.reference-p.axis_min
    p["capacity_distal"] = p.axis_max-p.reference
    out["primed_P7_S5_pool"] = p.rename(columns={"gene_id": "human_gene_id", "gene_symbol": "h_gene_symbol",
        "delta_DWUI": "delta", "PAS_number": "PAS", "coverage": "coverage", "supported_baseline": "baseline"})[keep].copy()
    s = state_values[state_values.study_id.eq("GSE150311") & state_values.setting.eq(old.SETTING) & state_values.library_key.str.endswith("7")].copy()
    s7 = s[s.state.eq("P")].set_index("gene_id")
    # P7 from exact N-P pairs fixes the otherwise arbitrary axis translation.
    offset = s7.value-q.reference
    span = q.axis_max-q.axis_min
    common_anchor = s7.index.intersection(q.index)
    assert s7.loc[common_anchor].available_pas.eq(q.loc[common_anchor].PAS_number).all()
    s["position"] = s.value-s.gene_id.map(offset)
    s["supported_span"] = s.gene_id.map(span)
    s["relative_baseline"] = s.position/s.supported_span
    s = s[s.gene_id.isin(q.index)]
    assert s.position.ge(-1e-12).all() and (s.position <= s.supported_span+1e-12).all()
    cp = control_pairs[control_pairs.study_id.eq("GSE150311") & control_pairs.setting.eq(old.SETTING) & control_pairs.library_key.str.endswith("7")]
    for label, contrast in [("primed_P7_S7", "P_S"), ("M7_S7", "M_S"), ("E7_S7", "E_S"), ("primed_P7_M7", "P_M")]:
        a, b = contrast.split("_")
        z = cp.loc[cp.contrast.eq(contrast), ["gene_id", "gene_symbol", "delta"]].copy()
        z = z.merge(s[s.state.eq(a)], on="gene_id", suffixes=("", "_a"), validate="one_to_one")
        z = z.merge(s[s.state.eq(b)], on="gene_id", suffixes=("", "_b"), validate="one_to_one")
        assert np.max(np.abs(z.position-z.position_b-z.delta)) < 1e-12
        assert z.available_pas.eq(z.available_pas_b).all()
        z["PAS"] = z.available_pas
        z["coverage"] = (z.pas_counts/z.n_mc+z.pas_counts_b/z.n_mc_b)/2
        z["baseline"] = z.relative_baseline_b
        z["capacity_proximal"] = z.position_b
        z["capacity_distal"] = z.supported_span-z.position_b
        out[label] = z.rename(columns={"gene_id": "human_gene_id", "gene_symbol": "h_gene_symbol"})[keep].copy()
        if contrast == "P_S":
            exact = reconstructed[reconstructed.node_a.eq("P") & reconstructed.node_b.eq("S") &
                reconstructed.library_a.str.endswith("7") & reconstructed.library_b.str.endswith("7")].set_index("gene_id")
            assert np.max(np.abs(out[label].set_index("human_gene_id").delta-exact.delta_DWUI)) < 1e-12
    n = out["pooled_hN_P7"]
    assert len(n) == 2575 and int(n.delta.le(-.02).sum()) == 662 and int(n.delta.ge(.02).sum()) == 422
    return out, len(reconstructed)


def human_recurrent(frozen, recurrence, threshold):
    out = {}
    f = frozen[frozen.species.eq("human") & frozen.threshold.eq(threshold)]
    for contrast in ["P_S", "P_M"]:
        z = f[f.contrast.eq(contrast)].copy()
        assert not z.duplicated(["gene_id", "unit"]).any()
        count = z.groupby("gene_id").unit.nunique()
        z = z[z.gene_id.isin(count[count.ge(2)].index)]
        z["observation_signature"] = z.unit.astype(str) + ":" + z.series_pattern.astype(str)
        cov = z.groupby("gene_id").agg(h_gene_symbol=("gene_symbol", "first"), PAS=("PAS_number", "median"),
            coverage=("APA_coverage", "median"), baseline=("baseline", "median"), GEX=("mean_GEX", "median"),
            capacity_proximal=("capacity_negative", lambda x: np.sort(x)[-2]),
            capacity_distal=("capacity_positive", lambda x: np.sort(x)[-2]))
        cov["observation_signature"] = z.groupby("gene_id").observation_signature.agg(lambda x: "|".join(sorted(x)))
        votes = z.assign(proximal=z.event.eq(-1), distal=z.event.eq(1)).groupby("gene_id")[["proximal", "distal"]].sum().ge(2)
        for sign, original in [("proximal", "negative"), ("distal", "positive")]:
            ids = set(recurrence[recurrence.species.eq("human") & recurrence.contrast.eq(contrast) & recurrence.threshold.eq(threshold) & recurrence.direction.eq(original)].gene_id)
            assert ids == set(votes[votes[sign]].index)
        assert not (votes.proximal & votes.distal).any()
        cov["proximal"] = votes.proximal
        cov["distal"] = votes.distal
        out["human_"+contrast+"_recurrent"] = cov.reset_index().rename(columns={"gene_id": "human_gene_id"})
    return out


def mouse(frozen, recurrence, vectors, axes, threshold):
    old_frames = old.mouse_frames(frozen, recurrence, vectors, axes, threshold)
    f = frozen[frozen.species.eq("mouse") & frozen.contrast.eq("P_M") & frozen.threshold.eq(threshold)]
    f = f[f.gene_id.isin(old_frames["P_M_proximal_recurrent"].mouse_gene_id)]
    dist = f.assign(selected=f.event.eq(1)).groupby("gene_id").selected.all()
    r = recurrence[recurrence.species.eq("mouse") & recurrence.contrast.eq("P_M") & recurrence.threshold.eq(threshold) & recurrence.direction.eq("positive")]
    assert set(r.gene_id) == set(dist[dist].index)
    v = vectors[vectors.species.eq("mouse") & vectors.setting.eq(old.SETTING) & vectors.route.isin(["CIP", "WOT"])]
    axis = axes[axes.species.eq("mouse") & axes.setting.eq(old.SETTING)].set_index("atlas_pas_id").WPP
    v = v.assign(axis_max=v.PAS_ids.map(lambda x: float(axis.loc[json.loads(x)].max())))
    v = v.assign(dist_pass=v.DWUI_PM.ge(threshold) & v.DWUI_PS.ge(threshold),
                 pos_capacity=v.axis_max-v[["M_DWUI", "S_DWUI"]].max(axis=1))
    out = {}
    for name, of in old_frames.items():
        z = of.rename(columns={"m_selected": "proximal", "m_capacity": "m_capacity_proximal"}).copy()
        z["distal"] = z.mouse_gene_id.map(dist).fillna(False)
        if name.startswith("P_M_PS"):
            z["distal"] &= z.mouse_gene_id.map(v.groupby("gene_id").dist_pass.all()).fillna(False)
            z["m_capacity_distal"] = z.mouse_gene_id.map(v.groupby("gene_id").pos_capacity.min())
            target = "P_M_PS_same_contributor_dual"
        else:
            z["m_capacity_distal"] = z.mouse_gene_id.map(f.groupby("gene_id").capacity_positive.min())
            target = "P_M_recurrent"
        assert not (z.proximal & z.distal).any()
        assert not (z.proximal & z.m_capacity_proximal.lt(threshold)).any()
        assert not (z.distal & z.m_capacity_distal.lt(threshold)).any()
        if threshold == .02:
            assert int(z.proximal.sum()) == (207 if "dual" in target else 372)
            assert int(z.distal.sum()) == (15 if "dual" in target else 48)
        out[target] = z
    return out


def main():
    for p in [HERE / "protocol.json", OLD_PATH, Path(__file__)]:
        INPUTS[str(p.relative_to(ROOT))] = digest(p)
    orth = read("cross_overlap/ortholog_matched_features.tsv.gz")
    assert len(orth) == 16092 and orth.human_gene_id.is_unique and orth.mouse_gene_id.is_unique
    assert orth.human_partner_count.eq(1).all() and orth.mouse_partner_count.eq(1).all()
    orth = orth[["human_gene_id", "mouse_gene_id"]]
    pairs = read("naive_state/library_pair_effects.tsv.gz")
    axes = read("state_vectors/PAS_axes.primary.tsv.gz", usecols=["atlas_pas_id", "species", "study_id", "setting", "WPP"])
    axes = old.complete_human_axes(pairs, axes, read("tandem/tandem_utr_eligibility.tsv.gz"))
    hf, reconstructed = human_base(pairs, axes, read("intermediate/library_state_values.tsv.gz"), read("intermediate/paired_library_effects.tsv.gz"))
    frozen = read("recurrence_controls/recurrence_frozen_covariates.tsv.gz")
    recurrence = read("recurrence/recurrence_members.tsv.gz")
    recurrence = recurrence[recurrence.metric.eq("WPP") & recurrence.series_matching & recurrence.variant.eq("primary")]
    vectors = read("matched_pms/matched_state_vectors.tsv.gz")
    models = ["coarse_joint_assay", "fine_joint_assay_mouse_expression"]
    # Human expression is an additional prespecified model, never a replacement.
    gex_path = ROOT / "notes/2026-10-10-naive-gex-recovery/all_gene_library_expression.tsv.gz"
    gex_used = gex_path.exists()
    if gex_used:
        INPUTS[str(gex_path.relative_to(ROOT))] = digest(gex_path)
        gex = pd.read_csv(gex_path, sep="\t", float_precision="round_trip")
        gex = gex[gex.qualified].rename(columns={"state": "node"})
        key = ["gene_id", "library_key", "node"]
        assert not gex.duplicated(key).any()
        expr = gex.set_index(key)["mean_log2CPM1"]
        endpoints = {
            "pooled_hN_P7": [("6", "N", .25), ("8", "N", .25), ("7", "P", .5)],
            "primed_P7_S5_pool": [("7", "P", .5), ("5", "S", .5)],
            "primed_P7_S7": [("7", "P", .5), ("7", "S", .5)],
            "M7_S7": [("7", "M", .5), ("7", "S", .5)],
            "E7_S7": [("7", "E", .5), ("7", "S", .5)],
            "primed_P7_M7": [("7", "P", .5), ("7", "M", .5)],
        }
        for name, frame in hf.items():
            values = []
            for number, node, weight in endpoints[name]:
                library = "GSE150311::gse150311_human_pluripotency_reprogramming_GSM454654"+number
                values.append(frame.human_gene_id.map(expr.xs((library, node), level=(1, 2))).to_numpy()*weight)
            frame["GEX"] = np.sum(values, axis=0)
        models += ["coarse_joint_assay_human_expression", "fine_joint_assay_mouse_expression_human_expression"]
    rows, direct, inventories, members, margins, draws = [], [], [], [], [], {}
    for threshold in [.02, .05]:
        humans = {**hf, **human_recurrent(frozen, recurrence, threshold)}
        mice = mouse(frozen, recurrence, vectors, axes, threshold)
        for block, names in LABELS.items():
            for target, mf in mice.items():
                f = orth.merge(mf, on="mouse_gene_id", validate="one_to_one")
                for name in names:
                    q = humans[name].copy()
                    if "delta" in q:
                        q[["proximal", "distal"]] = direction_labels(q.delta, threshold)
                    assert not (q.proximal & q.distal).any()
                    for direction in ["proximal", "distal"]:
                        assert not (q[direction] & q["capacity_"+direction].lt(threshold)).any()
                    renamed = q.rename(columns={c: name+"_"+c for c in q if c != "human_gene_id"})
                    f = f.merge(renamed, on="human_gene_id", validate="one_to_one")
                f = f.sort_values("human_gene_id").reset_index(drop=True)
                context = dict(block=block, mouse_target=target, threshold=threshold)
                inventories.append(context | dict(n_common_orthologs=len(f), mouse_proximal=int(f.proximal.sum()), mouse_distal=int(f.distal.sum()),
                    **{name+"_measurable": len(humans[name]) for name in names},
                    **{name+"_proximal": int(f[name+"_proximal"].sum()) for name in names},
                    **{name+"_distal": int(f[name+"_distal"].sum()) for name in names}))
                for cov in ["PAS", "coverage", "baseline"] + (["GEX"] if gex_used else []):
                    f["h_"+cov] = f[[name+"_"+cov for name in names]].median(axis=1)
                hcols = [name+"_"+direction for name in names for direction in ["proximal", "distal"]]
                h = f[hcols].to_numpy(bool)
                m = f[["proximal", "distal"]].to_numpy(bool)
                observed = h.astype(int).T @ m.astype(int)
                members.append(f.assign(**context))
                for model in models:
                    hs, hstrata = strata(f, "h", names, threshold, model)
                    ms, mstrata = strata(f, "m", names, threshold, model)
                    hp, mp = probabilities(h, hs), probabilities(m, ms)
                    expectation = hp.T @ mp
                    run_seed = old.seed("proximal_module_specificity", block, target, threshold, model)
                    null = simulate(h, m, hs, ms, B, np.random.default_rng(run_seed))
                    sd = null.std(axis=0, ddof=1)
                    assert np.all(np.abs(null.mean(axis=0)-expectation) < 6*sd/np.sqrt(B)+1e-10)
                    key = "|".join(map(str, [block, target, threshold, model]))
                    draws[key] = null
                    for side, lbl, ss, prob in [("human", h, hs, hp), ("mouse", m, ms, mp)]:
                        for j, ii in enumerate(ss):
                            signatures, counts = np.unique(lbl[ii], axis=0, return_counts=True)
                            margins.append(context | dict(model=model, side=side, stratum=j, genes=len(ii),
                                joint_pattern_counts=json.dumps({"".join(map(str, sig.astype(int))): int(count) for sig, count in zip(signatures, counts)})))
                    for i, hc in enumerate(hcols):
                        human_label, direction = hc.rsplit("_", 1)
                        for j, mouse_direction in enumerate(["proximal", "distal"]):
                            nul = null[:, i, j]
                            obs, exp = int(observed[i, j]), float(expectation[i, j])
                            informative = bool(nul.std() > 0)
                            tail = int((nul >= obs).sum())
                            rows.append(context | dict(model=model, human_label=human_label, human_direction=direction, mouse_direction=mouse_direction,
                                n_orthologs=len(f), human_selected=int(h[:, i].sum()), mouse_selected=int(m[:, j].sum()), observed=obs,
                                analytic_expected=exp, null_mean=float(nul.mean()), null_sd=float(nul.std(ddof=1)),
                                null_p025=float(np.quantile(nul, .025)), null_p975=float(np.quantile(nul, .975)),
                                fold=obs/exp if exp else np.nan, excess=obs-exp, tail_count=tail, p=(1+tail)/(B+1) if informative else 1.,
                                informative=informative, permutations=B, seed=run_seed,
                                human_movable_fraction=float(((hp[:, i]>0)&(hp[:, i]<1)).mean()),
                                mouse_movable_fraction=float(((mp[:, j]>0)&(mp[:, j]<1)).mean()),
                                human_selected_movable_fraction=float(((hp[:, i]>0)&(hp[:, i]<1))[h[:, i]].mean()) if h[:, i].any() else 0.,
                                mouse_selected_movable_fraction=float(((mp[:, j]>0)&(mp[:, j]<1))[m[:, j]].mean()) if m[:, j].any() else 0.))
                    for ci, control in enumerate(names[1:], start=1):
                        enaive, econtrol = expectation[0, 0], expectation[2*ci, 0]
                        if enaive and econtrol:
                            stat = observed[0, 0]/enaive-observed[2*ci, 0]/econtrol
                            nd = null[:, 0, 0]/enaive-null[:, 2*ci, 0]/econtrol
                            informative = bool(nd.std() > 1e-12)
                            tail = int((nd >= stat-1e-12).sum())
                        else:
                            stat, nd, informative, tail = np.nan, np.zeros(B), False, B
                        direct.append(context | dict(model=model, control=control, n_orthologs=len(f),
                            naive_observed=int(observed[0, 0]), control_observed=int(observed[2*ci, 0]),
                            naive_expected=float(enaive), control_expected=float(econtrol),
                            naive_fold=float(observed[0, 0]/enaive) if enaive else np.nan,
                            control_fold=float(observed[2*ci, 0]/econtrol) if econtrol else np.nan,
                            fold_difference=float(stat), null_mean=float(nd.mean()), null_sd=float(nd.std(ddof=1)),
                            null_p025=float(np.quantile(nd, .025)), null_p975=float(np.quantile(nd, .975)),
                            tail_count=tail, p=(1+tail)/(B+1) if informative else 1., informative=informative,
                            permutations=B, seed=run_seed))
                    print(key, "N", len(f), "prox overlaps", observed[::2, 0].tolist(), "expected", np.round(expectation[::2, 0], 3).tolist(), flush=True)
    result, contrast = pd.DataFrame(rows), pd.DataFrame(direct)
    contrast["q_specificity14_within_model_threshold"] = contrast.groupby(["model", "threshold"]).p.transform(old.bh)
    contrast["q_all_specificity_exploratory"] = old.bh(contrast.p)
    result["direction_control"] = result.block.eq("naive_direction") & result.human_direction.ne(result.mouse_direction)
    result["q_direction4_within_model_threshold"] = np.nan
    mask = result.direction_control
    result.loc[mask, "q_direction4_within_model_threshold"] = result[mask].groupby(["model", "threshold"]).p.transform(old.bh)
    result["q_all_overlaps_within_model_threshold"] = result.groupby(["model", "threshold"]).p.transform(old.bh)
    result["q_all_overlaps_exploratory"] = old.bh(result.p)
    assert contrast.groupby(["model", "threshold"]).size().eq(14).all()
    assert result[mask].groupby(["model", "threshold"]).size().eq(4).all()
    save(result, "overlap_results.tsv")
    save(contrast, "direct_specificity_results.tsv")
    save(pd.DataFrame(inventories), "measurability_inventory.tsv")
    save(pd.concat(members, ignore_index=True), "joint_membership_and_covariates.tsv.gz")
    save(pd.DataFrame(margins), "joint_pattern_stratum_margins.tsv.gz")
    np.savez_compressed(HERE / "joint_null_draws.npz", **draws)
    assert all(digest(ROOT / name) == value for name, value in INPUTS.items())
    outputs = ["overlap_results.tsv", "direct_specificity_results.tsv", "measurability_inventory.tsv", "joint_membership_and_covariates.tsv.gz", "joint_pattern_stratum_margins.tsv.gz", "joint_null_draws.npz"]
    manifest = dict(status="PASS", input_sha256=INPUTS, output_sha256={name: digest(HERE/name) for name in outputs},
        python_executable=sys.executable, python=sys.version, numpy=np.__version__, pandas=pd.__version__, human_pair_exact_reconstruction_rows=reconstructed,
        human_control_exact_reconstruction=True, fixed_naive2575_662prox422dist=True, mouse_rosters_verified=True,
        joint_label_patterns_preserved=True, finite_p=True, models=models, human_GEX_sensitivity_available=gex_used)
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


if __name__ == "__main__":
    main()
