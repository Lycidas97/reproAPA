"""Independent arithmetic replay from frozen membership and saved joint draws."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def file_hash(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def adjust(values):
    values = np.asarray(values, float)
    order = sorted(range(len(values)), key=lambda i: values[i])
    result = np.ones(len(values))
    minimum = 1.
    for rank in range(len(values), 0, -1):
        ii = order[rank-1]
        minimum = min(minimum, values[ii]*len(values)/rank)
        result[ii] = minimum
    return result


def probability(frame, columns, labels):
    keyed = {}
    for i, key in enumerate(frame[columns].itertuples(index=False, name=None)):
        keyed.setdefault(key, []).append(i)
    out = np.zeros(labels.shape)
    for indices in keyed.values():
        for j in range(labels.shape[1]):
            out[indices, j] = float(np.count_nonzero(labels[indices, j]))/len(indices)
    return out


def bin_covariate(values, number):
    out = pd.Series(-1, index=values.index)
    valid = values.notna()
    if valid.any():
        out.loc[valid] = pd.qcut(values[valid], number, labels=False, duplicates="drop").fillna(0)
    return out.astype(int)


def independent_strata(f, names, threshold, model, side):
    z = pd.DataFrame(index=f.index)
    if side == "h":
        for name in names:
            if name+"_observation_signature" in f and f[name+"_observation_signature"].notna().any():
                z[name+"_sig"] = f[name+"_observation_signature"]
            for direction in ["proximal", "distal"]:
                z[name+direction] = f[name+"_capacity_"+direction] >= threshold
    else:
        z["proximal"] = f.m_capacity_proximal >= threshold
        z["distal"] = f.m_capacity_distal >= threshold
    for cov, number in [("PAS", 3), ("coverage", 3), ("baseline", 5 if "fine" in model else 3)]:
        z[cov] = bin_covariate(f[side+"_"+cov], number)
    if side == "m" and "fine" in model or side == "h" and "human_expression" in model:
        z["GEX"] = bin_covariate(f[side+"_GEX"], 3)
    return z


def main():
    manifest = json.loads((HERE/"manifest.json").read_text())
    for name, expected in manifest["input_sha256"].items():
        assert file_hash(ROOT/name) == expected, name
    for name, expected in manifest["output_sha256"].items():
        assert file_hash(HERE/name) == expected, name
    protocol = json.loads((HERE/"protocol.json").read_text())
    result = pd.read_csv(HERE/"overlap_results.tsv", sep="\t", float_precision="round_trip")
    direct = pd.read_csv(HERE/"direct_specificity_results.tsv", sep="\t", float_precision="round_trip")
    members = pd.read_csv(HERE/"joint_membership_and_covariates.tsv.gz", sep="\t", float_precision="round_trip", low_memory=False)
    draws = np.load(HERE/"joint_null_draws.npz")
    verified = 0
    for context, q in result.groupby(["block", "mouse_target", "threshold", "model"], sort=False):
        block, target, threshold, model = context
        names = protocol["blocks"][block]
        f = members[members.block.eq(block) & members.mouse_target.eq(target) & members.threshold.eq(threshold)].reset_index(drop=True)
        assert f.human_gene_id.is_unique and f.mouse_gene_id.is_unique
        assert all(f[name+"_proximal"].notna().all() for name in names)
        cols = [name+"_"+direction for name in names for direction in ["proximal", "distal"]]
        h = f[cols].to_numpy(bool)
        m = f[["proximal", "distal"]].to_numpy(bool)
        hz = independent_strata(f, names, threshold, model, "h")
        mz = independent_strata(f, names, threshold, model, "m")
        hp = probability(hz, list(hz), h)
        mp = probability(mz, list(mz), m)
        expected = np.array([[sum(hp[:, i]*mp[:, j]) for j in range(2)] for i in range(len(cols))])
        null = draws["|".join(map(str, context))]
        assert len(null) == 20000
        for r in q.itertuples():
            i = cols.index(r.human_label+"_"+r.human_direction)
            j = ["proximal", "distal"].index(r.mouse_direction)
            observed = int(np.count_nonzero(h[:, i] & m[:, j]))
            assert observed == r.observed
            assert abs(expected[i, j]-r.analytic_expected) < 1e-10
            tail = int(np.count_nonzero(null[:, i, j] >= observed))
            assert tail == r.tail_count
            p = (1+tail)/20001 if null[:, i, j].std() > 0 else 1.
            assert abs(p-r.p) < 1e-12
            assert abs(null[:, i, j].mean()-r.null_mean) < 1e-12
            assert abs(((hp[:, i]>0)&(hp[:, i]<1)).mean()-r.human_movable_fraction) < 1e-12
            assert abs(((mp[:, j]>0)&(mp[:, j]<1)).mean()-r.mouse_movable_fraction) < 1e-12
            verified += 1
        d = direct[direct.block.eq(block) & direct.mouse_target.eq(target) & direct.threshold.eq(threshold) & direct.model.eq(model)]
        for r in d.itertuples():
            i = cols.index(r.control+"_proximal")
            enaive, econtrol = expected[0, 0], expected[i, 0]
            if enaive and econtrol:
                observed = np.count_nonzero(h[:, 0]&m[:, 0])/enaive-np.count_nonzero(h[:, i]&m[:, 0])/econtrol
                stat = null[:, 0, 0]/enaive-null[:, i, 0]/econtrol
                assert abs(observed-r.fold_difference) < 1e-10
                tail = int(np.count_nonzero(stat >= observed-1e-12))
                p = (1+tail)/20001 if stat.std() > 1e-12 else 1.
            else:
                p, tail = 1., 20000
            assert tail == r.tail_count and abs(p-r.p) < 1e-12
            verified += 1
    for _, q in direct.groupby(["model", "threshold"]):
        assert len(q) == 14
        np.testing.assert_allclose(q.q_specificity14_within_model_threshold, adjust(q.p), atol=1e-12)
    np.testing.assert_allclose(direct.q_all_specificity_exploratory, adjust(direct.p), atol=1e-12)
    for _, q in result.groupby(["model", "threshold"]):
        np.testing.assert_allclose(q.q_all_overlaps_within_model_threshold, adjust(q.p), atol=1e-12)
        d = q[q.direction_control]
        assert len(d) == 4
        np.testing.assert_allclose(d.q_direction4_within_model_threshold, adjust(d.p), atol=1e-12)
    np.testing.assert_allclose(result.q_all_overlaps_exploratory, adjust(result.p), atol=1e-12)
    pilot = pd.read_csv(HERE/"before-observation-signature-review/overlap_results.tsv", sep="\t", float_precision="round_trip")
    lhs = result[~result.block.eq("human_recurrent_controls")].reset_index(drop=True)
    rhs = pilot[~pilot.block.eq("human_recurrent_controls")].reset_index(drop=True)
    pd.testing.assert_frame_equal(lhs.drop(columns=["q_all_overlaps_within_model_threshold", "q_all_overlaps_exploratory"]),
                                  rhs.drop(columns=["q_all_overlaps_within_model_threshold", "q_all_overlaps_exploratory"]))
    pilot_draws = np.load(HERE/"before-observation-signature-review/joint_null_draws.npz")
    unchanged_draw_contexts = 0
    for key in draws.files:
        if not key.startswith("human_recurrent_controls|"):
            np.testing.assert_array_equal(draws[key], pilot_draws[key])
            unchanged_draw_contexts += 1
    report = dict(status="PASS", arithmetic_checks=verified, overlap_rows=len(result), direct_rows=len(direct),
        null_contexts=len(draws.files), input_hashes=len(manifest["input_sha256"]), output_hashes=len(manifest["output_sha256"]),
        analytic_expectation_independently_rebuilt=True, finite_p_and_BH_replayed=True,
        main_Liu_and_direction_null_unchanged_by_signature_repair=True,
        exactly_unchanged_draw_contexts=unchanged_draw_contexts,
        skipped=0, verifier_sha256=file_hash(Path(__file__)))
    (HERE/"validation.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
