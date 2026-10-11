"""Five synchronized metric tests on the frozen WOT independent-well cohort."""
from pathlib import Path
import hashlib
import importlib.util
import json
import platform
import time

import numpy as np
import pandas as pd
import scipy

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
METRICS = ["DWUI", "RWUI", "PDUI", "PPUI", "MPRO"]
INPUTS = {}


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def module(name, path):
    INPUTS[str(path)] = sha(path)
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


EXTRACT = module("frozen_metric_extract", PARENT / "extract.py")
CORE = module("frozen_well_inference", PARENT / "inference/analyze.py")


def read(path, **kwargs):
    INPUTS[str(path)] = sha(path)
    return pd.read_csv(path, sep="\t", **kwargs)


def save(frame, name):
    frame.to_csv(HERE / name, sep="\t", index=False, na_rep="NA", float_format="%.17g",
                 compression={"method": "gzip", "mtime": 0} if name.endswith(".gz") else None)


def statistics(difference, weights, metric_weights):
    """Apply synchronized signed-vector nulls before the nonlinear MPRO statistic."""
    difference = np.asarray(difference, float)
    weights = np.asarray(weights, float)
    assert difference.ndim == 2 and difference.shape[0] == len(weights)
    assert np.allclose(difference.sum(axis=1), 0, atol=1e-12)
    assert abs(weights.sum() - 1) < 1e-12
    scalar = difference @ np.asarray(metric_weights, float)
    weighted_scalar = scalar * weights[:, None]
    denom = np.sqrt(np.square(weighted_scalar).sum(axis=0))
    standardized = np.divide(weighted_scalar, denom, out=np.zeros_like(weighted_scalar), where=denom > 1e-14)
    signs, exact = CORE.signs(len(weights))
    observed_linear = np.abs(standardized.sum(axis=0))
    null_linear = np.abs(signs @ standardized)
    weighted_vector = difference * weights[:, None]
    mean_vector = weighted_vector.sum(axis=0)
    observed_mpro = np.ptp(mean_vector)
    null_vector = signs @ weighted_vector
    null_mpro = np.ptp(null_vector, axis=1)
    observed = np.r_[observed_linear, observed_mpro]
    extreme = np.r_[(null_linear >= observed_linear - 1e-10).sum(axis=0),
                    (null_mpro >= observed_mpro - 1e-10).sum()]
    p = extreme / len(signs) if exact else (extreme + 1) / (len(signs) + 1)
    effects = np.r_[weighted_scalar.sum(axis=0), observed_mpro]
    return effects, p, observed, dict(exact_sign_enumeration=bool(exact), sign_patterns=len(signs)), mean_vector


def extract_frozen():
    legacy = read(PARENT / "inference/gene_tests.tsv.gz")
    pairs = read(PARENT / "inference/qualified_well_pairs.tsv.gz")
    pairs = pairs[pairs.gene_id.isin(legacy.gene_id) & pairs.replicate_complete_date].copy()
    pairs = pairs.sort_values(["gene_id", "time", "library_key"]).reset_index(drop=True)
    assert len(legacy) == 268 and len(pairs) == 3566 and legacy.n_wells.sum() == 3566
    assert not pairs.duplicated(["gene_id", "library_key"]).any()
    assert pairs.replicate_status.eq("biological_independent").all() and pairs.condition.eq("2i").all()
    original = read(PARENT / "pair_effects.tsv.gz", low_memory=False)
    original = original[original.study_id.eq("GSE115943") & original.contrast.eq("P_M")].copy()
    selected = original.merge(pairs[["gene_id", "library_key"]], left_on=["gene_id", "library_key_a"],
                              right_on=["gene_id", "library_key"], validate="one_to_one")
    assert len(selected) == 3566 and selected.library_key_a.eq(selected.library_key_b).all()
    vectors, states, features, groups, axis, error = EXTRACT.study_vectors(CORE.STUDY, selected)
    old_axis = read(PARENT / "inference/supported50_PAS_axes.tsv.gz")
    rows, data, matrices, max_pair_error = [], {}, {}, 0.
    for gene, ix in groups.items():
        feature = features.iloc[ix]
        old = old_axis[old_axis.primary_gene_id.eq(gene)]
        assert feature.atlas_pas_id.tolist() == old.atlas_pas_id.tolist()
        assert np.allclose(feature.WPP, old.position, rtol=0, atol=1e-12)
        matrices[gene] = feature[["WPP", "RWUI", "PDUI", "PPUI"]].to_numpy()
    for row in pairs.itertuples():
        a = vectors["WOT", row.library_key, "P", row.gene_id]
        b = vectors["WOT", row.library_key, "M", row.gene_id]
        assert np.array_equal(a[1], np.array(json.loads(row.P_site_counts)))
        assert np.array_equal(b[1], np.array(json.loads(row.M_site_counts)))
        assert np.array_equal(a[2], np.array(json.loads(row.P_site_mc)))
        assert np.array_equal(b[2], np.array(json.loads(row.M_site_mc)))
        difference = a[0] - b[0]
        scalar = difference @ matrices[row.gene_id]
        max_pair_error = max(max_pair_error, abs(scalar[0] - row.equal_MC_WPP_delta))
        data[row.gene_id, row.library_key] = difference
        record = dict(gene_id=row.gene_id, gene_symbol=row.gene_symbol, library_key=row.library_key,
                      time=row.time, sample_unit=row.sample_unit, date_replicates=row.date_replicates,
                      target_usage=json.dumps(a[0].tolist()), reference_usage=json.dumps(b[0].tolist()),
                      delta_usage=json.dumps(difference.tolist()), n_PAS=len(difference))
        record.update({"delta_" + metric: float(effect) for metric, effect in zip(METRICS[:4], scalar)})
        rows.append(record)
    assert max_pair_error < 1e-8
    save(pd.DataFrame(rows), "frozen_well_vectors.tsv.gz")
    save(features, "frozen_supported50_PAS.tsv.gz")
    save(states, "frozen_library_state_metrics.tsv.gz")
    INPUTS.update({path: value["sha256"] for path, value in EXTRACT.INPUTS.items()})
    return legacy, pairs, data, matrices, dict(max_state_error=error, max_pair_DWUI_error=max_pair_error)


def summaries(tests, legacy):
    rows = []
    cohorts = [("all_268", tests), ("PAS_ge3", tests[tests.n_PAS.ge(3)])]
    for subset, frame in cohorts:
        for metric, group in frame.groupby("metric", sort=False):
            record = dict(subset=subset, metric=metric, genes=len(group),
                          significant_joint=int(group.q_BH_joint.lt(.05).sum()),
                          significant_metric=int(group.q_BH_metric.lt(.05).sum()))
            for threshold in [.10, .20] if metric == "MPRO" else [.02, .05]:
                label = str(threshold).replace(".", "p")
                record["effect_ge_" + label] = int(group.effect.abs().ge(threshold).sum())
                record["joint_q05_effect_ge_" + label] = int((group.q_BH_joint.lt(.05) & group.effect.abs().ge(threshold)).sum())
            rows.append(record)
    save(pd.DataFrame(rows), "detection_summary.tsv")
    overlap = []
    for subset, frame in cohorts:
        calls = {metric: set(frame.loc[frame.metric.eq(metric) & frame.q_BH_joint.lt(.05), "gene_id"])
                 for metric in METRICS}
        for i, first in enumerate(METRICS):
            for second in METRICS[i + 1:]:
                a, b = calls[first], calls[second]
                overlap.append(dict(subset=subset, first=first, second=second, both=len(a & b),
                                    first_only=len(a - b), second_only=len(b - a), union=len(a | b)))
    save(pd.DataFrame(overlap), "significant_set_overlap.tsv")
    wide = tests.pivot(index="gene_id", columns="metric", values="q_BH_joint").add_suffix("_q_BH_joint")
    save(wide.reset_index(), "joint_q_comparison.tsv")
    candidate = read(PARENT / "inference/original_candidate_statistical_support.tsv")
    candidate = candidate[candidate.contrast.eq("P_M")][["gene_id", "direction"]]
    linked = tests.merge(candidate, on="gene_id", how="inner", validate="many_to_one")
    signed_effect = linked.effect * np.where(linked.metric.eq("PPUI"), -1, 1)
    linked["original_direction_retained"] = pd.Series(signed_effect * linked.direction > 0, dtype="boolean")
    linked.loc[linked.metric.eq("MPRO"), "original_direction_retained"] = pd.NA
    save(linked, "original_PM_candidate_intersection.tsv.gz")
    cr = []
    for metric, frame in linked.groupby("metric", sort=False):
        record = dict(metric=metric, original_candidates=216, inference_covered_genes=frame.gene_id.nunique(),
                       significant_joint=int(frame.q_BH_joint.lt(.05).sum()),
                       significant_original_direction=int((frame.q_BH_joint.lt(.05) & frame.original_direction_retained).sum())
                       if metric != "MPRO" else np.nan,
                       direction_interpretation="unsigned; no inherited direction" if metric == "MPRO" else
                       "PPUI sign reversed for original proximal/distal label" if metric == "PPUI" else "proximal/distal sign")
        for threshold in [.10, .20] if metric == "MPRO" else [.02, .05]:
            label = str(threshold).replace(".", "p")
            supported = frame.q_BH_joint.lt(.05) & frame.effect.abs().ge(threshold)
            record["joint_q05_effect_ge_" + label] = int(supported.sum())
            record["joint_q05_original_direction_effect_ge_" + label] = int((supported & frame.original_direction_retained).sum()) if metric != "MPRO" else np.nan
        cr.append(record)
    save(pd.DataFrame(cr), "candidate_detection_summary.tsv")


def main():
    started = time.monotonic()
    # Existing inference outputs are immutable inputs to this supplement.
    protected = {str(path): sha(path) for path in (PARENT / "inference").glob("*")
                 if path.is_file() and path.suffix in [".gz", ".tsv", ".json"]}
    legacy, pairs, data, matrices, reconstruction = extract_frozen()
    rows, effect_errors, p_errors = [], [], []
    for row in legacy.sort_values("gene_id").itertuples():
        wells = pairs[pairs.gene_id.eq(row.gene_id)].sort_values(["time", "library_key"])
        assert len(wells) == row.n_wells and wells.time.nunique() == row.n_dates
        assert ";".join(wells.library_key) == row.libraries
        weights = 1 / row.n_dates / wells.date_replicates.to_numpy()
        difference = np.stack([data[row.gene_id, lib] for lib in wells.library_key])
        effects, p, stats, metadata, mean_vector = statistics(difference, weights, matrices[row.gene_id])
        p_errors.append(abs(p[0] - row.equal_MC_WPP_p))
        effect_errors.append(abs(effects[0] - row.equal_MC_WPP_estimate))
        if max(p_errors) >= 1e-8:
            (HERE / "numerical_conflict.json").write_text(json.dumps(dict(gene_id=row.gene_id,
                original_p=row.equal_MC_WPP_p, reconstructed_p=float(p[0]), maximum_p_error=max(p_errors)), indent=2) + "\n")
            raise AssertionError("DWUI original P recovery failed; see numerical_conflict.json")
        for i, metric in enumerate(METRICS):
            record = dict(gene_id=row.gene_id, gene_symbol=row.gene_symbol, metric=metric,
                          effect=float(effects[i]), p=float(p[i]), statistic=float(stats[i]),
                          n_PAS=row.n_PAS, n_wells=row.n_wells, n_dates=row.n_dates,
                          dates=row.dates, libraries=row.libraries, **metadata,
                          weighted_delta_usage=json.dumps(mean_vector.tolist()))
            for threshold in [.10, .20] if metric == "MPRO" else [.02, .05]:
                record["effect_ge_" + str(threshold).replace(".", "p")] = bool(abs(effects[i]) >= threshold)
            rows.append(record)
        if len(rows) % 250 == 0:
            print("TESTED", len(rows) // 5, "genes", flush=True)
    tests = pd.DataFrame(rows)
    assert len(tests) == 1340
    tests["q_BH_joint"] = CORE.bh(tests.p)
    tests["q_BH_metric"] = tests.groupby("metric").p.transform(CORE.bh)
    tests["joint_q05"] = tests.q_BH_joint.lt(.05)
    for threshold in [.02, .05, .10, .20]:
        label = str(threshold).replace(".", "p")
        flag = "effect_ge_" + label
        tests["joint_q05_effect_ge_" + label] = tests.joint_q05 & tests[flag].fillna(False).astype(bool)
    save(tests, "gene_metric_tests.tsv.gz")
    summaries(tests, legacy)
    assert protected == {path: sha(path) for path in protected}
    INPUTS[str(HERE / "PROTOCOL.zh.md")] = sha(HERE / "PROTOCOL.zh.md")
    INPUTS[str(Path(__file__))] = sha(__file__)
    audit = dict(status="PASS", frozen_genes=268, frozen_gene_well_records=3566,
                 metric_tests=len(tests), joint_BH_family=1340, per_metric_BH_family=268,
                 maximum_DWUI_p_error=max(p_errors), maximum_DWUI_effect_error=max(effect_errors),
                 reconstruction=reconstruction, synchronized_signs=True,
                 nonnegative_per_well_MPRO_tested=False, seed=CORE.SEED, random_patterns=CORE.B,
                 inputs_sha256=INPUTS, existing_inference_outputs_sha256=protected,
                 versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
                 elapsed_seconds=time.monotonic() - started)
    (HERE / "audit.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({key: value for key, value in audit.items() if key not in ["inputs_sha256", "existing_inference_outputs_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
