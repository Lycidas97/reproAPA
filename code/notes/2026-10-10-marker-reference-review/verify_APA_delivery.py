"""Independent raw-profile, contributor-weight and gradient verification.

Profile sampling uses the minimum SHA-256 gene key in every
study/setting/cohort/status stratum, without reading APA effect magnitudes.
"""
from pathlib import Path
import csv
import gzip
import hashlib
import json

import numpy as np
import pandas as pd

from verify_marker_delivery import SOURCE, HERE, STUDIES, records, sha


KEYS = ["setting", "version", "study_id", "route", "gene_id", "contrast"]


def read(name, **kwargs):
    return pd.read_csv(SOURCE / name, sep="\t", float_precision="round_trip", **kwargs)


def source_weights(frame):
    """Same-library date/unit/capture hierarchy, counted independently."""
    assert frame.library_key_a.eq(frame.library_key_b).all()
    assert frame.time_a.astype(str).eq(frame.time_b.astype(str)).all()
    assert frame.sample_unit_a.astype(str).eq(frame.sample_unit_b.astype(str)).all()
    if not len(frame):
        return np.zeros(0)
    dates = frame.groupby(KEYS, dropna=False).time_a.transform("nunique")
    units = frame.groupby(KEYS + ["time_a"], dropna=False).sample_unit_a.transform("nunique")
    captures = frame.groupby(KEYS + ["time_a", "sample_unit_a"], dropna=False).library_key.transform("size")
    return (1 / dates / units / captures).to_numpy()


def check_weights_and_identities():
    pairs = read("pair_effects.tsv.gz")
    cellkeys = KEYS + ["library_key"]
    assert not pairs.duplicated(cellkeys).any()
    full = source_weights(pairs)
    assert np.allclose(full, pairs.full_weight, rtol=0, atol=1e-14)
    libkeys = ["setting", "version", "study_id", "route", "gene_id", "library_key"]
    count = pairs.groupby(libkeys).contrast.transform("nunique")
    common = count.eq(4)
    assert np.array_equal(common, pairs.common_four)
    matched = pairs[common].copy()
    assert np.allclose(source_weights(matched), matched.common_weight, rtol=0, atol=1e-14)
    assert pairs.loc[~common, "common_weight"].isna().all()
    delta = matched.pivot(index=libkeys, columns="contrast", values="delta_DWUI")
    assert set(delta.columns) == {"D00", "D01", "D10", "D11"}
    identity = delta.D11 - delta.D00 - (delta.D01 - delta.D00) - (delta.D10 - delta.D00)
    assert np.abs(identity).max() < 1e-12
    assert matched.groupby(libkeys).common_weight.nunique().eq(1).all()
    baseline = read("baseline_pair_effects.tsv.gz")
    basekeys = ["study_id", "route", "gene_id", "library_key_a", "library_key_b"]
    baseline = baseline.set_index(basekeys)
    d00 = pairs[pairs.setting.eq("supported50") & pairs.contrast.eq("D00")]
    for version, group in d00.groupby("version"):
        original = group.set_index(basekeys)
        assert set(original.index) == set(baseline.index)
        reference = baseline.loc[original.index]
        assert np.allclose(original.delta_DWUI, reference.delta_DWUI, rtol=0, atol=1e-12)
        assert np.allclose(original.full_weight, reference.weight, rtol=0, atol=1e-14)
    bridge = read("bridge_pair_weights.tsv.gz")
    for background, contrasts in [("common_M_bridge", {"D00", "D01"}), ("common_P_bridge", {"D00", "D10"})]:
        selected = pairs[pairs.contrast.isin(contrasts)].copy()
        selected = selected[selected.groupby(libkeys).contrast.transform("nunique").eq(2)]
        actual = bridge[bridge.background.eq(background)].set_index(cellkeys)
        expected = selected.set_index(cellkeys)
        assert actual.index.is_unique and set(actual.index) == set(expected.index)
        expected["independent_weight"] = source_weights(selected)
        assert np.allclose(actual.loc[expected.index].shared_weight, expected.independent_weight, rtol=0, atol=1e-14)
    gradient = read("gradient_pair_effects.tsv.gz")
    g = gradient.assign(setting="supported50", version=gradient.panel, contrast="gradient")
    assert np.allclose(source_weights(g), g.weight, rtol=0, atol=1e-14)
    assert np.allclose(g.delta_high - g.delta_low, g.direct_M_low_minus_high, rtol=0, atol=1e-12)
    assert g.identity_abs_error.max() < 1e-12
    for row in gradient.to_dict("records"):
        p = np.array(json.loads(row["P_usage"]))
        positions = np.array(json.loads(row["position"]))
        for layer in ["low", "middle", "high", "whole"]:
            m = np.array(json.loads(row["M_" + layer + "_usage"]))
            assert abs((p - m) @ positions - row["delta_" + layer]) < 1e-12
    return dict(status="PASS", measurable_pairs=len(pairs), common_four_pairs=len(matched),
                common_four_library_gene_identities=len(delta), bridge_weight_rows=len(bridge),
                gradient_pairs=len(gradient), original_pair_replay_rows=len(d00))


def independent_cohorts(archive, classifications):
    groups = {}
    meta_by_study = {}
    for study in STUDIES:
        meta = records(archive, study + "_meta_json")
        meta_by_study[study] = meta
        endpoint = records(archive, study + "_endpoint_json")
        design = endpoint[endpoint.contrast.eq("P_M") & endpoint.design.eq("same_library")]
        oldkeys = set(zip(design.route, design.library_key_a))
        for (route, library), part in meta[meta.in_planned_route].groupby(["route", "library_key"]):
            old = part[part.original_selected & part.original_purity.ge(.8)]
            for node in ["P", "M"]:
                if (route, library) in oldkeys:
                    groups[(study, route, library, "original_" + node)] = old.index[old.original_node.eq(node)].to_numpy()
            classes = classifications[classifications.study_id.eq(study) & classifications.route.eq(route) & classifications.library_key.eq(library)]
            for version, cls in classes.groupby("version"):
                for label in ["Pstar", "Mstar"]:
                    ids = cls.loc[cls.expected_label.eq(label), "metacell_id"]
                    groups[(study, route, library, version + "_" + label)] = meta.index[meta.metacell_id.isin(ids)].to_numpy()
            old_m = old[old.original_node.eq("M")]
            for panel, version in [("core", "core_third_strict"), ("native", "native_third_strict")]:
                score = classes[classes.version.eq(version)].set_index("metacell_id").P_score.reindex(old_m.metacell_id).to_numpy()
                low, high = np.quantile(score, [1/3, 2/3], method="linear") if len(score) else [np.nan, np.nan]
                masks = {"low": score <= low, "middle": (score > low) & (score < high), "high": score >= high}
                for layer, mask in masks.items():
                    groups[(study, route, library, "gradient_" + panel + "_" + layer)] = old_m.index.to_numpy()[mask]
    return groups, meta_by_study


def check_axes_and_pre_gene_coverage(archive, groups, metas, classifications, coverage):
    """Check complete native axes and reasons before any gene-specific gate."""
    axes = read("primary_axis_features.tsv.gz").set_index(["study_id", "atlas_pas_id"])
    background = read("gene_background.tsv.gz").set_index(["study_id", "setting", "gene_id"])
    assert axes.index.is_unique and background.index.is_unique
    expected_axis_keys, expected_background_keys, old_libraries = set(), set(), set()
    for study in STUDIES:
        features = records(archive, study + "_features_json")
        for row in features.to_dict("records"):
            key = (study, row["atlas_pas_id"]); expected_axis_keys.add(key)
            actual = axes.loc[key]
            for field in ["primary_gene_id", "feature_index_1based", "chrom", "strand", "cut_site", "coordinate_eligible", "matched_source_count_50nt"]:
                assert actual[field] == row[field], (key, field)
            assert np.isclose(actual.position, row["position"], rtol=0, atol=1e-14, equal_nan=True)
        for setting in ["supported50", "expanded_primary"]:
            retained = features[features.coordinate_eligible]
            if setting == "supported50":
                retained = retained[retained.matched_source_count_50nt.ge(1)]
            retained_counts = retained.groupby("primary_gene_id").size()
            for gene, part in features.groupby("primary_gene_id"):
                key = (study, setting, gene); expected_background_keys.add(key)
                actual = background.loc[key]
                n = int(retained_counts.get(gene, 0)); coordinate = bool(part.coordinate_eligible.all())
                assert actual.n_full_primary_PAS == len(part) and actual.n_retained_PAS == n
                assert actual.coordinate_eligible == coordinate
                assert actual.status == ("measurable_axis_requalification_required" if n >= 2 else "not_measurable")
                reason = "" if n >= 2 else "unusable_primary_axis" if not coordinate else "fewer_than_two_retained_PAS"
                assert ("" if pd.isna(actual.reason) else actual.reason) == reason
        endpoint = records(archive, study + "_endpoint_json")
        design = endpoint[endpoint.contrast.eq("P_M") & endpoint.design.eq("same_library")]
        old_libraries.update((study, r.route, r.library_key_a) for r in design.itertuples())
    assert set(axes.index) == expected_axis_keys
    assert set(background.index) == expected_background_keys
    gradient = read("gradient_library_eligibility.tsv").set_index(["study_id", "route", "library_key", "panel"])
    assert gradient.index.is_unique and len(gradient) == 138 * 2
    expected_flags = {}
    for key, members in groups.items():
        study, route, library, cohort = key
        sufficient = len(members) >= 3 and metas[study].n_cells.iloc[members].sum() >= 30
        if not cohort.startswith("gradient_"):
            expected_flags[key] = sufficient
    for key, row in gradient.iterrows():
        study, route, library, panel = key
        version = "core_third_strict" if panel == "core" else "native_third_strict"
        cls = classifications[classifications.study_id.eq(study) & classifications.route.eq(route) & classifications.library_key.eq(library) & classifications.version.eq(version)]
        meta = metas[study]
        old_m = meta[meta.route.eq(route) & meta.library_key.eq(library) & meta.original_selected & meta.original_purity.ge(.8) & meta.original_node.eq("M")]
        scores = cls.set_index("metacell_id").P_score.reindex(old_m.metacell_id).to_numpy()
        cuts = np.quantile(scores, [1/3, 2/3], method="linear") if len(scores) else [np.nan, np.nan]
        assert np.allclose([row.q_low, row.q_high], cuts, rtol=0, atol=1e-12, equal_nan=True)
        assert row.n_original_M == len(old_m) and row.cells_original_M == old_m.n_cells.sum()
        why = "original_same_library_design_absent" if (study, route, library) not in old_libraries else ""
        if int(cls.P_varying.iloc[0]) < 2:
            why = why or "P_programme_unavailable"
        if len(old_m) < 3:
            why = why or "fewer_than_three_original_M"
        if not np.isfinite(cuts).all() or cuts[1] <= cuts[0]:
            why = why or "M_programme_quantiles_not_separated"
        for layer in ["low", "middle", "high"]:
            members = groups[(study, route, library, "gradient_" + panel + "_" + layer)]
            cells = meta.n_cells.iloc[members].sum()
            assert row["n_" + layer] == len(members) and row["cells_" + layer] == cells
            if len(members) < 3 or cells < 30:
                why = why or "M_layer_below_3MC_30cells"
        if not expected_flags.get((study, route, library, "original_P"), False):
            why = why or "original_P_below_3MC_30cells"
        assert ("" if pd.isna(row.reason) else row.reason) == why, key
        assert row.status == ("not_estimable" if why else "pre_gene_eligible"), key
        for layer in ["low", "middle", "high"]:
            expected_flags[(study, route, library, "gradient_" + panel + "_" + layer)] = not why
    actual_flags = coverage.set_index(["study_id", "route", "library_key", "cohort"]).pre_gene_eligible
    assert set(actual_flags.index) == set(expected_flags)
    assert all(actual_flags.loc[key] == value for key, value in expected_flags.items())
    eligibility = read("library_apa_eligibility.tsv")
    context_keys = ["study_id", "route", "library_key", "setting", "version", "contrast"]
    assert not eligibility.duplicated(context_keys).any() and len(eligibility) == 138 * 5 * 4
    for row in eligibility.to_dict("records"):
        lib = (row["study_id"], row["route"], row["library_key"])
        permitted = row["contrast"] == "D11" or lib in old_libraries
        a, b = lib + (row["target_cohort"],), lib + (row["reference_cohort"],)
        why = "original_same_library_design_absent" if not permitted else "" if expected_flags.get(a, False) and expected_flags.get(b, False) else "cohort_absent_or_below_pre_gene_3MC_30cells"
        assert ("" if pd.isna(row["reason"]) else row["reason"]) == why
        assert row["status"] == ("not_estimable" if why else "pre_gene_eligible")
        for side, group in [("a", a), ("b", b)]:
            members = groups.get(group, np.zeros(0, dtype=int))
            assert row["pre_gene_n_mc_" + side] == len(members)
            assert row["pre_gene_n_cells_" + side] == metas[lib[0]].n_cells.iloc[members].sum()
    return dict(status="PASS", full_primary_axis_rows=len(axes), complete_gene_background_rows=len(background),
                library_comparison_contexts=len(eligibility), gradient_library_contexts=len(gradient))


def select_pair_qualification_samples():
    """Check the full eligibility grid; sample reasons independently of effects."""
    libraries = read("library_apa_eligibility.tsv")
    contexts = libraries[libraries.status.eq("pre_gene_eligible")]
    fields = ["study_id", "route", "library_key", "setting", "version", "contrast"]
    background = read("gene_background.tsv.gz")
    usable = background[background.n_retained_PAS.ge(2)]
    hashes = {gene: int.from_bytes(hashlib.sha256(gene.encode()).digest()[:16], "big") for gene in usable.gene_id.unique()}
    expected = {}
    allowed_genes = {}
    for (study, setting), part in usable.groupby(["study_id", "setting"]):
        genes = set(part.gene_id); allowed_genes[(study, setting)] = genes
    for row in contexts.to_dict("records"):
        key = tuple(row[f] for f in fields)
        genes = allowed_genes[(row["study_id"], row["setting"])]
        expected[key] = (len(genes), sum(hashes[gene] for gene in genes))
    seen = {key: [0, 0] for key in expected}
    chosen = {}
    total = 0
    with gzip.open(SOURCE / "pair_eligibility.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            key = tuple(row[f] for f in fields); gene = row["gene_id"]
            assert key in seen and gene in allowed_genes[(row["study_id"], row["setting"])]
            seen[key][0] += 1; seen[key][1] += hashes[gene]; total += 1
            stratum = tuple(row[f] for f in ["study_id", "setting", "version", "contrast", "status", "reason"])
            signature = hashlib.sha256(("|".join(key) + "|" + gene).encode()).hexdigest()
            if stratum not in chosen or signature < chosen[stratum][0]:
                chosen[stratum] = (signature, row)
    assert all(tuple(seen[key]) == values for key, values in expected.items())
    return [row for _, row in chosen.values()], dict(status="PASS", eligibility_rows=total,
        complete_pre_gene_contexts=len(expected), raw_pair_qualification_samples=len(chosen),
        completeness_check="count and SHA256-derived 128-bit sum for every complete context gene multiset")


def check_profiles(archive, groups, metas, pair_samples):
    chosen = {}
    needed_ids = {row[side + "_group"] for row in pair_samples for side in ["target", "reference"]}
    needed = {}
    total = 0
    with gzip.open(SOURCE / "gene_group_profiles.tsv.gz", "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            total += 1
            stratum = (row["study_id"], row["setting"], row["cohort"], row["status"])
            signature = hashlib.sha256(row["group_id"].encode()).hexdigest()
            if stratum not in chosen or signature < chosen[stratum][0]:
                chosen[stratum] = (signature, row)
            if row["group_id"] in needed_ids:
                assert row["group_id"] not in needed
                needed[row["group_id"]] = row
    assert set(needed) == needed_ids
    rows_to_check = {row["group_id"]: row for _, row in chosen.values()}
    rows_to_check.update(needed)
    checks = []
    for study in STUDIES:
        meta = metas[study]
        features = records(archive, study + "_features_json")
        features["native_column"] = np.arange(len(features))
        data = archive[study + "_pas_data"]
        indices = archive[study + "_pas_indices"]
        pointers = archive[study + "_pas_indptr"]
        for row in [v for v in rows_to_check.values() if v["study_id"] == study]:
            sites = features[features.primary_gene_id.eq(row["gene_id"]) & features.coordinate_eligible]
            if row["setting"] == "supported50":
                sites = sites[sites.matched_source_count_50nt.ge(1)]
            assert json.loads(row["PAS_ids"]) == sites.atlas_pas_id.tolist()
            assert np.allclose(json.loads(row["position"]), sites.position, rtol=0, atol=1e-14)
            selected = groups[(study, row["route"], row["library_key"], row["cohort"])]
            full = np.zeros((len(meta), len(sites)), dtype=data.dtype)
            for j, column in enumerate(sites.native_column):
                start, end = pointers[column:column+2]
                full[indices[start:end], j] = data[start:end]
            counts = full[selected]
            qualified = selected[counts.sum(axis=1) >= 20]
            counts = full[qualified]
            assert int(row["n_pre_gene_mc"]) == len(selected)
            assert int(row["n_pre_gene_cells"]) == meta.n_cells.iloc[selected].sum()
            assert int(row["n_mc"]) == len(qualified)
            assert int(row["n_cells"]) == meta.n_cells.iloc[qualified].sum()
            assert json.loads(row["qualified_MC_ids"]) == meta.metacell_id.iloc[qualified].tolist()
            assert np.array_equal(json.loads(row["qualified_MC_pas_depths"]), counts.sum(axis=1))
            assert np.array_equal(json.loads(row["pooled_molecules"]), counts.sum(axis=0))
            assert np.array_equal(json.loads(row["nonzero_MC"]), (counts > 0).sum(axis=0))
            expected_reason = "insufficient_gene_qualified_metacells" if len(qualified) < 3 else "insufficient_gene_qualified_cells" if meta.n_cells.iloc[qualified].sum() < 30 else ""
            assert row["reason"] == expected_reason
            assert row["status"] == ("not_measurable" if expected_reason else "available")
            error = 0.
            if len(qualified):
                usage = (counts / counts.sum(axis=1)[:, None]).mean(axis=0)
                error = float(np.max(np.abs(usage - np.array(json.loads(row["usage"])))))
                assert error < 1e-12
            else:
                assert row["usage"] == "NA"
            checks.append(dict(study_id=study, group_id=row["group_id"], setting=row["setting"],
                               cohort=row["cohort"], status=row["status"], n_mc=len(qualified), usage_abs_error=error))
        del data, indices, pointers
    for row in pair_samples:
        a, b = needed[row["target_group"]], needed[row["reference_group"]]
        assert a["PAS_ids"] == b["PAS_ids"] and a["position"] == b["position"]
        ma, mb = np.array(json.loads(a["pooled_molecules"])), np.array(json.loads(b["pooled_molecules"]))
        na, nb = np.array(json.loads(a["nonzero_MC"])), np.array(json.loads(b["nonzero_MC"]))
        joint = int(((ma + mb >= 10) & (na + nb >= 2)).sum())
        reason = "insufficient_gene_qualified_metacells" if min(int(a["n_mc"]), int(b["n_mc"])) < 3 else "insufficient_gene_qualified_cells" if min(int(a["n_cells"]), int(b["n_cells"])) < 30 else "fewer_than_two_joint_supported_PAS" if joint < 2 else ""
        assert row["reason"] == reason and row["status"] == ("not_measurable" if reason else "eligible")
        assert int(row["n_joint_PAS"]) == joint
        for side, profile in [("a", a), ("b", b)]:
            assert int(row["n_mc_" + side]) == int(profile["n_mc"])
            assert int(row["n_cells_" + side]) == int(profile["n_cells"])
    pd.DataFrame(checks).to_csv(HERE / "raw_profile_checks.tsv", sep="\t", index=False)
    return dict(status="PASS", source_profile_rows=total, independently_reconstructed_profiles=len(checks),
                fixed_sampling="minimum SHA256 group key per study/setting/cohort/status and pair key per study/setting/version/contrast/status/reason",
                usage_max_abs_difference=max((r["usage_abs_error"] for r in checks), default=0.))


def main():
    acceptance = json.loads((SOURCE / "analysis_validation.json").read_text())
    assert acceptance["status"] == "COMPLETE_AWAITING_INDEPENDENT_QA"
    classification = json.loads((HERE / "delivery_validation.json").read_text())
    assert classification["classification_check"]["status"] == "PASS"
    archive = np.load(SOURCE / "inputs.npz", allow_pickle=False)
    assert sha(SOURCE / "inputs.npz") == classification["inputs_sha256"]
    expected = pd.read_csv(HERE / "expected_marker_classification.tsv.gz", sep="\t", float_precision="round_trip")
    groups, metas = independent_cohorts(archive, expected)
    coverage = read("cohort_coverage.tsv")
    assert len(coverage) == len(groups)
    for row in coverage.to_dict("records"):
        key = (row["study_id"], row["route"], row["library_key"], row["cohort"])
        indices = groups[key]; meta = metas[row["study_id"]]
        assert json.loads(row["candidate_MC_ids"]) == meta.metacell_id.iloc[indices].tolist()
        assert row["n_mc"] == len(indices) and row["n_cells"] == meta.n_cells.iloc[indices].sum()
    print("COHORT_MEMBERSHIP_PASS", len(groups), flush=True)
    qualifications = check_axes_and_pre_gene_coverage(archive, groups, metas, expected, coverage)
    print("AXES_AND_PRE_GENE_ELIGIBILITY_PASS", json.dumps(qualifications), flush=True)
    weight_checks = check_weights_and_identities()
    print("WEIGHT_AND_IDENTITIES_PASS", json.dumps(weight_checks), flush=True)
    pair_samples, eligibility_checks = select_pair_qualification_samples()
    print("COMPLETE_GENE_ELIGIBILITY_PASS", json.dumps(eligibility_checks), flush=True)
    profiles = check_profiles(archive, groups, metas, pair_samples)
    report = dict(status="PASS", marker_classification=classification["classification_check"],
                  cohort_memberships=len(groups), complete_axes_and_pre_gene_eligibility=qualifications,
                  contributor_weights=weight_checks, complete_gene_eligibility=eligibility_checks, raw_profile_checks=profiles,
                  source_analyses=acceptance, verifier_sha256=sha(__file__))
    (HERE / "apa_delivery_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
