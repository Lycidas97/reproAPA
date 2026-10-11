"""Independently rebuild the locked marker classifications from integer GEX.

This verifier does not import production analysis code. APA and effect checks
are added only against explicit delivered schemas, never through its helpers.
"""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "notes/2026-10-10-marker-reference-sensitivity"
HERE = Path(__file__).resolve().parent
LOCK = "1fc9ff8f14ee74302e032b52d4a6a68b6e89463b73f5bc457fa5ece266c75040"
STUDIES = {
    "GSE103221": ("mouse", {"CIP"}),
    "GSE115943": ("mouse", {"WOT"}),
    "GSE118258": ("human", {"GSE118258"}),
    "GSE150311": ("human", {"GSE150311"}),
    "GSE178325": ("human", {"178325_0618", "178325_0809"}),
    "GSE221739": ("human", {"GSE221739"}),
}
P_CORE = ["DPPA2", "DPPA4", "SALL4", "TDGF1"]
P_NATIVE = {
    "human": ["DPPA4", "DNMT3B", "ESRG", "L1TD1", "TDGF1"],
    "mouse": ["DPPA4", "ESRRB", "ZFP42", "DPPA5A", "TDGF1"],
}
ECM = ["COL1A1", "COL3A1", "DCN", "LUM", "FN1"]
VERSIONS = {
    "core_third_strict": ("core", 1 / 3, 2 / 3, True),
    "core_quartile_strict": ("core", 1 / 4, 3 / 4, True),
    "core_third_relaxed": ("core", 1 / 3, 2 / 3, False),
    "native_third_strict": ("native", 1 / 3, 2 / 3, True),
}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def records(archive, key):
    return pd.DataFrame(json.loads(str(archive[key])))


def marker_score(raw):
    logs = np.log1p(raw)
    spread = logs.std(axis=0, ddof=0)
    varying = spread > 1e-12
    standardized = np.zeros(logs.shape)
    standardized[:, varying] = (
        logs[:, varying] - logs[:, varying].mean(axis=0)
    ) / spread[varying]
    return standardized.mean(axis=1), int(varying.sum())


def reconstruct_classes(archive):
    outputs = []
    input_rows = []
    panels = json.loads(str(archive["panels_json"]))
    assert panels == dict(P_shared=P_CORE, P_legacy=P_NATIVE, ECM=ECM)
    frozen = records(archive, "frozen_states_json").set_index("metacell_id")
    assert frozen.index.is_unique
    for study, (species, routes) in STUDIES.items():
        meta = records(archive, study + "_meta_json")
        selected = records(archive, study + "_selected_json")
        mapping = {name: j for j, name in enumerate(archive[study + "_markers"])}
        counts = archive[study + "_marker_counts"]
        totals = archive[study + "_total_gex"]
        assert meta.metacell_id.is_unique and selected.metacell_id.is_unique
        assert np.issubdtype(counts.dtype, np.integer)
        assert np.all(counts >= 0) and np.all(totals > 0)
        assert np.array_equal(meta.matrix_row, np.arange(len(meta)))
        assert counts.shape == (len(meta), len(mapping))
        assert np.array_equal(meta.total_gex, totals)
        assert np.array_equal(meta.original_selected, meta.metacell_id.isin(selected.metacell_id))
        assert np.array_equal(meta.in_planned_route, meta.route.isin(routes))
        assert set(selected.metacell_id).issubset(meta.metacell_id)
        old = meta.set_index("metacell_id").loc[selected.metacell_id]
        assert np.array_equal(old.original_node, selected.node)
        assert np.allclose(old.original_purity, selected.purity, atol=1e-14, rtol=0)
        assert np.array_equal(old.route, selected.route)
        candidates = meta.loc[meta.in_planned_route].copy()
        labels = frozen.reindex(candidates.metacell_id)
        assert labels.state_group.notna().all()
        assert np.array_equal(labels.library_key, candidates.library_key)
        input_rows.append(dict(
            study_id=study, species=species, native_metacells=len(meta),
            planned_metacells=len(candidates),
            unselected_planned_metacells=int((~candidates.original_selected).sum()),
            old_selected_low_purity_planned_metacells=int(candidates.original_purity.lt(.8).sum()),
            full_frozen_low_purity_planned_metacells=int(labels.state_purity.lt(.8).sum()),
            planned_libraries=int(candidates.library_key.nunique()),
        ))
        for (route, library), part in candidates.groupby(["route", "library_key"], sort=True):
            indices = part.matrix_row.to_numpy()
            cp = counts[indices] / totals[indices, None] * 10000
            ecp = cp[:, [mapping[gene] for gene in ECM]]
            escore, evariable = marker_score(ecp)
            edetected = (ecp >= 1).sum(axis=1)
            for version, (panel, qlo, qhi, opposite_gate) in VERSIONS.items():
                pgenes = P_CORE if panel == "core" else P_NATIVE[species]
                pcp = cp[:, [mapping[gene] for gene in pgenes]]
                pscore, pvariable = marker_score(pcp)
                pdetected = (pcp >= 1).sum(axis=1)
                plo, phi = np.quantile(pscore, [qlo, qhi], method="linear")
                elo, ehi = np.quantile(escore, [qlo, qhi], method="linear")
                classifiable = pvariable >= 2 and evariable >= 2 and phi > plo and ehi > elo
                relative_p = (pscore >= phi) & (escore <= elo)
                relative_m = (escore >= ehi) & (pscore <= plo)
                pabsolute = (pdetected >= 2) & ((edetected <= 1) if opposite_gate else True)
                mabsolute = (edetected >= 2) & ((pdetected <= 1) if opposite_gate else True)
                pstar = relative_p & pabsolute & classifiable
                mstar = relative_m & mabsolute & classifiable
                assert not np.any(pstar & mstar)
                for i, record in enumerate(part.to_dict("records")):
                    outputs.append(dict(
                        study_id=study, species=species, route=route, library_key=library,
                        metacell_id=record["metacell_id"], version=version,
                        expected_label="Pstar" if pstar[i] else "Mstar" if mstar[i] else "unassigned",
                        classifiable=classifiable,
                        P_score=float(pscore[i]), ECM_score=float(escore[i]),
                        P_q_low=float(plo), P_q_high=float(phi),
                        ECM_q_low=float(elo), ECM_q_high=float(ehi),
                        P_varying=pvariable, ECM_varying=evariable,
                        P_detected=int(pdetected[i]), ECM_detected=int(edetected[i]),
                        relative_P=bool(relative_p[i]), relative_M=bool(relative_m[i]),
                        absolute_P=bool(pabsolute[i]), absolute_M=bool(mabsolute[i]),
                        original_selected=bool(record["original_selected"]),
                        original_node=record["original_node"],
                    ))
    return pd.DataFrame(outputs), pd.DataFrame(input_rows)


def reconstruct_axes(archive):
    audit = []
    for study in STUDIES:
        features = records(archive, study + "_features_json")
        assert features.feature_index_1based.is_unique
        # Map the signed coordinate extrema across the entire primary catalog,
        # before any supported50 choice. This is independent of APA responses.
        signed = features.cut_site.to_numpy() * np.where(features.strand.eq("+"), 1, -1)
        columns = pd.DataFrame(dict(
            gene=features.primary_gene_id, signed=signed, cut=features.cut_site,
            chrom=features.chrom, strand=features.strand,
        ))
        bounds = columns.groupby("gene").agg(
            left=("signed", "min"), right=("signed", "max"),
            cuts=("cut", "nunique"), chromosomes=("chrom", "nunique"),
            strands=("strand", "nunique"),
        )
        eligibility = bounds.cuts.ge(2) & bounds.chromosomes.eq(1) & bounds.strands.eq(1)
        expected_eligibility = features.primary_gene_id.map(eligibility).to_numpy()
        assert np.array_equal(expected_eligibility, features.coordinate_eligible)
        low = features.primary_gene_id.map(bounds.left).to_numpy()
        high = features.primary_gene_id.map(bounds.right).to_numpy()
        positions = np.full(len(features), np.nan)
        varied = high > low
        positions[varied] = (signed[varied] - low[varied]) / (high[varied] - low[varied])
        stored = features.position.to_numpy(float)
        assert np.array_equal(np.isnan(positions), np.isnan(stored))
        error = float(np.max(np.abs(positions[varied] - stored[varied])))
        assert error < 1e-14
        audit.append(dict(
            study_id=study, primary_PAS=len(features),
            native_gene_axes=len(bounds), coordinate_eligible_axes=int(eligibility.sum()),
            position_max_abs_difference=error,
        ))
    return audit


def check_classifications(expected, archive):
    path = SOURCE / "metacell_classification.tsv.gz"
    if not path.exists():
        return dict(status="awaiting_production_output")
    delivered = pd.read_csv(path, sep="\t", float_precision="round_trip")
    keys = ["study_id", "route", "library_key", "metacell_id", "version"]
    actual = delivered.set_index(keys)
    wanted = expected.set_index(keys)
    assert actual.index.is_unique and wanted.index.is_unique
    assert set(actual.index) == set(wanted.index), "Classification candidate keys changed"
    actual = actual.loc[wanted.index]
    normalized = actual.classification.replace({"unassigned/mixed": "unassigned"})
    assert np.array_equal(normalized, wanted.expected_label), "Marker membership mismatch"
    assert np.array_equal(actual.classifiable, wanted.classifiable)
    for column, reference in [
        ("P_expressed_markers", "P_detected"),
        ("ECM_expressed_markers", "ECM_detected"),
        ("P_variable_markers", "P_varying"),
        ("ECM_variable_markers", "ECM_varying"),
    ]:
        assert np.array_equal(actual[column], wanted[reference]), column
    numeric = ["P_score", "ECM_score"]
    errors = {}
    for column in numeric:
        delta = np.abs(actual[column] - wanted[column])
        errors[column] = float(delta.max())
        assert errors[column] < 1e-12, column
    # Production deliberately reports NA thresholds when variability is
    # insufficient; those unavailable thresholds are not numerical zeroes.
    classified = wanted.classifiable
    for column in ["P_q_low", "P_q_high", "ECM_q_low", "ECM_q_high"]:
        delta = np.abs(actual.loc[classified, column] - wanted.loc[classified, column])
        errors[column] = float(delta.max())
        assert errors[column] < 1e-12, column
    assert np.array_equal(actual.original_selected, wanted.original_selected)
    old_actual = actual.original_node.fillna("<unselected>")
    old_wanted = wanted.original_node.fillna("<unselected>")
    assert np.array_equal(old_actual, old_wanted)
    frozen = records(archive, "frozen_states_json").set_index("metacell_id")
    for study in STUDIES:
        frame = records(archive, study + "_meta_json").set_index("metacell_id")
        rows = actual.xs(study, level="study_id", drop_level=False)
        ids = rows.index.get_level_values("metacell_id")
        assert np.array_equal(rows.frozen_state_group, frozen.loc[ids, "state_group"])
        native_column = "coarse_state" if "coarse_state" in frame else "working_comparison_group"
        assert np.array_equal(rows.native_original_label, frame.loc[ids, native_column].fillna("unknown"))
        assert np.array_equal(rows.n_cells, frame.loc[ids, "n_cells"])
    library = pd.read_csv(SOURCE / "library_classification.tsv", sep="\t", float_precision="round_trip")
    libkeys = ["study_id", "route", "library_key", "version"]
    wanted_libraries = expected.groupby(libkeys)
    assert not library.duplicated(libkeys).any()
    assert set(library.set_index(libkeys).index) == set(wanted_libraries.groups)
    lookup = library.set_index(libkeys)
    for key, group in wanted_libraries:
        row = lookup.loc[key]
        assert row.n_candidate_mc == len(group)
        assert row.n_Pstar == group.expected_label.eq("Pstar").sum()
        assert row.n_Mstar == group.expected_label.eq("Mstar").sum()
        assert row.n_unassigned == group.expected_label.eq("unassigned").sum()
        assert row.n_not_in_original_selected == (~group.original_selected).sum()
    return dict(
        status="PASS", classification_rows=len(actual),
        library_version_rows=len(library), numerical_max_abs_differences=errors,
    )


def main():
    assert sha(SOURCE / "PROTOCOL.zh.md") == LOCK, "Locked protocol changed"
    archive = np.load(SOURCE / "inputs.npz", allow_pickle=False)
    expected, candidate_summary = reconstruct_classes(archive)
    axes = reconstruct_axes(archive)
    classification_check = check_classifications(expected, archive)
    expected.to_csv(HERE / "expected_marker_classification.tsv.gz", sep="\t", index=False)
    candidate_summary.to_csv(HERE / "candidate_pool_check.tsv", sep="\t", index=False)
    report = dict(
        phase="input_and_marker_classification_reconstructed",
        complete_delivery_verified=False,
        expected_classification_rows=len(expected),
        candidate_pool=candidate_summary.to_dict("records"),
        full_primary_axes=axes,
        classification_check=classification_check,
        protocol_sha256=LOCK,
        inputs_sha256=sha(SOURCE / "inputs.npz"),
        verifier_sha256=sha(Path(__file__)),
        production_APA_checked=False,
    )
    (HERE / "delivery_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
