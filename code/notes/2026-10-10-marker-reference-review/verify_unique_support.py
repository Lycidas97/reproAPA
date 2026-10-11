"""Independently reconstruct all contributor sets after exact unique-MC losses."""
import json

import numpy as np
import pandas as pd

from verify_APA_delivery import HERE, SOURCE, KEYS, read, source_weights
from verify_marker_delivery import sha


BACKGROUNDS = {
    "full": {"D00", "D01", "D10", "D11"},
    "common_four": {"D00", "D01", "D10", "D11"},
    "common_M_bridge": {"D00", "D01"},
    "common_P_bridge": {"D00", "D10"},
}
PAIR_KEYS = KEYS + ["library_key"]
LIB_KEYS = ["setting", "version", "study_id", "route", "gene_id", "library_key"]


def background_pairs(frame, background):
    selected = frame[frame.contrast.isin(BACKGROUNDS[background])]
    if background != "full":
        selected = selected[selected.groupby(LIB_KEYS).contrast.transform("nunique").eq(len(BACKGROUNDS[background]))]
    return selected


def verify_intent():
    # Losing D01 must remove its library from shared backgrounds, while its
    # unrelated D00/D10 bridge keeps the original library and normalized weights.
    rows = [dict(setting="supported50", version="test", study_id="study", route="route", gene_id="gene",
        contrast=contrast, library_key=library, library_key_a=library, library_key_b=library,
        time_a=date, time_b=date, sample_unit_a=library, sample_unit_b=library)
        for library, date in [("A", "D1"), ("B", "D2")] for contrast in BACKGROUNDS["full"]]
    frame = pd.DataFrame(rows)
    retained = frame[~(frame.library_key.eq("A") & frame.contrast.eq("D01"))]
    assert set(background_pairs(retained, "common_four").library_key) == {"B"}
    assert set(background_pairs(retained, "common_M_bridge").library_key) == {"B"}
    assert set(background_pairs(retained, "common_P_bridge").library_key) == {"A", "B"}
    assert np.all(source_weights(background_pairs(retained, "common_four")) == 1)
    assert np.all(source_weights(background_pairs(retained, "common_P_bridge")) == .5)


def verify_unique_support():
    verify_intent()
    assert sha(SOURCE / "QA-ADDENDUM.zh.md") == "f0b3be44713f72e8845c57384c978e5ef6a0f87f3679c22054df411fb1b2836b"
    assert sha(HERE / "shared_MC_pairs_lost.tsv") == "f9fb32fd29a90af67ff19daf5ef57a29274196e50c8a55567997e3c6ffa0524d"
    columns = PAIR_KEYS + ["library_key_a", "library_key_b", "time_a", "time_b", "sample_unit_a", "sample_unit_b", "delta_DWUI"]
    frozen = read("pair_effects.tsv.gz", usecols=columns)
    assert not frozen.duplicated(PAIR_KEYS).any()
    lost = pd.read_csv(HERE / "shared_MC_pairs_lost.tsv", sep="\t").set_index(PAIR_KEYS)
    assert lost.index.is_unique and len(lost) == 15
    assert len(lost.reset_index()[["study_id", "route", "gene_id", "library_key"]].drop_duplicates()) == 7
    original_index = pd.MultiIndex.from_frame(frozen[PAIR_KEYS])
    removed = original_index.isin(lost.index)
    assert removed.sum() == 15 and set(original_index[removed]) == set(lost.index)
    retained = frozen.loc[~removed].copy()
    actual = read("unique_support_pair_weights.tsv.gz")
    assert set(actual.background) == set(BACKGROUNDS)
    assert not actual.duplicated(["background"] + PAIR_KEYS).any()
    expected_routes = {}
    expected_route_deltas = {}
    changes = []
    totals = []
    for background in BACKGROUNDS:
        old = background_pairs(frozen, background)
        expected = background_pairs(retained, background).copy()
        expected["independent_weight"] = source_weights(expected)
        observed = actual[actual.background.eq(background)].set_index(PAIR_KEYS)
        indexed = expected.set_index(PAIR_KEYS)
        assert len(observed) == len(indexed) and indexed.index.isin(observed.index).all(), background
        aligned = observed.reindex(indexed.index)
        assert np.allclose(aligned.weight, indexed.independent_weight, rtol=0, atol=1e-14), background
        for column in ["time_a", "time_b", "sample_unit_a", "sample_unit_b"]:
            if column in aligned:
                assert aligned[column].astype(str).eq(indexed[column].astype(str)).all(), (background, column)
        before_keys = set(pd.MultiIndex.from_frame(old[PAIR_KEYS]))
        after_keys = set(indexed.index)
        dropped = before_keys - after_keys
        assert not after_keys - before_keys
        allowed_libs = set(pd.MultiIndex.from_frame(lost.reset_index()[LIB_KEYS]))
        for pairkey in dropped:
            record = dict(zip(PAIR_KEYS, pairkey))
            assert tuple(record[k] for k in LIB_KEYS) in allowed_libs, (background, pairkey)
            assert background != "common_P_bridge"
            if background == "full":
                assert pairkey in set(lost.index)
            changes.append(dict(record, background=background))
        for key, part in expected.groupby(KEYS, sort=False):
            expected_routes[(background,) + key] = dict(zip(part.library_key, part.independent_weight))
            expected_route_deltas[(background,) + key] = float(part.delta_DWUI @ part.independent_weight)
        affected_genes = len({tuple(dict(zip(PAIR_KEYS, key))[k] for k in ["setting", "version", "study_id", "route", "gene_id"]) for key in dropped})
        totals.append(dict(background=background, retained_pair_weight_rows=len(expected), dropped_contributor_rows=len(dropped),
                           affected_setting_version_route_gene_sets=affected_genes))
    route = read("unique_support_route_effects.tsv.gz")
    bridge = read("unique_support_bridge_route_effects.tsv.gz")
    delivered = pd.concat([route, bridge], ignore_index=True)
    route_keys = ["background"] + KEYS
    assert not delivered.duplicated(route_keys).any() and len(delivered) == len(expected_routes)
    maximum_delta_error = 0.
    for row in delivered.to_dict("records"):
        key = tuple(row[k] for k in route_keys)
        expected = expected_routes[key]
        libs, weights = json.loads(row["contributor_libraries"]), json.loads(row["weights"])
        assert len(libs) == len(set(libs)) == row["n_libraries"] == len(expected)
        assert set(libs) == set(expected)
        assert np.allclose(weights, [expected[lib] for lib in libs], rtol=0, atol=1e-14), key
        assert abs(sum(weights) - 1) < 1e-12
        error = abs(row["delta_DWUI"] - expected_route_deltas[key])
        assert error < 1e-12, key
        maximum_delta_error = max(maximum_delta_error, error)
    frozen_routes = pd.concat([read("route_effects.tsv.gz", usecols=route_keys + ["delta_DWUI", "n_libraries"]),
                               read("bridge_route_effects.tsv.gz", usecols=route_keys + ["delta_DWUI", "n_libraries"])], ignore_index=True).set_index(route_keys)
    new_routes = delivered.set_index(route_keys)
    assert frozen_routes.index.is_unique and new_routes.index.is_unique
    allowed_changed_routes = {tuple(row[k] for k in route_keys) for row in changes}
    removed_routes = set(frozen_routes.index) - set(new_routes.index)
    assert removed_routes.issubset(allowed_changed_routes)
    assert not set(new_routes.index) - set(frozen_routes.index)
    differences = []
    common_keys = new_routes.index.intersection(frozen_routes.index)
    old_common, new_common = frozen_routes.loc[common_keys], new_routes.loc[common_keys]
    delta_diff = new_common.delta_DWUI - old_common.delta_DWUI
    libs_diff = new_common.n_libraries - old_common.n_libraries
    affected = (delta_diff.abs() > 1e-12) | libs_diff.ne(0)
    for key in common_keys[affected]:
        assert key in allowed_changed_routes and key[0] == "full", key
        differences.append(dict(zip(route_keys, key), status="retained_with_contributor_change",
            old_delta=float(old_common.loc[key].delta_DWUI), new_delta=float(new_common.loc[key].delta_DWUI),
            old_n_libraries=int(old_common.loc[key].n_libraries), new_n_libraries=int(new_common.loc[key].n_libraries)))
    for key in removed_routes:
        differences.append(dict(zip(route_keys, key), status="not_estimable_after_loss",
            old_delta=float(frozen_routes.loc[key].delta_DWUI), new_delta=np.nan,
            old_n_libraries=int(frozen_routes.loc[key].n_libraries), new_n_libraries=0))
    for background in ["common_four", "common_M_bridge", "common_P_bridge"]:
        old = frozen_routes.xs(background, level="background")
        new = new_routes.xs(background, level="background")
        assert set(old.index) == set(new.index)
        assert np.allclose(new.delta_DWUI, old.reindex(new.index).delta_DWUI, rtol=0, atol=1e-12)
        assert new.n_libraries.eq(old.reindex(new.index).n_libraries).all()
    pd.DataFrame(changes).to_csv(HERE / "unique_support_dropped_contributors.tsv", sep="\t", index=False)
    pd.DataFrame(totals).to_csv(HERE / "unique_support_background_check.tsv", sep="\t", index=False)
    pd.DataFrame(differences).to_csv(HERE / "unique_support_route_changes.tsv", sep="\t", index=False)
    report = dict(status="PASS", qualification_rule="distinct physical source MC union", loss_pair_rows=15,
        distinct_gene_library_pairs=7, original_pairs=len(frozen), retained_pairs=len(retained),
        background_checks=totals, delivered_route_contributor_sets=len(delivered),
        independently_reconstructed_route_delta_max_error=maximum_delta_error,
        changed_or_lost_route_effect_rows=len(differences), matched_backgrounds_complete_N_libraries_weights_deltas_unchanged=True,
        verification_scope="exact losses, all four backgrounds, contributors, hierarchy weights and route means; threshold summaries checked separately by root",
        verifier_sha256=sha(__file__))
    (HERE / "unique_support_validation.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    print(json.dumps(verify_unique_support(), indent=2))
