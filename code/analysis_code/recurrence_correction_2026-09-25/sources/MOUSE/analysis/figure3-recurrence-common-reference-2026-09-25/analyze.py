from __future__ import annotations

import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OUT = PROJECT / "outputs" / HERE.name
OUT.mkdir(parents=True, exist_ok=True)

STATE_ANALYSIS = PROJECT / "analysis/figure3-apa-gex-state-analysis-2026-09-19"
STATE_OUTPUT = PROJECT / "outputs/figure3-apa-gex-state-analysis-2026-09-19"
RECURRENCE_OUTPUT = PROJECT / "outputs/figure3-baseline-wpp-null-2026-09-23"
GEX_OUTPUT = PROJECT / "outputs/figure3-symmetric-apa-gex-2026-09-19"
OLD_AUDIT = PROJECT / "outputs/figure3-recurrence-strengthening-2026-09-25"

PATHS = {
    "recurrent": RECURRENCE_OUTPUT / "recurrent_gene_sets.tsv.gz",
    "shared_states": STATE_ANALYSIS / "shared/state_values.tsv.gz",
    "shared_regions": STATE_ANALYSIS / "shared/PAS_regions.tsv.gz",
    "conditional": STATE_ANALYSIS / "shared/conditional_sensitivity.tsv.gz",
    "triangle_studies": STATE_OUTPUT / "common_window_triangle_studies.tsv",
    "triangle": STATE_OUTPUT / "common_window_triangle.tsv.gz",
    "gex_pairs": GEX_OUTPUT / "shared_gene_pairs.tsv.gz",
    "old_D_four": OLD_AUDIT / "D_four_class_nulls.tsv.gz",
    "old_D_exchangeability": OLD_AUDIT / "D_exchangeability.tsv",
    "old_D_baseline": OLD_AUDIT / "D_shared_baseline_study_summary.tsv",
    "sampling_units": STATE_ANALYSIS / "raw/library_sampling_units.tsv",
}

AXES = {"mouse": "GSE103221__GSE115943", "human": "GSE118258__GSE178325"}
STUDIES = {"mouse": ("GSE103221", "GSE115943"), "human": ("GSE118258", "GSE178325")}
PROCESS_ROUTES = {
    "mouse": {"GSE103221": "CIP", "GSE115943": "WOT"},
    "human": {"GSE118258": "GSE118258", "GSE178325": "178325_0618"},
}
SOURCE_TIMES = {"D0", "source", "Stage_I_day_0"}
PRIMARY_THRESHOLD = 0.02
N_PERMUTATIONS = 10_000
SEED = 20260925


def read(name: str) -> pd.DataFrame:
    return pd.read_csv(PATHS[name], sep="\t", low_memory=False)


def save(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(OUT / name, sep="\t", index=False, float_format="%.10g")


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values: pd.Series) -> pd.Series:
    x = values.to_numpy(float)
    order = np.argsort(x)
    ranked = x[order]
    adjusted = np.minimum.accumulate((ranked * len(x) / np.arange(1, len(x) + 1))[::-1])[::-1]
    out = np.empty(len(x), float)
    out[order] = np.minimum(adjusted, 1.0)
    return pd.Series(out, index=values.index)


def markdown_table(frame: pd.DataFrame) -> str:
    values = frame.copy()
    for column in values.columns:
        values[column] = values[column].map(
            lambda value: "" if pd.isna(value) else (f"{value:.4g}" if isinstance(value, float) else str(value))
        )
    header = "| " + " | ".join(map(str, values.columns)) + " |"
    divider = "| " + " | ".join(["---"] * len(values.columns)) + " |"
    rows = ["| " + " | ".join(row) + " |" for row in values.astype(str).to_numpy()]
    return "\n".join([header, divider, *rows])


def frozen_sets() -> dict[tuple[str, str], pd.DataFrame]:
    data = read("recurrent")
    frozen = data[
        data.model.eq("baseline_q5_capacity")
        & data.omitted_study.eq("none")
        & data.threshold.eq(PRIMARY_THRESHOLD)
    ].copy()
    expected = {
        ("mouse", "P_S"): 488,
        ("human", "P_S"): 510,
        ("mouse", "P_M"): 216,
        ("human", "P_M"): 170,
    }
    result = {}
    for key, expected_n in expected.items():
        part = frozen[frozen.species.eq(key[0]) & frozen.contrast.eq(key[1])].copy()
        assert part.gene_id.nunique() == expected_n, (key, part.gene_id.nunique(), expected_n)
        assert not part.gene_id.duplicated().any()
        part["direction_sign"] = part.direction.map({"positive": 1, "negative": -1})
        assert part.direction_sign.notna().all()
        result[key] = part
    return result


def weighted_window(group: pd.DataFrame, value: str = "APA") -> tuple[float, dict[str, float]]:
    """Equal dates, equal sampling units within dates, equal libraries within units."""
    if group.empty:
        return np.nan, {}
    keys = ["time", "sample_unit", "library_key"]
    assert not group.duplicated(keys).any(), group[keys].loc[group.duplicated(keys, keep=False)]
    dates = sorted(group.time.astype(str).unique())
    weights: dict[str, float] = {}
    for date in dates:
        at_date = group[group.time.astype(str).eq(date)]
        units = sorted(at_date.sample_unit.astype(str).unique())
        for unit in units:
            at_unit = at_date[at_date.sample_unit.astype(str).eq(unit)]
            libraries = sorted(at_unit.library_key.astype(str).unique())
            weight = 1.0 / len(dates) / len(units) / len(libraries)
            for library in libraries:
                weights[library] = weight
    lookup = group.set_index("library_key")[value].astype(float)
    assert set(weights) == set(lookup.index.astype(str))
    estimate = float(sum(weights[str(library)] * lookup.loc[library] for library in lookup.index))
    assert np.isclose(sum(weights.values()), 1.0)
    return estimate, weights


def aggregate_fixed_reference(
    states: pd.DataFrame,
    frozen: dict[tuple[str, str], pd.DataFrame],
    variant: str,
    source_kind: str,
) -> pd.DataFrame:
    """Build route-state effects with one fixed true-source estimate per gene × route."""
    frozen_ids = pd.concat(
        [part[["species", "gene_id", "gene_symbol"]] for (species, contrast), part in frozen.items() if contrast == "P_S"],
        ignore_index=True,
    ).drop_duplicates(["species", "gene_id"])
    data = states.merge(
        frozen_ids[["species", "gene_id"]], on=["species", "gene_id"], how="inner", validate="many_to_one"
    )
    data = data[~data.route.astype(str).str.contains("reference")].copy()
    keys = ["species", "study_id", "route", "setting", "gene_id", "gene_symbol"]
    rows: list[dict] = []
    for identity, group in data.groupby(keys, sort=False):
        source = group[group.state.eq("S") & group.time.astype(str).isin(SOURCE_TIMES)]
        if source.empty:
            continue
        s_value, s_weights = weighted_window(source)
        source_signature = hashlib.sha256(
            json.dumps(sorted(s_weights.items()), separators=(",", ":")).encode()
        ).hexdigest()
        for state in "SREP":
            target = source if state == "S" else group[group.state.eq(state)]
            target_value, target_weights = weighted_window(target)
            if not np.isfinite(target_value):
                continue
            row = dict(zip(keys, identity))
            row.update(
                {
                    "variant": variant,
                    "source_kind": source_kind,
                    "state": state,
                    "source_WPP": s_value,
                    "target_WPP": target_value,
                    "delta_WPP": target_value - s_value,
                    "source_libraries": ";".join(sorted(s_weights)),
                    "source_weights_json": json.dumps(s_weights, sort_keys=True, separators=(",", ":")),
                    "source_signature": source_signature,
                    "target_libraries": ";".join(sorted(target_weights)),
                    "target_weights_json": json.dumps(target_weights, sort_keys=True, separators=(",", ":")),
                    "n_source_libraries": len(s_weights),
                    "n_target_libraries": len(target_weights),
                    "n_source_dates": source.time.nunique(),
                    "n_target_dates": target.time.nunique(),
                }
            )
            rows.append(row)
    result = pd.DataFrame(rows)
    assert result.groupby(keys).source_signature.nunique().max() == 1
    assert result[result.state.eq("S")].delta_WPP.abs().max() < 1e-12
    return result


def native_state_values() -> pd.DataFrame:
    parts = []
    for path in sorted((STATE_ANALYSIS / "raw").glob("gse*/apa_gex_states.tsv.gz")):
        data = pd.read_csv(path, sep="\t", low_memory=False)
        if data.study_id.iloc[0] not in {study for studies in STUDIES.values() for study in studies}:
            continue
        parts.append(data[data.setting.eq("supported50_p08")])
    return pd.concat(parts, ignore_index=True)


def reconstruct_mouse_membership_states(
    frozen: dict[tuple[str, str], pd.DataFrame], exclusions: dict[str, tuple[str, str]], variant: str
) -> pd.DataFrame:
    """Recreate mouse shared-axis library WPP after optional local-state exclusions."""
    axis = AXES["mouse"]
    regions = read("shared_regions")
    regions = regions[regions.axis_id.eq(axis)].copy()
    wanted = set(frozen[("mouse", "P_S")].gene_id)
    regions = regions[regions.gene_id.isin(wanted)]
    symbols = (
        read("shared_states")
        .loc[lambda x: x.axis_id.eq(axis) & x.gene_id.isin(wanted), ["gene_id", "gene_symbol"]]
        .drop_duplicates("gene_id")
        .set_index("gene_id")
        .gene_symbol.to_dict()
    )
    sampling = read("sampling_units").set_index("library_key")
    route_info = {
        "CIP": ("GSE103221", "gse103221_pluripotency_reprogramming_core"),
        "WOT": ("GSE115943", "gse115943_pluripotency_reprogramming_core"),
    }
    rows: list[dict] = []
    for route, (study_id, study_folder) in route_info.items():
        profile = PROJECT / f"data/apa-majority-shared-v1/profiles/{study_folder}"
        source = PROJECT / f"data/metacell-pas-v7-v1/studies/{study_folder}"
        meta = pd.read_csv(profile / "selected_metacells.tsv.gz", sep="\t")
        original_meta = pd.read_csv(source / "metacells.tsv.gz", sep="\t")
        assert np.array_equal(original_meta.metacell_id.iloc[meta.matrix_row], meta.metacell_id)
        local = original_meta[["metacell_id", "local_state_id"]]
        meta = meta.merge(local, on="metacell_id", how="left", validate="one_to_one")
        assert meta.local_state_id.notna().all()
        keep = np.ones(len(meta), dtype=bool)
        if route in exclusions:
            state, local_id = exclusions[route]
            keep &= ~(meta.node.eq(state).to_numpy() & meta.local_state_id.eq(local_id).to_numpy())
        features = pd.read_csv(
            source / "features.tsv.gz", sep="\t", usecols=["feature_index_1based", "atlas_pas_id"]
        )
        feature_index = features.set_index(features.atlas_pas_id.astype(str)).feature_index_1based.astype(int) - 1
        study_regions = regions[regions.study_id.eq(study_id)]
        used_atlas_ids = study_regions.atlas_pas_id.astype(str)
        used_atlas_ids = used_atlas_ids[used_atlas_ids.isin(feature_index.index)].drop_duplicates()
        raw_columns = np.sort(feature_index.loc[used_atlas_ids].to_numpy(int))
        raw_to_local = {int(column): index for index, column in enumerate(raw_columns)}
        matrix = sparse.load_npz(source / "metacell_by_pas_counts.npz").tocsr()[meta.matrix_row.to_numpy()][keep]
        matrix = matrix[:, raw_columns].tocsr()
        meta = meta.loc[keep].reset_index(drop=True)
        for gene_id, mapping in study_regions.groupby("gene_id", sort=False):
            mapping = mapping[mapping.atlas_pas_id.astype(str).isin(feature_index.index)]
            region_ids = sorted(mapping.region.unique())
            if len(region_ids) < 2:
                continue
            columns = []
            axes = []
            for region_id in region_ids:
                block = mapping[mapping.region.eq(region_id)]
                raw = feature_index.loc[block.atlas_pas_id.astype(str)].to_numpy(int)
                columns.append(np.asarray([raw_to_local[int(column)] for column in raw], dtype=int))
                axes.append(float(block.axis.iloc[0]))
            region_counts = np.column_stack(
                [np.asarray(matrix[:, column].sum(axis=1)).ravel() for column in columns]
            )
            total = region_counts.sum(axis=1)
            apa = np.divide(
                region_counts @ np.asarray(axes), total, out=np.full(len(meta), np.nan), where=total > 0
            )
            eligible = meta.purity.ge(0.8).to_numpy() & meta.node.isin(list("SREPM")).to_numpy() & (total >= 20)
            gene_symbol = symbols.get(gene_id, gene_id)
            for (library, state), indices0 in meta[eligible].groupby(["library_key", "node"]).groups.items():
                indices = np.asarray(list(indices0), dtype=int)
                if len(indices) < 3 or meta.iloc[indices].n_cells.sum() < 30:
                    continue
                route_name = str(meta.iloc[indices].route.iloc[0])
                if "reference" in route_name:
                    continue
                library = str(library)
                rows.append(
                    {
                        "axis_id": axis,
                        "species": "mouse",
                        "study_id": study_id,
                        "route": route_name,
                        "library_key": library,
                        "time": str(meta.iloc[indices].experimental_time.iloc[0]),
                        "sample_unit": str(sampling.loc[library, "sampling_unit"]),
                        "replicate_status": str(sampling.loc[library, "replicate_status"]),
                        "setting": "shared_p08",
                        "state": str(state),
                        "gene_id": gene_id,
                        "gene_symbol": gene_symbol,
                        "APA": float(np.nanmean(apa[indices])),
                        "n_mc": len(indices),
                        "n_cells": int(meta.iloc[indices].n_cells.sum()),
                        "pas_counts": int(total[indices].sum()),
                        "variant": variant,
                    }
                )
    return pd.DataFrame(rows)


def validate_membership_reconstruction(rebuilt: pd.DataFrame, frozen: dict[tuple[str, str], pd.DataFrame]) -> dict:
    original = read("shared_states")
    keys = ["axis_id", "species", "study_id", "route", "library_key", "time", "setting", "state", "gene_id"]
    original = original[
        original.axis_id.eq(AXES["mouse"])
        & original.setting.eq("shared_p08")
        & original.gene_id.isin(frozen[("mouse", "P_S")].gene_id)
    ][keys + ["APA"]]
    check = original.merge(rebuilt[keys + ["APA"]], on=keys, how="outer", suffixes=("_original", "_rebuilt"), indicator=True)
    max_error = float((check.APA_original - check.APA_rebuilt).abs().max())
    result = {
        "original_rows": len(original),
        "rebuilt_rows": len(rebuilt),
        "key_mismatches": int(check._merge.ne("both").sum()),
        "max_abs_WPP_error": max_error,
    }
    assert result["key_mismatches"] == 0, result
    assert max_error < 1e-10, result
    return result


def endpoint_table(
    effects: pd.DataFrame,
    genes: pd.DataFrame,
    species: str,
    endpoint_threshold: float,
    eligibility_source: str,
) -> pd.DataFrame:
    endpoint = effects[effects.species.eq(species) & effects.state.eq("P")].merge(
        genes[["gene_id", "direction_sign"]], on="gene_id", how="inner", validate="many_to_one"
    )
    endpoint["aligned_delta"] = endpoint.direction_sign * endpoint.delta_WPP
    endpoint["route_pass"] = endpoint.aligned_delta.ge(endpoint_threshold)
    rows = []
    for gene in genes.itertuples(index=False):
        group = endpoint[endpoint.gene_id.eq(gene.gene_id)]
        per_study = {}
        for study in STUDIES[species]:
            block = group[group.study_id.eq(study)]
            per_study[study] = {
                "measurable": len(block) > 0,
                "pass": len(block) > 0 and bool(block.route_pass.all()),
                "routes": ";".join(sorted(block.route.astype(str).unique())),
                "n_routes": int(block.route.nunique()),
                "minimum_aligned_delta": float(block.aligned_delta.min()) if len(block) else np.nan,
            }
        measurable = all(per_study[study]["measurable"] for study in STUDIES[species])
        retained = measurable and all(per_study[study]["pass"] for study in STUDIES[species])
        row = {
            "species": species,
            "gene_id": gene.gene_id,
            "gene_symbol": gene.gene_symbol,
            "direction": gene.direction,
            "direction_sign": gene.direction_sign,
            "endpoint_threshold": endpoint_threshold,
            "eligibility_source": eligibility_source,
            "endpoint_measurable": measurable,
            "endpoint_retained": retained,
        }
        for index, study in enumerate(STUDIES[species], 1):
            for key, value in per_study[study].items():
                row[f"study_{index}_{key}"] = value
        rows.append(row)
    return pd.DataFrame(rows)


def signed_class(value: float, threshold: float) -> str:
    if not np.isfinite(value):
        return "unavailable"
    if value >= threshold:
        return "endpoint_aligned"
    if value <= -threshold:
        return "reverse"
    return "within_threshold"


def process_tables(
    effects: pd.DataFrame,
    genes: pd.DataFrame,
    species: str,
    endpoints: pd.DataFrame,
    intermediate_threshold: float,
    analysis_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    endpoint_map = endpoints.set_index("gene_id").endpoint_retained.to_dict()
    states = ("R", "E") if species == "mouse" else ("R",)
    rows = []
    for gene in genes.itertuples(index=False):
        for state in states:
            record = {
                "analysis": analysis_label,
                "species": species,
                "gene_id": gene.gene_id,
                "gene_symbol": gene.gene_symbol,
                "direction": gene.direction,
                "direction_sign": gene.direction_sign,
                "state": state,
                "intermediate_threshold": intermediate_threshold,
                "endpoint_retained": bool(endpoint_map.get(gene.gene_id, False)),
            }
            route_classes = []
            route_eligible = []
            for index, study in enumerate(STUDIES[species], 1):
                route = PROCESS_ROUTES[species][study]
                block = effects[
                    effects.species.eq(species)
                    & effects.study_id.eq(study)
                    & effects.route.eq(route)
                    & effects.gene_id.eq(gene.gene_id)
                ]
                by_state = block.set_index("state").delta_WPP.to_dict()
                same_series = all(name in by_state for name in ("S", state, "P"))
                aligned_value = gene.direction_sign * by_state.get(state, np.nan)
                category = signed_class(aligned_value, intermediate_threshold) if same_series else "unavailable"
                record.update(
                    {
                        f"study_{index}": study,
                        f"route_{index}": route,
                        f"same_series_eligible_{index}": same_series,
                        f"aligned_delta_{index}": aligned_value,
                        f"category_{index}": category,
                        f"P_S_measurable_{index}": "P" in by_state and "S" in by_state,
                    }
                )
                route_classes.append(category)
                route_eligible.append(same_series)
            pair_measurable = all(route_eligible)
            if not record["endpoint_retained"]:
                pair_class = "endpoint_not_retained"
            elif not pair_measurable:
                pair_class = "unavailable"
            elif route_classes[0] == route_classes[1]:
                pair_class = route_classes[0]
            elif set(route_classes) == {"endpoint_aligned", "reverse"}:
                pair_class = "direction_difference"
            else:
                pair_class = "mixed_or_threshold_undetermined"
            record["pair_measurable"] = pair_measurable
            record["pair_class"] = pair_class
            rows.append(record)
    pairs = pd.DataFrame(rows)
    profiles = []
    for gene in genes.itertuples(index=False):
        block = pairs[pairs.gene_id.eq(gene.gene_id)].set_index("state")
        classes = block.pair_class.to_dict()
        measurable = block.pair_measurable.to_dict()
        if not bool(endpoint_map.get(gene.gene_id, False)):
            profile = "endpoint_not_retained"
        elif not all(measurable.get(state, False) for state in states):
            profile = "incomplete_same_series_process"
        elif species == "human":
            profile = classes["R"]
        elif classes["R"] == classes["E"] == "endpoint_aligned":
            profile = "both_intermediates_endpoint_aligned"
        elif "direction_difference" in classes.values():
            profile = "at_least_one_direction_difference"
        elif set(classes.values()).issubset({"endpoint_aligned", "reverse"}) and "reverse" in classes.values():
            profile = "consistent_signed_configuration_with_reverse"
        else:
            profile = "mixed_or_threshold_undetermined"
        profiles.append(
            {
                "analysis": analysis_label,
                "species": species,
                "gene_id": gene.gene_id,
                "gene_symbol": gene.gene_symbol,
                "profile_class": profile,
            }
        )
    return pairs, pd.DataFrame(profiles)


def module_A(frozen: dict[tuple[str, str], pd.DataFrame]) -> dict[str, pd.DataFrame | dict]:
    shared = read("shared_states")
    variant_effects: dict[str, pd.DataFrame] = {}
    for setting, label in [
        ("shared_p08", "primary_0.02_0.02"),
        ("shared_p09", "purity_0.9"),
        ("shared_depth40", "depth_40"),
    ]:
        state_rows = shared[
            shared.setting.eq(setting) & shared.axis_id.isin(AXES.values())
        ].copy()
        variant_effects[label] = aggregate_fixed_reference(state_rows, frozen, label, "shared_axis")

    native = native_state_values()
    variant_effects["native_PAS_coordinates"] = aggregate_fixed_reference(
        native, frozen, "native_PAS_coordinates", "native_axis"
    )

    rebuilt_original = reconstruct_mouse_membership_states(frozen, {}, "membership_original_rebuilt")
    reconstruction = validate_membership_reconstruction(rebuilt_original, frozen)
    membership_variants = {
        "exclude_CIP_C012": {"CIP": ("E", "C012")},
        "exclude_WOT_L00": {"WOT": ("R", "WOT::common_to_2i::L00")},
        "exclude_CIP_C012_and_WOT_L00": {
            "CIP": ("E", "C012"),
            "WOT": ("R", "WOT::common_to_2i::L00"),
        },
    }
    membership_states = [rebuilt_original]
    for label, exclusions in membership_variants.items():
        rebuilt = reconstruct_mouse_membership_states(frozen, exclusions, label)
        membership_states.append(rebuilt)
        variant_effects[label] = aggregate_fixed_reference(rebuilt, frozen, label, "shared_axis_member_exclusion")
    membership_state_table = pd.concat(membership_states, ignore_index=True)

    primary_effects = variant_effects["primary_0.02_0.02"]
    endpoints_parts = []
    pairs_parts = []
    profiles_parts = []

    # Corrected primary: endpoint 0.02, intermediate 0.02.
    primary_endpoints: dict[str, pd.DataFrame] = {}
    for species in ("mouse", "human"):
        genes = frozen[(species, "P_S")]
        endpoints = endpoint_table(primary_effects, genes, species, 0.02, "recomputed_0.02")
        endpoints["analysis"] = "primary_0.02_0.02"
        primary_endpoints[species] = endpoints
        pairs, profiles = process_tables(
            primary_effects, genes, species, endpoints, 0.02, "primary_0.02_0.02"
        )
        endpoints_parts.append(endpoints)
        pairs_parts.append(pairs)
        profiles_parts.append(profiles)

        # Clean classification sensitivity: freeze the primary endpoint eligibility.
        fixed = endpoints.copy()
        fixed["analysis"] = "fixed_primary_eligibility_intermediate_0.05"
        fixed["eligibility_source"] = "frozen_from_primary_0.02"
        pairs, profiles = process_tables(
            primary_effects,
            genes,
            species,
            fixed,
            0.05,
            "fixed_primary_eligibility_intermediate_0.05",
        )
        endpoints_parts.append(fixed)
        pairs_parts.append(pairs)
        profiles_parts.append(profiles)

        # Strong endpoint subset: independently require endpoint 0.05.
        strong = endpoint_table(primary_effects, genes, species, 0.05, "recomputed_0.05")
        strong["analysis"] = "strong_endpoint_0.05_intermediate_0.05"
        pairs, profiles = process_tables(
            primary_effects,
            genes,
            species,
            strong,
            0.05,
            "strong_endpoint_0.05_intermediate_0.05",
        )
        endpoints_parts.append(strong)
        pairs_parts.append(pairs)
        profiles_parts.append(profiles)

    # Measurement / definition sensitivity at the corrected 0.02 / 0.02 rule.
    sensitivity_labels = [
        "purity_0.9",
        "depth_40",
        "native_PAS_coordinates",
        "exclude_CIP_C012",
        "exclude_WOT_L00",
        "exclude_CIP_C012_and_WOT_L00",
    ]
    for label in sensitivity_labels:
        effects = variant_effects[label]
        species_values = ("mouse",) if label.startswith("exclude_") else ("mouse", "human")
        for species in species_values:
            genes = frozen[(species, "P_S")]
            endpoints = endpoint_table(effects, genes, species, 0.02, "recomputed_0.02")
            endpoints["analysis"] = label
            pairs, profiles = process_tables(effects, genes, species, endpoints, 0.02, label)
            endpoints_parts.append(endpoints)
            pairs_parts.append(pairs)
            profiles_parts.append(profiles)

    endpoints = pd.concat(endpoints_parts, ignore_index=True)
    pairs = pd.concat(pairs_parts, ignore_index=True)
    profiles = pd.concat(profiles_parts, ignore_index=True)
    effects = pd.concat(variant_effects.values(), ignore_index=True)

    primary_pairs = pairs[pairs.analysis.eq("primary_0.02_0.02")].copy()
    transition_rows = []
    for label in [
        "fixed_primary_eligibility_intermediate_0.05",
        "strong_endpoint_0.05_intermediate_0.05",
        *sensitivity_labels,
    ]:
        sensitivity = pairs[pairs.analysis.eq(label)].copy()
        if sensitivity.empty:
            continue
        merged = primary_pairs.merge(
            sensitivity,
            on=["species", "gene_id", "state"],
            how="inner",
            suffixes=("_primary", "_sensitivity"),
            validate="one_to_one",
        )
        for row in merged.itertuples(index=False):
            if not row.endpoint_retained_primary:
                outcome = "primary_endpoint_not_eligible"
            elif not row.pair_measurable_sensitivity:
                outcome = "unavailable"
            elif not row.endpoint_retained_sensitivity:
                outcome = "endpoint_no_longer_passes"
            elif not row.pair_measurable_primary:
                outcome = "primary_process_unavailable"
            elif row.pair_class_primary != row.pair_class_sensitivity:
                outcome = "category_changed"
            else:
                outcome = "stable"
            transition_rows.append(
                {
                    "variant": label,
                    "species": row.species,
                    "gene_id": row.gene_id,
                    "state": row.state,
                    "primary_endpoint_retained": row.endpoint_retained_primary,
                    "sensitivity_endpoint_retained": row.endpoint_retained_sensitivity,
                    "primary_pair_measurable": row.pair_measurable_primary,
                    "sensitivity_pair_measurable": row.pair_measurable_sensitivity,
                    "primary_class": row.pair_class_primary,
                    "sensitivity_class": row.pair_class_sensitivity,
                    "outcome": outcome,
                }
            )
    transitions = pd.DataFrame(transition_rows)
    transition_summary = (
        transitions.groupby(["variant", "species", "state", "outcome"], dropna=False)
        .size()
        .rename("n_gene_state_rows")
        .reset_index()
    )

    flow_rows = []
    for analysis, block in endpoints.groupby("analysis", sort=False):
        for species, group in block.groupby("species"):
            total = len(group)
            flow_rows.extend(
                [
                    {"analysis": analysis, "species": species, "step": "frozen", "n_genes": total},
                    {
                        "analysis": analysis,
                        "species": species,
                        "step": "endpoint_measurable",
                        "n_genes": int(group.endpoint_measurable.sum()),
                    },
                    {
                        "analysis": analysis,
                        "species": species,
                        "step": "endpoint_retained",
                        "n_genes": int(group.endpoint_retained.sum()),
                    },
                ]
            )
        process = pairs[pairs.analysis.eq(analysis)]
        for (species, state), group in process.groupby(["species", "state"]):
            flow_rows.append(
                {
                    "analysis": analysis,
                    "species": species,
                    "state": state,
                    "step": "same_series_process_measurable",
                    "n_genes": int(group.pair_measurable.sum()),
                }
            )
    flow = pd.DataFrame(flow_rows)

    guan_background = pairs[
        pairs.species.eq("human")
        & pairs.state.eq("R")
        & pairs.route_2.eq("178325_0618")
        & pairs["aligned_delta_2"].notna()
        & ~pairs.same_series_eligible_2
    ].copy()

    return {
        "state_effects": effects,
        "membership_states": membership_state_table,
        "endpoints": endpoints,
        "pairs": pairs,
        "profiles": profiles,
        "transitions": transitions,
        "transition_summary": transition_summary,
        "flow": flow,
        "guan_background": guan_background,
        "membership_reconstruction": reconstruction,
    }


def b_configuration(aligned_ps: float, aligned_ms: float, pm: float, threshold: float) -> str:
    if not all(np.isfinite([aligned_ps, aligned_ms, pm])):
        return "incomparable"
    if aligned_ps < threshold:
        return "common_window_endpoint_not_retained"
    if aligned_ms >= threshold and abs(pm) < threshold:
        return "P_M_common_shift_PM_within_threshold"
    if aligned_ms >= threshold and abs(pm) >= threshold:
        return "P_M_common_shift_with_PM_difference"
    if aligned_ms <= -threshold:
        return "P_M_opposite_relative_to_S"
    if abs(aligned_ms) < threshold:
        return "P_changed_M_within_threshold"
    return "other"


def threshold_direction(value: float, threshold: float = 0.02) -> str:
    if not np.isfinite(value):
        return "unavailable"
    if value >= threshold:
        return "positive"
    if value <= -threshold:
        return "negative"
    return "within_threshold"


def module_B(frozen: dict[tuple[str, str], pd.DataFrame]) -> dict[str, pd.DataFrame]:
    ps = frozen[("mouse", "P_S")][["gene_id", "gene_symbol", "direction", "direction_sign"]]
    pm = frozen[("mouse", "P_M")][["gene_id", "direction_sign"]].rename(
        columns={"direction_sign": "pm_direction_sign"}
    )
    ps_ids, pm_ids = set(ps.gene_id), set(pm.gene_id)
    union = ps_ids | pm_ids
    axis = AXES["mouse"]

    triangle = read("triangle")
    triangle = triangle[
        triangle.axis_id.eq(axis) & triangle.setting.eq("shared_p08") & triangle.gene_id.isin(union)
    ].copy()
    assert triangle.APA_delta_triangle_error.abs().max() <= 1e-8
    assert triangle.GEX_delta_triangle_error.abs().max() <= 1e-8
    assert triangle.GEX_all_delta_triangle_error.abs().max() <= 1e-8

    effects = read("triangle_studies")
    effects = effects[
        effects.axis_id.eq(axis) & effects.setting.eq("shared_p08") & effects.gene_id.isin(union)
    ].copy()
    effects = effects.merge(ps, on=["gene_id", "gene_symbol"], how="left", validate="many_to_one")
    effects = effects.merge(pm, on="gene_id", how="left", validate="many_to_one")
    effects["cohort"] = np.where(
        effects.gene_id.isin(ps_ids), "P_S_recurrent", np.where(effects.gene_id.isin(pm_ids), "P_M_only", "neither")
    )
    effects["aligned_PS"] = effects.direction_sign * effects.APA_delta_PS
    effects["aligned_MS"] = effects.direction_sign * effects.APA_delta_MS
    effects["configuration"] = [
        b_configuration(a, b, c, 0.02) if cohort == "P_S_recurrent" else "not_applicable_P_M_only"
        for a, b, c, cohort in zip(effects.aligned_PS, effects.aligned_MS, effects.APA_delta, effects.cohort)
    ]
    effects["PM_direction"] = [threshold_direction(value) for value in effects.APA_delta]

    comparison_rows = []
    for (cohort, gene_id), group in effects.groupby(["cohort", "gene_id"]):
        by = group.set_index("study_id")
        if not all(study in by.index for study in STUDIES["mouse"]):
            continue
        first = by.loc[STUDIES["mouse"][0]]
        second = by.loc[STUDIES["mouse"][1]]
        if isinstance(first, pd.DataFrame) or isinstance(second, pd.DataFrame):
            raise AssertionError((gene_id, "duplicate common-window study summary"))
        pm_directional = first.PM_direction in {"positive", "negative"} and second.PM_direction in {
            "positive",
            "negative",
        }
        same_pm_direction = pm_directional and first.PM_direction == second.PM_direction
        endpoint_retained_both = (
            cohort == "P_S_recurrent"
            and first.configuration != "common_window_endpoint_not_retained"
            and second.configuration != "common_window_endpoint_not_retained"
        )
        coarse_comparable = endpoint_retained_both
        same_coarse = coarse_comparable and first.configuration == second.configuration
        comparison_rows.append(
            {
                "cohort": cohort,
                "gene_id": gene_id,
                "gene_symbol": first.gene_symbol,
                "configuration_study_1": first.configuration,
                "configuration_study_2": second.configuration,
                "common_window_endpoint_retained_both": endpoint_retained_both,
                "coarse_configuration_comparable": coarse_comparable,
                "same_coarse_configuration": same_coarse if coarse_comparable else pd.NA,
                "PM_direction_study_1": first.PM_direction,
                "PM_direction_study_2": second.PM_direction,
                "same_PM_category": first.PM_direction == second.PM_direction,
                "PM_direction_comparable": pm_directional,
                "same_PM_direction": same_pm_direction,
                "same_coarse_and_PM_direction": (
                    same_coarse and same_pm_direction
                ) if coarse_comparable else pd.NA,
            }
        )
    configurations = pd.DataFrame(comparison_rows)

    conditional = read("conditional")
    conditional = conditional[
        conditional.axis_id.eq(axis)
        & conditional.contrast.eq("P_M")
        & conditional.study_id.isin(STUDIES["mouse"])
        & conditional.gene_id.isin(union)
    ].copy()
    expected_rank = {"cycle": 4, "cycle_and_gex": 5, "cycle_gex_depth": 6}
    conditional["expected_rank"] = conditional.model.map(expected_rank)
    conditional["diagnostic_pass"] = (
        conditional.adjusted.notna()
        & conditional["rank"].eq(conditional.expected_rank)
        & conditional["condition"].lt(1e6)
    )
    key = ["axis_id", "species", "study_id", "route", "gene_id", "gene_symbol", "library_key", "time", "sample_unit"]
    model_count = conditional[conditional.diagnostic_pass].groupby(key).model.nunique()
    common_keys = set(model_count[model_count.eq(len(expected_rank))].index)
    conditional["common_across_models"] = [tuple(row) in common_keys for row in conditional[key].itertuples(index=False, name=None)]

    summary_rows = []
    for background, use in [
        ("model_specific", conditional.diagnostic_pass),
        ("common_across_models", conditional.diagnostic_pass & conditional.common_across_models),
    ]:
        data = conditional[use]
        for identity, group in data.groupby(
            ["axis_id", "species", "study_id", "route", "gene_id", "gene_symbol", "model"], sort=False
        ):
            unadjusted, weights_before = weighted_window(group, "unadjusted")
            adjusted, weights_after = weighted_window(group, "adjusted")
            assert weights_before == weights_after
            row = dict(
                zip(["axis_id", "species", "study_id", "route", "gene_id", "gene_symbol", "model"], identity)
            )
            row.update(
                {
                    "background": background,
                    "unadjusted": unadjusted,
                    "adjusted": adjusted,
                    "same_weight_signature": hashlib.sha256(
                        json.dumps(sorted(weights_before.items()), separators=(",", ":")).encode()
                    ).hexdigest(),
                    "n_libraries": len(weights_before),
                    "direction_before": threshold_direction(unadjusted),
                    "direction_after": threshold_direction(adjusted),
                    "direction_preserved": threshold_direction(unadjusted) == threshold_direction(adjusted),
                    "max_condition": float(group["condition"].max()),
                    "minimum_rank": int(group["rank"].min()),
                }
            )
            summary_rows.append(row)
    conditional_summary = pd.DataFrame(summary_rows)
    conditional_summary["in_P_S_frozen"] = conditional_summary.gene_id.isin(ps_ids)
    conditional_summary["in_P_M_frozen"] = conditional_summary.gene_id.isin(pm_ids)

    conditional_cross = []
    for (background, model, gene_id), group in conditional_summary.groupby(["background", "model", "gene_id"]):
        by = group.set_index("study_id")
        if not all(study in by.index for study in STUDIES["mouse"]):
            continue
        a = by.loc[STUDIES["mouse"][0]]
        b = by.loc[STUDIES["mouse"][1]]
        if isinstance(a, pd.DataFrame) or isinstance(b, pd.DataFrame):
            raise AssertionError((background, model, gene_id, "duplicate conditional summary"))
        unadjusted_category_1 = threshold_direction(a.unadjusted)
        unadjusted_category_2 = threshold_direction(b.unadjusted)
        adjusted_category_1 = threshold_direction(a.adjusted)
        adjusted_category_2 = threshold_direction(b.adjusted)
        unadjusted_direction_comparable = (
            unadjusted_category_1 in {"positive", "negative"}
            and unadjusted_category_2 in {"positive", "negative"}
        )
        adjusted_direction_comparable = (
            adjusted_category_1 in {"positive", "negative"}
            and adjusted_category_2 in {"positive", "negative"}
        )
        conditional_cross.append(
            {
                "background": background,
                "model": model,
                "gene_id": gene_id,
                "gene_symbol": a.gene_symbol,
                "unadjusted_study_1": a.unadjusted,
                "unadjusted_study_2": b.unadjusted,
                "adjusted_study_1": a.adjusted,
                "adjusted_study_2": b.adjusted,
                "unadjusted_category_study_1": unadjusted_category_1,
                "unadjusted_category_study_2": unadjusted_category_2,
                "adjusted_category_study_1": adjusted_category_1,
                "adjusted_category_study_2": adjusted_category_2,
                "unadjusted_cross_study_same_category": unadjusted_category_1 == unadjusted_category_2,
                "adjusted_cross_study_same_category": adjusted_category_1 == adjusted_category_2,
                "unadjusted_direction_comparable": unadjusted_direction_comparable,
                "adjusted_direction_comparable": adjusted_direction_comparable,
                "unadjusted_cross_study_same_direction": (
                    unadjusted_direction_comparable and unadjusted_category_1 == unadjusted_category_2
                ),
                "adjusted_cross_study_same_direction": (
                    adjusted_direction_comparable and adjusted_category_1 == adjusted_category_2
                ),
                "in_P_S_frozen": gene_id in ps_ids,
                "in_P_M_frozen": gene_id in pm_ids,
            }
        )
    conditional_cross = pd.DataFrame(conditional_cross)

    diagnostics = (
        conditional.groupby(["model"], dropna=False)
        .agg(
            rows=("gene_id", "size"),
            diagnostic_pass=("diagnostic_pass", "sum"),
            common_background_rows=("common_across_models", "sum"),
            genes=("gene_id", "nunique"),
        )
        .reset_index()
    )
    return {
        "triangle_design": triangle,
        "triangle_effects": effects,
        "configurations": configurations,
        "conditional_rows": conditional,
        "conditional_summary": conditional_summary,
        "conditional_cross": conditional_cross,
        "conditional_diagnostics": diagnostics,
    }


def effect_sign(low: np.ndarray, high: np.ndarray, threshold: float) -> np.ndarray:
    return np.where(low >= threshold, 1, np.where(high <= -threshold, -1, 0))


def four_counts(x: np.ndarray, y: np.ndarray, g: np.ndarray) -> np.ndarray:
    product = x * y
    return np.array(
        [
            np.count_nonzero((product > 0) & (g == 1)),
            np.count_nonzero((product > 0) & (g == -1)),
            np.count_nonzero((product < 0) & (g == 1)),
            np.count_nonzero((product < 0) & (g == -1)),
        ],
        dtype=int,
    )


def qbin(values: pd.Series, bins: int) -> pd.Series:
    return pd.qcut(values.rank(method="average"), bins, labels=False, duplicates="drop").astype("Int64")


def guan_pattern(routes: str) -> str:
    values = set(str(routes).split(";"))
    has_0618 = "178325_0618" in values
    has_0809 = "178325_0809" in values
    if has_0618 and has_0809:
        return "both_measurable"
    if has_0618:
        return "0618_only"
    if has_0809:
        return "0809_only"
    return "neither_prespecified_series"


def prepare_D_pairs() -> pd.DataFrame:
    pairs = read("gex_pairs")
    pairs = pairs[
        pairs.setting.eq("shared_p08")
        & pairs.contrast.eq("P_S")
        & pairs.axis_id.isin(AXES.values())
    ].copy()
    baseline = read("old_D_baseline")
    for side in (1, 2):
        renamed = baseline.rename(
            columns={
                "study_id": f"study_id_{side}",
                "reference_WPP": f"reference_WPP_{side}",
                "reference_WPP_min": f"reference_WPP_min_{side}",
                "reference_WPP_max": f"reference_WPP_max_{side}",
                "baseline_routes": f"baseline_routes_{side}",
                "route_names": f"baseline_route_names_{side}",
            }
        )
        keep = [
            "axis_id",
            "species",
            f"study_id_{side}",
            "setting",
            "gene_id",
            f"reference_WPP_{side}",
            f"reference_WPP_min_{side}",
            f"reference_WPP_max_{side}",
            f"baseline_routes_{side}",
            f"baseline_route_names_{side}",
        ]
        pairs = pairs.merge(
            renamed[keep],
            on=["axis_id", "species", f"study_id_{side}", "setting", "gene_id"],
            how="left",
            validate="many_to_one",
        )
    assert pairs.filter(regex="^reference_WPP").notna().all().all()
    human = pairs.species.eq("human")
    pairs["guan_series_pattern"] = "not_applicable_mouse"
    pairs.loc[human, "guan_series_pattern"] = pairs.loc[human, "routes_2"].map(guan_pattern)
    assert set(pairs.loc[human, "guan_series_pattern"]) <= {
        "0618_only",
        "0809_only",
        "both_measurable",
    }
    return pairs


def permute_human_D(
    data: pd.DataFrame, threshold: float, null_model: str, baseline_bins: int
) -> tuple[pd.DataFrame, pd.DataFrame]:
    w = data[data.species.eq("human")].copy().reset_index(drop=True)
    x = effect_sign(w.alo_1.to_numpy(), w.ahi_1.to_numpy(), threshold)
    y = effect_sign(w.alo_2.to_numpy(), w.ahi_2.to_numpy(), threshold)
    gx = effect_sign(w.glo_1.to_numpy(), w.ghi_1.to_numpy(), 0.25)
    gy = effect_sign(w.glo_2.to_numpy(), w.ghi_2.to_numpy(), 0.25)
    gx = np.where(w.meas_1.to_numpy(bool), gx, 0)
    gy = np.where(w.meas_2.to_numpy(bool), gy, 0)
    g = gx * gy
    observed = four_counts(x, y, g)
    rows, diagnostics = [], []
    labels_by_study = {1: x, 2: y}
    other_by_study = {1: y, 2: x}
    for shuffled in (1, 2):
        label = labels_by_study[shuffled]
        suffix = str(shuffled)
        w[f"depth_bin_{suffix}"] = qbin(w[f"depth_{suffix}"], 3)
        w[f"baseline_bin_{suffix}"] = qbin(w[f"reference_WPP_{suffix}"], baseline_bins)
        w[f"gex_response_{suffix}"] = pd.cut(
            w[f"GEX_{suffix}"], [-np.inf, -0.25, 0.25, np.inf], labels=["down", "flat", "up"]
        ).astype(str)
        w[f"pas_cap_{suffix}"] = np.minimum(w[f"PAS_{suffix}"].astype(int), 4)
        can_negative = w[f"reference_WPP_min_{suffix}"].ge(threshold - 1e-9)
        can_positive = w[f"reference_WPP_max_{suffix}"].le(1 - threshold + 1e-9)
        assert not ((label == 1) & ~can_positive.to_numpy()).any()
        assert not ((label == -1) & ~can_negative.to_numpy()).any()
        pieces = [
            w[f"gex_response_{suffix}"],
            w[f"depth_bin_{suffix}"].astype(str),
            w[f"pas_cap_{suffix}"].astype(str),
            w[f"baseline_bin_{suffix}"].astype(str),
            can_negative.astype(str),
            can_positive.astype(str),
            w["guan_series_pattern"].astype(str),
        ]
        if null_model == "direction_specific":
            pieces.append(pd.Series(label != 0, index=w.index).astype(str))
        stratum = pd.concat(pieces, axis=1).astype(str).agg("|".join, axis=1)
        unique = sorted(stratum.unique())
        groups = [np.flatnonzero(stratum.to_numpy() == value) for value in unique]
        for value, indices in zip(unique, groups):
            labels = label[indices]
            counts = {int(v): int(np.count_nonzero(labels == v)) for v in (-1, 0, 1)}
            probabilities = np.asarray(list(counts.values()), float) / len(indices)
            diagnostics.append(
                {
                    "species": "human",
                    "threshold": threshold,
                    "null_model": null_model,
                    "baseline_bins": baseline_bins,
                    "shuffled_study": shuffled,
                    "study_id": w[f"study_id_{suffix}"].iloc[0],
                    "stratum": value,
                    "guan_series_in_stratum": value.split("|")[-2 if null_model == "direction_specific" else -1],
                    "n_genes": len(indices),
                    "negative": counts[-1],
                    "non_event": counts[0],
                    "positive": counts[1],
                    "mixed_labels": len(set(labels)) > 1,
                    "expected_labels_changed": len(indices) * (1 - np.sum(probabilities**2)),
                }
            )
        seed_key = f"human|{threshold}|{null_model}|{baseline_bins}|{shuffled}|guan_series"
        seed = SEED + int(hashlib.sha256(seed_key.encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng(seed)
        null = np.empty((N_PERMUTATIONS, 4), dtype=int)
        for iteration in range(N_PERMUTATIONS):
            permuted = label.copy()
            for indices in groups:
                permuted[indices] = rng.permutation(permuted[indices])
            if shuffled == 1:
                null[iteration] = four_counts(permuted, other_by_study[shuffled], g)
            else:
                null[iteration] = four_counts(other_by_study[shuffled], permuted, g)
        if null_model == "direction_specific":
            assert np.all(null[:, 0] + null[:, 2] == observed[0] + observed[2])
            assert np.all(null[:, 1] + null[:, 3] == observed[1] + observed[3])
        shuffleable = int(sum(len(indices) for indices in groups if len(set(label[indices])) > 1))
        for index, cell in enumerate("ABCD"):
            mean = float(null[:, index].mean())
            rows.append(
                {
                    "species": "human",
                    "axis_id": AXES["human"],
                    "threshold": threshold,
                    "null_model": null_model,
                    "baseline_bins": baseline_bins,
                    "shuffled_study": shuffled,
                    "study_id": w[f"study_id_{suffix}"].iloc[0],
                    "cell": cell,
                    "observed": int(observed[index]),
                    "null_mean": mean,
                    "enrichment": observed[index] / mean if mean > 0 else np.nan,
                    "null_p025": float(np.quantile(null[:, index], 0.025)),
                    "null_p975": float(np.quantile(null[:, index], 0.975)),
                    "p_enrich": (1 + int(np.count_nonzero(null[:, index] >= observed[index]))) / (N_PERMUTATIONS + 1),
                    "p_deplete": (1 + int(np.count_nonzero(null[:, index] <= observed[index]))) / (N_PERMUTATIONS + 1),
                    "comparable_genes": len(w),
                    "shuffleable_genes": shuffleable,
                    "permutations": N_PERMUTATIONS,
                    "seed": seed,
                    "source": "rerun_with_Guan_series_stratum",
                }
            )
    return pd.DataFrame(rows), pd.DataFrame(diagnostics)


def summarize_D_tests(four: pd.DataFrame) -> pd.DataFrame:
    tests = (
        four[four.cell.eq("B")]
        .sort_values("shuffled_study")
        .groupby(["species", "axis_id", "threshold", "null_model", "baseline_bins"], dropna=False)
        .agg(
            observed=("observed", "first"),
            p_both=("p_enrich", "max"),
            enrichment_min=("enrichment", "min"),
            p_shuffle_study_1=("p_enrich", "first"),
            p_shuffle_study_2=("p_enrich", "last"),
            null_mean_max=("null_mean", "max"),
            comparable_genes=("comparable_genes", "first"),
            shuffleable_genes_min=("shuffleable_genes", "min"),
        )
        .reset_index()
    )
    tests["q_family"] = np.nan
    for bins, group in tests.groupby("baseline_bins"):
        assert len(group) == 8, (bins, len(group))
        tests.loc[group.index, "q_family"] = bh(group.p_both)
    tests["family_definition"] = np.where(
        tests.baseline_bins.eq(5),
        "primary 8 tests: 2 species x 2 thresholds x 2 nulls",
        "q3 sensitivity 8 tests: 2 species x 2 thresholds x 2 nulls",
    )
    return tests


def module_D() -> dict[str, pd.DataFrame]:
    pairs = prepare_D_pairs()
    human_four_parts, human_diag_parts = [], []
    for bins in (5, 3):
        for threshold in (0.02, 0.05):
            for null_model in ("response_identity", "direction_specific"):
                four, diagnostics = permute_human_D(pairs, threshold, null_model, bins)
                human_four_parts.append(four)
                human_diag_parts.append(diagnostics)
    human_four = pd.concat(human_four_parts, ignore_index=True)
    human_diagnostics = pd.concat(human_diag_parts, ignore_index=True)

    old_four = read("old_D_four")
    mouse_four = old_four[old_four.species.eq("mouse")].copy()
    mouse_four["source"] = "reused_unchanged_mouse_permutation"
    combined = pd.concat([mouse_four, human_four], ignore_index=True, sort=False)
    tests = summarize_D_tests(combined)

    old_tests = summarize_D_tests(old_four)
    comparison = old_tests.merge(
        tests,
        on=["species", "axis_id", "threshold", "null_model", "baseline_bins"],
        suffixes=("_old", "_corrected"),
        validate="one_to_one",
    )
    assert comparison.loc[comparison.species.eq("mouse"), "p_both_old"].equals(
        comparison.loc[comparison.species.eq("mouse"), "p_both_corrected"]
    )

    exchange_summary = (
        human_diagnostics.groupby(
            ["threshold", "null_model", "baseline_bins", "shuffled_study", "study_id"], dropna=False
        )
        .agg(
            strata=("stratum", "size"),
            mixed_strata=("mixed_labels", "sum"),
            genes_in_mixed_strata=("n_genes", lambda values: int(values[human_diagnostics.loc[values.index, "mixed_labels"]].sum())),
            expected_labels_changed=("expected_labels_changed", "sum"),
        )
        .reset_index()
    )
    assert exchange_summary.mixed_strata.gt(0).all()
    assert exchange_summary.expected_labels_changed.gt(0).all()

    pattern_counts = (
        pairs[pairs.species.eq("human")]
        .groupby("guan_series_pattern")
        .size()
        .rename("n_genes")
        .reset_index()
    )
    return {
        "pairs": pairs,
        "four": combined,
        "human_exchangeability": human_diagnostics,
        "exchange_summary": exchange_summary,
        "tests": tests,
        "old_vs_corrected": comparison,
        "guan_pattern_counts": pattern_counts,
    }


def write_report(
    frozen: dict[tuple[str, str], pd.DataFrame],
    module_a: dict,
    module_b: dict,
    module_d: dict,
    elapsed: float,
) -> None:
    endpoints = module_a["endpoints"]
    pairs = module_a["pairs"]
    primary_endpoint = endpoints[endpoints.analysis.eq("primary_0.02_0.02")]
    primary_pairs = pairs[pairs.analysis.eq("primary_0.02_0.02")]
    endpoint_summary = (
        primary_endpoint.groupby("species")
        .agg(
            frozen=("gene_id", "size"),
            endpoint_measurable=("endpoint_measurable", "sum"),
            endpoint_retained=("endpoint_retained", "sum"),
        )
        .reset_index()
    )
    process_summary = (
        primary_pairs.groupby(["species", "state", "pair_class"], dropna=False)
        .size()
        .rename("n_genes")
        .reset_index()
    )
    threshold_summary = (
        endpoints[
            endpoints.analysis.isin(
                [
                    "primary_0.02_0.02",
                    "fixed_primary_eligibility_intermediate_0.05",
                    "strong_endpoint_0.05_intermediate_0.05",
                ]
            )
        ]
        .groupby(["analysis", "species"])
        .agg(endpoint_measurable=("endpoint_measurable", "sum"), endpoint_retained=("endpoint_retained", "sum"))
        .reset_index()
    )
    sensitivity = module_a["transition_summary"]
    config = module_b["configurations"]
    config_summary = (
        config.groupby("cohort")
        .agg(
            genes=("gene_id", "size"),
            coarse_configuration_comparable=("coarse_configuration_comparable", "sum"),
            same_coarse_configuration=("same_coarse_configuration", "sum"),
            PM_direction_comparable=("PM_direction_comparable", "sum"),
            same_PM_direction=("same_PM_direction", "sum"),
            same_both=("same_coarse_and_PM_direction", "sum"),
        )
        .reset_index()
    )
    d_primary = module_d["tests"][module_d["tests"].baseline_bins.eq(5)][
        ["species", "threshold", "null_model", "observed", "p_both", "q_family", "enrichment_min", "shuffleable_genes_min"]
    ]
    exchange = module_d["exchange_summary"]
    text = f"""# Corrected fixed-source-reference recurrence audit

## Scope and estimand

This audit preserved the frozen recurrent sets and their original directions ({len(frozen[("mouse", "P_S")])} mouse and {len(frozen[("human", "P_S")])} human P–S genes). For every gene × route × measurement setting, one fixed true-source S estimate was formed by equal weighting of source dates, sampling units within dates and libraries within units. R, E and P retained their observed target windows and each subtracted this same S estimate.

The endpoint and intermediate thresholds were represented independently. `primary_0.02_0.02` is the corrected main analysis. `fixed_primary_eligibility_intermediate_0.05` keeps the main endpoint-eligible genes fixed and changes only the intermediate classification threshold. `strong_endpoint_0.05_intermediate_0.05` is a different, stronger endpoint subset.

## A. Corrected endpoint and process summaries

{markdown_table(endpoint_summary)}

{markdown_table(threshold_summary)}

The corrected same-series process categories are:

{markdown_table(process_summary)}

Human Guan0618 genes with an R estimate but without a measurable same-series P–S contrast were retained in `A_Guan0618_R_background.tsv.gz`; they were not counted as process disagreement and do not enter the same-series consistency denominator.

The member-exclusion reconstruction reproduced the original mouse shared-axis library WPP table before exclusions: `{json.dumps(module_a['membership_reconstruction'], sort_keys=True)}`.

Sensitivity outcomes explicitly distinguish missing measurement, endpoint failure and category change:

{markdown_table(sensitivity)}

## B. Common-window and conditional-model audit

The exact P–S = P–M + M–S triangle was retained. Cross-study coarse configuration identity and P–M direction identity were recorded separately:

{markdown_table(config_summary)}

Conditional models were admitted only when the adjusted coefficient was finite, the design matrix had the expected full rank and the condition number was below 1e6. Within each model, adjusted and unadjusted effects used the identical library set and weights. Cross-model comparisons used only library rows passing all three models.

## D. Guan-series-stratified human permutations

Human permutations were re-run after adding `0618_only`, `0809_only` or `both_measurable` to every matching stratum. Mouse raw permutation results were reused because their inputs and rules were unchanged. BH values below were recomputed over the complete pre-specified eight-test family.

{markdown_table(d_primary)}

Exchangeability after the additional human stratum remained non-zero for every threshold, null model, binning sensitivity and shuffle direction:

{markdown_table(exchange)}

## Interpretation boundary

These analyses audit the frozen candidates under a corrected comparison estimand. They do not reselect recurrent genes and do not replace the original Figure 3A recurrence test. Classification stability is calculated among common measurable genes; loss of measurement is reported as coverage loss rather than direction disagreement.

Runtime: {elapsed:.1f} seconds.
"""
    (OUT / "report.md").write_text(text)


def main() -> None:
    started = time.time()
    frozen = frozen_sets()
    status_rows = []

    module_a = module_A(frozen)
    save(module_a["state_effects"], "A_fixed_reference_state_effects.tsv.gz")
    save(module_a["membership_states"], "A_membership_recomputed_state_values.tsv.gz")
    save(module_a["endpoints"], "A_endpoint_eligibility.tsv")
    save(module_a["pairs"], "A_same_series_process_classes.tsv.gz")
    save(module_a["profiles"], "A_process_profiles.tsv")
    save(module_a["transitions"], "A_sensitivity_gene_state_outcomes.tsv.gz")
    save(module_a["transition_summary"], "A_sensitivity_outcome_summary.tsv")
    save(module_a["flow"], "A_coverage_flow.tsv")
    save(module_a["guan_background"], "A_Guan0618_R_background.tsv.gz")
    status_rows.append({"module": "A", "status": "completed", **module_a["membership_reconstruction"]})

    module_b = module_B(frozen)
    save(module_b["triangle_design"], "B_common_window_triangle_library_rows.tsv.gz")
    save(module_b["triangle_effects"], "B_common_window_study_effects.tsv.gz")
    save(module_b["configurations"], "B_cross_study_configuration_identity.tsv")
    save(module_b["conditional_rows"], "B_conditional_model_library_diagnostics.tsv.gz")
    save(module_b["conditional_summary"], "B_conditional_model_weighted_effects.tsv.gz")
    save(module_b["conditional_cross"], "B_conditional_model_cross_study.tsv")
    save(module_b["conditional_diagnostics"], "B_conditional_model_diagnostic_summary.tsv")
    status_rows.append({"module": "B", "status": "completed"})

    module_d = module_D()
    save(module_d["pairs"], "D_permutation_input_with_Guan_pattern.tsv.gz")
    save(module_d["four"], "D_four_class_nulls.tsv.gz")
    save(module_d["human_exchangeability"], "D_human_exchangeability.tsv.gz")
    save(module_d["exchange_summary"], "D_human_exchangeability_summary.tsv")
    save(module_d["tests"], "D_prespecified_eight_test_family.tsv")
    save(module_d["old_vs_corrected"], "D_old_vs_Guan_stratified.tsv")
    save(module_d["guan_pattern_counts"], "D_Guan_series_pattern_counts.tsv")
    status_rows.append({"module": "D", "status": "completed"})

    elapsed = time.time() - started
    write_report(frozen, module_a, module_b, module_d, elapsed)
    save(pd.DataFrame(status_rows), "analysis_status.tsv")

    inputs = {name: {"path": str(path), "sha256": digest(path)} for name, path in PATHS.items()}
    for path in sorted((STATE_ANALYSIS / "raw").glob("gse*/apa_gex_states.tsv.gz")):
        inputs[f"native_state_{path.parent.name}"] = {"path": str(path), "sha256": digest(path)}
    for study_folder in [
        "gse103221_pluripotency_reprogramming_core",
        "gse115943_pluripotency_reprogramming_core",
    ]:
        membership_paths = {
            "selected_metacells": PROJECT / f"data/apa-majority-shared-v1/profiles/{study_folder}/selected_metacells.tsv.gz",
            "metacells": PROJECT / f"data/metacell-pas-v7-v1/studies/{study_folder}/metacells.tsv.gz",
            "features": PROJECT / f"data/metacell-pas-v7-v1/studies/{study_folder}/features.tsv.gz",
            "counts": PROJECT / f"data/metacell-pas-v7-v1/studies/{study_folder}/metacell_by_pas_counts.npz",
        }
        for role, path in membership_paths.items():
            inputs[f"membership_{study_folder}_{role}"] = {"path": str(path), "sha256": digest(path)}
    (OUT / "input_manifest.json").write_text(json.dumps(inputs, indent=2, sort_keys=True))
    (OUT / "run_config.json").write_text(
        json.dumps(
            {
                "primary_endpoint_threshold": 0.02,
                "primary_intermediate_threshold": 0.02,
                "classification_sensitivity_intermediate_threshold": 0.05,
                "strong_endpoint_threshold": 0.05,
                "permutations": N_PERMUTATIONS,
                "seed": SEED,
                "fixed_source_weighting": "equal dates; equal sample units within dates; equal libraries within units",
                "runtime_seconds": elapsed,
            },
            indent=2,
            sort_keys=True,
        )
    )
    verification = {
        "frozen_sets": {f"{species}_{contrast}": len(frame) for (species, contrast), frame in frozen.items()},
        "fixed_source_signature_max_per_gene_route": int(
            module_a["state_effects"]
            .groupby(["variant", "species", "study_id", "route", "setting", "gene_id"])
            .source_signature.nunique()
            .max()
        ),
        "membership_reconstruction": module_a["membership_reconstruction"],
        "human_exchangeability_all_nonzero": bool(
            module_d["exchange_summary"].mixed_strata.gt(0).all()
            and module_d["exchange_summary"].expected_labels_changed.gt(0).all()
        ),
        "eight_tests_per_BH_family": {
            str(bins): int(len(group)) for bins, group in module_d["tests"].groupby("baseline_bins")
        },
        "completed": True,
    }
    assert verification["fixed_source_signature_max_per_gene_route"] == 1
    assert all(value == 8 for value in verification["eight_tests_per_BH_family"].values())
    (OUT / "verification.json").write_text(json.dumps(verification, indent=2, sort_keys=True))
    print(json.dumps(verification, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
