"""Descriptive common-PAS comparison with read-only original metacell matrices."""
from pathlib import Path
import json

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

from expanded_measurements import HERE, OUT, BASE, ENDPOINT, save
from matched_analysis import digest, read

PROFILE_BASE = Path("/root/nfsdata/apadetect/analysis_projects/mouse_time_aware_integration_pilot_v1/results_apa_majority_shared_v1/profiles")
ROUTES = {"gse103221_pluripotency_reprogramming_core": "CIP",
          "gse115943_pluripotency_reprogramming_core": "WOT"}


def endpoint_mean(records, role):
    subset = [r for r in records if r["role"] == role]
    times = sorted({r["time"] for r in subset})
    means = [np.mean([r["proportions"] for r in subset if r["time"] == time], axis=0) for time in times]
    return np.mean(means, axis=0)


def run():
    frozen = read(HERE / "inputs/frozen_candidates.tsv")
    symbols = frozen.drop_duplicates("gene_id").set_index("gene_id").gene_symbol.to_dict()
    external = read(OUT / "expanded_site_mapping.tsv.gz")
    external = external[external.gene_id.isin(symbols)]
    hash_inputs = {"vector-protocol.json": digest(HERE / "vector-protocol.json"),
        "vector_results/expanded_site_mapping.tsv.gz": digest(OUT / "expanded_site_mapping.tsv.gz")}
    rows, site_rows, native_checks = [], [], []
    for study, route in ROUTES.items():
        directory = BASE / study
        feature_path = directory / "primary_features.tsv.gz"
        meta_path = directory / "selected_metacells.tsv.gz"
        matrix_path = directory / "primary_mc.npz"
        rows_path = directory / "matrix_rows.tsv.gz"
        expected_path = PROFILE_BASE / study / "supported50_library_positions.tsv.gz"
        for path in [feature_path, meta_path, matrix_path, rows_path, expected_path]:
            hash_inputs[str(path)] = digest(path)
        features = read(feature_path)
        meta = read(meta_path)
        matrix = sparse.load_npz(matrix_path).tocsr()
        matrix_rows = read(rows_path)
        assert matrix.shape == (len(matrix_rows), len(features))
        assert meta.matrix_row.is_unique
        assert matrix_rows.metacell_id.iloc[meta.matrix_row].to_list() == meta.metacell_id.to_list()
        assert meta.purity.ge(.8).all()
        expected = read(expected_path)
        expected = expected[expected.route.eq(route) & expected.gene_id.isin(symbols)].set_index(["node", "library_key", "gene_id"]).position_index
        all_groups = [(key, group) for key, group in meta[meta.route.eq(route)].groupby(["node", "library_key", "experimental_time"])
                      if len(group) >= 3 and group.n_cells.sum() >= 30]
        feature_groups = {g: v for g, v in features.groupby("primary_gene_id")}
        external_groups = {g: v for g, v in external.groupby("gene_id")}
        for gene in symbols:
            row = {"study": study, "route": route, "gene_id": gene, "gene_symbol": symbols[gene], "status": "no_reference_or_external_identity"}
            if gene not in feature_groups or gene not in external_groups:
                rows.append(row)
                continue
            all_features = feature_groups[gene]
            f = all_features[all_features.supported50].copy()
            if not len(f) or all_features.relative.isna().any():
                row["status"] = "no_valid_native_supported_axis"
                rows.append(row)
                continue
            counts = matrix[:, f.index.to_numpy()].toarray()
            total = counts.sum(axis=1)
            good = total >= 20
            qualified_route = np.concatenate([group.matrix_row.to_numpy(int) for _, group in all_groups])
            contributing = qualified_route[good[qualified_route]]
            route_counts = counts[contributing]
            if ((route_counts.sum(axis=0) >= 10) & ((route_counts > 0).sum(axis=0) >= 2)).sum() < 2:
                row["status"] = "original_route_pas_support_unavailable"
                rows.append(row)
                continue
            records = []
            for (node, library, time), group in all_groups:
                ix = group.matrix_row.to_numpy(int)
                ix = ix[good[ix]]
                if len(ix) < 3:
                    continue
                p = (counts[ix] / total[ix, None]).mean(axis=0)
                native = p @ f.relative.to_numpy()
                key = (node, library, gene)
                assert key in expected.index and np.isfinite(expected.loc[key]), (study, key)
                error = abs(native - expected.loc[key])
                assert error <= 1e-10, (study, key, error)
                native_checks.append({"study": study, "gene_id": gene, "node": node,
                    "library_key": library, "absolute_error": error})
                if node == "P" or (node == "S" and time == "D0"):
                    records.append({"role": "target" if node == "P" else "source",
                        "time": time, "library": library, "proportions": p, "qualified_metacells": len(ix)})
            if not any(r["role"] == "source" for r in records) or not any(r["role"] == "target" for r in records):
                row["status"] = "no_qualified_original_source_or_target"
                rows.append(row)
                continue
            ext = external_groups[gene]
            pairs = []
            for item in ext.itertuples():
                if (item.coordinate_base_sensitive or item.mapped_chrom != f.chrom.iloc[0]
                    or item.orientation != f.strand.iloc[0] or item.source_strand != f.strand.iloc[0]):
                    continue
                d = (f.cut_site - item.mapped_cut).abs()
                a = (f.cut_site - item.alternate_cut).abs()
                if (d.min() <= 50 and a.min() <= 50 and d.eq(d.min()).sum() == 1
                    and a.eq(a.min()).sum() == 1 and d.idxmin() == a.idxmin()):
                    pairs.append({"COORDS": item.COORDS, "feature_index": d.idxmin()})
            pair = pd.DataFrame(pairs, columns=["COORDS", "feature_index"])
            pair = pair[~pair.feature_index.duplicated(False)]
            row["n_common_pas"] = len(pair)
            if len(pair) < 2:
                row["status"] = "fewer_than_two_common_pas"
                rows.append(row)
                continue
            ordered = ext.set_index("COORDS").loc[pair.COORDS]
            positions = features.loc[pair.feature_index, "cut_site"].to_numpy() * (1 if f.strand.iloc[0] == "+" else -1)
            positions = positions - positions.min()
            span = np.ptp(positions)
            ext_counts = ordered[ENDPOINT].to_numpy().T
            ext_total = ext_counts.sum(axis=1)
            external_fraction = ext_total / ext[ENDPOINT].sum().to_numpy()
            feature_indices = f.index.get_indexer(pair.feature_index)
            assert (feature_indices >= 0).all()
            original_fraction = min(r["proportions"][feature_indices].sum() for r in records)
            row.update(shared_span_bp=span, min_external_common_fraction=external_fraction.min(),
                min_original_common_fraction=original_fraction,
                min_external_common_count=ext_total.min(),
                n_original_source_libraries=sum(r["role"] == "source" for r in records),
                n_original_target_libraries=sum(r["role"] == "target" for r in records))
            overlapping = ordered.pos.sort_values().diff().min() <= 80
            supported = (ext_counts.sum(axis=0) >= 10).sum() >= 2
            if span < 100 or overlapping or not supported or ext_total.min() < 20:
                row["status"] = "external_common_site_coverage_or_windows_fail"
                rows.append(row)
                continue
            if external_fraction.min() < .8 or original_fraction < .8:
                row["status"] = "common_pas_mass_below_80percent"
                rows.append(row)
                continue
            source = endpoint_mean(records, "source")[feature_indices]
            target = endpoint_mean(records, "target")[feature_indices]
            source, target = source / source.sum(), target / target.sum()
            ext_p = ext_counts / ext_total[:, None]
            ext_source, ext_target = ext_p[:2].mean(axis=0), ext_p[2:].mean(axis=0)
            original_delta, external_delta = target - source, ext_target - ext_source
            norm = np.linalg.norm(original_delta) * np.linalg.norm(external_delta)
            row.update(status="comparable", original_tvd=np.abs(original_delta).sum()/2,
                external_tvd=np.abs(external_delta).sum()/2,
                original_delta_axis=(original_delta @ positions)/span,
                external_delta_axis=(external_delta @ positions)/span,
                original_delta_bp=original_delta @ positions,
                external_delta_bp=external_delta @ positions,
                pas_change_cosine=(original_delta @ external_delta)/norm if norm else np.nan,
                same_position_direction=bool((original_delta @ positions)*(external_delta @ positions) > 0))
            rows.append(row)
            for j, coord in enumerate(pair.COORDS):
                site_rows.append({"study": study, "gene_id": gene, "gene_symbol": symbols[gene],
                    "COORDS": coord, "atlas_pas_id": features.loc[pair.feature_index.iloc[j], "atlas_pas_id"],
                    "local_cut_bp": positions[j], "original_source_usage": source[j],
                    "original_target_usage": target[j], "external_source_usage": ext_source[j],
                    "external_target_usage": ext_target[j], "original_usage_change": original_delta[j],
                    "external_usage_change": external_delta[j]})
        del matrix
        print('Compared original', study, flush=True)
    all_rows = pd.DataFrame(rows)
    candidates = frozen.merge(all_rows, on="gene_id", suffixes=("", "_reference"), validate="many_to_many")
    assert len(candidates) == 618 * 2
    save(candidates, "original_comparison_candidates.tsv")
    save(pd.DataFrame(site_rows), "original_comparison_pas.tsv")
    save(pd.DataFrame(native_checks), "original_native_reconstruction_checks.tsv.gz")
    summary = []
    for (name, study), group in candidates.groupby(["set_name", "study"]):
        measured = group[group.status.eq("comparable")]
        rho = spearmanr(measured.original_delta_axis, measured.external_delta_axis).statistic if len(measured) >= 3 else np.nan
        summary.append({"set_name": name, "study": study, "frozen_n": len(group), "comparable_n": len(measured),
            "same_position_direction_n": int(measured.same_position_direction.sum()),
            "descriptive_spearman_common_axis": rho, "median_pas_change_cosine": measured.pas_change_cosine.median(),
            "median_original_tvd": measured.original_tvd.median(), "median_external_tvd": measured.external_tvd.median()})
    save(pd.DataFrame(summary), "original_comparison_summary.tsv")
    receipt = {"status": "completed", "script_sha256": digest(Path(__file__)), "inputs_sha256": hash_inputs,
        "native_library_gene_measurements_reconstructed": len(native_checks),
        "maximum_native_absolute_error": max(r["absolute_error"] for r in native_checks),
        "all_frozen_member_study_records": len(candidates), "scope": "Descriptive common-PAS endpoint correspondence, not per-gene validation significance"}
    (OUT / "comparison_run.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    run()
