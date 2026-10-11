"""Export target/reference WPP for the frozen endpoint recount.

Run inside the frozen remote endpoint-counting project after all library recounts
pass.  It mirrors estimate.py but retains both endpoint WPP values so the
Figure 3 reference-WPP matched recurrence null can be rerun.
"""
import argparse

from common import *


def collapse_wpp(pairs):
    keys = ["species", "study_id", "route", "view", "mode", "gene_id"]
    values = ["target_WPP", "reference_WPP"]
    x = pairs.groupby(
        keys + ["time_a", "sample_unit_a", "time_b", "sample_unit_b"]
    )[values].mean()
    x = x.groupby(keys + ["time_a", "sample_unit_a", "time_b"])[values].mean()
    x = x.groupby(keys + ["time_a", "sample_unit_a"])[values].mean()
    dates = x.groupby(keys + ["time_a"])[values].mean().reset_index()
    effects = dates.groupby(keys)[values].mean().reset_index()
    effects["delta"] = effects.target_WPP - effects.reference_WPP
    return effects, dates


def export(excluded_library=None):
    output = HERE / ("provisional_results" if excluded_library else "results")
    output.mkdir(exist_ok=True)
    all_pairs = []
    excluded_routes = []
    for study, folder_name in STUDIES.items():
        dest = HERE / "data" / folder_name
        features = read(dest / "features.tsv.gz")
        genes = pd.Index(features.primary_gene_id.unique())
        meta = read(dest / "metacells.tsv.gz")
        if excluded_library is not None and excluded_library in set(meta.library_key):
            excluded_routes.extend(
                (study, route)
                for route in sorted(
                    meta.loc[meta.library_key.eq(excluded_library), "route"].unique()
                )
            )
            meta = meta[meta.library_key.ne(excluded_library)].copy()
        allowed = json.loads((dest / "gex_gene_ids.json").read_text())
        states = {}
        for key, frame in meta.groupby("library_key", sort=False):
            folder = dest / "libraries" / key
            frame = frame.sort_values("local_row")
            matrices = {}
            good = []
            for mode in MODES:
                audit = json.loads((folder / "recount/audit.json").read_text())
                assert audit["status"] == "PASS"
                matrix = sparse.load_npz(folder / "recount" / (mode + "_mc.npz")).tocsr()
                assert matrix.shape == (len(frame), len(features))
                matrices[mode] = matrix
                stats, valid = library_stats(matrix, features, frame, allowed[key], mode)
                states["reapplied", mode, key] = stats
                good.append(valid)
            common = np.logical_and.reduce(good)
            for mode, matrix in matrices.items():
                states["common_contributors", mode, key] = library_stats(
                    matrix, features, frame, allowed[key], mode, common
                )[0]

        libmeta = meta.groupby("library_key", sort=False).first().reset_index()
        rows = []
        for route, route_meta in libmeta.groupby("route"):
            sources = route_meta[route_meta.node.eq("S")]
            targets = route_meta[route_meta.node.eq("P")]
            if excluded_library is not None and (not len(sources) or not len(targets)):
                continue
            assert len(sources) and len(targets)
            for target in targets.itertuples():
                for source in sources.itertuples():
                    for view in ["reapplied", "common_contributors"]:
                        comparisons = {
                            mode: paired_stats(
                                states[view, mode, target.library_key],
                                states[view, mode, source.library_key],
                                features,
                            )
                            for mode in MODES
                        }
                        common_valid = np.logical_and.reduce(
                            [valid for _, valid in comparisons.values()]
                        )
                        for mode, (_, valid) in comparisons.items():
                            if view == "common_contributors":
                                valid = common_valid
                            target_values = states[view, mode, target.library_key]["APA"]
                            source_values = states[view, mode, source.library_key]["APA"]
                            frame = pd.DataFrame(
                                {
                                    "gene_id": genes[valid],
                                    "target_WPP": target_values[valid],
                                    "reference_WPP": source_values[valid],
                                }
                            )
                            frame = frame.assign(
                                species=target.species,
                                study_id=study,
                                route=route,
                                view=view,
                                mode=mode,
                                library_key_a=target.library_key,
                                library_key_b=source.library_key,
                                time_a=target.experimental_time,
                                time_b=source.experimental_time,
                                sample_unit_a=target.sample_unit,
                                sample_unit_b=source.sample_unit,
                            )
                            rows.append(frame)
        pairs = pd.concat(rows, ignore_index=True)
        effects, _ = collapse_wpp(pairs)
        eligible = (
            effects[effects.view.eq("reapplied")]
            .groupby(["route", "gene_id"])["mode"]
            .nunique()
        )
        eligible = eligible[eligible.eq(len(MODES))].reset_index()[["route", "gene_id"]]
        fixed_pairs = pairs[pairs.view.eq("reapplied")].merge(
            eligible, on=["route", "gene_id"], validate="many_to_one"
        )
        fixed_pairs["view"] = "fixed_genes"
        all_pairs.append(pd.concat([pairs, fixed_pairs], ignore_index=True))

    pairs = pd.concat(all_pairs, ignore_index=True)
    effects, dates = collapse_wpp(pairs)
    reference = read(HERE / "inputs/qualified_route_effects.tsv.gz")
    reference = reference[
        reference.setting.eq("supported50_p08")
        & reference.design.eq("true_source")
        & reference.contrast.eq("P_S")
    ]
    if excluded_routes:
        excluded_index = pd.MultiIndex.from_tuples(
            excluded_routes, names=["study_id", "route"]
        )
        reference_index = pd.MultiIndex.from_frame(reference[["study_id", "route"]])
        reference = reference[~reference_index.isin(excluded_index)]
    keys = ["species", "study_id", "route", "gene_id"]
    rebuilt = effects[
        effects.view.eq("reapplied") & effects["mode"].eq("original")
    ].set_index(keys)
    expected = reference.set_index(keys)
    assert rebuilt.index.equals(expected.index)
    error = float((rebuilt.delta - expected.delta).abs().max())
    assert error < 1e-8, error
    save(effects, output / "route_wpp.tsv.gz")
    save(dates, output / "date_wpp.tsv.gz")
    save(pairs, output / "pair_wpp.tsv.gz")
    dump(
        {
            "status": "PROVISIONAL_EXCLUDED_LIBRARY" if excluded_library else "PASS",
            "excluded_library": excluded_library,
            "excluded_routes": excluded_routes,
            "route_rows": len(effects),
            "date_rows": len(dates),
            "pair_rows": len(pairs),
            "original_effect_max_error": error,
            "fixed_axis": True,
            "modes": MODES,
        },
        output / "route_wpp_audit.json",
    )
    print(
        "ROUTE_WPP_EXPORT_PROVISIONAL"
        if excluded_library
        else "ROUTE_WPP_EXPORT_PASS",
        len(effects),
        error,
        flush=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--exclude-library")
    args = parser.parse_args()
    export(args.exclude_library)
