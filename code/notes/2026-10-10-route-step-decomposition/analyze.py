"""Finite-background accounting of matched state contrasts; no species p value."""
import ast
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / "notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables"
SETTING = "supported50_p08"
ROUTES = ("GSE118258", "CIP", "WOT")
LABELS = {"GSE118258": "Xing2020", "CIP": "CIP", "WOT": "WOT/2i"}
STEPS = ("PS", "MS", "PM")
THRESHOLDS = (.02, .05)
SOURCES = {}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path, **kwargs):
    SOURCES[str(path.relative_to(ROOT))] = digest(path)
    return pd.read_csv(path, sep="\t", float_precision="round_trip", **kwargs)


def save(frame, name):
    compression = {"method": "gzip", "mtime": 0} if name.endswith(".gz") else None
    pd.DataFrame(frame).to_csv(HERE / name, sep="\t", index=False, na_rep="NA",
                              float_format="%.17g", compression=compression)


def load_weights():
    path = ROOT / "analysis_code/pas_rerun/polyapipe_effects.py"
    SOURCES[str(path.relative_to(ROOT))] = digest(path)
    tree = ast.parse(path.read_text())
    functions = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "weights"]
    assert len(functions) == 1
    namespace = dict(np=np, pd=pd)
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), "exec"), namespace)
    return namespace["weights"]


def validate_matched(frame, triplets, axes):
    """Check addition independently of gene admission, then replay primary weights."""
    checks = []
    assert not frame.duplicated(["species", "route", "setting", "gene_id"]).any()
    assert set(frame.route) == set(ROUTES)
    for metric in ["DWUI", "RWUI", "PDUI", "PPUI"]:
        error = float(np.abs(frame[metric+"_PS"] - frame[metric+"_MS"] - frame[metric+"_PM"]).max())
        assert error < 1e-12
        checks.append(dict(check=metric+"_additivity_all_settings", n=len(frame), max_abs_error=error))
    vector_error = 0.
    for row in frame.itertuples():
        states = {s: np.asarray(json.loads(getattr(row, s+"_usage"))) for s in "PMS"}
        assert all(abs(v.sum()-1) < 1e-12 for v in states.values())
        assert all(len(v) == row.PAS_number for v in states.values())
        for step, a, b in [("PS", "P", "S"), ("MS", "M", "S"), ("PM", "P", "M")]:
            vector_error = max(vector_error, float(np.abs(
                np.asarray(json.loads(getattr(row, "delta_usage_"+step))) - states[a]+states[b]).max()))
            vector_error = max(vector_error, abs(getattr(row, "DWUI_"+step) - getattr(row,a+"_DWUI") + getattr(row,b+"_DWUI")))
    assert vector_error < 1e-12
    checks.append(dict(check="full_vectors_and_state_indices_all_settings", n=len(frame), max_abs_error=vector_error))
    keys = ["species", "study_id", "route", "setting", "gene_id"]
    primary = frame[frame.setting.eq(SETTING)].set_index(keys)
    triplets = triplets[triplets.setting.eq(SETTING)]
    assert set(triplets.set_index(keys).index) == set(primary.index)
    assert triplets.design.eq("true_source").all()
    assert triplets.contrast.eq("P_S").all()
    axis = {(s,g): part.set_index("atlas_pas_id").WPP for (s,g),part in axes.groupby(["study_id", "primary_gene_id"],sort=False)}
    weights = load_weights()
    error = 0.
    for key, part in triplets.groupby(keys, sort=False):
        row = primary.loc[key]
        w = weights(part)
        assert np.allclose(w, part.triplet_weight, rtol=0, atol=1e-12)
        assert len(part) == row.n_triplets
        assert part.library_key_a.nunique() == row.n_target_libraries
        assert part.library_key_b.nunique() == row.n_source_libraries
        assert part.time_a.nunique() == row.n_target_dates
        assert part.MS_joint_PAS.ge(2).all()
        assert part.PAS_ids.eq(part.M_PAS_ids).all()
        assert part.target_usage.eq(part.PM_P_usage).all()
        ids = json.loads(row.PAS_ids)
        position = axis[(key[1],key[-1])].reindex(ids).to_numpy(float)
        assert np.isfinite(position).all()
        for state, column in [("P","target_usage"),("M","M_usage"),("S","reference_usage")]:
            estimate = w @ np.stack(part[column].map(json.loads))
            error = max(error, float(np.abs(estimate-np.asarray(json.loads(row[state+"_usage"]))).max()),
                        abs(float(estimate@position) - row[state+"_DWUI"]))
    assert error < 1e-12
    checks.append(dict(check="primary_original_triplets_weights_and_PAS_projection", n=len(primary), max_abs_error=error))
    return checks


def direction_summary(values, threshold):
    """Every measured effect enters exactly one of three signed contributions."""
    values = np.asarray(values, float)
    assert len(values) and np.isfinite(values).all()
    masks = {"proximal": values <= -threshold, "subthreshold": np.abs(values) < threshold,
             "distal": values >= threshold}
    assert np.all(sum(masks.values()) == 1)
    rows = []
    for direction, mask in masks.items():
        q = values[mask]
        rows.append(dict(direction=direction, denominator=len(values), count=len(q),
            fraction=float(mask.mean()), conditional_signed_mean=float(q.mean()) if len(q) else np.nan,
            conditional_abs_mean=float(np.abs(q).mean()) if len(q) else np.nan,
            signed_contribution=float(q.sum()/len(values)),
            abs_contribution=float(np.abs(q).sum()/len(values))))
    assert abs(sum(x["signed_contribution"] for x in rows)-values.mean()) < 1e-12
    return rows


def combine_direction_rows(route_rows):
    """Equal route mass, including the zero mass of an absent direction class."""
    result = []
    assert len(set(r["denominator"] for rows in route_rows for r in rows)) == 1
    for i in range(3):
        q = [rows[i] for rows in route_rows]
        assert len(set(x["direction"] for x in q)) == 1
        p = float(np.mean([x["fraction"] for x in q]))
        signed = float(np.mean([x["signed_contribution"] for x in q]))
        absolute = float(np.mean([x["abs_contribution"] for x in q]))
        result.append(dict(direction=q[0]["direction"], denominator=q[0]["denominator"],
            count=float(np.mean([x["count"] for x in q])), fraction=p,
            conditional_signed_mean=signed/p if p else np.nan,
            conditional_abs_mean=absolute/p if p else np.nan,
            signed_contribution=signed, abs_contribution=absolute))
    return result


def symmetric_difference(reference, comparator):
    p0, p1 = reference["fraction"], comparator["fraction"]
    mu0, mu1 = reference["conditional_signed_mean"], comparator["conditional_signed_mean"]
    delta = comparator["signed_contribution"]-reference["signed_contribution"]
    identified = bool(np.isfinite(mu0) and np.isfinite(mu1))
    number = (p1-p0)*(mu1+mu0)/2 if identified else np.nan
    amplitude = (mu1-mu0)*(p1+p0)/2 if identified else np.nan
    if identified:
        assert abs(number+amplitude-delta) < 1e-12
    same_n = reference["denominator"] == comparator["denominator"]
    return dict(reference_count=reference["count"], comparator_count=comparator["count"],
        reference_denominator=reference["denominator"], comparator_denominator=comparator["denominator"],
        reference_fraction=p0, comparator_fraction=p1, delta_fraction=p1-p0,
        equal_denominator_count_difference=comparator["count"]-reference["count"] if same_n else np.nan,
        reference_conditional_mean=mu0, comparator_conditional_mean=mu1,
        delta_conditional_mean=mu1-mu0 if identified else np.nan,
        reference_signed_contribution=reference["signed_contribution"],
        comparator_signed_contribution=comparator["signed_contribution"],
        difference_signed_contribution=delta, number_contribution=number,
        amplitude_contribution=amplitude, symmetric_allocation_identified=identified)


def ortholog_background(frame, orthologs, routes):
    """Join frozen identities only; inherited ortholog response fields are ignored."""
    out = orthologs[["human_gene_id", "mouse_gene_id"]].copy()
    for route in routes:
        q = frame[frame.route.eq(route)]
        column = "human_gene_id" if route == "GSE118258" else "mouse_gene_id"
        fields = ["gene_id", "gene_symbol", "study_id", "PAS_number", "n_triplets", "n_target_libraries",
                  "n_source_libraries", "n_target_dates"] + [s+"_DWUI" for s in "PMS"] + ["DWUI_"+s for s in STEPS]
        out = out.merge(q[fields].rename(columns={"gene_id":column} | {c:route+"__"+c for c in fields if c!="gene_id"}),
                        on=column, validate="one_to_one")
    assert out.human_gene_id.is_unique and out.mouse_gene_id.is_unique
    assert len(out)
    return out.sort_values("human_gene_id").reset_index(drop=True)


def wide_to_routes(wide, routes):
    out = {}
    for route in routes:
        q = wide[[c for c in wide if c.startswith(route+"__")]].copy()
        q.columns = [c.split("__",1)[1] for c in q]
        q["gene_id"] = wide.human_gene_id if route == "GSE118258" else wide.mouse_gene_id
        q["route"] = route
        out[route] = q
    return out


def add_summary_background(frames, setting, background, summaries, directions, differences, allocations, common):
    """A background can contain disjoint native genes; its scope is explicit."""
    means = {}
    detail = {}
    for route, frame in frames.items():
        assert frame.gene_id.is_unique
        means[route] = {step: float(frame["DWUI_"+step].mean()) for step in STEPS}
        assert abs(means[route]["PS"]-means[route]["MS"]-means[route]["PM"]) < 1e-12
        context = dict(setting=setting, background=background, route=route, route_label=LABELS[route],
                       equal_route_group=False, denominator=len(frame), identical_cross_route_gene_background=common)
        summaries.append(context | {"mean_"+step:means[route][step] for step in STEPS} |
                         {"mean_"+s:float(frame[s+"_DWUI"].mean()) for s in "PMS"})
        for threshold in THRESHOLDS:
            for step in STEPS:
                detail[(route,threshold,step)] = direction_summary(frame["DWUI_"+step],threshold)
                directions.extend(context | dict(threshold=threshold,contrast=step) | r for r in detail[(route,threshold,step)])
    if "CIP" in frames and "WOT" in frames and common:
        route = "CIP_WOT_equal"
        means[route] = {s: (means["CIP"][s]+means["WOT"][s])/2 for s in STEPS}
        context = dict(setting=setting,background=background,route=route,route_label="CIP + WOT equal-route",
                       equal_route_group=True,denominator=len(frames["CIP"]),identical_cross_route_gene_background=True)
        summaries.append(context | {"mean_"+s:means[route][s] for s in STEPS} |
                         {"mean_"+s:float((frames["CIP"][s+"_DWUI"].mean()+frames["WOT"][s+"_DWUI"].mean())/2) for s in "PMS"})
        for threshold in THRESHOLDS:
            for step in STEPS:
                detail[(route,threshold,step)] = combine_direction_rows([detail[(r,threshold,step)] for r in ["CIP","WOT"]])
                directions.extend(context | dict(threshold=threshold,contrast=step) | r for r in detail[(route,threshold,step)])
    if "GSE118258" not in frames:
        return
    for route in means:
        if route == "GSE118258":
            continue
        delta = {s: means[route][s]-means["GSE118258"][s] for s in STEPS}
        assert abs(delta["PS"]-delta["MS"]-delta["PM"]) < 1e-12
        stable = abs(delta["PS"]) >= .005
        mass = abs(delta["MS"])+abs(delta["PM"])
        context = dict(setting=setting,background=background,reference_route="GSE118258",comparator_route=route,
                       identical_cross_route_gene_background=common)
        differences.append(context | {"difference_"+s:delta[s] for s in STEPS} |
            dict(step_contrast_PM_minus_MS=delta["PM"]-delta["MS"],PS_gap_ratio_stable=stable,
                MS_signed_fraction_of_PS_gap=delta["MS"]/delta["PS"] if stable else np.nan,
                PM_signed_fraction_of_PS_gap=delta["PM"]/delta["PS"] if stable else np.nan,
                MS_signed_fraction_of_step_abs_mass=delta["MS"]/mass if mass else np.nan,
                PM_signed_fraction_of_step_abs_mass=delta["PM"]/mass if mass else np.nan))
        for threshold in THRESHOLDS:
            for step in STEPS:
                refs = detail[("GSE118258",threshold,step)]
                comps = detail[(route,threshold,step)]
                values = [symmetric_difference(a,b) for a,b in zip(refs,comps)]
                assert abs(sum(x["difference_signed_contribution"] for x in values)-delta[step]) < 1e-12
                allocations.extend(context | dict(threshold=threshold,contrast=step,direction=a["direction"],
                    fraction_of_PS_gap=x["difference_signed_contribution"]/delta["PS"] if stable else np.nan,
                    PS_gap_ratio_stable=stable) | x for a,x in zip(refs,values))


def null_audit(wide, covariates):
    """Demonstrate why an identity null does not test an all-background mean."""
    frame = wide.copy()
    q = covariates[covariates.threshold.eq(.02)]
    columns = []
    edges = {}
    for route, unit, species in [("GSE118258","GSE118258","human"),("CIP","GSE103221","mouse"),("WOT","GSE115943","mouse")]:
        gene = "human_gene_id" if species == "human" else "mouse_gene_id"
        native = q[q.unit.eq(unit)].copy()
        native["log_GEX"] = np.log1p(native.mean_GEX)
        native["log_APA"] = np.log1p(native.APA_coverage)
        cov = native.groupby("gene_id").agg(mean_log_GEX=("log_GEX","mean"),mean_log_APA=("log_APA","mean"),native_PAS=("PAS_number","mean"))
        frame = frame.merge(cov.add_prefix(route+"__"),left_on=gene,right_index=True,validate="one_to_one")
    for species, gene in [("human","human_gene_id"),("mouse","mouse_gene_id")]:
        native = q[q.species.eq(species)]
        sig = native.groupby("gene_id").apply(lambda a: ";".join(
            f"{r.contrast}:{r.unit}:{r.series_pattern}" for r in a.sort_values(["contrast","unit"]).itertuples()),include_groups=False)
        col = species+"_actual_availability"
        frame[col] = frame[gene].map(sig)
        assert frame[col].notna().all()
        columns.append(col)
    for route in ROUTES:
        for feature in ["S_DWUI","M_DWUI","mean_log_GEX","mean_log_APA","native_PAS","n_target_libraries"]:
            col = route+"__"+feature
            values = frame[col].to_numpy(float)
            assert np.isfinite(values).all()
            inner = np.unique(np.quantile(values,[0,.5,1]))[1:-1]
            frame[col+"_bin"] = np.searchsorted(inner,values,side="left")
            columns.append(col+"_bin")
            edges[col] = inner.tolist()
    cells = list(frame.groupby(columns,sort=True,dropna=False).indices.values())
    frame["stratum"] = -1
    for i, cell in enumerate(cells):
        frame.loc[cell,"stratum"] = i
    values = {r: frame[[r+"__DWUI_"+s for s in STEPS]].to_numpy(float) for r in ROUTES}
    def statistic(arrays):
        h = arrays["GSE118258"].mean(0)
        m = (arrays["CIP"].mean(0)+arrays["WOT"].mean(0))/2
        delta = m-h
        return delta[0], delta[2]-delta[1]
    observed = statistic(values)
    rng = np.random.default_rng(20261010)
    draws = []
    for i in range(100):
        shuffled = {r:v.copy() for r,v in values.items()}
        for r in ROUTES:
            for cell in cells:
                shuffled[r][cell] = values[r][rng.permutation(cell)]
        ps, interaction = statistic(shuffled)
        draws.append(dict(draw=i,difference_PS=ps,step_contrast_PM_minus_MS=interaction))
    error = max(abs(r["difference_PS"]-observed[0]) for r in draws)
    error = max(error,max(abs(r["step_contrast_PM_minus_MS"]-observed[1]) for r in draws))
    assert error < 1e-12
    save(frame,"null_audit_background.tsv.gz")
    save(draws,"null_audit_draws.tsv")
    assignments = []
    for hroute in ROUTES:
        others = [r for r in ROUTES if r != hroute]
        h = values[hroute].mean(0)
        m = (values[others[0]].mean(0)+values[others[1]].mean(0))/2
        assignments.append(dict(single_route=hroute,two_route_group=";".join(others),
            observed_species_assignment=hroute=="GSE118258",difference_PS=m[0]-h[0],
            step_contrast_PM_minus_MS=(m[2]-h[2])-(m[1]-h[1]),
            valid_biological_permutation=False))
    save(assignments,"route_label_resolution_audit.tsv")
    result = dict(status="complete",primary_ortholog_pairs=len(frame),strata=len(cells),
        singleton_strata=sum(len(c)==1 for c in cells),movable_genes=sum(len(c) for c in cells if len(c)>1),
        joint_matching_columns=columns,matching_edges=edges,permutations=100,
        observed_difference_PS=observed[0],observed_step_contrast_PM_minus_MS=observed[1],
        max_abs_permutation_change=error,informative=False,p=None,q=None,
        reason="An all-background within-route complete-vector identity permutation preserves every route/step marginal sum, with or without matching.",
        route_label_assignments=3,minimum_one_sided_resolution_if_exchangeable=1/3,
        route_labels_exchangeable=False,species_inference_identifiable=False)
    (HERE/"null_identifiability_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    return dict(check="matched_complete_vector_permutation_invariance",n=100,max_abs_error=error)


def main():
    for path in [HERE/"protocol.json",HERE/"analyze.py",ROOT/"analysis_code/pas_rerun/downstream_cohorts.py"]:
        SOURCES[str(path.relative_to(ROOT))] = digest(path)
    matched = read(BASE/"matched_pms/matched_state_vectors.tsv.gz")
    keys = ["species","study_id","route","setting","gene_id"]
    fields = keys+["design","contrast","PAS_ids","M_PAS_ids","target_usage","PM_P_usage","reference_usage","M_usage",
                  "MS_joint_PAS","triplet_weight","library_key_a","library_key_b","time_a","time_b","sample_unit_a","sample_unit_b"]
    triplets = read(BASE/"matched_pms/weighted_triplets.tsv.gz",usecols=fields)
    axes = read(BASE/"state_vectors/PAS_axes.primary.tsv.gz",usecols=["study_id","primary_gene_id","atlas_pas_id","WPP"])
    checks = validate_matched(matched,triplets,axes)
    del triplets,axes
    orthologs = read(BASE/"cross_overlap/ortholog_matched_features.tsv.gz",usecols=["human_gene_id","mouse_gene_id","homology_type","is_high_confidence","human_partner_count","mouse_partner_count"])
    assert orthologs.human_gene_id.is_unique and orthologs.mouse_gene_id.is_unique
    assert orthologs.homology_type.eq("ortholog_one2one").all()
    assert orthologs.is_high_confidence.eq(1).all()
    assert orthologs.human_partner_count.eq(1).all() and orthologs.mouse_partner_count.eq(1).all()
    endpoints = read(BASE/"state_vectors/route_effects.primary.tsv.gz",usecols=["species","study_id","route","contrast","design","gene_id"])
    availability = endpoints.groupby(["species","study_id","route","contrast","design"]).gene_id.nunique().reset_index(name="genes")
    save(availability,"native_endpoint_availability.tsv")
    summaries,directions,differences,allocations,coverage,allwide = [],[],[],[],[],[]
    primary_wide = None
    for setting, frame in matched.groupby("setting",sort=True):
        allroutes = {r:frame[frame.route.eq(r)].copy() for r in ROUTES}
        add_summary_background(allroutes,setting,"route_all_matched",summaries,directions,differences,allocations,False)
        mouse_common = set(allroutes["CIP"].gene_id)&set(allroutes["WOT"].gene_id)
        add_summary_background({r:d[d.gene_id.isin(mouse_common)].copy() for r,d in allroutes.items() if r in ["CIP","WOT"]},
                               setting,"mouse_two_route_common",summaries,directions,differences,allocations,True)
        for background, routes in [("three_route_ortholog_common",ROUTES),("Xing_CIP_ortholog_common",("GSE118258","CIP")),
                                   ("Xing_WOT_ortholog_common",("GSE118258","WOT"))]:
            wide = ortholog_background(frame,orthologs,routes)
            add_summary_background(wide_to_routes(wide,routes),setting,background,summaries,directions,differences,allocations,True)
            coverage.append(dict(setting=setting,background=background,ortholog_pairs=len(wide)))
            if background == "three_route_ortholog_common":
                allwide.append(wide.assign(setting=setting))
                if setting == SETTING:
                    primary_wide = wide
        coverage.extend(dict(setting=setting,background="route_all_matched",route=r,measured_genes=len(d)) for r,d in allroutes.items())
        coverage.append(dict(setting=setting,background="mouse_two_route_common",measured_genes=len(mouse_common)))
    save(matched[[c for c in matched if not c.endswith("usage") and "MPRO" not in c]],"all_matched_gene_state_effects.tsv.gz")
    save(pd.concat(allwide,ignore_index=True),"three_route_ortholog_gene_effects.tsv.gz")
    save(coverage,"background_coverage.tsv")
    save(summaries,"route_step_summary.tsv")
    save(directions,"route_step_direction_contributions.tsv")
    save(differences,"route_step_differences.tsv")
    save(allocations,"direction_number_amplitude_differences.tsv")
    covariates = read(BASE/"recurrence_controls/recurrence_frozen_covariates.tsv.gz")
    checks.append(null_audit(primary_wide,covariates))
    save(checks,"verification.tsv")
    (HERE/"input_sha256.json").write_text(json.dumps(SOURCES,indent=2,sort_keys=True)+"\n")
    (HERE/"validation.json").write_text(json.dumps(dict(status="complete",matched_rows=len(matched),
        primary_matched_rows=int(matched.setting.eq(SETTING).sum()),primary_three_route_ortholog_pairs=len(primary_wide),
        route_step_rows=len(summaries),direction_contribution_rows=len(directions),number_amplitude_rows=len(allocations),
        all_arithmetic_checks_passed=True,biological_species_test_performed=False,manuscript_changed=False),indent=2)+"\n")
    print(json.dumps(dict(status="complete",primary_three_route_ortholog_pairs=len(primary_wide),
        primary_summary=pd.DataFrame(summaries).query("setting == @SETTING and background == 'three_route_ortholog_common'").to_dict("records")),indent=2))


if __name__ == "__main__":
    main()
