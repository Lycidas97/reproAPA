"""Independent recomputation from original matched effects and ortholog IDs."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT/"notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables"
ROUTES = ["GSE118258","CIP","WOT"]
STEPS = ["PS","MS","PM"]


def load(path):
    return pd.read_csv(path,sep="\t",float_precision="round_trip")


def near(a,b):
    if pd.isna(a) or pd.isna(b):
        assert pd.isna(a) and pd.isna(b),(a,b)
    else:
        assert abs(float(a)-float(b)) < 1e-12,(a,b)


def main():
    source = json.loads((HERE/"input_sha256.json").read_text())
    for path,expected in source.items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == expected,path
    native = load(BASE/"matched_pms/matched_state_vectors.tsv.gz")
    orth = load(BASE/"cross_overlap/ortholog_matched_features.tsv.gz")
    summary = load(HERE/"route_step_summary.tsv")
    directional = load(HERE/"route_step_direction_contributions.tsv")
    allocations = load(HERE/"direction_number_amplitude_differences.tsv")
    gaps = load(HERE/"route_step_differences.tsv")
    exported = load(HERE/"three_route_ortholog_gene_effects.tsv.gz")
    expected = {}
    export_checks = 0
    for setting in sorted(native.setting.unique()):
        d = {r:native[native.setting.eq(setting)&native.route.eq(r)].set_index("gene_id") for r in ROUTES}
        bg = {}
        bg["route_all_matched"] = {r:q for r,q in d.items()}
        common_mouse = sorted(set(d["CIP"].index)&set(d["WOT"].index))
        bg["mouse_two_route_common"] = {r:d[r].loc[common_mouse] for r in ["CIP","WOT"]}
        for name,routes in [("three_route_ortholog_common",ROUTES),("Xing_CIP_ortholog_common",["GSE118258","CIP"]),
                            ("Xing_WOT_ortholog_common",["GSE118258","WOT"])]:
            keep = pd.Series(True,index=orth.index)
            for route in routes:
                col="human_gene_id" if route == "GSE118258" else "mouse_gene_id"
                keep &= orth[col].isin(d[route].index)
            ids = orth[keep].sort_values("human_gene_id")
            bg[name] = {r:d[r].loc[ids.human_gene_id if r=="GSE118258" else ids.mouse_gene_id] for r in routes}
            if name == "three_route_ortholog_common":
                q=exported[exported.setting.eq(setting)].sort_values("human_gene_id")
                assert q.human_gene_id.tolist() == ids.human_gene_id.tolist()
                assert q.mouse_gene_id.tolist() == ids.mouse_gene_id.tolist()
                for r in routes:
                    for field in ["P_DWUI","M_DWUI","S_DWUI","DWUI_PS","DWUI_MS","DWUI_PM","n_target_libraries","n_triplets"]:
                        assert np.allclose(q[r+"__"+field],bg[name][r][field],atol=1e-12,rtol=0)
                        export_checks += len(q)
        for background,frames in bg.items():
            for route,q in frames.items():
                expected[(setting,background,route)] = {s:q["DWUI_"+s].to_numpy() for s in STEPS} | {
                    s:q[s+"_DWUI"].to_numpy() for s in "PMS"}
            if background != "route_all_matched" and "CIP" in frames and "WOT" in frames:
                # Mixture rows retain each route's half-mass, including directions.
                expected[(setting,background,"CIP_WOT_equal")] = {
                    s:np.concatenate([expected[(setting,background,r)][s] for r in ["CIP","WOT"]]) for s in STEPS+list("PMS")}
    scalar_checks=0
    for row in summary.itertuples():
        v=expected[(row.setting,row.background,row.route)]
        n=len(v["PS"])/(2 if row.equal_route_group else 1)
        near(row.denominator,n)
        for field in STEPS+list("PMS"):
            near(getattr(row,"mean_"+field),v[field].mean())
            scalar_checks+=1
        near(row.mean_PS,row.mean_MS+row.mean_PM)
    direction_checks=0
    for row in directional.itertuples():
        v=expected[(row.setting,row.background,row.route)][row.contrast]
        mask={"proximal":v<=-row.threshold,"subthreshold":np.abs(v)<row.threshold,"distal":v>=row.threshold}[row.direction]
        q=v[mask]
        divisor=2 if row.equal_route_group else 1
        near(row.count,len(q)/divisor)
        near(row.denominator,len(v)/divisor)
        near(row.fraction,mask.mean())
        near(row.conditional_signed_mean,q.mean() if len(q) else np.nan)
        near(row.conditional_abs_mean,np.abs(q).mean() if len(q) else np.nan)
        near(row.signed_contribution,q.sum()/len(v))
        near(row.abs_contribution,np.abs(q).sum()/len(v))
        direction_checks+=7
    direction_index=directional.set_index(["setting","background","route","threshold","contrast","direction"])
    allocation_checks=0
    for row in allocations.itertuples():
        prefix=(row.setting,row.background)
        suffix=(row.threshold,row.contrast,row.direction)
        h=direction_index.loc[prefix+(row.reference_route,)+suffix]
        m=direction_index.loc[prefix+(row.comparator_route,)+suffix]
        number=(m.fraction-h.fraction)*(m.conditional_signed_mean+h.conditional_signed_mean)/2
        amplitude=(m.conditional_signed_mean-h.conditional_signed_mean)*(m.fraction+h.fraction)/2
        near(row.number_contribution,number)
        near(row.amplitude_contribution,amplitude)
        near(row.difference_signed_contribution,m.signed_contribution-h.signed_contribution)
        if np.isfinite(number) and np.isfinite(amplitude):
            near(row.difference_signed_contribution,number+amplitude)
        assert row.symmetric_allocation_identified == bool(np.isfinite(number) and np.isfinite(amplitude))
        allocation_checks+=4
    gap_checks=0
    for row in gaps.itertuples():
        h=expected[(row.setting,row.background,row.reference_route)]
        m=expected[(row.setting,row.background,row.comparator_route)]
        delta={s:m[s].mean()-h[s].mean() for s in STEPS}
        for s in STEPS:
            near(getattr(row,"difference_"+s),delta[s])
        near(row.difference_PS,row.difference_MS+row.difference_PM)
        near(row.step_contrast_PM_minus_MS,delta["PM"]-delta["MS"])
        stable=abs(delta["PS"])>=.005
        assert row.PS_gap_ratio_stable == stable
        for s in ["MS","PM"]:
            near(getattr(row,s+"_signed_fraction_of_PS_gap"),delta[s]/delta["PS"] if stable else np.nan)
        gap_checks+=7
    null=json.loads((HERE/"null_identifiability_audit.json").read_text())
    assert null["informative"] is False and null["p"] is None and null["q"] is None
    assert null["route_labels_exchangeable"] is False and null["species_inference_identifiable"] is False
    draws=load(HERE/"null_audit_draws.tsv")
    near(draws.difference_PS.max(),draws.difference_PS.min())
    near(draws.step_contrast_PM_minus_MS.max(),draws.step_contrast_PM_minus_MS.min())
    near(null["observed_difference_PS"],draws.difference_PS.iloc[0])
    near(null["observed_step_contrast_PM_minus_MS"],draws.step_contrast_PM_minus_MS.iloc[0])
    assignment=load(HERE/"route_label_resolution_audit.tsv")
    assert len(assignment)==3 and assignment.single_route.nunique()==3
    assert (~assignment.valid_biological_permutation).all()
    near(null["minimum_one_sided_resolution_if_exchangeable"],1/3)
    result=dict(status="passed",source_hashes=len(source),export_scalar_checks=export_checks,
        summary_scalar_checks=scalar_checks,direction_scalar_checks=direction_checks,
        allocation_scalar_checks=allocation_checks,gap_scalar_checks=gap_checks,
        all_primary_and_sensitivity_backgrounds_reconstructed=True,constant_null_correctly_noninformative=True,
        no_species_biological_p_value=True,skipped=0)
    (HERE/"independent_validation.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
