"""Compare independent replay against root and available production tables; no null."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
TASK=HERE.parents[1]
ANALYSIS=TASK/"analysis"
read=lambda p:pd.read_csv(p,sep="\t",float_precision="round_trip")
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def compare(a,b,keys,columns):
    aa=a.set_index(keys).sort_index();bb=b.set_index(keys).sort_index()
    assert aa.index.equals(bb.index)
    errors={}
    for c in columns:
        if pd.api.types.is_bool_dtype(aa[c]):
            assert aa[c].equals(bb[c]),c
        elif pd.api.types.is_numeric_dtype(aa[c]):
            np.testing.assert_allclose(aa[c],bb[c],atol=1e-12,rtol=0,equal_nan=True)
            errors[c]=float(np.max(abs(aa[c]-bb[c]))) if len(aa) else 0
        else:assert aa[c].equals(bb[c]),c
    return dict(rows=len(aa),columns=len(columns),max_numeric_difference=max(errors.values(),default=0),numeric_differences=errors)


def main():
    paths=[HERE/"paired_labels_independent.tsv",HERE/"study_events_independent.tsv",TASK/"root_observed_labels.tsv",TASK/"root_observed_3x3.tsv",TASK/"protocol.json"]
    independent=read(paths[0]);root=read(paths[2]);keys=["definition","threshold","human_gene_id","mouse_gene_id"]
    results={"root_labels":compare(independent,root,keys,["human_direction","mouse_direction"])}
    counts=independent.groupby(["definition","threshold","human_direction","mouse_direction"]).size()
    grid=pd.MultiIndex.from_product([independent.definition.unique(),sorted(independent.threshold.unique()),["proximal","distal","nonrecurrent"],["proximal","distal","nonrecurrent"]],names=["definition","threshold","human_direction","mouse_direction"])
    expected=counts.reindex(grid,fill_value=0).rename("observed").reset_index()
    results["root_3x3"]=compare(expected,read(paths[3]),grid.names,["observed"])
    production_members=ANALYSIS/"measured_ortholog_members.tsv.gz"
    production_events=ANALYSIS/"study_events.tsv.gz"
    if production_members.exists():
        production=read(production_members)
        labels=read(HERE/"recurrence_labels_independent.tsv")
        protocol=json.loads(paths[4].read_text())
        independent_full=[]
        for definition,source in protocol["definitions"].items():
            q=independent[independent.definition.eq(definition)].copy()
            for sp in ["human","mouse"]:
                z=labels[labels.species.eq(sp)&labels.contrast.eq(source[sp])]
                columns=["proximal_support","distal_support","measured_studies"]
                z=z[["gene_id","threshold"]+columns].rename(columns={"gene_id":sp+"_gene_id",**{c:sp+"_"+c for c in columns}})
                q=q.merge(z,on=[sp+"_gene_id","threshold"],validate="many_to_one")
                q[sp+"_opposing_study_allowed"]=(q[sp+"_proximal_support"].gt(0)&q[sp+"_distal_support"].gt(0))
            independent_full.append(q)
        independent_full=pd.concat(independent_full,ignore_index=True)
        columns=[sp+"_"+c for sp in ["human","mouse"] for c in ["direction","proximal_support","distal_support","measured_studies","opposing_study_allowed"]]
        results["production_labels"]=compare(independent_full,production,keys,columns)
        paths.append(production_members)
        core=read(HERE/"original_core27_flow_independent.tsv")[["human_gene_id","mouse_gene_id","human_gene_symbol","mouse_gene_symbol"]]
        contexts=independent[["definition","threshold"]].drop_duplicates()
        expected_core=contexts.merge(core,how="cross").merge(independent_full,on=keys,how="left",validate="one_to_one")
        expected_core["measurable"]=expected_core.human_direction.notna()
        expected_core["status"]=np.where(~expected_core.measurable,"not_estimable_background_loss",np.where(expected_core.human_direction.eq("proximal")&expected_core.mouse_direction.eq("proximal"),"retained_shared_proximal","measurable_changed_class"))
        flow_path=ANALYSIS/"original_core27_flow.tsv"
        if flow_path.exists():
            results["production_core27_flow"]=compare(expected_core,read(flow_path),keys,["human_direction","mouse_direction","measurable","status","human_gene_symbol","mouse_gene_symbol"])
            paths.append(flow_path)
        shared_path=ANALYSIS/"shared_proximal_primary_members.tsv"
        if shared_path.exists():
            expected_shared=independent_full[independent_full.threshold.eq(.02)&independent_full.human_direction.eq("proximal")&independent_full.mouse_direction.eq("proximal")].copy()
            expected_shared["original_core27"]=expected_shared.human_gene_id.isin(core.human_gene_id)
            results["production_primary_shared_members"]=compare(expected_shared,read(shared_path),keys,columns+["original_core27"])
            paths.append(shared_path)
    if production_events.exists():
        e=read(paths[1]);protocol=json.loads(paths[4].read_text());parts=[]
        for definition,source in protocol["definitions"].items():
            q=e[[r.contrast==source[r.species] for r in e.itertuples(index=False)]].copy();q["definition"]=definition;parts.append(q)
        expected_events=pd.concat(parts,ignore_index=True)
        production=read(production_events)
        ek=["definition","threshold","species","unit","gene_id"]
        results["production_events"]=compare(expected_events,production,ek,["event","series_pattern","effect_mean","effect_min","effect_max","baseline","capacity_negative","capacity_positive"])
        paths.append(production_events)
    qa=dict(status="PRODUCTION_EVENTS_AND_ROOT_LABELS_PASS" if "production_labels" in results and "production_events" in results else "ROOT_LABELS_PASS_PRODUCTION_PENDING",comparisons=results,new_null_computed=False,input_sha256={str(p.relative_to(TASK)):sha(p) for p in paths})
    (HERE/"comparison_QA.json").write_text(json.dumps(qa,indent=2)+"\n")
    print(json.dumps(qa,indent=2))


if __name__=="__main__":main()
