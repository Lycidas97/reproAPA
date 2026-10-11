"""Independent library-to-study events replay; no import from production analysis."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
PARENT=HERE.parent
TASK=HERE.parents[1]
MARKER=ROOT/"notes/2026-10-10-marker-reference-sensitivity"
OLD=ROOT/"notes/2026-10-10-ortholog-threshold-scan"
read=lambda p,**kw:pd.read_csv(p,sep="\t",float_precision="round_trip",**kw)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
SIDS={"GSE103221":"mouse","GSE115943":"mouse","GSE118258":"human","GSE150311":"human","GSE221739":"human"}
checked=0


def check(condition,n=1):
    global checked
    assert condition
    checked+=n


def weights_independent(q):
    # Equal child groups at each successive sampling level, then equal captures.
    cols=["time_a","sample_unit_a","time_b","sample_unit_b"]
    values=[tuple(v) for v in q[cols].itertuples(index=False,name=None)]
    result=[]
    for v in values:
        denominator=sum(x==v for x in values)
        for depth in range(1,5):
            denominator*=len({x[:depth] for x in values if x[:depth-1]==v[:depth-1]})
        result.append(1/denominator)
    return np.asarray(result)


def event(deltas,t):
    # Mean does not define a study vote; all series must cross the threshold.
    return -1 if max(deltas)<=-t else 1 if min(deltas)>=t else 0


def class_from_events(values):
    near=sum(v==-1 for v in values);far=sum(v==1 for v in values)
    check(not(near>=2 and far>=2))
    return "proximal" if near>=2 else "distal" if far>=2 else "nonrecurrent",near,far


def save(d,name):
    pd.DataFrame(d).to_csv(HERE/name,sep="\t",index=False,na_rep="NA",float_format="%.17g")


def main():
    global checked
    files=[TASK/"PROTOCOL.zh.md",TASK/"protocol.json",TASK/"protocol_freeze_sha256.json",PARENT/"qualification_members.tsv",MARKER/"unique_support_pair_weights.tsv.gz",MARKER/"pair_effects.tsv.gz",MARKER/"unique_support_bridge_route_effects.tsv.gz",OLD/"measured_ortholog_members.tsv.gz"]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in files}
    (HERE/"input_sha256.json").write_text(json.dumps(hashes,indent=2)+"\n")
    proto=json.loads((TASK/"protocol.json").read_text());thresholds=proto["thresholds"]
    qualification=read(PARENT/"qualification_members.tsv")
    pairs=qualification[qualification.old_background.eq("P_M")&qualification.version.eq("core_third_strict")&qualification.setting.eq("supported50")&qualification.both_species_ge2].copy()
    check(len(pairs)==728 and pairs.human_gene_id.nunique()==728 and pairs.mouse_gene_id.nunique()==728)
    wanted={"human":set(pairs.human_gene_id),"mouse":set(pairs.mouse_gene_id)}
    w=read(MARKER/"unique_support_pair_weights.tsv.gz")
    w=w[w.version.eq("core_third_strict")&w.setting.eq("supported50")&w.background.eq("common_M_bridge")]
    w=w[[g in wanted[SIDS[sid]] for sid,g in zip(w.study_id,w.gene_id)]].copy()
    key=["study_id","route","gene_id","library_key","contrast"]
    check(not w.duplicated(key).any())
    needed=set(w[key].itertuples(index=False,name=None))
    columns=key+["setting","version","PAS_ids","position","target_usage","reference_usage","delta_DWUI"]
    chunks=[]
    for q in read(MARKER/"pair_effects.tsv.gz",usecols=columns,chunksize=100000):
        q=q[q.version.eq("core_third_strict")&q.setting.eq("supported50")&q.contrast.isin(["D00","D01"])]
        chunks.append(q[[v in needed for v in q[key].itertuples(index=False,name=None)]])
    library=pd.concat(chunks,ignore_index=True).merge(w[key+["weight","time_a","sample_unit_a","time_b","sample_unit_b"]],on=key,validate="one_to_one")
    check(len(library)==len(needed)==len(w))
    for c in ["position","target_usage","reference_usage"]:library[c]=library[c].map(lambda v:np.asarray(json.loads(v),float))
    library_delta=np.array([(r.target_usage-r.reference_usage)@r.position for r in library.itertuples(index=False)])
    max_library_error=float(np.max(abs(library_delta-library.delta_DWUI)))
    check(max_library_error<1e-12,len(library));library["recomputed_delta"]=library_delta
    routes=[];max_weight_error=0
    for (sid,route,gene),q in library.groupby(["study_id","route","gene_id"]):
        a=q[q.contrast.eq("D00")].sort_values("library_key");b=q[q.contrast.eq("D01")].sort_values("library_key")
        check(a.library_key.tolist()==b.library_key.tolist())
        check(np.array_equal(a.weight,b.weight))
        check(all(np.array_equal(x,y) for x,y in zip(a.target_usage,b.target_usage)),len(a))
        expected=weights_independent(a)
        err=float(np.max(abs(expected-a.weight.to_numpy())));max_weight_error=max(max_weight_error,err)
        check(err<1e-12 and abs(expected.sum()-1)<1e-12,len(a))
        for contrast,p in [("D00",a),("D01",b)]:
            check(p.PAS_ids.nunique()==1 and all(np.array_equal(x,p.position.iloc[0]) for x in p.position))
            delta=float(np.dot(expected,p.recomputed_delta))
            baseline=float(np.dot(expected,[u@pos for u,pos in zip(p.reference_usage,p.position)]))
            routes.append(dict(species=SIDS[sid],study_id=sid,route=route,gene_id=gene,contrast=contrast,delta_DWUI=delta,baseline=baseline,n_libraries=len(p),contributor_libraries=json.dumps(p.library_key.tolist()),weights=json.dumps(expected.tolist()),PAS_ids=p.PAS_ids.iloc[0]))
    route=pd.DataFrame(routes);save(route,"route_effects_independent.tsv")
    accepted=read(MARKER/"unique_support_bridge_route_effects.tsv.gz",usecols=["study_id","route","gene_id","contrast","version","setting","background","delta_DWUI"])
    accepted=accepted[accepted.version.eq("core_third_strict")&accepted.setting.eq("supported50")&accepted.background.eq("common_M_bridge")]
    comparison=route.merge(accepted,on=["study_id","route","gene_id","contrast"],suffixes=("_independent","_accepted"),validate="one_to_one")
    max_route_error=float(np.max(abs(comparison.delta_DWUI_independent-comparison.delta_DWUI_accepted)))
    check(len(comparison)==len(route) and max_route_error<1e-12,len(route))
    events=[]
    for (sp,sid,gene,contrast),q in route.groupby(["species","study_id","gene_id","contrast"]):
        check(q.route.nunique()==len(q))
        for t in thresholds:
            events.append(dict(species=sp,unit=sid,gene_id=gene,contrast=contrast,threshold=t,series_pattern=";".join(sorted(q.route)),n_routes=len(q),effect_mean=q.delta_DWUI.mean(),effect_min=q.delta_DWUI.min(),effect_max=q.delta_DWUI.max(),baseline=q.baseline.median(),capacity_negative=q.baseline.min(),capacity_positive=(1-q.baseline).min(),event=event(q.delta_DWUI.tolist(),t)))
    events=pd.DataFrame(events);save(events,"study_events_independent.tsv")
    labels=[]
    for (sp,gene,contrast,t),q in events.groupby(["species","gene_id","contrast","threshold"]):
        check(len(q)>=2)
        direction,near,far=class_from_events(q.event.tolist())
        labels.append(dict(species=sp,gene_id=gene,contrast=contrast,threshold=t,measured_studies=len(q),direction=direction,proximal_support=near,distal_support=far,measured_pattern=";".join(q.unit+":"+q.series_pattern)))
    labels=pd.DataFrame(labels);save(labels,"recurrence_labels_independent.tsv")
    memberships=[]
    for definition,co in proto["definitions"].items():
        for t in thresholds:
            z=pairs[["human_gene_id","mouse_gene_id"]].copy()
            for sp in ["human","mouse"]:
                roster=labels[labels.species.eq(sp)&labels.contrast.eq(co[sp])&labels.threshold.eq(t)].set_index("gene_id")
                z[sp+"_direction"]=z[sp+"_gene_id"].map(roster.direction)
            check(z[["human_direction","mouse_direction"]].notna().all().all())
            z["definition"]=definition;z["threshold"]=t;memberships.append(z)
    memberships=pd.concat(memberships,ignore_index=True);save(memberships,"paired_labels_independent.tsv")
    old=read(OLD/"measured_ortholog_members.tsv.gz")
    old=old[old.contrast.eq("P_M")&old.threshold.eq(.02)&old.human_direction.eq("proximal")&old.mouse_direction.eq("proximal")]
    check(len(old)==27)
    flow=old[["human_gene_id","mouse_gene_id","human_gene_symbol","mouse_gene_symbol"]].copy()
    flow["joint_measurable"]=flow.human_gene_id.isin(pairs.human_gene_id)
    check(flow.joint_measurable.sum()==25)
    check(set(flow[~flow.joint_measurable].human_gene_symbol)=={"GNPNAT1","KLC1"})
    for definition in proto["definitions"]:
        z=memberships[memberships.definition.eq(definition)&memberships.threshold.eq(.02)].set_index("human_gene_id")
        for sp in ["human","mouse"]:flow[definition+"_"+sp+"_direction"]=flow.human_gene_id.map(z[sp+"_direction"])
        flow[definition+"_shared_proximal"]=flow.joint_measurable&flow[definition+"_human_direction"].eq("proximal")&flow[definition+"_mouse_direction"].eq("proximal")
    save(flow,"original_core27_flow_independent.tsv")
    shared=memberships[memberships.human_direction.eq("proximal")&memberships.mouse_direction.eq("proximal")].copy()
    shared["member_of_original_core27"]=shared.human_gene_id.isin(old.human_gene_id)
    save(shared,"shared_members_independent.tsv")
    # Intent cases: same mean can conceal a non-voting series; third opposition allowed.
    check(event([-.04,-.01],.02)==0);check(event([-.04,-.02],.02)==-1)
    check(class_from_events([-1,-1,1])[0]=="proximal")
    for p,h in hashes.items():check(sha(ROOT/p)==h)
    qa=dict(status="INDEPENDENT_EVENTS_REBUILT_PRODUCTION_COMPARISON_PENDING",checks=checked,failed=0,library_effects=len(library),route_effects=len(route),study_events=len(events),background_pairs=728,original_core27_measurable=25,original_core27_unmeasurable=["GNPNAT1","KLC1"],max_library_delta_error=max_library_error,max_route_delta_error=max_route_error,max_sampling_weight_error=max_weight_error,new_null_computed=False,production_comparison_complete=False)
    (HERE/"QA.json").write_text(json.dumps(qa,indent=2)+"\n")
    print(json.dumps(qa,indent=2))
    print(flow.groupby("joint_measurable")[[d+"_shared_proximal" for d in proto["definitions"]]].sum().to_string())
    print(shared.groupby(["definition","threshold"]).size().to_string())


if __name__=="__main__":main()
