"""Only read qualification/contributor fields; never inspect new effect or p columns."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MARKER=ROOT/"notes/2026-10-10-marker-reference-sensitivity"
OLD=ROOT/"notes/2026-10-10-ortholog-threshold-scan"
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
read=lambda p,cols:pd.read_csv(p,sep="\t",usecols=cols,float_precision="round_trip")
STUDIES={"GSE103221":"mouse","GSE115943":"mouse","GSE118258":"human","GSE150311":"human","GSE178325":"human","GSE221739":"human"}


def main():
    sources=[MARKER/"unique_support_bridge_route_effects.tsv.gz",MARKER/"bridge_route_effects.tsv.gz",OLD/"fixed_ortholog_background.tsv",OLD/"study_events.tsv.gz",MARKER/"unique_support_pair_weights.tsv.gz",MARKER/"inputs.npz",OLD/"measured_ortholog_members.tsv.gz"]
    hashes={str(p.relative_to(ROOT)):sha(p) for p in sources}
    (HERE/"input_sha256.json").write_text(json.dumps(hashes,indent=2)+"\n")
    cols=["study_id","route","gene_id","contrast","background","n_libraries","contributor_libraries","weights","PAS_ids","version","setting"]
    q=read(sources[0],cols);q=q[q.background.eq("common_M_bridge")]
    before=read(sources[1],cols);before=before[before.background.eq("common_M_bridge")]
    key=["study_id","route","gene_id","version","setting"]
    pairs=q.pivot(index=key,columns="contrast",values="n_libraries")
    assert set(pairs)=={"D00","D01"} and pairs.notna().all().all() and pairs.D00.eq(pairs.D01).all()
    d0=q[q.contrast.eq("D00")].set_index(key).sort_index();d1=q[q.contrast.eq("D01")].set_index(key).sort_index()
    assert d0.index.equals(d1.index)
    for c in ["contributor_libraries","weights","PAS_ids"]:assert d0[c].equals(d1[c]),c
    assert q[cols].sort_values(key+["contrast"]).reset_index(drop=True).equals(before[cols].sort_values(key+["contrast"]).reset_index(drop=True))
    opportunity=d0.reset_index().copy()
    opportunity["species"]=opportunity.study_id.map(STUDIES)
    assert opportunity.species.notna().all()
    original=read(sources[2],["contrast","human_gene_id","mouse_gene_id"])
    assert original.groupby("contrast").size().to_dict()=={"P_M":751,"P_S":1656}
    old_pm=set(original[original.contrast.eq("P_M")].human_gene_id)
    old_ps=set(original[original.contrast.eq("P_S")].human_gene_id)
    old_labels=read(sources[6],["contrast","threshold","human_gene_id","human_direction","mouse_direction"])
    core=set(old_labels[old_labels.contrast.eq("P_M")&old_labels.threshold.eq(.02)&old_labels.human_direction.eq("proximal")&old_labels.mouse_direction.eq("proximal")].human_gene_id)
    assert len(core)==27
    events=read(sources[3],["species","contrast","threshold","unit","gene_id","series_pattern"])
    events=events[events.threshold.eq(.02)&events.contrast.eq("P_M")]
    patterns={(r.species,r.gene_id):[] for r in events.itertuples(index=False)}
    for r in events.itertuples(index=False):patterns[(r.species,r.gene_id)].append((r.unit,r.series_pattern))
    rows=[];routes=[]
    for (version,setting),part in opportunity.groupby(["version","setting"]):
        for (sp,sid,route),f in part.groupby(["species","study_id","route"]):
            routes.append(dict(version=version,setting=setting,species=sp,study_id=sid,route=route,matched_route_genes=len(f),library_union=len({lib for v in f.contributor_libraries for lib in json.loads(v)}),max_libraries_per_gene=int(f.n_libraries.max())))
        signatures={}
        for (sp,gene),f in part.groupby(["species","gene_id"]):
            signature=[]
            for sid,t in f.groupby("study_id"):
                signature.append((sid,";".join(sorted(t.route))))
            signatures[(sp,gene)]=tuple(sorted(signature))
        for r in original.itertuples(index=False):
            out=dict(version=version,setting=setting,old_background=r.contrast,human_gene_id=r.human_gene_id,mouse_gene_id=r.mouse_gene_id,inside_old_PM=r.human_gene_id in old_pm,inside_old_PS=r.human_gene_id in old_ps,original_PM_core27=r.human_gene_id in core)
            for sp,gene in [("human",r.human_gene_id),("mouse",r.mouse_gene_id)]:
                new=signatures.get((sp,gene),())
                old=tuple(sorted(patterns.get((sp,gene),[])))
                out[sp+"_new_measured_studies"]=len(new)
                out[sp+"_new_measured_pattern"]=";".join(u+":"+p for u,p in new)
                out[sp+"_old_PM_measured_pattern"]=";".join(u+":"+p for u,p in old)
                out[sp+"_old_PM_opportunity_intact"]=bool(old) and new==old
            out["both_species_ge2"]=min(out["human_new_measured_studies"],out["mouse_new_measured_studies"])>=2
            out["human_ge2"]=out["human_new_measured_studies"]>=2
            out["all_old_PM_opportunities_intact"]=out["human_old_PM_opportunity_intact"] and out["mouse_old_PM_opportunity_intact"]
            rows.append(out)
    members=pd.DataFrame(rows)
    members.to_csv(HERE/"qualification_members.tsv",sep="\t",index=False)
    pd.DataFrame(routes).to_csv(HERE/"route_opportunity.tsv",sep="\t",index=False)
    summary=members.groupby(["version","setting","old_background"]).agg(fixed_pairs=("human_gene_id","size"),human_ge2=("human_ge2","sum"),both_species_ge2=("both_species_ge2","sum"),all_old_PM_opportunities_intact=("all_old_PM_opportunities_intact","sum"),inside_old_PM=("inside_old_PM","sum"),inside_old_PS=("inside_old_PS","sum")).reset_index()
    summary.to_csv(HERE/"qualification_summary.tsv",sep="\t",index=False)
    z=np.load(MARKER/"inputs.npz",allow_pickle=False)
    gex=[k for k in z.files if "gex" in k.lower()]
    assert all(k.endswith("_total_gex") for k in gex)
    for p,h in hashes.items():assert sha(ROOT/p)==h,p
    qa=dict(status="QUALIFICATION_ONLY_PASS_NO_NEW_SHARING_TEST",qualified_route_contexts=len(opportunity),matched_D00_D01_contributors_identical=True,unique_support_bridge_qualification_unchanged=True,old_background_sizes={"P_M":751,"P_S":1656},full_target_gene_GEX_available=False,loaded_GEX_fields=gex,new_effect_columns_read=False,new_direction_labels_computed=False,new_null_or_p_values_computed=False)
    (HERE/"QA.json").write_text(json.dumps(qa,indent=2)+"\n")
    print(summary.to_string(index=False))
    print(pd.DataFrame(routes).query("version=='core_third_strict' and setting=='supported50'").to_string(index=False))


if __name__=="__main__":main()
