"""Recover exact group requests from sealed source profiles before scalar GEX review."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
TASK=HERE.parents[1]
MARKER=ROOT/"notes/2026-10-10-marker-reference-sensitivity"
ANALYSIS=TASK/"analysis"
read=lambda p,**k:pd.read_csv(p,sep="\t",float_precision="round_trip",**k)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    qualification=read(HERE.parent/"qualification_members.tsv")
    bg=qualification[qualification.old_background.eq("P_M")&qualification.version.eq("core_third_strict")&qualification.setting.eq("supported50")&qualification.both_species_ge2]
    wanted=set(bg.human_gene_id)|set(bg.mouse_gene_id)
    w=read(MARKER/"unique_support_pair_weights.tsv.gz")
    w=w[w.version.eq("core_third_strict")&w.setting.eq("supported50")&w.background.eq("common_M_bridge")&w.gene_id.isin(wanted)]
    key=["study_id","route","gene_id","library_key","contrast"]
    keys=set(w[key].itertuples(index=False,name=None))
    pieces=[]
    columns=key+["setting","version","target_group","reference_group"]
    for c in read(MARKER/"pair_effects.tsv.gz",usecols=columns,chunksize=100000):
        c=c[c.version.eq("core_third_strict")&c.setting.eq("supported50")]
        pieces.append(c[[v in keys for v in c[key].itertuples(index=False,name=None)]])
    pairs=pd.concat(pieces,ignore_index=True)
    assert len(pairs)==len(keys)==30470 and not pairs.duplicated(key).any()
    pairs.to_csv(HERE/"library_group_ids_independent.tsv",sep="\t",index=False)
    ids=set(pairs.target_group)|set(pairs.reference_group)
    cols=["species","study_id","route","library_key","setting","gene_id","group_id","cohort","status","n_mc","n_cells","qualified_MC_ids","qualified_MC_pas_depths","PAS_ids","position","pooled_molecules","usage"]
    pieces=[]
    for c in read(MARKER/"gene_group_profiles.tsv.gz",usecols=cols,chunksize=100000):pieces.append(c[c.group_id.isin(ids)])
    groups=pd.concat(pieces,ignore_index=True).sort_values("group_id")
    assert len(groups)==len(ids)==45705 and groups.group_id.is_unique and groups.status.eq("available").all()
    producer=read(ANALYSIS/"qualified_groups.tsv.gz",usecols=cols).sort_values("group_id")
    assert groups[cols].reset_index(drop=True).equals(producer[cols].reset_index(drop=True))
    meta=read(MARKER/"metacell_metadata.tsv.gz",usecols=["study_id","metacell_id","total_gex","n_cells","library_key"]).set_index(["study_id","metacell_id"])
    derived=[]
    for r in groups.itertuples(index=False):
        mids=json.loads(r.qualified_MC_ids);depths=json.loads(r.qualified_MC_pas_depths)
        assert len(set(mids))==len(mids)==r.n_mc and len(depths)==r.n_mc and min(depths)>=20
        assert sum(depths)==sum(json.loads(r.pooled_molecules))
        q=meta.loc[[(r.study_id,mid) for mid in mids]]
        assert q.n_cells.sum()==r.n_cells and q.library_key.eq(r.library_key).all()
        derived.append(dict(group_id=r.group_id,all_gene_GEX_sum_independent=int(q.total_gex.sum()),APA_coverage_group=sum(depths)/r.n_mc))
    derived=pd.DataFrame(derived);groups=groups.merge(derived,on="group_id",validate="one_to_one")
    groups.to_csv(HERE/"qualified_groups_independent.tsv.gz",sep="\t",index=False,float_format="%.17g",compression={"method":"gzip","mtime":0})
    request=read(ANALYSIS/"gex_requests.tsv.gz")
    rcols=request.columns.tolist()
    assert groups[rcols].sort_values("group_id").reset_index(drop=True).equals(request.sort_values("group_id").reset_index(drop=True))
    inputs=[MARKER/"pair_effects.tsv.gz",MARKER/"gene_group_profiles.tsv.gz",MARKER/"unique_support_pair_weights.tsv.gz",MARKER/"metacell_metadata.tsv.gz",ANALYSIS/"qualified_groups.tsv.gz",ANALYSIS/"gex_requests.tsv.gz"]
    (HERE/"group_input_sha256.json").write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in inputs},indent=2)+"\n")
    qa=dict(status="GROUP_REQUEST_REPLAY_PASS",groups=len(groups),library_pairs=len(pairs),all_requests_equal_original_profiles=True,group_denominators_independently_recovered=True,remote_GEX_executed_here=False)
    (HERE/"group_request_QA.json").write_text(json.dumps(qa,indent=2)+"\n")
    print(json.dumps(qa,indent=2))


if __name__=="__main__":main()
