"""Separate observed-S diagnostic on seven-route jointly measurable genes."""
import hashlib
import json
from pathlib import Path
import pandas as pd
H=Path(__file__).resolve().parent;R=H.parents[2]
V=['exact_Cycle','exact_Adjusted_Cycle','cohort_Cycle','cohort_Adjusted_Cycle']

def save(f,n):
    f.to_csv(H/n,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if n.endswith('.gz') else None)

def main():
    protocol=R/'notes/2026-10-10-cycle-source-attribution/PROTOCOL.zh.md'
    assert hashlib.sha256(protocol.read_bytes()).hexdigest()==json.loads((H/'authorization.json').read_text())['protocol_sha256']
    source=pd.read_csv(H/'route_gene_cycles.tsv.gz',sep='\t',float_precision='round_trip');source=source[source.contrast.eq('P_S')].copy()
    path=H.parent/'sm-core/observed_PS_route_effects.tsv.gz'
    obs=pd.read_csv(path,sep='\t',float_precision='round_trip').rename(columns={'delta_Cycle_exact':'exact_Cycle_delta','delta_Adjusted_Cycle_exact':'exact_Adjusted_Cycle_delta','delta_Cycle_fullcohort':'cohort_Cycle_delta','delta_Adjusted_Cycle_fullcohort':'cohort_Adjusted_Cycle_delta'})
    assert set(obs.study_id)=={'GSE150311','GSE221739'}
    source['documented_source']=True;source['design']='true_source';obs['documented_source']=False;obs['design']='observed_S';obs['contrast']='P_S'
    cols=['species','study_id','route','contrast','design','gene_id','documented_source','delta_DWUI']+[v+'_delta' for v in V]
    assert set(cols).issubset(obs.columns)
    original=pd.read_csv(H/'common_P_S_orthologs.tsv',sep='\t')
    source=pd.concat([q[q.gene_id.isin(original[sp+'_gene_id'])] for sp,q in source.groupby('species')],ignore_index=True)
    f=pd.concat([source[cols],obs[cols]],ignore_index=True)
    assert not f.duplicated(['species','route','gene_id']).any()
    common=original.copy()
    for key,q in obs.groupby(['species','route']):common=common[common[key[0]+'_gene_id'].isin(q.gene_id)]
    assert len(common)>0
    rows=[]
    for key,q in f.groupby(['species','study_id','route','documented_source','design']):
        sp,st,rt,doc,design=key
        for label,z in [('frozen1499_available',q),('seven_route_common',q[q.gene_id.isin(common[sp+'_gene_id'])])]:
            if label=='seven_route_common':assert len(z)==len(common)
            for v in V:
                rows.append(dict(species=sp,study_id=st,route=rt,contrast='P_S',documented_source=doc,design=design,background=label,genes=len(z),cycle_variant=v,mean_DWUI=float(z.delta_DWUI.mean()),median_DWUI=float(z.delta_DWUI.median()),mean_Cycle_delta=float(z[v+'_delta'].mean()),proximal=int(z.delta_DWUI.le(-.02).sum()),distal=int(z.delta_DWUI.ge(.02).sum())))
    summary=pd.DataFrame(rows)
    save(common,'seven_route_PS_common_orthologs.tsv');save(f,'source_and_observed_PS_gene_cycles.tsv.gz');save(summary,'source_and_observed_PS_summary.tsv')
    validation=json.loads((H/'validation.json').read_text());validation.update(observed_source_augmentation='PASS',seven_route_PS_common=len(common),observed_source_isolated_from_primary_model=True,observed_comparisons_rebuilt_from_integer_PAS_with_matched_GEX_cycle_contributors=True)
    (H/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    inputs=json.loads((H/'input_sha256.json').read_text());inputs[str(path.relative_to(R))]=hashlib.sha256(path.read_bytes()).hexdigest();(H/'input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
    print('Seven-route common orthologs:',len(common));print(summary[summary.cycle_variant.eq('exact_Cycle')].to_string(index=False))
if __name__=='__main__':main()
