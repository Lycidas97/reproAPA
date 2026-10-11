"""Independent bridge effects from original pair usage and source hierarchy."""
from pathlib import Path
import hashlib,json

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
MARKER=ROOT/'notes/2026-10-10-marker-reference-sensitivity'
PAIRKEY=['study_id','route','library_key','gene_id','version','setting','contrast']
ROUTEKEY=['study_id','route','gene_id','contrast']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def hierarchy_weights(frame):
    # A balanced hierarchy is the intent; this recursive implementation does
    # not call the production product-of-cardinalities weight function.
    result=pd.Series(0.,index=frame.index)
    def walk(part,columns,total):
        if not columns:
            result.loc[part.index]=total/len(part)
            return
        groups=list(part.groupby(columns[0],dropna=False,sort=True))
        for _,child in groups:walk(child,columns[1:],total/len(groups))
    walk(frame,['time_a','sample_unit_a','time_b','sample_unit_b'],1.)
    assert abs(result.sum()-1)<2e-12
    return result.to_numpy()

def read(p,cols=None):return pd.read_csv(p,sep='\t',float_precision='round_trip',usecols=cols,low_memory=False)
def main():
    protocol=json.loads((HERE.parent/'protocol.json').read_text())
    ids=read(HERE/'independent_background728.tsv')
    genes=set(ids.human_gene_id)|set(ids.mouse_gene_id)
    weight_file=MARKER/'unique_support_pair_weights.tsv.gz'
    weights=read(weight_file)
    weights=weights[weights.version.eq(protocol['marker_version']) & weights.setting.eq(protocol['setting']) & weights.background.eq(protocol['bridge_background']) & weights.gene_id.isin(genes)].copy()
    assert not weights.duplicated(PAIRKEY).any()
    fields=PAIRKEY+['PAS_ids','position','target_usage','reference_usage','target_group','reference_group']
    pair_file=MARKER/'pair_effects.tsv.gz'
    pairs=read(pair_file,fields)
    pairs=pairs[pairs.version.eq(protocol['marker_version']) & pairs.setting.eq(protocol['setting']) & pairs.gene_id.isin(genes)]
    assert not pairs.duplicated(PAIRKEY).any()
    data=weights.merge(pairs,on=PAIRKEY,how='left',validate='one_to_one')
    assert data.PAS_ids.notna().all() and data.position.notna().all()
    rows=[];max_weight_error=0.;max_delta_error=0.
    for key,g in data.groupby(ROUTEKEY,sort=True):
        g=g.sort_values('library_key')
        w=hierarchy_weights(g);max_weight_error=max(max_weight_error,float(abs(w-g.weight).max()))
        assert np.allclose(w,g.weight,atol=2e-12,rtol=0)
        pas=g.PAS_ids.map(json.loads).tolist();pos=g.position.map(json.loads).tolist()
        assert all(v==pas[0] for v in pas) and all(v==pos[0] for v in pos)
        position=np.asarray(pos[0],float)
        target=np.asarray(g.target_usage.map(json.loads).tolist(),float)
        ref=np.asarray(g.reference_usage.map(json.loads).tolist(),float)
        assert np.allclose(target.sum(1),1,atol=2e-12,rtol=0) and np.allclose(ref.sum(1),1,atol=2e-12,rtol=0)
        pair_delta=(target-ref)@position
        max_delta_error=max(max_delta_error,float(abs(pair_delta-g.delta_DWUI).max()))
        assert np.allclose(pair_delta,g.delta_DWUI,atol=2e-12,rtol=0)
        row=dict(zip(ROUTEKEY,key));row.update(n_libraries=len(g),
            contributor_libraries=json.dumps(g.library_key.tolist(),separators=(',',':')),
            weights=json.dumps(w.tolist(),separators=(',',':')),
            PAS_ids=json.dumps(pas[0],separators=(',',':')),
            position=json.dumps(pos[0],separators=(',',':')),
            target_usage=json.dumps((w@target).tolist(),separators=(',',':')),
            reference_usage=json.dumps((w@ref).tolist(),separators=(',',':')),
            delta_DWUI=float(w@pair_delta),baseline=float((w@ref)@position),
            retained_position_min=float(position.min()),retained_position_max=float(position.max()),
            PAS_number=len(position))
        rows.append(row)
    out=pd.DataFrame(rows)
    before=read(MARKER/'unique_support_bridge_route_effects.tsv.gz')
    before=before[before.version.eq(protocol['marker_version']) & before.setting.eq(protocol['setting']) & before.background.eq(protocol['bridge_background']) & before.gene_id.isin(genes)]
    a=out.set_index(ROUTEKEY).sort_index();b=before.set_index(ROUTEKEY).sort_index()
    assert a.index.equals(b.index)
    assert np.allclose(a.delta_DWUI,b.delta_DWUI,atol=2e-12,rtol=0)
    for column in ['contributor_libraries','PAS_ids','target_usage','reference_usage','weights']:
        for x,y in zip(a[column],b[column]):
            xx,yy=json.loads(x),json.loads(y)
            if column in ['contributor_libraries','PAS_ids']:assert xx==yy
            else:assert np.allclose(xx,yy,atol=2e-12,rtol=0)
    out.to_csv(HERE/'independent_marker_routes.tsv.gz',sep='\t',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
    data.to_csv(HERE/'independent_marker_pairs.tsv.gz',sep='\t',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
    inputs=[HERE.parent/'PROTOCOL.zh.md',HERE.parent/'protocol.json',HERE/'independent_background728.tsv',weight_file,pair_file,MARKER/'unique_support_bridge_route_effects.tsv.gz']
    qa={'status':'PASS','n_pair_contributions':len(data),'n_route_effects':len(out),
        'max_weight_abs_error':max_weight_error,'max_pair_delta_abs_error':max_delta_error,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
        'actual_gene_GEX_audited':False,'new_direction_or_null_results_audited':False}
    (HERE/'route_reconstruction_QA.json').write_text(json.dumps(qa,indent=2)+'\n')
    print(json.dumps(qa,indent=2))

if __name__=='__main__':main()
