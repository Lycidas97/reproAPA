"""Join newly qualified APA contributors to unchanged all-gene expression."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from downstream_parallel import ordered_map
from downstream_runtime import Context, Study, STUDIES, SETTINGS, KEY
from polyapipe_effects import weights


def reconstruct_setting(ctx,item):
    sid,setting=item
    all_pairs,all_routes,all_axes,all_states,checks=[],[],[],[],[]
    study=Study(ctx,sid)
    study.load_gex()
    pairs=study.core('pair_effects',setting)
    routes=study.core('route_effects',setting)
    if pairs.empty:
        checks.append(dict(study_id=sid,setting=setting,pairs=0,status='no_qualified_endpoint_pairs'))
        return (all_pairs,all_routes,all_axes,all_states,checks)
    assert not pairs.duplicated(KEY+['library_key_a','library_key_b']).any()
    added=[];states=[]
    for gene,part in pairs.groupby('gene_id',sort=False):
        axis=study.axis(gene,setting);data=study.gene_counts(axis)
        axis['setting']=setting;all_axes.append(axis)
        cache={}
        for rec in part.to_dict('records'):
            akey=(rec['route'],rec['library_key_a'],'P')
            bkey=(rec['route'],rec['library_key_b'],rec['contrast'][-1])
            for role in [akey,bkey]:
                if role in cache:continue
                s=study.state(data,*role,setting)
                assert s['endpoint_eligible']
                s.update(study.expression_summary(gene,s['rows']))
                all_rows=np.flatnonzero(study.meta.route.eq(role[0])&study.meta.library_key.eq(role[1])&
                    study.meta.node.eq(role[2])&study.meta.purity.ge(SETTINGS[setting][1]))
                s['GEX_all']=study.expression_summary(gene,all_rows)['GEX']
                cache[role]=s
                states.append(dict(species=rec['species'],study_id=sid,setting=setting,route=role[0],
                    library_key=role[1],state=role[2],gene_id=gene,gene_symbol=rec['gene_symbol'],
                    n_mc=s['n_mc'],n_cells=s['n_cells'],pas_counts=int(s['site_counts'].sum()),
                    GEX_CPM=s['GEX_CPM'],GEX=s['GEX'],GEX_all=s['GEX_all'],GEX_count=s['GEX_count'],
                    APA=float(s['usage']@axis.position),
                    site_counts=json.dumps(s['site_counts'].tolist()),site_mc=json.dumps(s['site_mc'].tolist()),
                    usage=json.dumps(s['usage'].tolist()),PAS_ids=json.dumps(axis.atlas_pas_id.tolist()),
                    metacell_ids=json.dumps(study.meta.metacell_id.iloc[s['rows']].tolist())))
            a,b=cache[akey],cache[bkey]
            assert np.allclose(a['usage'],json.loads(rec['target_usage']),rtol=0,atol=1e-12)
            assert np.allclose(b['usage'],json.loads(rec['reference_usage']),rtol=0,atol=1e-12)
            assert np.array_equal(a['site_counts'],json.loads(rec['molecules_a']))
            assert np.array_equal(b['site_counts'],json.loads(rec['molecules_b']))
            assert rec['PAS_ids']==json.dumps(axis.atlas_pas_id.tolist())
            rec.update(mean_GEX=(a['GEX_CPM']+b['GEX_CPM'])/2,
                GEX_CPM_a=a['GEX_CPM'],GEX_CPM_b=b['GEX_CPM'],
                GEX_log2FC=a['GEX']-b['GEX'],GEX_delta=a['GEX']-b['GEX'],
                GEX_all_delta=a['GEX_all']-b['GEX_all'],GEX_count_a=a['GEX_count'],GEX_count_b=b['GEX_count'],
                GEX_measurable=bool(a['GEX_count']+b['GEX_count']>=20 and max(a['GEX_CPM'],b['GEX_CPM'])>=1),
                APA_coverage=(a['site_counts'].sum()/a['n_mc']+b['site_counts'].sum()/b['n_mc'])/2,
                n_joint_PAS=int(((a['site_counts']+b['site_counts']>=10)&(a['site_mc']+b['site_mc']>=2)).sum()),
                delta_WPP=rec['delta_DWUI'],APA_delta=rec['delta_DWUI'],
                site_mc_a=json.dumps(a['site_mc'].tolist()),site_mc_b=json.dumps(b['site_mc'].tolist()))
            assert rec['n_joint_PAS']>=2
            added.append(rec)
    pairs=pd.DataFrame(added)
    joined=[]
    route_index={tuple('' if pd.isna(r[k]) else r[k] for k in KEY):r
                 for r in routes.to_dict('records')}
    assert len(route_index)==len(routes)
    for keys,part in pairs.groupby(KEY,dropna=False,sort=False):
        r=route_index[tuple('' if pd.isna(value) else value for value in keys)].copy();w=weights(part)
        assert np.allclose(w,part.weight,rtol=0,atol=1e-12)
        a=w@np.stack(part.target_usage.map(json.loads));b=w@np.stack(part.reference_usage.map(json.loads))
        assert np.allclose(a,json.loads(r['target_usage']),rtol=0,atol=1e-12)
        assert np.allclose(b,json.loads(r['reference_usage']),rtol=0,atol=1e-12)
        assert abs(np.ptp(a-b)-r['MPRO'])<1e-12
        r.update(mean_GEX=float(w@part.mean_GEX),APA_coverage=float(w@part.APA_coverage),
                 GEX_log2FC=float(w@part.GEX_log2FC),GEX_delta=float(w@part.GEX_delta),
                 GEX_all_delta=float(w@part.GEX_all_delta),GEX_measurable=bool(part.GEX_measurable.all()),delta_WPP=r['delta_DWUI'],
                 reference_WPP=r['reference_DWUI'],target_WPP=r['target_DWUI'],
                 n_target_libraries=part.library_key_a.nunique(),n_source_libraries=part.library_key_b.nunique())
        joined.append(r)
    all_pairs.append(pairs);all_routes.append(pd.DataFrame(joined));all_states.append(pd.DataFrame(states))
    checks.append(dict(study_id=sid,setting=setting,pairs=len(pairs),routes=len(joined),status='matrix_and_vector_recovered'))
    print(sid,setting,len(pairs),len(joined),flush=True)
    return (all_pairs,all_routes,all_axes,all_states,checks)


def run(ctx):
    all_pairs,all_routes,all_axes,all_states,checks=[],[],[],[],[]
    for result in ordered_map(ctx,reconstruct_setting,[(sid,setting) for sid in STUDIES for setting in SETTINGS]):
        for target,values in zip([all_pairs,all_routes,all_axes,all_states,checks],result):target.extend(values)
    assert all_pairs, 'No measurable endpoints across the seven-study release'
    for name,data in [('pair_effects',all_pairs),('route_effects',all_routes),('PAS_axes',all_axes),('library_states',all_states)]:
        frame=pd.concat(data,ignore_index=True)
        ctx.table(frame,name+'.tsv.gz')
        ctx.table(frame[frame.setting.eq('supported50_p08')],name+'.primary.tsv.gz')
    ctx.table(checks,'reconstruction.tsv')
    ctx.finish(dict(integer_states_recovered=True,complete_usage_vectors_recovered=True,
        fresh_gene_qualification=True,contributor_weights_recovered=True,GEX_requeried_on_new_contributors=True),
        metric_alias='WPP is a compatibility alias for DWUI, not a separate measure')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'state_vectors',a.out))
