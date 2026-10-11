"""Independent group denominators and pair/route/study assay covariates."""
from pathlib import Path
import hashlib,json

import numpy as np
import pandas as pd

from reconstruct_marker_routes import hierarchy_weights

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2]
ANALYSIS=HERE.parent/'analysis';MARKER=ROOT/'notes/2026-10-10-marker-reference-sensitivity'
PAIRKEY=['study_id','route','library_key','gene_id','version','setting','contrast']
ROUTEKEY=['study_id','route','gene_id','contrast']

def read(p):return pd.read_csv(p,sep='\t',float_precision='round_trip',low_memory=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(f,n):f.to_csv(HERE/n,sep='\t',index=False,float_format='%.17g',compression={'method':'gzip','mtime':0})
def close(a,b,tol=2e-10):assert np.allclose(a,b,atol=tol,rtol=0,equal_nan=True)
def compare(a,b,key,columns):
    aa=a.set_index(key).sort_index();bb=b.set_index(key).sort_index()
    assert aa.index.equals(bb.index) and aa.index.is_unique
    for c in columns:
        if pd.api.types.is_numeric_dtype(aa[c]):close(aa[c],bb[c])
        else:assert aa[c].fillna('').equals(bb[c].fillna('')),c
    return len(aa)

def main():
    groups=read(ANALYSIS/'qualified_groups.tsv.gz')
    gex=read(ANALYSIS/'group_GEX_aggregates.tsv.gz')
    requests=read(ANALYSIS/'gex_requests.tsv.gz')
    assert groups.group_id.is_unique and gex.group_id.is_unique and requests.group_id.is_unique
    assert set(groups.group_id)==set(gex.group_id)==set(requests.group_id)
    m=groups.merge(gex,on=['species','study_id','group_id','gene_id','n_mc','n_cells','qualified_MC_ids'],validate='one_to_one')
    assert len(m)==len(groups)==45705
    integer_columns=['gene_count_sum','all_gene_GEX_sum','GEX_gene_column']
    for c in integer_columns:assert pd.api.types.is_integer_dtype(m[c]) and m[c].ge(0).all()
    assert m.all_gene_GEX_sum.gt(0).all() and m.gene_count_sum.le(m.all_gene_GEX_sum).all()
    close(m.GEX_CPM,m.gene_count_sum/m.all_gene_GEX_sum*1e6)
    cfg=json.loads((ANALYSIS/'server_source_config.json').read_text())
    manifest=json.loads((ANALYSIS/'GEX_source_manifest.json').read_text())
    assert manifest['request_sha256']==sha(ANALYSIS/'gex_requests.tsv.gz')
    assert manifest['output_sha256']==sha(ANALYSIS/'group_GEX_aggregates.tsv.gz')
    for p,x in manifest['sources'].items():assert cfg['expected_hashes'][p]==x['sha256']
    # These frozen total-GEX vectors and MC identities predate the new extraction.
    # They independently verify all denominators and cells without downloading
    # full matrices or trusting the fresh server's metadata rows alone.
    old=np.load(MARKER/'inputs.npz',allow_pickle=False)
    ap=[]
    for sid,g in m.groupby('study_id',sort=True):
        meta=pd.DataFrame(json.loads(str(old[sid+'_meta_json'])))
        total=old[sid+'_total_gex'];assert len(total)==len(meta) and meta.metacell_id.is_unique
        assert np.issubdtype(total.dtype,np.integer) and (total>=0).all()
        lookup={v:(int(total[i]),int(meta.n_cells.iloc[i])) for i,v in enumerate(meta.metacell_id)}
        for r in g.itertuples():
            ids=json.loads(r.qualified_MC_ids);depths=json.loads(r.qualified_MC_pas_depths)
            assert len(ids)==len(set(ids))==r.n_mc and len(depths)==r.n_mc and min(depths)>=20
            assert all(v in lookup for v in ids)
            den=sum(lookup[v][0] for v in ids);cells=sum(lookup[v][1] for v in ids)
            assert den==r.all_gene_GEX_sum and cells==r.n_cells
            molecules=json.loads(r.pooled_molecules)
            assert sum(molecules)==sum(depths)
            ap.append((r.group_id,float(sum(molecules)/r.n_mc)))
    m['independent_APA_coverage']=m.group_id.map(dict(ap))
    assert m.independent_APA_coverage.notna().all()
    if (HERE/'source_GEX_spotcheck.json').exists():
        spot=json.loads((HERE/'source_GEX_spotcheck.json').read_text());assert spot['status']=='PASS' and len(spot['results'])==15
        fresh=pd.DataFrame(spot['results']).set_index('group_id');cached=gex.set_index('group_id').loc[fresh.index]
        for c in ['gene_id','study_id','n_mc','n_cells','gene_count_sum','all_gene_GEX_sum']:
            assert fresh[c].equals(cached[c]),c
        close(fresh.GEX_CPM,cached.GEX_CPM)
        assert np.array_equal(fresh.gene_column,cached.GEX_gene_column)
        for p,h in spot['source_sha256'].items():assert manifest['sources'][p]['sha256']==h
        spot_pass=True
    else:spot_pass=False
    pairs=read(HERE/'independent_marker_pairs.tsv.gz')
    lookup=m.set_index('group_id')
    for c,s in [('mean_GEX','GEX_CPM'),('APA_coverage','independent_APA_coverage')]:
        av=lookup[s].reindex(pairs.target_group).to_numpy()
        bv=lookup[s].reindex(pairs.reference_group).to_numpy()
        assert np.isfinite(av).all() and np.isfinite(bv).all()
        pairs[c]=(av+bv)/2
    pairs['baseline']=[np.dot(json.loads(r.reference_usage),json.loads(r.position)) for r in pairs.itertuples()]
    pairs['PAS_number']=pairs.PAS_ids.map(lambda v:len(json.loads(v)))
    compare(pairs,read(ANALYSIS/'pair_covariates.tsv.gz'),PAIRKEY,['mean_GEX','APA_coverage','baseline','PAS_number','delta_DWUI','weight'])
    base=read(HERE/'independent_marker_routes.tsv.gz').set_index(ROUTEKEY)
    rows=[]
    for k,g in pairs.groupby(ROUTEKEY,sort=True):
        g=g.sort_values('library_key');w=hierarchy_weights(g)
        r=dict(zip(ROUTEKEY,k));r['species']='mouse' if k[0] in ['GSE103221','GSE115943'] else 'human'
        for c in ['mean_GEX','APA_coverage','baseline']:r[c]=float(w@g[c])
        for c in ['delta_DWUI','PAS_number','retained_position_min','retained_position_max']:r[c]=base.loc[k,c]
        rows.append(r)
    route=pd.DataFrame(rows)
    compare(route,read(ANALYSIS/'route_covariates.tsv.gz'),ROUTEKEY,['species','mean_GEX','APA_coverage','baseline','PAS_number','delta_DWUI','retained_position_min','retained_position_max'])
    studies=[]
    for key,g in route.groupby(['species','contrast','study_id','gene_id'],sort=True):
        assert g.route.is_unique
        row=dict(zip(['species','source_contrast','unit','gene_id'],key))
        row.update(contrast='P_M',metric='WPP',n_routes=len(g),series_pattern=';'.join(sorted(g.route)),
                   effect_mean=float(g.delta_DWUI.mean()),effect_min=float(g.delta_DWUI.min()),effect_max=float(g.delta_DWUI.max()),
                   baseline=float(g.baseline.median()),baseline_min=float(g.baseline.min()),baseline_max=float(g.baseline.max()),
                   capacity_negative=float(g.baseline.min()),capacity_positive=float(1-g.baseline.max()))
        for c in ['mean_GEX','APA_coverage','PAS_number']:row[c]=float(g[c].median())
        row['supported_catalogue_PAS_number']=row['PAS_number'];studies.append(row)
    study=pd.DataFrame(studies)
    compare(study,read(ANALYSIS/'study_covariates.tsv.gz'),['species','source_contrast','unit','gene_id'],
            [c for c in study if c not in ['species','source_contrast','unit','gene_id']])
    save(study,'independent_study_covariates.tsv.gz');save(route,'independent_route_covariates.tsv.gz')
    inputs=[HERE.parent/'PROTOCOL.zh.md',HERE.parent/'protocol.json',ANALYSIS/'qualified_groups.tsv.gz',ANALYSIS/'group_GEX_aggregates.tsv.gz',
            ANALYSIS/'GEX_source_manifest.json',MARKER/'inputs.npz',HERE/'independent_marker_pairs.tsv.gz',HERE/'independent_marker_routes.tsv.gz']
    qa={'status':'PASS' if spot_pass else 'LOCAL_COVARIATES_PASS_SOURCE_SPOTCHECK_PENDING',
        'all_group_denominators_and_cells_checked':len(m),'pair_covariates_checked':len(pairs),'route_covariates_checked':len(route),
        'study_covariates_checked':len(study),'source_matrix_fixed15_scalar_replay_pass':spot_pass,
        'actual_CPM_not_old_proxy':True,'new_null_or_BH_audited':False,
        'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs}}
    (HERE/'covariate_QA.json').write_text(json.dumps(qa,indent=2)+'\n');print(json.dumps(qa,indent=2))

if __name__=='__main__':main()
