"""Restore frozen P/M/S PAS vectors and aggregate on identical contributors."""
from pathlib import Path
import importlib.util
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OLD = ROOT / 'notes/2026-09-27-source-state-connection'
MET = ROOT / 'notes/2026-09-30-multimetric-statistics'
KEY = ['species', 'study_id', 'route', 'setting', 'gene_id']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(frame, name):
    frame.to_csv(HERE/name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
                 compression={'method':'gzip', 'mtime':0} if name.endswith('.gz') else None)


def main():
    ext = module('frozen_metric_extract', MET/'extract.py')
    original = module('frozen_source_analysis', OLD/'analyze.py')
    protected = [p for folder in ['main','supplement','Fig'] for p in (ROOT/folder).rglob('*') if p.is_file()]
    protected += [ROOT/'main.tex', ROOT/'supplement.tex']
    before = {str(p.relative_to(ROOT)):sha(p) for p in protected}
    (HERE/'protected-sha256.json').write_text(json.dumps(before,indent=2)+'\n')
    t = pd.read_csv(OLD/'matched_library_source_triplets.tsv.gz',sep='\t')
    t = t[t.setting.isin(['supported50_p08','supported50_p08_mc5'])].copy()
    profiles = {}
    for path in ext.PROFILES.glob('*/selected_metacells.tsv.gz'):
        sid = pd.read_csv(path,sep='\t',nrows=1).study_id.iloc[0]
        profiles[sid] = path.parent.name
    all_states, axes, audits = [], [], []
    for study, part in t[t.setting.eq('supported50_p08')].groupby('study_id'):
        cache = HERE/f'{study}_library_vectors.tsv.gz'
        axis_path = HERE/f'{study}_PAS_axis.tsv.gz'
        audit_path = HERE/f'{study}_extraction_audit.json'
        if cache.exists() and axis_path.exists() and audit_path.exists():
            audit = json.loads(audit_path.read_text())
            assert audit['cache_sha256'] == sha(cache) and audit['axis_sha256'] == sha(axis_path)
            assert audit['helper_sha256'] == sha(MET/'extract.py')
            states = pd.read_csv(cache,sep='\t'); axis = pd.read_csv(axis_path,sep='\t')
        else:
            cols = ['species','study_id','route','gene_id','gene_symbol','library_key_a','library_key_b']
            ps = part[cols].drop_duplicates().assign(contrast='P_S',design='true_source')
            pm = ps.drop(columns='library_key_b').drop_duplicates().assign(contrast='P_M',design='same_library')
            pm['library_key_b'] = pm.library_key_a
            pairs = pd.concat([ps,pm],ignore_index=True)
            vectors, states, features, groups, axis, error = ext.study_vectors(profiles[study],pairs)
            ids, weights, usage = [], [], []
            for row in states.itertuples():
                ix = groups[row.gene_id]
                ids.append(json.dumps(features.atlas_pas_id.iloc[ix].tolist()))
                weights.append(json.dumps(features.WPP.iloc[ix].tolist()))
                usage.append(json.dumps(vectors[(row.route,row.library_key,row.state,row.gene_id)][0].tolist()))
            states['PAS_ids'],states['position_weights'],states['usage'] = ids,weights,usage
            save(states,cache.name);save(axis,axis_path.name)
            audit = {'study_id':study,'library_state_gene_records':len(states),'max_DWUI_recovery_error':error,
                     'helper_sha256':sha(MET/'extract.py'),'triplet_input_sha256':sha(OLD/'matched_library_source_triplets.tsv.gz'),
                     'cache_sha256':sha(cache),'axis_sha256':sha(axis_path),'inputs':ext.INPUTS.copy()}
            audit_path.write_text(json.dumps(audit,indent=2)+'\n')
        all_states.append(states);axes.append(axis);audits.append(audit)
        print('VECTORS_READY',study,len(states),flush=True)
    states = pd.concat(all_states,ignore_index=True)
    lookup = states.set_index(['study_id','route','library_key','state','gene_id'])
    assert lookup.index.is_unique
    store = {}
    for key,row in lookup.iterrows():
        store[key]=(np.array(json.loads(row.usage)),np.array(json.loads(row.position_weights)),row.PAS_ids)
    records, weighted = [], []
    for setting, part in t.groupby('setting'):
        weighted_part = ext.pair_weights(part.assign(contrast='TRI',design='common_P_M_S'))
        weighted.append(weighted_part)
        for keys, g in weighted_part.groupby(KEY,sort=False):
            rec = dict(zip(KEY,keys)); rec['gene_symbol']=g.gene_symbol.iloc[0]
            pms = {}
            for state,side in [('P','a'),('M','a'),('S','b')]:
                vals=[]
                for row in g.itertuples():
                    vec,w,pids=store[(row.study_id,row.route,getattr(row,'library_key_'+side),state,row.gene_id)]
                    vals.append(vec)
                    assert pids == store[(row.study_id,row.route,row.library_key_a,'P',row.gene_id)][2]
                pms[state]=np.stack(vals)
                averaged=g.weight.to_numpy()@pms[state]
                assert abs(averaged.sum()-1)<1e-12
                rec[state+'_usage']=json.dumps(averaged.tolist())
                rec[state+'_DWUI']=float(averaged@w)
                rec[state+'_min_mc']=int(g[state+'_n_mc'].min())
                assert abs(rec[state+'_DWUI']-g.weight.to_numpy()@g[state+'_WPP'])<1e-8
            rec['PAS_ids'],rec['position_weights'],rec['PAS_number']=pids,json.dumps(w.tolist()),len(w)
            rec.update(n_triplets=len(g),n_target_libraries=g.library_key_a.nunique(),
                       n_source_libraries=g.library_key_b.nunique(),n_target_dates=g.time_a.nunique())
            for contrast,a,b in [('PS','P','S'),('PM','P','M'),('MS','M','S')]:
                differences=pms[a]-pms[b]
                delta=g.weight.to_numpy()@differences
                rec['delta_usage_'+contrast]=json.dumps(delta.tolist())
                rec['DWUI_'+contrast]=float(delta@w)
                rec['MPRO_'+contrast]=float(np.ptp(delta))
                rec['mean_triplet_MPRO_'+contrast]=float(g.weight.to_numpy()@np.ptp(differences,axis=1))
                assert rec['MPRO_'+contrast]<=rec['mean_triplet_MPRO_'+contrast]+1e-12
            assert abs(rec['DWUI_PS']-rec['DWUI_PM']-rec['DWUI_MS'])<1e-12
            assert np.allclose(np.array(json.loads(rec['delta_usage_PS'])),
                               np.array(json.loads(rec['delta_usage_PM']))+np.array(json.loads(rec['delta_usage_MS'])),atol=1e-12)
            records.append(rec)
    result=pd.DataFrame(records)
    old=pd.read_csv(OLD/'matched_state_summaries.tsv',sep='\t')
    old=old[old.setting.isin(result.setting.unique())]
    cmp=result.merge(old,on=KEY,how='outer',validate='one_to_one',indicator=True,suffixes=('','_old'))
    assert cmp._merge.eq('both').all()
    errors={co:float((cmp['DWUI_'+co]-cmp['delta_'+co]).abs().max()) for co in ['PS','PM','MS']}
    assert max(errors.values())<1e-8,errors
    # Independent scalar aggregation checks the hierarchical weights, not just linearity.
    check=original.aggregate(t,['P_WPP','M_WPP','S_WPP']).merge(result,on=KEY,validate='one_to_one')
    for state in ['P','M','S']:
        assert np.max(abs(check[state+'_WPP']-check[state+'_DWUI']))<1e-8
    save(result,'matched_state_vectors.tsv.gz')
    save(pd.concat(weighted,ignore_index=True),'weighted_triplets.tsv.gz')
    after={str(p.relative_to(ROOT)):sha(p) for p in protected}
    assert before==after
    report={'status':'PASS','studies':audits,'matched_records':len(result),'triplets':len(t),
            'DWUI_recovery_errors':errors,'all_probability_and_vector_identities_passed':True,
            'nonlinear_MPRO_aggregation_checked':True,'formal_files_unchanged':True,
            'script_sha256':sha(Path(__file__)),'protocol_sha256':sha(HERE/'PROTOCOL.zh.md')}
    (HERE/'extraction_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('COMPLETE',len(result),len(t),errors,flush=True)


if __name__=='__main__':
    main()
