"""Compare independently recomputed v2 and v3 results without mixing their inputs."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from downstream_runtime import Context, KEY
from downstream_cohorts import CONTRASTS
from polyapipe_prepare import digest


def effect_comparison(primary, control):
    keys=['setting']+KEY
    metrics=['delta_DWUI','delta_RWUI','delta_PDUI','delta_PPUI','MPRO']
    assert not primary.duplicated(keys).any() and not control.duplicated(keys).any()
    joined=primary[keys+metrics].merge(control[keys+metrics],on=keys,how='outer',
        suffixes=('_v2','_v3'),indicator='availability',validate='one_to_one')
    expanded=joined[joined.setting.eq('expanded_primary_p08')]
    assert expanded.availability.eq('both').all(), 'Reference-independent eligibility changed'
    for metric in metrics:
        np.testing.assert_allclose(expanded[metric+'_v2'],expanded[metric+'_v3'],rtol=0,atol=1e-12,
                                   err_msg='Reference-independent effects changed')
    rows=[]
    groups=['species','study_id','route','contrast','design','setting']
    for key, frame in joined.groupby(groups,dropna=False,sort=True):
        info=dict(zip(groups,key));both=frame[frame.availability.eq('both')]
        for metric in metrics:
            a=both[metric+'_v2'].to_numpy(float);b=both[metric+'_v3'].to_numpy(float)
            assert np.isfinite(a).all() and np.isfinite(b).all()
            nonzero=(a!=0)&(b!=0);varying=len(a)>=3 and np.ptp(a)>0 and np.ptp(b)>0
            for threshold in ([.1,.2] if metric=='MPRO' else [.02,.05]):
                rows.append(dict(info,metric=metric,threshold=threshold,measured_v2=int(frame.availability.ne('right_only').sum()),
                    measured_v3=int(frame.availability.ne('left_only').sum()),common=len(both),
                    unmeasured_in_v2=int(frame.availability.eq('right_only').sum()),unmeasured_in_v3=int(frame.availability.eq('left_only').sum()),
                    pearson=float(np.corrcoef(a,b)[0,1]) if varying else np.nan,
                    spearman=float(spearmanr(a,b).statistic) if varying else np.nan,
                    median_absolute_difference=float(np.median(abs(a-b))) if len(a) else np.nan,
                    direction_denominator=int(nonzero.sum()) if metric!='MPRO' else 0,
                    direction_agreement=float((np.sign(a[nonzero])==np.sign(b[nonzero])).mean()) if metric!='MPRO' and nonzero.any() else np.nan,
                    threshold_agreement=float(((abs(a)>=threshold)==(abs(b)>=threshold)).mean()) if len(a) else np.nan))
    return joined,pd.DataFrame(rows)


def member_comparison(primary, control, keys, groups, planned=None):
    a=primary[keys].drop_duplicates().assign(selected_v2=True)
    b=control[keys].drop_duplicates().assign(selected_v3=True)
    union=a.merge(b,on=keys,how='outer',validate='one_to_one').fillna({'selected_v2':False,'selected_v3':False})
    rows=[]
    for key, frame in union.groupby(groups,dropna=False,sort=True):
        key=key if isinstance(key,tuple) else (key,)
        n2=int(frame.selected_v2.sum());n3=int(frame.selected_v3.sum());both=int((frame.selected_v2&frame.selected_v3).sum())
        rows.append(dict(zip(groups,key),v2_members=n2,v3_members=n3,shared_members=both,
            v2_only=n2-both,v3_only=n3-both,union=len(frame),jaccard=both/len(frame),
            fraction_v2_retained=both/n2 if n2 else np.nan,fraction_v3_in_primary=both/n3 if n3 else np.nan))
    counts=pd.DataFrame(rows,columns=groups+['v2_members','v3_members','shared_members','v2_only','v3_only','union','jaccard','fraction_v2_retained','fraction_v3_in_primary'])
    if planned is not None:
        design=planned[groups].drop_duplicates()
        counts=design.merge(counts,on=groups,how='outer',validate='one_to_one',indicator=True)
        assert not counts._merge.eq('right_only').any(), 'Selected members outside the planned comparison'
        for col in ['v2_members','v3_members','shared_members','v2_only','v3_only','union']:
            counts[col]=counts[col].fillna(0).astype(int)
        counts=counts.drop(columns='_merge')
        # Empty/empty is a measured zero membership count, not perfect agreement.
    return union,counts


def run(args):
    ctx=Context(args.primary_config,'reference_comparison',args.out)
    control=json.loads(ctx.register(args.control_config).read_text())
    assert ctx.config['reference_variant']=='v2_primary' and control['reference_variant']=='v3_ge62_comparison'
    assert ctx.run_id!=control['run_id'] and ctx.config['input_fingerprint']!=control['input_fingerprint']
    configs={'v2':ctx.config,'v3':control}
    def read(variant,relative):
        config=configs[variant];root=Path(config['analysis_root']);node=Path(relative).parts[0]
        receipt=json.loads(ctx.register(root/node/'result.json').read_text())
        assert receipt['status']=='complete' and all(receipt['checks'].values())
        assert receipt['run_id']==config['run_id'] and receipt['input_fingerprint']==config['input_fingerprint']
        hashes={a['path']:a['sha256'] for a in receipt['artifacts']}
        path=ctx.register(root/relative)
        assert digest(path)==hashes[str(Path(relative).relative_to(node))]
        return ctx.read(path)
    native={v:Path(c['native_run']) for v,c in configs.items()};support=[]
    for name in sorted((native['v2']/'studies').iterdir()):
        if not name.is_dir():continue
        a=ctx.read(name/'features.tsv.gz',low_memory=False)
        b=ctx.read(native['v3']/'studies'/name.name/'features.tsv.gz',low_memory=False)
        same=['atlas_pas_id','feature_index_1based','chrom','strand','cut_site','primary_gene_id','keep_for_apa_primary']
        pd.testing.assert_frame_equal(a[same],b[same])
        matrix_a=ctx.register(name/'metacell_by_pas_counts.npz')
        matrix_b=ctx.register(native['v3']/'studies'/name.name/'metacell_by_pas_counts.npz')
        assert digest(matrix_a)==digest(matrix_b)
        for window in [25,50,100]:
            primary=a.keep_for_apa_primary.astype(str).str.lower().eq('true')
            for layer,mask in [('all',np.ones(len(a),bool)),('primary',primary)]:
                x=a[f'matched_source_count_{window}nt'].ge(1)&mask;y=b[f'matched_source_count_{window}nt'].ge(1)&mask
                support.append(dict(study_id=name.name,species=a.species.iloc[0],window_nt=window,layer=layer,
                    denominator=int(np.sum(mask)),supported_v2=int(x.sum()),supported_v3=int(y.sum()),
                    supported_both=int((x&y).sum()),v2_only=int((x&~y).sum()),v3_only=int((~x&y).sum())))
    ctx.table(support,'PAS_support_comparison.tsv')
    effects,summary=effect_comparison(read('v2','state_vectors/route_effects.tsv.gz'),read('v3','state_vectors/route_effects.tsv.gz'))
    ctx.table(effects,'endpoint_effects.tsv.gz');ctx.table(summary,'endpoint_agreement.tsv')
    specifications=[('trajectory','profiles/gene_membership_k3.tsv.gz',['species','gene_id'],['species']),
        ('recurrence','recurrence/recurrence_members.tsv.gz',['metric','species','contrast','threshold','series_matching','variant','gene_id','direction'],
         ['metric','species','contrast','threshold','series_matching','variant']),
        ('cross_shared','cross_overlap/new_shared_members.tsv.gz',['layer','threshold','contrast','human_gene_id'],['layer','threshold','contrast']),
        ('figure5_union','union_pool/figure5_union_claims.tsv',['species','contrast','gene_id'],['species','contrast'])]
    for name,relative,keys,groups in specifications:
        if name=='trajectory':planned=pd.DataFrame({'species':['human','mouse']})
        elif name=='recurrence':
            planned=pd.concat([read(v,'recurrence/recurrence_summary.tsv')[groups] for v in configs])
        elif name=='cross_shared':
            planned=pd.concat([read(v,'cross_overlap/overlap_threshold_layer_sensitivity.tsv')[groups] for v in configs])
        else:
            planned=pd.DataFrame([(s,c) for s,design in CONTRASTS.items() for c in design],columns=groups)
        members,counts=member_comparison(read('v2',relative),read('v3',relative),keys,groups,planned)
        ctx.table(members,name+'_members.tsv.gz');ctx.table(counts,name+'_summary.tsv')
    ctx.table([dict(reference=v,run_id=c['run_id'],input_fingerprint=c['input_fingerprint'],reference_variant=c['reference_variant'])
        for v,c in configs.items()],'reference_runs.tsv')
    ctx.finish(dict(two_independent_runs_verified=True,identical_native_axes_and_counts=True,
        expanded_primary_effects_identical=True,unmeasured_effects_not_zero=True,all_selected_members_retained=True),
        comparison_scope='Annotation-reference sensitivity; not PAS caller accuracy; both run fingerprints recorded',
        primary_run=ctx.run_id,control_run=control['run_id'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--primary-config',required=True);p.add_argument('--control-config',required=True);p.add_argument('--out',required=True)
    run(p.parse_args())
