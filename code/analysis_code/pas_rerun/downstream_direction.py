"""V2-only distal/proximal contrasts, separate from sharing enrichment."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from downstream_runtime import Context, bh

SOURCE='__direction_math__'
HUMAN=['GSE118258','GSE178325']
MOUSE=['GSE103221','GSE115943']


def correction_families(synchronized,pair):
    a=synchronized.copy();b=pair.copy()
    assert set(a.version)==set(b.version)=={'v2'}
    assert len(a)==48 and len(b)==32,'All planned v2 configurations must be present'
    absolute=a.test.isin(['distal_minus_proximal_excess_rate','centered_direction_variation_between_pairs'])
    assert absolute.sum()==24
    a.loc[absolute,'q_absolute_24']=bh(a.loc[absolute,'p_primary'])
    a.loc[~absolute,'q_relative_24']=bh(a.loc[~absolute,'p_primary'])
    a['q_synchronized_48']=bh(a.p_primary)
    b['q_pair_sensitivity_32']=bh(b.p_primary)
    q=bh(np.r_[a.p_primary,b.p_primary]);a['q_complete_80']=q[:48];b['q_complete_80']=q[48:]
    a['null_model']='study_vector_synchronized';b['null_model']='original_pair_specific'
    return a,b


def family_registry(a,b):
    rows=[]
    for origin,frame in [('synchronized',a),('pair_specific',b)]:
        for i,r in frame.iterrows():
            key='|'.join(str(r[c]) for c in ['version','metric','threshold','pair','scope','test','null_model'])
            for col in [c for c in frame if c.startswith('q_')]:
                if pd.notna(r[col]):rows.append(dict(test_id=key,source=origin,source_row=i+1,
                    correction_family=col,p=float(r.p_primary),q=float(r[col])))
    result=pd.DataFrame(rows)
    assert result.test_id.nunique()==80 and not result.duplicated(['test_id','correction_family']).any()
    return result


def run(ctx):
    assert ctx.config['reference_variant']=='v2_primary','v3 is excluded by the current manuscript scope'
    folder=ctx.analysis/'cross_sharing'
    for name in ['sharing_strata.tsv.gz','sharing_backgrounds.tsv.gz','sharing_summary.tsv','null_draws.npz','result.json']:
        ctx.register(folder/name)
    names=['sha','read','prepare_vectors','vector_shuffle','p_record','log_relative',
           'analyze_configuration','original_pair_sensitivity']
    num=ctx.functions(SOURCE,names,dict(HUMAN=HUMAN,MOUSE=MOUSE,B=10000,SOURCES={}))
    strata=num.read(folder,'sharing_strata.tsv.gz');background=num.read(folder,'sharing_backgrounds.tsv.gz')
    summary=num.read(folder,'sharing_summary.tsv')
    original=num.original_pair_sensitivity(folder,summary,'v2')
    rows=[];components=[];mobility=[];cells=[];archives={}
    for metric in ['WPP','RWUI']:
        for threshold in [.02,.05]:
            r,c,m,s,bank=num.analyze_configuration(strata,background,summary,'v2',metric,threshold)
            rows+=r;components+=c;mobility+=m;cells+=s
            archives.update({f'v2__{metric}__{threshold}__{key}':value for key,value in bank.items()})
    a,b=correction_families(pd.DataFrame(rows),pd.DataFrame(original))
    ctx.table(a,'direction_contrast_tests.tsv');ctx.table(b,'original_pair_direction_sensitivity.tsv')
    ctx.table(components,'direction_components.tsv');ctx.table(mobility,'joint_null_mobility.tsv')
    ctx.table(cells,'joint_permutation_strata.tsv.gz');ctx.table(family_registry(a,b),'correction-family-membership.tsv')
    ctx.table(summary[summary.metric.isin(['WPP','RWUI'])],'v2_original_direction_decomposition.tsv')
    np.savez_compressed(ctx.out/'synchronized_null_draws.npz',**archives)
    ctx.finish(dict(v2_only=True,complete_48_plus_32_tests=True,study_vectors_permuted_intact=True,
        eligibility_and_stratum_margins_fixed=True,ortholog_equal_weighted=True,separate_from_q24_q56_enrichment=True),
        correction_families=dict(absolute=24,relative=24,synchronized=48,pair_specific=32,complete=80),
        permutations=10000,scope='Fixed sampled studies and current route-specific jointly measurable orthologs')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'direction_difference',a.out))
