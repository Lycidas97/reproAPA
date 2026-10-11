"""Fresh endpoint inventories, matched triplets and regulatory candidate pools."""
import argparse
import json
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from downstream_runtime import Context,Study,STUDIES,KEY,SETTINGS
from polyapipe_effects import weights

HUMAN={'GSE118258':['GSE118258'],'GSE150311':['GSE150311'],'GSE178325':['178325_0618'],'GSE221739':['GSE221739']}
BOTH={'GSE115943':['WOT'],'GSE103221':['CIP','OSK']}
CIP={'GSE115943':['WOT'],'GSE103221':['CIP']}
CONTRASTS={'human':{'SP':HUMAN,'SM':HUMAN,'MP':HUMAN},
    'mouse':{'SR':BOTH,'RE':BOTH,'EP':CIP,'SP':CIP,'MP':CIP,
             'SM':{**BOTH,'GSE171127':['171127_GETM','171127_OSKM']}}}


def endpoints(ctx):
    p=ctx.read(ctx.analysis/'state_vectors/pair_effects.tsv.gz')
    r=ctx.read(ctx.analysis/'state_vectors/route_effects.tsv.gz')
    assert p.loc[p.contrast.eq('P_M'),'library_key_a'].eq(p.loc[p.contrast.eq('P_M'),'library_key_b']).all()
    audit=[]
    for key,part in p.groupby(KEY+['setting'],dropna=False):
        expected=weights(part);assert np.allclose(expected,part.weight,atol=1e-12,rtol=0)
        audit.append(dict(zip(KEY+['setting'],key),n_pairs=len(part),n_source_libraries=part.library_key_b.nunique(),
            n_target_libraries=part.library_key_a.nunique(),weight_sum=expected.sum()))
    ctx.table(p,'qualified_pairs.tsv.gz');ctx.table(r,'route_effects.tsv.gz');ctx.table(audit,'contributor_audit.tsv.gz')
    ctx.finish(dict(physical_PM_pairing=True,fresh_gene_qualification=True,source_target_capture_weights=True))


def metrics(ctx):
    p=ctx.read(ctx.analysis/'state_vectors/pair_effects.tsv.gz')
    r=ctx.read(ctx.analysis/'state_vectors/route_effects.tsv.gz')
    records=[];correlations=[]
    names=['DWUI','RWUI','PDUI','PPUI','MPRO'];group=['species','study_id','route','setting','contrast','design']
    for level,frame in [('pair',p),('route',r)]:
        for key,f in frame.groupby(group,dropna=False):
            base=dict(zip(group,key),level=level)
            for metric in names:
                x=f['MPRO' if metric=='MPRO' else 'delta_'+metric]
                for threshold in [.1,.2] if metric=='MPRO' else [.02,.05]:
                    records.append(base|dict(metric=metric,threshold=threshold,n=len(f),genes=f.gene_id.nunique(),
                        median=x.median(),mean=x.mean(),positive=int(x.ge(threshold).sum()),
                        negative=0 if metric=='MPRO' else int(x.le(-threshold).sum()),above_threshold=int(x.abs().ge(threshold).sum())))
            for first,second in combinations(names,2):
                a=f['MPRO' if first=='MPRO' else 'delta_'+first];b=f['MPRO' if second=='MPRO' else 'delta_'+second]
                correlations.append(base|dict(first=first,second=second,n=len(f),spearman_r=spearmanr(a,b).statistic if len(f)>1 and a.nunique()>1 and b.nunique()>1 else np.nan))
    ctx.table(records,'metric_summaries.tsv');ctx.table(correlations,'metric_correlations.tsv')
    ctx.table(r,'route_effects.tsv.gz');ctx.table(p,'pair_effects.tsv.gz')
    ctx.finish(dict(full_usage_vectors=True,MPRO_after_vector_weighting=True,five_metrics_all_settings=True))


def dwui_pool(ctx):
    nodes=[]
    for sid in STUDIES:
        study=Study(ctx,sid);f=study.core('node_positions')
        if len(f):nodes.append(f)
    data=pd.concat(nodes,ignore_index=True)
    math=ctx.functions('notes/2026-09-14-integration-reliability/evaluate_reaggregated_apa_effects.py',
        ['study_votes','majority','mean_available'])
    claims=[];votes=[];effects=[];candidate_rows=[]
    for sp,design in CONTRASTS.items():
        source=data[data.species.eq(sp)];genes=pd.Index(sorted(source.gene_id.unique()))
        symbols=source.drop_duplicates('gene_id').set_index('gene_id').gene_symbol.to_dict();records=[]
        for contrast,studies in design.items():
            signed=[];all_deltas=[];ncomplete=np.zeros(len(genes),int)
            for sid,routes in studies.items():
                route_deltas=[]
                for route in routes:
                    f=source[source.study_id.eq(sid)&source.route.eq(route)].set_index(['gene_id','node'])
                    a=f.position_index.reindex(pd.MultiIndex.from_product([genes,[contrast[0]]])).to_numpy()
                    b=f.position_index.reindex(pd.MultiIndex.from_product([genes,[contrast[1]]])).to_numpy()
                    delta=b-a;route_deltas.append(delta);all_deltas.append(delta)
                    effects.extend(dict(species=sp,study_id=sid,route=route,gene_id=g,contrast=contrast,delta=d)
                                   for g,d in zip(genes,delta))
                vote,status=math.study_votes(np.stack(route_deltas),.01);signed.append(vote)
                ncomplete+=np.isfinite(route_deltas).all(axis=0)
                votes.extend(dict(species=sp,study_id=sid,gene_id=g,contrast=contrast,vote=int(v),status=s)
                             for g,v,s in zip(genes,vote,status))
            direction,pos,neg,need=math.majority(np.stack(signed));means=math.mean_available(np.stack(all_deltas))
            for i,gene in enumerate(genes):
                records.append(dict(claim_id=f'APA:{sp}:{gene}:{contrast}',species=sp,gene_id=gene,gene_symbol=symbols[gene],
                    contrast=contrast,delta_definition=contrast[1]+' minus '+contrast[0],baseline_amplitude=.01,
                    baseline_supported=bool(direction[i]),original_full_direction=int(direction[i]),
                    positive_votes=int(pos[i]),negative_votes=int(neg[i]),required_votes=need,planned_studies=len(studies),
                    complete_studies=int(ncomplete[i]),baseline_full_mean_delta=means[i],admission='new_DWUI_fixed_majority_pool'))
        selected={r['gene_id'] for r in records if r['baseline_supported']}
        selected_records=[r|dict(fixed_candidate_denominator=len(selected)) for r in records if r['gene_id'] in selected]
        claims.extend(selected_records)
        for gene in sorted(selected):
            supported=[r['contrast'] for r in selected_records if r['gene_id']==gene and r['baseline_supported']]
            candidate_rows.append(dict(species=sp,gene_id=gene,gene_symbol=symbols[gene],supported_contrasts=';'.join(supported),
                cohort='new_DWUI_fixed_majority_pool'))
    ctx.table(claims,'claim_index.tsv',columns=['claim_id','species','gene_id','gene_symbol','contrast','baseline_supported','original_full_direction'])
    ctx.table(candidate_rows,'candidates.tsv',columns=['species','gene_id','gene_symbol','supported_contrasts','cohort'])
    ctx.table(votes,'all_gene_study_votes.tsv.gz');ctx.table(effects,'all_gene_route_effects.tsv.gz')
    ctx.finish(dict(all_current_profile_genes_considered=True,fixed_planned_study_denominators=True,
        missing_and_conflicting_routes_abstain=True,any_supported_contrast_admits_gene=True,
        all_planned_claims_retained_for_admitted_gene=True),DWUI_threshold=.01)


def union_pool(ctx):
    dwui=ctx.read(ctx.analysis/'dwui_pool/claim_index.tsv')
    rec=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz')
    rec=rec[rec.series_matching&rec.metric.eq('MPRO')&rec.threshold.eq(.1)]
    fields=['species','gene_id','gene_symbol','contrast','baseline_supported','admission']
    claims=dwui[fields].to_dict('records');lookup={(r['species'],r['gene_id'],r['contrast']):r for r in claims};admissions=[]
    for r in rec.to_dict('records'):
        co={'P_S':'SP','P_M':'MP'}[r['contrast']];key=(r['species'],r['gene_id'],co)
        admissions.append(dict(species=r['species'],gene_id=r['gene_id'],gene_symbol=r['gene_symbol'],contrast=co,
            admission='MPRO_recurrence_0p1',already_DWUI_gene_contrast=key in lookup))
        if key in lookup:lookup[key]['admission']+=';MPRO_recurrence_0p1'
        else:
            row=dict(species=r['species'],gene_id=r['gene_id'],gene_symbol=r['gene_symbol'],contrast=co,
                baseline_supported=False,admission='MPRO_recurrence_0p1');claims.append(row);lookup[key]=row
    ctx.table(claims,'figure5_union_claims.tsv',columns=fields)
    ctx.table(admissions,'figure5_MPRO_admission.tsv',columns=fields+['already_DWUI_gene_contrast'])
    ctx.table(pd.DataFrame(claims,columns=fields).groupby('species').agg(genes=('gene_id','nunique'),gene_contrasts=('gene_id','size')).reset_index(),'union_summary.tsv')
    ctx.finish(dict(new_DWUI_pool=True,new_MPRO_recurrence=True,contrast_specific_admission=True,no_fixed_old_gene_counts=True))


def complete_triplets(ps,pm):
    keys=['species','study_id','route','setting','gene_id','gene_symbol','library_key_a']
    renamed={'reference_usage':'M_usage','PAS_ids':'M_PAS_ids','n_mc_b':'M_n_mc',
        'molecules_b':'M_molecules','site_mc_b':'M_site_mc','target_usage':'PM_P_usage'}
    out=ps.merge(pm[keys+list(renamed)].rename(columns=renamed),on=keys,validate='many_to_one')
    assert out.PAS_ids.eq(out.M_PAS_ids).all()
    assert out.target_usage.eq(out.PM_P_usage).all()
    if len(out):
        counts=np.stack(out.molecules_b.map(json.loads))+np.stack(out.M_molecules.map(json.loads))
        metacells=np.stack(out.site_mc_b.map(json.loads))+np.stack(out.M_site_mc.map(json.loads))
        # P--S and P--M eligibility does not imply support for the third contrast.
        out['MS_joint_PAS']=((counts>=10)&(metacells>=2)).sum(1)
    else:out['MS_joint_PAS']=pd.Series(dtype=int)
    return out[out.MS_joint_PAS.ge(2)].copy()


def matched(ctx):
    pairs=ctx.read(ctx.analysis/'state_vectors/pair_effects.tsv.gz')
    axes=ctx.read(ctx.analysis/'state_vectors/PAS_axes.tsv.gz')
    lookup={(s,g,t):a for (s,g,t),a in axes.groupby(['study_id','primary_gene_id','setting'],sort=False)}
    rows=[];triplets=[]
    for setting in SETTINGS:
        q=pairs[pairs.setting.eq(setting)];ps=q[q.contrast.eq('P_S')&q.design.eq('true_source')]
        pm=q[q.contrast.eq('P_M')]
        # PAS-vector lengths differ between genes, so qualify each gene separately.
        matched=pd.concat([complete_triplets(part,pm[pm.study_id.eq(sid)&pm.gene_id.eq(gene)])
            for (sid,gene),part in ps.groupby(['study_id','gene_id'])],
            ignore_index=True) if len(ps) else ps.copy()
        for key,f in matched.groupby(['species','study_id','route','setting','gene_id'],sort=False):
            w=weights(f);axis=lookup[key[1],key[4],setting];linear=axis[['WPP','RWUI','PDUI','PPUI']].to_numpy()
            usage={'P':np.stack(f.target_usage.map(json.loads)),'M':np.stack(f.M_usage.map(json.loads)),'S':np.stack(f.reference_usage.map(json.loads))}
            means={s:w@v for s,v in usage.items()};record=dict(zip(['species','study_id','route','setting','gene_id'],key))
            record.update(gene_symbol=f.gene_symbol.iloc[0],PAS_ids=f.PAS_ids.iloc[0],PAS_number=len(axis),n_triplets=len(f),
                n_target_libraries=f.library_key_a.nunique(),n_source_libraries=f.library_key_b.nunique(),n_target_dates=f.time_a.nunique())
            for state,value in means.items():record[state+'_usage']=json.dumps(value.tolist());record[state+'_DWUI']=float(value@linear[:,0])
            deltas={}
            for co,a,b in [('PS','P','S'),('PM','P','M'),('MS','M','S')]:
                delta=means[a]-means[b];deltas[co]=delta;record['delta_usage_'+co]=json.dumps(delta.tolist())
                for j,metric in enumerate(['DWUI','RWUI','PDUI','PPUI']):record[metric+'_'+co]=float(delta@linear[:,j])
                record['MPRO_'+co]=float(np.ptp(delta));record['mean_triplet_MPRO_'+co]=float(w@np.ptp(usage[a]-usage[b],axis=1))
                assert record['MPRO_'+co]<=record['mean_triplet_MPRO_'+co]+1e-12
            assert np.allclose(deltas['PS'],deltas['PM']+deltas['MS'],atol=1e-12,rtol=0)
            rows.append(record);f=f.copy();f['triplet_weight']=w;triplets.append(f)
    ctx.table(rows,'matched_state_vectors.tsv.gz',columns=['species','study_id','route','setting','gene_id','n_triplets'])
    ctx.table(pd.concat(triplets) if triplets else pairs.iloc[:0],'weighted_triplets.tsv.gz')
    recurrent=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz');recurrent=recurrent[recurrent.series_matching]
    inventory=[]
    for sp in ['mouse','human']:
        for co in ['P_S','P_M']:
            for wt,mt in [(.02,.1),(.05,.2)]:
                f=recurrent[recurrent.species.eq(sp)&recurrent.contrast.eq(co)]
                a=set(f[f.metric.eq('WPP')&f.threshold.eq(wt)].gene_id);b=set(f[f.metric.eq('MPRO')&f.threshold.eq(mt)].gene_id)
                for gene in sorted(a|b):inventory.append(dict(species=sp,contrast=co,DWUI_threshold=wt,MPRO_threshold=mt,gene_id=gene,
                    cohort='shared' if gene in a&b else 'DWUI_only' if gene in a else 'MPRO_only'))
    ctx.table(inventory,'metric_membership.tsv.gz',columns=['species','contrast','DWUI_threshold','MPRO_threshold','gene_id','cohort'])
    tests=ctx.read(ctx.analysis/'well_inference/gene_metric_tests.tsv.gz')
    wide=tests.pivot(index='gene_id',columns='metric',values='q_BH_joint')
    wide['DWUI_only']=wide.DWUI.lt(.05)&wide.MPRO.ge(.05);wide['MPRO_only']=wide.MPRO.lt(.05)&wide.DWUI.ge(.05)
    ctx.table(wide.reset_index(),'well_significance_overlap.tsv')
    ctx.finish(dict(new_complete_PMS_triplets=True,same_contributors_all_contrasts=True,
        full_vector_additivity=True,MPRO_nonadditivity_preserved=True,new_metric_membership=True,
        all_three_contrasts_joint_PAS_qualified=True,all_four_settings=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--node',choices=['endpoints','metrics','dwui_pool','union_pool','matched_pms'],required=True)
    a=p.parse_args();{'endpoints':endpoints,'metrics':metrics,'dwui_pool':dwui_pool,'union_pool':union_pool,'matched_pms':matched}[a.node](Context(a.config,a.node,a.out))
