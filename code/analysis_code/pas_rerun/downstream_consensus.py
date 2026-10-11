"""Six-rule recurrence and joint direction consensus on fresh recount vectors."""
import argparse
import itertools

import numpy as np
import pandas as pd

from downstream_runtime import Context,bh,seed
from downstream_recount import MODES
from downstream_recurrence import events_from_routes

SOURCE='notes/2026-09-30-multimetric-statistics/nulls/analyze.py'


def consensus(matrix):
    positive=(matrix==1).sum(1);negative=(matrix==-1).sum(1)
    recurrent=np.maximum(positive,negative)>=2;direction=positive>=negative
    all_recurrent=recurrent.all(1);same=direction.all(1)|(~direction).all(1)
    return all_recurrent&same,all_recurrent,direction


def joint_draw(matrix,groups,rng):
    shuffled=matrix.copy()
    for unit,cells in groups.items():
        for ix in cells:shuffled[ix,unit,:]=matrix[rng.permutation(ix),unit,:]
    return shuffled


def run(ctx):
    route=ctx.read(ctx.analysis/'recount/route_effects.tsv.gz')
    math=ctx.functions(SOURCE,['qbin','collapse','stratify','shuffle','summarize'],dict(THRESHOLDS={'WPP':[.02,.05]}))
    summary=[];members=[];events=[];diagnostics=[];archives={};joint_summary=[];joint_members=[]
    schema=['species','contrast','unit','gene_id','gene_symbol','n_routes','mean_GEX','APA_coverage','PAS_number',
        'supported_catalogue_PAS_number','metric','effect_mean','effect_min','effect_max','baseline','baseline_min','baseline_max',
        'capacity_negative','capacity_positive','capacity_unsigned','threshold','event','eligible_units','series_pattern','view','mode']
    for view,mode in itertools.product(['fixed_genes','common_contributors'],MODES):
        r=route[route.view.eq(view)&route['mode'].eq(mode)].copy();r['delta_WPP']=r.delta
        e=events_from_routes(r,math) if len(r) else pd.DataFrame(columns=schema)
        e['view']=view;e['mode']=mode;events.append(e)
        for sp,t in itertools.product(['human','mouse'],[.02,.05]):
            key=dict(view=view,mode=mode,species=sp,contrast='P_S',threshold=t);q=e[e.species.eq(sp)&e.threshold.eq(t)]
            genes=pd.Index(sorted(q.gene_id.unique()));rng=np.random.default_rng(seed('six_rule',*key.values()))
            pos=np.zeros((10000,len(genes)),np.int8);neg=np.zeros_like(pos);op=np.zeros(len(genes),int);om=np.zeros(len(genes),int)
            for unit,f in q.groupby('unit'):
                f=math.stratify(f,'WPP',t);f['stratum']=f.stratum+'|'+f.series_pattern
                draw,_,diag=math.shuffle(f,rng,10000);ix=genes.get_indexer(f.gene_id);pos[:,ix]+=draw==1;neg[:,ix]+=draw==-1
                op[ix]+=f.event.eq(1).to_numpy();om[ix]+=f.event.eq(-1).to_numpy();diagnostics.append(diag.assign(unit=unit,**key))
            observed=np.maximum(op,om)>=2;null=(np.maximum(pos,neg)>=2).sum(1)
            summary.append(dict(key,**math.summarize(observed.sum(),null),genes=len(genes),studies=q.unit.nunique(),
                status='tested' if q.unit.nunique()>=2 else 'fewer_than_two_studies'))
            for i in np.flatnonzero(observed):members.append(dict(key,gene_id=genes[i],direction='positive' if op[i]>=om[i] else 'negative'))
            archives['single|'+'|'.join(map(str,key.values()))]=null
    all_events=pd.concat(events,ignore_index=True);ctx.table(all_events,'study_events.tsv.gz')
    joint_strata=[]
    for view,sp,t in itertools.product(['fixed_genes','common_contributors'],['human','mouse'],[.02,.05]):
        key=dict(view=view,species=sp,contrast='P_S',threshold=t);q=all_events[all_events.view.eq(view)&all_events.species.eq(sp)&all_events.threshold.eq(t)]
        genes=pd.Index(sorted(q.gene_id.unique()));units=sorted(q.unit.unique());mat=np.full((len(genes),len(units),6),9,np.int8);groups={}
        for u,unit in enumerate(units):
            part=q[q.unit.eq(unit)];sets=[]
            for k,mode in enumerate(MODES):
                f=part[part['mode'].eq(mode)].copy()
                if not len(f):sets.append(pd.DataFrame(columns=['gene_id','stratum_'+mode,'series_'+mode]));continue
                f=math.stratify(f,'WPP',t);mat[genes.get_indexer(f.gene_id),u,k]=f.event
                sets.append(f[['gene_id','stratum','series_pattern']].rename(columns={'stratum':'stratum_'+mode,'series_pattern':'series_'+mode}))
            merged=sets[0]
            for frame in sets[1:]:merged=merged.merge(frame,on='gene_id',validate='one_to_one')
            # These two views require the same study-gene measurability in every mode.
            assert len(merged)==part.gene_id.nunique(),('incomplete_six_mode_vector',view,sp,unit,t)
            strata=merged.drop(columns='gene_id').astype(str).agg('||'.join,axis=1) if len(merged) else pd.Series(dtype=str)
            cells=[]
            for label,ix in merged.groupby(strata).indices.items():
                ids=genes.get_indexer(merged.gene_id.iloc[ix]);cells.append(ids)
                joint_strata.append(dict(key,unit=unit,stratum=label,n=len(ids),distinct_vectors=len(np.unique(mat[ids,u,:],axis=0))))
            groups[u]=cells
        observed,all_recurrent,directions=consensus(mat);rng=np.random.default_rng(seed('joint_consensus',*key.values()));null=np.empty(10000,int)
        # One mode vector is moved intact; avoid a 10000 x gene x study x mode allocation.
        for draw in range(len(null)):
            null[draw]=consensus(joint_draw(mat,groups,rng))[0].sum()
        joint_summary.append(dict(key,**math.summarize(observed.sum(),null),genes=len(genes),studies=len(units),
            observed_all_mode_recurrent_any_direction=int(all_recurrent.sum()),status='tested' if len(units)>=2 else 'fewer_than_two_studies'))
        for i in np.flatnonzero(observed):joint_members.append(dict(key,gene_id=genes[i],direction='positive' if directions[i,0] else 'negative'))
        archives['joint|'+'|'.join(map(str,key.values()))]=null
    s=pd.DataFrame(summary);s['q_four_per_view_mode']=s.groupby(['view','mode']).empirical_p.transform(bh)
    j=pd.DataFrame(joint_summary);j['q_four_per_view']=j.groupby('view').empirical_p.transform(bh)
    ctx.table(s,'six_rule_recurrence_summary.tsv');ctx.table(j,'joint_consensus_summary.tsv')
    ctx.table(members,'six_rule_members.tsv.gz',columns=['view','mode','species','contrast','threshold','gene_id','direction'])
    ctx.table(joint_members,'joint_consensus_members.tsv.gz',columns=['view','species','contrast','threshold','gene_id','direction'])
    ctx.table(pd.concat(diagnostics,ignore_index=True) if diagnostics else [],'mode_strata.tsv.gz',columns=['view','mode','species','threshold','unit'])
    ctx.table(joint_strata,'joint_strata.tsv.gz',columns=['view','species','contrast','threshold','unit','stratum','n','distinct_vectors'])
    np.savez_compressed(ctx.out/'null_draws.npz',**archives)
    ctx.finish(dict(new_all_gene_counting_rule_events=True,full_six_mode_vectors_shuffled_together=True,
        exact_series_strata=True,new_counting_covariates=True,four_tests_per_view=True,empty_positions_P1=True),P_M_scope='not recounted')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'consensus',a.out))
