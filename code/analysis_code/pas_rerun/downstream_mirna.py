"""External miRNA groups on each run's newly selected lost windows.

External expression defines groups; these tests do not establish binding/activity.
All eight primary and thirty combined control slots remain present when empty.
"""
import argparse
from collections import defaultdict
import itertools
import pickle
import xml.etree.ElementTree as ET
import zipfile
import numpy as np
import pandas as pd
from downstream_runtime import Context,bh

SOURCE='__mirna_math__'
EXTERNAL=['human_mirna_annotations.tsv','mouse_mirna_annotations.tsv','human_S2.xlsx','GSM307156.soft','GSM307157.soft']
MARKERS={'barrier':'GGCAGTG','ESCC':'AAGTGCT'}
TESTS=[('marker_barrier_minus_ESCC','barrier','ESCC'),('external_somatic_minus_PSC','somatic','PSC')]
SITE_COLS=['gene_id','seed','chrom','strand','start0','end0']
WINDOW_COLS=['species','contrast','gene_id','start0','end0','length','control_start0','control_end0','group',
    'n_seeds','lost_sites','control_sites','lost_density','control_density','difference','random_control_density','random_control_sd']


def seed_sets(expression):
    return dict(barrier={MARKERS['barrier']},ESCC={MARKERS['ESCC']},
        somatic=set(expression.loc[expression.expression_group.eq('somatic'),'seed']),
        PSC=set(expression.loc[expression.expression_group.eq('PSC'),'seed']))


def grouped_difference(frame,left,right,nleft,nright,value='difference'):
    if not nleft or not nright:return pd.Series(dtype=float,name='difference')
    w=frame.pivot(index='gene_id',columns='group',values=value).reindex(columns=[left,right])
    return (w[left]/nleft-w[right]/nright).dropna()


def influence(values):
    n=len(values)
    return pd.DataFrame(dict(gene_id=values.index,difference=values.to_numpy(),n_genes=n,
        leave_one_out_mean=(values.sum()-values.to_numpy())/(n-1) if n>1 else np.full(n,np.nan)))


def primary(ctx,num,rng):
    annotations={s:pd.read_csv(ctx.resource('mirna_'+s+'_mirna_annotations.tsv'),sep='\t') for s in ['human','mouse']}
    expressions={'human':num.classify(num.human_expression(annotations['human'])),
                 'mouse':num.classify(num.mouse_expression(annotations['mouse']))}
    expression=pd.concat([d.assign(species=s) for s,d in expressions.items()],ignore_index=True)
    ctx.table(expression,'external_seed_expression.tsv')
    with ctx.register(ctx.analysis/'motif_background/window_context.pickle').open('rb') as f:contexts,changed,_=pickle.load(f)
    lost=changed[changed.max_effect.le(-.1)].copy()
    windows=ctx.read(ctx.analysis/'motif_background/matched_windows.tsv.gz')
    windows=windows[windows.category.eq('miRNA')].merge(lost,
        on=['species','contrast','gene_id','start0','end0','length'],validate='one_to_one',suffixes=('','_atomic'))
    hits={};audit=[]
    for species in ['human','mouse']:
        raw=ctx.read(ctx.analysis/f'site_scan/full-cohort-v2/data/{species}/site_effects.tsv.gz',
            usecols=['gene_id','category','chrom','strand','start0','end0','feature_groups'])
        raw=raw[raw.category.eq('miRNA')].drop_duplicates()
        raw['feature_groups']=raw.feature_groups.str.split(';');raw=raw.explode('feature_groups').drop_duplicates()
        assert raw.feature_groups.str.fullmatch('[ACGT]{7}').all()
        physical=num.collapse_sites(raw.rename(columns={'feature_groups':'seed'})) if len(raw) else pd.DataFrame(columns=SITE_COLS)
        hits[species]=dict(tuple(physical.groupby('gene_id')));ctx.table(physical,f'{species}_physical_seed_sites.tsv.gz')
        e=expressions[species]
        audit.append(dict(species=species,scan_rows_unique=len(raw),physical_seed_sites=len(physical),tested_seed_universe=physical.seed.nunique(),
            external_seeds=len(e),somatic_seeds=int(e.expression_group.eq('somatic').sum()),PSC_seeds=int(e.expression_group.eq('PSC').sum())))
    ctx.table(audit,'seed_universe.tsv');sets={s:seed_sets(d) for s,d in expressions.items()}
    rows=[];full=[]
    for q in lost.itertuples():
        for label,seeds in sets[q.species].items():
            full.append(dict(species=q.species,contrast=q.contrast,gene_id=q.gene_id,start0=q.start0,end0=q.end0,length=q.length,
                group=label,n_seeds=len(seeds),lost_sites=num.count(hits[q.species],q.gene_id,seeds,q.start0,q.end0),effect20=q.effect20,effect50=q.effect50))
    ctx.table(full,'all_lost_window_inventory.tsv',columns=['species','contrast','gene_id','start0','end0','length','group','n_seeds','lost_sites','effect20','effect50'])
    for q in windows.itertuples():
        c=num.candidates(contexts[q.species,q.contrast,q.gene_id],q.length,q.GC,q.relative_position)
        c=c[(c[:,2]<=.05)&(c[:,3]<=.2)];assert len(c)
        starts=rng.choice(c[:,0].astype(int),size=100,replace=True)
        for label,seeds in sets[q.species].items():
            n=num.count(hits[q.species],q.gene_id,seeds,q.start0,q.end0)
            nc=num.count(hits[q.species],q.gene_id,seeds,q.control_start0,q.control_end0)
            rn=np.array([num.count(hits[q.species],q.gene_id,seeds,a,a+q.length) for a in starts])
            rows.append(dict(species=q.species,contrast=q.contrast,gene_id=q.gene_id,start0=q.start0,end0=q.end0,length=q.length,
                control_start0=q.control_start0,control_end0=q.control_end0,group=label,n_seeds=len(seeds),lost_sites=n,control_sites=nc,
                lost_density=1000*n/q.length,control_density=1000*nc/q.length,difference=1000*(n-nc)/q.length,
                random_control_density=1000*rn.mean()/q.length,random_control_sd=1000*rn.std()/q.length))
    w=pd.DataFrame(rows,columns=WINDOW_COLS);ctx.table(w,'matched_family_windows.tsv')
    genes=w.groupby(['species','contrast','gene_id','group'],as_index=False).agg(difference=('difference','mean'),
        lost_density=('lost_density','mean'),control_density=('control_density','mean'),random_control_density=('random_control_density','mean'),n_windows=('start0','size'))
    ctx.table(genes,'gene_family_densities.tsv');results=[];random_control=[];attrition=[]
    for species in ['human','mouse']:
        for co in ['SP','MP']:
            g=genes[genes.species.eq(species)&genes.contrast.eq(co)].copy()
            for test,left,right in TESTS:
                nl,nr=len(sets[species][left]),len(sets[species][right])
                r=dict(species=species,contrast=co,test=test,n_left_seeds=nl,n_right_seeds=nr)
                r.update(num.infer(grouped_difference(g,left,right,nl,nr)))
                if not nl or not nr:r['status']='missing_external_group'
                r['external_scope']=expressions[species].scope.iloc[0];results.append(r)
            g['random_diff']=g.lost_density-g.random_control_density
            for test,left,right in TESTS:
                v=grouped_difference(g,left,right,len(sets[species][left]),len(sets[species][right]),'random_diff')
                random_control.append(dict(species=species,contrast=co,test=test,n_genes=len(v),mean_random_background_difference=v.mean()))
            a=lost[lost.species.eq(species)&lost.contrast.eq(co)];b=windows[windows.species.eq(species)&windows.contrast.eq(co)]
            attrition.append(dict(species=species,contrast=co,all_lost_windows=len(a),all_lost_genes=a.gene_id.nunique(),
                matched_windows=len(b),matched_genes=b.gene_id.nunique(),bp_all=a.length.sum(),bp_matched=b.length.sum()))
    result=pd.DataFrame(results);assert len(result)==8;result['q_BH8']=bh(result.p)
    ctx.table(result,'primary_tests.tsv');ctx.table(random_control,'random_control_sensitivity.tsv');ctx.table(attrition,'attrition.tsv')
    return dict(annotations=annotations,expressions=expressions,sets=sets,contexts=contexts,windows=windows,hits=hits,genes=genes,primary=result)


def mapping_sensitivity(ctx,num,state):
    ann=state['annotations']['mouse'];prefixes=defaultdict(list)
    for r in ann.itertuples():
        seq=r.mature_sequence.replace('U','T');prefixes[seq[:18]].append((len(seq),r.seed_dna))
    rows=[];audit=[]
    for acc,label in [('GSM307156','PSC'),('GSM307157','somatic')]:
        seeds=defaultdict(int);active=False;total=0
        with ctx.resource('mirna_'+acc+'.soft').open() as f:
            for line in f:
                if line.startswith('!sample_table_begin'):active=True;continue
                if line.startswith('!sample_table_end'):break
                if not active or line.startswith('SEQUENCE'):continue
                seq,ct=line.strip().split('\t');n=int(ct);total+=n;possible=set()
                for length,seed in prefixes.get(seq[:18],[]):
                    if any(seq[j:].startswith('TCGTAT') for j in range(max(18,length-2),length+3)):possible.add(seed)
                if len(possible)==1:seeds[next(iter(possible))]+=n
        mapped=sum(seeds.values());assert mapped>0 and total>0
        audit.append(dict(accession=acc,state=label,total=total,mapped=mapped,mapped_fraction=mapped/total))
        for seed in sorted(set(ann.seed_dna)):rows.append(dict(seed=seed,state=label,count=seeds[seed],mapped_RPM=1e6*seeds[seed]/mapped,total_RPM=1e6*seeds[seed]/total))
    d=pd.DataFrame(rows);ctx.table(audit,'isomir_mapping_sensitivity_audit.tsv');ctx.table(d,'isomir_expression_sensitivity.tsv')
    original=state['expressions']['mouse'];rows=[];all_sets=[];windows=state['windows']
    num.infer.__globals__['RNG']=np.random.default_rng(20261004)
    for setting,frame,field in [('exact_totalRPM',original,None),('isomir_mappedRPM',d,'mapped_RPM'),('isomir_totalRPM',d,'total_RPM')]:
        e=frame.pivot(index='seed',columns='state',values=field).reset_index() if field else frame.rename(
            columns={'PSC':'main_PSC','somatic':'main_somatic','PSC_total_RPM':'PSC','somatic_total_RPM':'somatic'})[['seed','PSC','somatic']]
        e=num.classify(e);e['setting']=setting;all_sets.append(e)
        left=set(e.loc[e.expression_group.eq('somatic'),'seed']);right=set(e.loc[e.expression_group.eq('PSC'),'seed'])
        for co in ['SP','MP']:
            rr=[]
            if left and right:
                for q in windows[windows.species.eq('mouse')&windows.contrast.eq(co)].itertuples():
                    v=[1000*(num.count(state['hits']['mouse'],q.gene_id,s,q.start0,q.end0)-num.count(state['hits']['mouse'],q.gene_id,s,q.control_start0,q.control_end0))/q.length/len(s) for s in [left,right]]
                    rr.append(dict(gene_id=q.gene_id,difference=v[0]-v[1]))
            g=pd.DataFrame(rr,columns=['gene_id','difference']).groupby('gene_id').difference.mean();r=num.infer(g)
            if not left or not right:r['status']='missing_external_group'
            r.update(setting=setting,contrast=co,n_somatic_seeds=len(left),n_PSC_seeds=len(right));rows.append(r)
    result=pd.DataFrame(rows);assert len(result)==6;result['q_sensitivity_BH6']=bh(result.p)
    ctx.table(result,'expression_mapping_sensitivity_tests.tsv');ctx.table(pd.concat(all_sets),'expression_mapping_sensitivity_sets.tsv')
    g=state['genes'];g=g[g.species.eq('mouse')&g.contrast.eq('SP')]
    s=state['sets']['mouse'];v=grouped_difference(g,'somatic','PSC',len(s['somatic']),len(s['PSC']))
    ctx.table(influence(v),'mouse_SP_influence.tsv')
    return result


def position_tests(ctx,num,state):
    n=10000;rng=np.random.default_rng(2026100401);tests=[];mobility=[];arrays={}
    for sp in ['human','mouse']:
        sets=state['sets'][sp]
        for co in ['SP','MP']:
            w=state['windows'];w=w[w.species.eq(sp)&w.contrast.eq(co)];ng=w.gene_id.nunique()
            uniforms={g:rng.random(n) for g in w.gene_id.unique()};counts=w.groupby('gene_id').size()
            stores={test:dict(observed=0.,coupled=np.zeros(n),independent=np.zeros(n)) for test,_,_ in TESTS}
            for q in w.itertuples():
                c=num.candidates(state['contexts'][sp,co,q.gene_id],q.length,q.GC,q.relative_position)
                starts=c[(c[:,2]<=.05)&(c[:,3]<=.2),0].astype(int);assert len(starts)>0 and np.all(np.diff(starts)>=0)
                h=state['hits'][sp].get(q.gene_id,pd.DataFrame(columns=SITE_COLS));dens={}
                for group,seeds in sets.items():
                    a=h.loc[h.seed.isin(seeds),'start0'].to_numpy();b=h.loc[h.seed.isin(seeds),'end0'].to_numpy()
                    dens[group]=1000*((a[None,:]>=starts[:,None])&(b[None,:]<=starts[:,None]+q.length)).sum(1)/q.length/max(1,len(seeds))
                ix=np.minimum((uniforms[q.gene_id]*len(starts)).astype(int),len(starts)-1);ix2=rng.integers(0,len(starts),size=n)
                weight=1/counts[q.gene_id]/ng
                for test,left,right in TESTS:
                    if not sets[left] or not sets[right]:continue
                    values=dens[left]-dens[right]
                    obs=1000*(num.count(state['hits'][sp],q.gene_id,sets[left],q.start0,q.end0)/len(sets[left])-num.count(state['hits'][sp],q.gene_id,sets[right],q.start0,q.end0)/len(sets[right]))/q.length
                    stores[test]['observed']+=weight*obs;stores[test]['coupled']+=weight*values[ix];stores[test]['independent']+=weight*values[ix2]
                    mobility.append(dict(species=sp,contrast=co,gene_id=q.gene_id,start0=q.start0,end0=q.end0,test=test,
                        n_feasible_positions=len(starts),n_unique_null_densities=len(np.unique(values)),null_density_sd=values.std(),gene_window_weight=weight))
            for test,left,right in TESTS:
                for model in ['coupled','independent']:
                    values=stores[test][model];observed=stores[test]['observed'];mean=values.mean();lo,hi=np.quantile(values,[.025,.975])
                    status='conditional_position_randomization';p=(1+np.sum(abs(values-mean)>=abs(observed-mean)-1e-12))/(n+1)
                    if ng<10:p=1.;status='descriptive_fewer_than_10_genes'
                    elif values.std()<1e-12:p=1.;status='frozen_position_null'
                    if not ng or not sets[left] or not sets[right]:
                        p=1.;status='no_genes' if not ng else 'missing_external_group';observed=mean=lo=hi=np.nan
                        values=np.full(n,np.nan)
                    tests.append(dict(species=sp,contrast=co,test=test,model=model,n_genes=ng,n_windows=len(w),
                        observed_density_contrast=observed,null_mean=mean,null_sd=values.std(),effect=observed-mean,
                        null_low=lo,null_high=hi,draws=n,p=float(p),status=status))
                    arrays[f'{sp}_{co}_{test}_{model}']=values
    out=pd.DataFrame(tests);assert len(out)==16
    for model,col in [('coupled','q_position_primary_BH8'),('independent','q_position_sensitivity_BH8')]:
        select=out.model.eq(model);out.loc[select,col]=bh(out.loc[select,'p'])
    out['q_position_complete_BH16']=bh(out.p);ctx.table(out,'position_null_tests.tsv')
    ctx.table(mobility,'position_null_mobility.tsv',columns=['species','contrast','gene_id','start0','end0','test','n_feasible_positions','n_unique_null_densities','null_density_sd','gene_window_weight'])
    np.savez_compressed(ctx.out/'position_null_draws.npz',**arrays)
    return out


def run(ctx):
    assert ctx.config['reference_variant']=='v2_primary'
    files=[ctx.resource('mirna_'+name) for name in EXTERNAL]
    assert len({p.parent for p in files})==1 and {p.name for p in files}==set(EXTERNAL)
    rng=np.random.default_rng(20261004)
    num=ctx.functions(SOURCE,['xlsx_rows','human_expression','mouse_expression','classify','collapse_sites','count','infer','candidates'],
        dict(defaultdict=defaultdict,ET=ET,zipfile=zipfile,itertools=itertools,HERE=ctx.out,EXT=files[0].parent,RNG=rng))
    state=primary(ctx,num,rng);sensitivity=mapping_sensitivity(ctx,num,state);position=position_tests(ctx,num,state)
    a=state['primary'].assign(family='gene_sign',model='primary')
    b=sensitivity.assign(species='mouse',test='external_somatic_minus_PSC',family='expression_sensitivity',model=lambda x:x.setting)
    cols=['species','contrast','test','family','model','p']
    total=pd.concat([a[cols],b[cols],position.assign(family='position')[cols]],ignore_index=True)
    assert len(total)==30;total['q_complete_BH30']=bh(total.p);ctx.table(total,'complete_test_family.tsv')
    ctx.finish(dict(current_lost_windows=True,external_groups_not_binding_evidence=True,exact_seed_coordinates=True,
        complete_8_primary_6_mapping_16_position_tests=True,gene_equal_weighted=True,unavailable_not_zero=True,
        no_fixed_17_gene_or_2920bp_requirement=True,full_seed_inventory_not_single_highlight=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'external_mirna',a.out))
