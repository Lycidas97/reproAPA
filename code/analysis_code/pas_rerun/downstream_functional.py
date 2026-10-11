"""Rebuild expression/testability/recurrence universes and exact functional nulls."""
import argparse
import csv
import gzip
from collections import defaultdict
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import logsumexp
from scipy.stats import hypergeom

from downstream_runtime import Context,Study,bh

SOURCE='github/reproapa-analysis/revisions/20260930/external/2026-09-26_reproapa-go-analysis/scripts/'
DESIGN={('mouse','P_S'):{'states':('S','P'),'studies':{'GSE103221':['CIP'],'GSE115943':['WOT']}},
    ('mouse','P_M'):{'states':('M','P'),'studies':{'GSE103221':['CIP'],'GSE115943':['WOT']}},
    ('human','P_S'):{'states':('S','P'),'studies':{'GSE118258':['GSE118258'],'GSE178325':['178325_0618','178325_0809']}},
    ('human','P_M'):{'states':('M','P'),'studies':{'GSE118258':['GSE118258'],'GSE150311':['GSE150311'],'GSE221739':['GSE221739']}}}


def annotations(ctx):
    root=ctx.resource('go_obo').parent
    for species in ['mouse','human']:
        assert ctx.resource('go_'+species)==root/f'{species}_go_ids.tsv'
    base=ctx.functions(SOURCE+'run_go.py',['ontology','load_go_genes','covariate_table','add_bins'],
        dict(SRC=root,defaultdict=defaultdict,lru_cache=lru_cache,ROOT_TERM='GO:0008150'))
    terms,aliases,ancestors=base.ontology()
    result={'GO':{sp:base.load_go_genes(sp,terms,aliases,ancestors) for sp in ['mouse','human']},'KEGG':{}}
    names={'GO':terms,'KEGG':{}}
    kroot=ctx.resource('kegg_human_gene_info').parent
    config={}
    for sp,org,title,prefix,suffix in [('human','hsa','Homo_sapiens','ENSG','Homo sapiens (human)'),
                                      ('mouse','mmu','Mus_musculus','ENSMUSG','Mus musculus (house mouse)')]:
        config[sp]=dict(organism=org,gene_info=title+'.gene_info.gz',ensembl_prefix=prefix,species_suffix=' - '+suffix)
        assert ctx.resource('kegg_'+sp+'_gene_info')==kroot/(title+'.gene_info.gz')
        for kind in ['ncbi_conv','gene_pathway','pathways']:
            assert ctx.resource('kegg_'+org+'_'+kind)==kroot/f'{org}_{kind}.tsv'
    kegg=ctx.functions(SOURCE+'run_kegg_nested.py',['parse_gene_info','parse_two_column','load_kegg'],
        dict(SRC=kroot,SPECIES_CONFIG=config,gzip=gzip,csv=csv))
    for sp in ['mouse','human']:
        result['KEGG'][sp],knames,mapping=kegg.load_kegg(sp);names['KEGG'].update(knames)
        ctx.table(mapping,f'{sp}_kegg_ensembl_mapping.tsv.gz')
    return result,names,base


def expression_universe(ctx):
    cache={}
    for sid in sorted({s for d in DESIGN.values() for s in d['studies']}):
        study=Study(ctx,sid);study.load_gex();rows=[]
        for (route,lib,node),part in study.meta[study.meta.purity.ge(.8)].groupby(['route','library_key','node'],sort=False):
            if len(part)<3 or part.n_cells.sum()<30:continue
            ix=part.index.to_numpy();c=study.gex_counts[ix].tocsr();total=study.gex_total[ix]
            scaled=c.multiply(1e6/total[:,None]).tocsr();scaled.data=np.log2(1+scaled.data)
            mean=np.asarray(scaled.mean(0)).ravel();summed=np.asarray(c.sum(0)).ravel();good=summed>0
            f=study.gex_genes.loc[good,['gene_id','gene_symbol']].copy()
            f['GEX']=mean[good];f['GEX_CPM']=summed[good]/total.sum()*1e6;f['GEX_count']=summed[good]
            f['route']=route;f['library_key']=lib;f['state']=node;rows.append(f)
        cache[sid]=pd.concat(rows,ignore_index=True)
        ctx.table(cache[sid],f'full_gex/{sid}.tsv.gz')
    nested=ctx.functions(SOURCE+'run_go_nested.py',['build_g0','add_expression_strata','term_family'],
        dict(DESIGN=DESIGN,read_full_gex=lambda s:cache[s].copy(),MIN_TERM=10,MAX_TERM=500,defaultdict=defaultdict))
    g0=nested.add_expression_strata(nested.build_g0())
    ctx.table(g0,'G0_expression_covariates.tsv.gz')
    return g0,nested


def functional(ctx):
    annotation,names,base=annotations(ctx);g0,nested=expression_universe(ctx)
    events=ctx.read(ctx.analysis/'recurrence/study_event_covariates.tsv.gz')
    events=events[events.metric.eq('WPP')&events.threshold.eq(.02)].copy();events['reference_WPP']=events.baseline
    recurrent=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz')
    recurrent=recurrent[recurrent.metric.eq('WPP')&recurrent.threshold.eq(.02)&recurrent.series_matching]
    math=ctx.functions('notes/2026-10-02-null-model-controls/functional.py',
        ['component','distribution','tail','rejection_probability','setup'],dict(hypergeom=hypergeom,logsumexp=logsumexp,lru_cache=lru_cache))
    membership=[];attrition=[];frozen=[];audits=[];power=[];mobility=[]
    for database in ['GO','KEGG']:
        records=[]
        for species,contrast in DESIGN:
            expr=g0[g0.species.eq(species)&g0.contrast.eq(contrast)].copy()
            apa=base.add_bins(base.covariate_table(events,species,contrast));apa=apa[apa.eligible_units.ge(2)]
            rec=recurrent[recurrent.species.eq(species)&recurrent.contrast.eq(contrast)]
            g0_ids=set(expr.gene_id);raw=set(apa.gene_id);g1=raw&g0_ids;rawg2=set(rec.gene_id)
            # New PAS genes absent from expression are an explicit attrition, not grounds to import old GEX.
            g2=rawg2&g1
            for gene in sorted(raw-g0_ids):attrition.append(dict(species=species,contrast=contrast,gene_id=gene,
                recurrent=gene in rawg2,reason='PAS_testable_but_outside_expression_G0'))
            assert rawg2<=raw
            apa=apa[apa.gene_id.isin(g1)]
            if database=='GO':
                direction=rec.set_index('gene_id').direction.to_dict()
                for gene in sorted(g0_ids):membership.append(dict(species=species,contrast=contrast,gene_id=gene,
                    in_G0=True,in_G1_raw=gene in raw,in_G1_nested=gene in g1,in_G2=gene in g2,recurrent_direction=direction.get(gene,'')))
            targets={'all':g2,'proximal':set(rec[rec.direction.eq('negative')].gene_id)&g2,
                     'distal':set(rec[rec.direction.eq('positive')].gene_id)&g2}
            comparisons=[('G2_vs_G0',target,ids,expr,'expression_stratum') for target,ids in targets.items()]
            comparisons+=[('G1_vs_G0','all',g1,expr,'expression_stratum')]
            comparisons+=[('G2_vs_G1',target,ids,apa,'stratum_full') for target,ids in targets.items()]
            ann=annotation[database][species]
            for comparison,target,ids,bg,stratum in comparisons:
                aid='|'.join([species,contrast,comparison,target]);before=len(bg)
                bg=bg[bg.gene_id.map(lambda g:bool(ann.get(g)))].reset_index(drop=True)
                bg['target']=bg.gene_id.isin(ids);bg['stratum']=bg[stratum]
                family=nested.term_family(set(bg.gene_id),ann);models=math.setup(bg,comparison) if len(bg) else {}
                frozen.append(bg[['gene_id','target','stratum']].assign(database=database,analysis_id=aid))
                audits.append(dict(database=database,analysis_id=aid,background_before_annotation=before,N=len(bg),
                    target_before_annotation=len(ids),n=int(bg.target.sum()),terms=len(family),
                    status='tested' if family else 'no_terms_in_10_500_background_family'))
                ordered=sorted(family,key=lambda term:(len(family[term]),term))
                selected={ordered[i] for i in np.round(np.linspace(0,len(ordered)-1,min(5,len(ordered)))).astype(int)} if ordered else set()
                index={g:i for i,g in enumerate(bg.gene_id)}
                for model,(codes,N,n) in models.items():
                    mobility.append(dict(database=database,analysis_id=aid,model=model,n_genes=len(bg),target=int(bg.target.sum()),
                        strata=len(N),frozen_gene_fraction=N[(n==0)|(n==N)].sum()/N.sum(),expected_changed_fraction=(2*n*(1-n/N)).sum()/N.sum()))
                for term in ordered:
                    ii=np.array([index[g] for g in family[term]],int);hits=int(bg.target.iloc[ii].sum())
                    tests={};params={};Ntotal=len(bg);ntotal=int(bg.target.sum());Ktotal=len(ii)
                    p_hyper=float(hypergeom.sf(hits-1,Ntotal,Ktotal,ntotal))
                    for model,(codes,N,n) in models.items():
                        K=np.bincount(codes[ii],minlength=len(N));parameters=list(zip(N,K,n));params[model]=parameters
                        dist=math.distribution(parameters);tests[model]=dist;pvalue=math.tail(dist,hits);mean=float((n*K/N).sum())
                        records.append(dict(database=database,analysis_id=aid,species=species,contrast=contrast,
                            comparison=comparison,target=target,term_id=term,term_name=names[database][term]['name'],
                            model=model,N=Ntotal,n=ntotal,K=Ktotal,hits=hits,expected=mean,p=pvalue,p_hyper=p_hyper,
                            fold_enrichment=hits/(ntotal*Ktotal/Ntotal) if ntotal else np.nan,
                            matched_fold=hits/mean if mean else np.nan,excess=hits-mean))
                    if term in selected:
                        for odds in [1,2,4,8]:
                            generated=math.distribution(params['L3_full'],odds)
                            for model,test in tests.items():
                                for label,alpha in [('nominal',.05),('bonferroni_terms',.05/len(family))]:
                                    power.append(dict(database=database,analysis_id=aid,term_id=term,model=model,odds_ratio=odds,
                                        threshold=label,alpha=alpha,rejection_probability=math.rejection_probability(test,generated,alpha)))
                print('functional',database,aid,'terms',len(family),'target',int(bg.target.sum()),flush=True)
        out=pd.DataFrame(records,columns=['database','analysis_id','species','contrast','comparison','target','term_id','term_name',
            'model','N','n','K','hits','expected','p','p_hyper','fold_enrichment','matched_fold','excess'])
        assert sum(r['database']==database for r in audits)==28
        out['q_comparison']=out.groupby(['analysis_id','model']).p.transform(bh)
        out['q_database_model']=out.groupby('model').p.transform(bh);out['q_all_new_positions']=bh(out.p)
        out['q_hyper_family']=out.groupby(['analysis_id','model']).p_hyper.transform(bh)
        out['q_hyper_grid']=out.groupby('model').p_hyper.transform(bh)
        ctx.table(out,f'{database.lower()}_ladder.tsv.gz')
        ctx.table(out[out.model.eq('L3_full')],f'{database.lower()}_exact.tsv.gz')
    ctx.table(membership,'nested_universe_membership.tsv.gz')
    ctx.table(attrition,'nested_universe_attrition.tsv',columns=['species','contrast','gene_id','recurrent','reason'])
    ctx.table(pd.concat(frozen),'functional_frozen_strata.tsv.gz');ctx.table(audits,'functional_families.tsv')
    ctx.table(power,'functional_power.tsv.gz',columns=['database','analysis_id','term_id','model','odds_ratio','threshold','alpha','rejection_probability'])
    ctx.table(mobility,'functional_mobility.tsv',columns=['database','analysis_id','model','n_genes','target','strata','frozen_gene_fraction','expected_changed_fraction'])
    ctx.finish(dict(G0_from_frozen_GEX_counts=True,G1_G2_from_new_PAS=True,all_28_families_per_database=True,
        all_new_term_families=True,exact_stratified_convolution=True,nested_controls_and_exact_power=True),
        pandas=pd.__version__,binning='rank(method=first), pandas qcut; runtime recorded; no inherited old strata')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();functional(Context(a.config,'functional',a.out))
