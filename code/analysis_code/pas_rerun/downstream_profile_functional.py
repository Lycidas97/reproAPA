"""Full-SREP amplitudes and functional families for the new trajectory cohort."""
from collections import defaultdict

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import hypergeom

from downstream_functional import annotations
from downstream_runtime import bh

SOURCE='github/reproapa-analysis/modules/07_apa_dynamics/code/figure2/analyze_figure2_cde_dynamics_v1.py'
FIT={'mouse':['WOT','CIP'],'human':['GSE118258']}


def rank_effect(x,mask,pvalue=True):
    x=np.asarray(x,float);mask=np.asarray(mask,bool);n=int(mask.sum());m=len(mask)-n
    assert n and m and np.isfinite(x).all()
    u=float(stats.rankdata(x)[mask].sum()-n*(n+1)/2)
    p=(1. if np.ptp(x)==0 else float(stats.mannwhitneyu(x[mask],x[~mask],alternative='two-sided',method='asymptotic').pvalue)) if pvalue else np.nan
    return 2*u/(n*m)-1,p


def term_sets(annotation,genes):
    out=defaultdict(set)
    for gene in genes:
        for term in annotation.get(gene,set()):out[term].add(gene)
    return {t:g for t,g in out.items() if 10<=len(g)<=500}


def choose_themes(ora):
    sets={(r.species,r.shared_term):set(r.background_genes.split(';')) for r in ora.drop_duplicates(['species','shared_term']).itertuples()}
    candidates=ora[ora.hits.ge(5)&ora.fold_enrichment.ge(1.5)].sort_values(['p_value','term_id'],kind='stable').drop_duplicates('shared_term')
    chosen=[];audit=[]
    for row in candidates.itertuples():
        reason='';related=[]
        if not all((sp,row.shared_term) in sets for sp in FIT):reason='not eligible in both species'
        if not reason:
            for prior in chosen:
                for sp in FIT:
                    a,b=sets[sp,row.shared_term],sets[sp,prior['shared_term']]
                    if len(a&b)/len(a|b)>.5 or len(a&b)/min(len(a),len(b))>.8:related.append(prior['shared_term'])
            if related:reason='redundant membership'
        if not reason and len(chosen)>=6:reason='below six-theme display cutoff'
        if not reason:chosen.append(dict(shared_term=row.shared_term,database=row.database,term_name=row.term_name,
            order=len(chosen),label=row.term_name,selection_species=row.species,selection_cluster=row.cluster,
            selection_p=row.p_value,selection_q=row.q_value))
        audit.append(dict(shared_term=row.shared_term,term_name=row.term_name,selected=not reason,reason=reason,
            related=';'.join(sorted(set(related))),selection_p=row.p_value))
    return pd.DataFrame(chosen,columns=['shared_term','database','term_name','order','label','selection_species','selection_cluster','selection_p','selection_q']),pd.DataFrame(audit,columns=['shared_term','term_name','selected','reason','related','selection_p'])


def functional_ranges(ctx,nodes,members):
    math=ctx.functions(SOURCE,['amplitudes','matched_effect'],dict(stats=stats,rank_effect=rank_effect))
    annotation,names,_=annotations(ctx)
    n=nodes[nodes.depth.eq(20)].copy();n['layer']=n.setting.map({'supported50_p08':'supported50','expanded_primary_p08':'expanded_primary'})
    amp=math.amplitudes(n)
    amp=amp.merge(members,on=['species','gene_id'],how='left',validate='many_to_one');amp['selected']=amp.cluster.notna()
    # Structural full-SREP coverage is established before range/correlation selection.
    backgrounds={}
    for sp,routes in FIT.items():
        f=amp[amp.species.eq(sp)&amp.layer.eq('supported50')&amp.route.isin(routes)]
        count=f.groupby('gene_id').route.nunique();backgrounds[sp]=set(count[count.eq(len(routes))].index)
        assert set(members[members.species.eq(sp)].gene_id)<=backgrounds[sp]
    amp=amp[amp.apply(lambda r:r.gene_id in backgrounds[r.species],axis=1)]
    ctx.table(amp,'background_amplitudes.tsv.gz');selected=amp[amp.selected]
    ctx.table(selected,'panel_C_points.tsv.gz')
    summary=[]
    for key,d in selected.groupby(['layer','species','route']):
        summary.append(dict(zip(['layer','species','route'],key))|dict(n=len(d),median_apa_range=d.apa_range.median(),
            median_gex_range=d.gex_range.median(),spearman_r=stats.spearmanr(d.apa_range,d.gex_range).statistic if d.apa_range.nunique()>1 and d.gex_range.nunique()>1 else np.nan))
    ctx.table(summary,'panel_C_summary.tsv')
    coverage=[]
    for key,d in n.groupby(['species','layer','route']):
        sp,layer,route=key;nc=int((selected.species.eq(sp)&selected.layer.eq(layer)&selected.route.eq(route)).sum())
        coverage.append(dict(species=sp,layer=layer,route=route,observed_nodes=''.join(s for s in 'SREPM' if s in set(d.node)),
            complete_SREP_selected_genes=nc,status='included' if nc else 'no_complete_matched_SREP'))
    ctx.table(coverage,'study_coverage.tsv')
    records=[];audit=[];mapping=[]
    for sp in FIT:
        for database in ['GO','KEGG']:
            ann=annotation[database][sp];bg={g for g in backgrounds[sp] if ann.get(g)};family=term_sets(ann,bg)
            for gene in sorted(backgrounds[sp]):mapping.append(dict(species=sp,database=database,gene_id=gene,annotated=gene in bg))
            for cluster in ['C1','C2','C3']:
                full=set(members[members.species.eq(sp)&members.cluster.eq(cluster)].gene_id);fg=full&bg
                audit.append(dict(species=sp,database=database,cluster=cluster,evaluable_background=len(backgrounds[sp]),
                    annotated_background=len(bg),foreground=len(full),annotated_foreground=len(fg),terms_tested=len(family)))
                for term,genes in sorted(family.items()):
                    hits=fg&genes;expected=len(fg)*len(genes)/len(bg)
                    records.append(dict(species=sp,database=database,cluster=cluster,term_id=term,
                        shared_term=term if database=='GO' else term[-5:],term_name=names[database][term]['name'],hits=len(hits),
                        foreground_size=len(fg),background_hits=len(genes),background_size=len(bg),expected_hits=expected,
                        fold_enrichment=len(hits)/expected if expected else np.nan,p_value=float(hypergeom.sf(len(hits)-1,len(bg),len(genes),len(fg))),
                        hit_genes=';'.join(sorted(hits)),background_genes=';'.join(sorted(genes))))
    ora=pd.DataFrame(records,columns=['species','database','cluster','term_id','shared_term','term_name','hits','foreground_size',
        'background_hits','background_size','expected_hits','fold_enrichment','p_value','hit_genes','background_genes'])
    ora['q_value']=ora.groupby(['species','database','cluster']).p_value.transform(bh)
    ora['q_all_clusters']=ora.groupby(['species','database']).p_value.transform(bh)
    ctx.table(ora,'enrichment_all_tests.tsv.gz');ctx.table(audit,'enrichment_background_audit.tsv');ctx.table(mapping,'gene_annotation_coverage.tsv')
    themes,audit=choose_themes(ora);ctx.table(themes,'themes.tsv');ctx.table(audit,'theme_selection_audit.tsv')
    ctx.table(ora.merge(themes[['shared_term','order','label']],on='shared_term',validate='many_to_one'),'panel_D_values.tsv')
    rows=[];audit=[]
    for (sp,layer,route),raw in amp.groupby(['species','layer','route']):
        for database in ['GO','KEGG']:
            ann=annotation[database][sp];d=raw[raw.gene_id.map(lambda g:bool(ann.get(g)))].sort_values('gene_id').reset_index(drop=True)
            bins=[pd.qcut(d[col],3,labels=False,duplicates='drop').fillna(0).astype(int).to_numpy() if len(d) else np.array([],int)
                for col in ['baseline_gex','min_pas_umi','min_qualified_metacells']]
            strata=bins[0]*9+bins[1]*3+bins[2]
            audit.append(dict(species=sp,layer=layer,route=route,database=database,complete_matched_genes=len(raw),
                annotated_matched_genes=len(d),selected_genes=int(d.selected.sum()),strata=len(np.unique(strata)),
                APA_flat_within_all_libraries=int(d.max_within_library_apa_range.le(1e-12).sum())))
            for term,genes in sorted(term_sets(ann,set(d.gene_id)).items()):
                mask=d.gene_id.isin(genes).to_numpy()
                if (~mask).sum()<10:continue
                for measure,col in [('APA','apa_range'),('GEX','gex_range')]:
                    x=d[col].to_numpy();effect,p=rank_effect(x,mask);adjusted,retained=math.matched_effect(x,mask,strata)
                    variable=d.max_within_library_apa_range.gt(1e-12).to_numpy();vmask=mask[variable]
                    nonflat=rank_effect(x[variable],vmask)[0] if vmask.sum()>=5 and (~vmask).sum()>=10 else np.nan
                    rows.append(dict(species=sp,layer=layer,route=route,database=database,term_id=term,shared_term=term if database=='GO' else term[-5:],
                        term_name=names[database][term]['name'],measure=measure,n_set=int(mask.sum()),n_background=len(d),rank_effect=effect,p_value=p,
                        matched_rank_effect=adjusted,matched_retained_fraction=retained,nonflat_APA_background_effect=nonflat,member_genes=';'.join(sorted(genes))))
    effects=pd.DataFrame(rows,columns=['species','layer','route','database','term_id','shared_term','term_name','measure','n_set','n_background',
        'rank_effect','p_value','matched_rank_effect','matched_retained_fraction','nonflat_APA_background_effect','member_genes'])
    effects['q_value']=effects.groupby(['layer','route','database','measure']).p_value.transform(bh)
    effects['q_both_databases']=effects.groupby(['layer','route','measure']).p_value.transform(bh)
    ctx.table(effects,'rank_enrichment_all_tests.tsv.gz');ctx.table(audit,'rank_background_audit.tsv')
    ctx.table(effects.merge(themes[['shared_term','order','label']],on='shared_term',validate='many_to_one'),'panel_E_values.tsv')
