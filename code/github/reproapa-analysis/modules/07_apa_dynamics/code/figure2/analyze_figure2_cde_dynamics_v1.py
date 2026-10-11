"""Full-SREP amplitudes, nonredundant ORA themes, and paired-layer rank effects."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy import stats
from enrich_figure2_cde_v1 import bh

P=Path(__file__).resolve().parents[1];OUT=P/'outputs/figure2-cde-dynamics-v1'
DATA=P/'data/figure2-dynamic-background-v1'
LAYERS=['supported50','expanded_primary'];ROUTES=['WOT','CIP','GSE118258']

def save(d,name):d.to_csv(OUT/name,sep='\t',index=False,compression='gzip' if name.endswith('.gz') else None)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def amplitudes(nodes):
    rows=[]
    for (sp,layer,route,gene),d in nodes[nodes.node.isin(list('SREP'))].groupby(['species','layer','route','gene_id']):
        assert d.node.is_unique
        d=d.set_index('node').reindex(list('SREP'))
        if d[['position_index','gex_log2_cpm1']].isna().any().any():continue
        # Equal effects must remain ties after TSV round trips; 1e-12 is far below
        # biological resolution and removes machine-precision pseudo-ordering.
        row=dict(species=sp,layer=layer,route=route,gene_id=gene,apa_range=float(np.round(np.ptp(d.position_index),12)),gex_range=float(np.round(np.ptp(d.gex_log2_cpm1),12)),baseline_gex=d.loc['S','gex_log2_cpm1'])
        if 'mean_pas_umi' in d:
            row.update(min_pas_umi=d.mean_pas_umi.min(),min_qualified_metacells=d.qualified_metacells.min(),max_within_library_apa_range=d.within_library_apa_range.max())
        rows.append(row)
    return pd.DataFrame(rows)

def themes(ora):
    sets={(r.species,r.shared_term):set(r.background_genes.split(';')) for r in ora.drop_duplicates(['species','shared_term']).itertuples()}
    candidates=ora[ora.hits.ge(5)&ora.fold_enrichment.ge(1.5)].sort_values(['p_value','term_id'],kind='stable').drop_duplicates('shared_term')
    chosen=[];audit=[]
    for row in candidates.itertuples():
        reason='';related=[]
        if not all((sp,row.shared_term) in sets for sp in ['mouse','human']):reason='no eligible term in both species'
        if not reason:
            for prior in chosen:
                for sp in ['mouse','human']:
                    a,b=sets[sp,row.shared_term],sets[sp,prior['shared_term']]
                    if len(a&b)/len(a|b)>.5 or len(a&b)/min(len(a),len(b))>.8:related.append(prior['shared_term'])
            if related:reason='redundant membership'
        if not reason and len(chosen)>=6:reason='below six-theme display cutoff'
        selected=not reason
        if selected:
            chosen.append(dict(shared_term=row.shared_term,database=row.database,term_name=row.term_name,order=len(chosen),selection_species=row.species,selection_cluster=row.cluster,selection_p=row.p_value,selection_q=row.q_value))
        audit.append(dict(shared_term=row.shared_term,term_name=row.term_name,selected=selected,reason=reason,related=';'.join(sorted(set(related))),selection_p=row.p_value))
    labels={'GO:0001501':'Skeletal development','GO:0022613':'RNP complex\nbiogenesis','GO:0040013':'Negative regulation\nof locomotion','GO:0030518':'Steroid receptor\nsignaling','GO:0051170':'Nuclear import','04068':'FoxO signaling'}
    selected=pd.DataFrame(chosen);selected['label']=selected.shared_term.map(labels);assert selected.label.notna().all()
    save(selected,'themes.tsv');save(pd.DataFrame(audit),'theme_selection_audit.tsv')
    d=ora.merge(selected[['shared_term','order','label']],on='shared_term',validate='many_to_one')
    d['log2_fold']=np.log2(d.fold_enrichment.replace(0,np.nan));d['display']=d.hits.ge(3)&d.fold_enrichment.gt(1)
    save(d,'panel_D_values.tsv');return selected

def rank_effect(x,mask,pvalue=True):
    n=int(mask.sum());m=len(mask)-n
    u=float(stats.rankdata(x)[mask].sum()-n*(n+1)/2)
    return 2*u/(n*m)-1,float(stats.mannwhitneyu(x[mask],x[~mask],alternative='two-sided',method='asymptotic').pvalue) if pvalue else np.nan

def matched_effect(x,mask,strata):
    values=[];weights=[]
    for k in np.unique(strata):
        use=strata==k;a=mask[use]
        if a.sum() and (~a).sum()>=5:
            effect,_=rank_effect(x[use],a,pvalue=False);values.append(effect);weights.append(int(a.sum()))
    return (float(np.average(values,weights=weights)),sum(weights)/mask.sum()) if weights else (np.nan,0.)

def main():
    OUT.mkdir(exist_ok=True)
    contract=dict(nodes='SREP exactly; no imputation; M excluded',C='fixed selected genes; gene by study range=max-min; all observations retained',D='complete existing ORA; hits>=5 and fold>=1.5 candidates; rank by minimum nominal p across species/groups; reject Jaccard>0.5 or overlap/min_size>0.8 in either species; six representatives; not confirmed pathways',E='all GO BP/KEGG sets with 10–500 measured genes, database-annotated matched background; 2*AUC-1, ties=0.5; no effect/significance selection',sensitivity='within baseline S GEX tertile × minimum node PAS UMI tertile × minimum node qualified MC tertile; average stratum effects weighted by function gene counts; require >=5 controls per stratum',BH='per layer/route/database/measure across all tested sets; additional combined-database correction',selection='E themes selected using D, not an independent validation; expanded layer retains same groups and themes')
    (OUT/'analysis_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    contract['rank_precision']='Amplitudes rounded to 12 decimal places before ranking, so machine precision residuals do not split ties.'
    (OUT/'analysis_contract.json').write_text(json.dumps(contract,indent=2)+'\n')
    selected=pd.read_csv(P/'data/figure2-cde-v1/matched_node_values.tsv.gz',sep='\t')
    c=amplitudes(selected);membership=selected[['species','gene_id','cluster','gene_symbol']].drop_duplicates()
    c=c.merge(membership,on=['species','gene_id'],validate='many_to_one');save(c,'panel_C_points.tsv.gz')
    style=json.loads((P/'sources/local/apa-publication-layout-style-v1.json').read_text());coverage=[]
    for layer in LAYERS:
        for r in style['routes']:
            d=selected[selected.layer.eq(layer)&selected.route.eq(r['route'])]
            n=int((c.layer.eq(layer)&c.route.eq(r['route'])).sum())
            coverage.append(dict(layer=layer,species=r['species'],route=r['route'],label=r['label'],observed_nodes=''.join(s for s in 'SREPM' if s in set(d.node)),complete_SREP_selected_genes=n,status='Included in C/E' if n else 'A/B only: no complete matched SREP'))
    save(pd.DataFrame(coverage),'study_coverage.tsv')
    assert set(c.route)==set(ROUTES)
    summaries=[]
    for (layer,sp,route),d in c.groupby(['layer','species','route']):
        summaries.append(dict(layer=layer,species=sp,route=route,n=len(d),median_apa_range=d.apa_range.median(),median_gex_range=d.gex_range.median(),spearman_r=stats.spearmanr(d.apa_range,d.gex_range).statistic))
    save(pd.DataFrame(summaries),'panel_C_summary.tsv')
    ora=pd.read_csv(P/'outputs/figure2-analysis-v1/enrichment_all_tests.tsv.gz',sep='\t')
    ora['shared_term']=ora.term_id.where(ora.database.eq('GO_BP'),ora.term_id.str[-5:]);chosen=themes(ora)
    nd=pd.read_csv(DATA/'matched_background_nodes.tsv.gz',sep='\t');amp=amplitudes(nd)
    bg=pd.concat([pd.read_csv(P/'data/figure2-cde-v1'/f'{sp}_evaluable_background.tsv',sep='\t') for sp in ['mouse','human']])
    amp=amp.merge(bg[['species','gene_id','selected','cluster']],on=['species','gene_id'],validate='many_to_one')
    cmp=c.merge(amp,on=['species','layer','route','gene_id'],suffixes=('_c','_bg'),validate='one_to_one')
    assert len(cmp)==len(c)
    for col in ['apa_range','gex_range']:assert np.allclose(cmp[col+'_c'],cmp[col+'_bg'],atol=1e-12,rtol=0)
    save(amp,'background_amplitudes.tsv.gz')
    annot=pd.read_csv(P/'outputs/figure2-analysis-v1/gene_annotation_coverage.tsv',sep='\t')
    records=[];bgaudit=[]
    for (sp,layer,route),raw in amp.groupby(['species','layer','route']):
        for database in ['GO_BP','KEGG']:
            ann=set(annot.loc[annot.species.eq(sp)&annot.database.eq(database)&annot.annotated,'gene_id'])
            d=raw[raw.gene_id.isin(ann)].sort_values('gene_id').reset_index(drop=True);ids=set(d.gene_id)
            bins=[]
            for col in ['baseline_gex','min_pas_umi','min_qualified_metacells']:
                bins.append(pd.qcut(d[col],3,labels=False,duplicates='drop').fillna(0).astype(int).to_numpy())
            strata=bins[0]*9+bins[1]*3+bins[2]
            bgaudit.append(dict(species=sp,layer=layer,route=route,database=database,complete_matched_genes=len(raw),annotated_matched_genes=len(d),selected_genes=int(d.selected.sum()),strata=len(np.unique(strata)),APA_flat_within_all_libraries=int(d.max_within_library_apa_range.le(1e-12).sum())))
            terms=ora[ora.species.eq(sp)&ora.database.eq(database)].drop_duplicates('term_id')
            for term in terms.itertuples():
                members=set(term.background_genes.split(';'))&ids;mask=d.gene_id.isin(members).to_numpy()
                if not 10<=mask.sum()<=500 or (~mask).sum()<10:continue
                for measure,col in [('APA','apa_range'),('GEX','gex_range')]:
                    x=d[col].to_numpy();effect,p=rank_effect(x,mask);adjusted,retained=matched_effect(x,mask,strata)
                    variable=d.max_within_library_apa_range.gt(1e-12).to_numpy();vmask=mask[variable]
                    nonflat=rank_effect(x[variable],vmask)[0] if vmask.sum()>=5 and (~vmask).sum()>=10 else np.nan
                    records.append(dict(species=sp,layer=layer,route=route,database=database,term_id=term.term_id,shared_term=term.shared_term,term_name=term.term_name,measure=measure,n_set=int(mask.sum()),n_background=len(d),rank_effect=effect,p_value=p,matched_rank_effect=adjusted,matched_retained_fraction=retained,nonflat_APA_background_effect=nonflat,member_genes=';'.join(sorted(members))))
        print('ranked',sp,layer,route,flush=True)
    effects=pd.DataFrame(records)
    effects['q_value']=effects.groupby(['layer','route','database','measure']).p_value.transform(bh)
    effects['q_both_databases']=effects.groupby(['layer','route','measure']).p_value.transform(bh)
    save(effects,'rank_enrichment_all_tests.tsv.gz');save(pd.DataFrame(bgaudit),'rank_background_audit.tsv')
    e=effects.merge(chosen[['shared_term','order','label']],on='shared_term',validate='many_to_one');save(e,'panel_E_values.tsv')
    inputs=[P/'data/figure2-cde-v1/matched_node_values.tsv.gz',DATA/'matched_background_nodes.tsv.gz',DATA/'manifest.json',P/'outputs/figure2-analysis-v1/enrichment_all_tests.tsv.gz',P/'outputs/figure2-analysis-v1/gene_annotation_coverage.tsv']
    (OUT/'analysis_manifest.json').write_text(json.dumps(dict(status='PASS',source_sha256={str(p):sha(p) for p in inputs},script_sha256=sha(Path(__file__)),C_points=len(c),rank_tests=len(effects),themes=len(chosen),old_selected_amplitudes_reproduced=len(cmp),contract=contract),indent=2)+'\n')
    print(pd.DataFrame(summaries).to_string(index=False));print(e[e.layer.eq('supported50')][['label','route','measure','rank_effect','matched_rank_effect','q_value']].to_string(index=False))

if __name__=='__main__':main()
