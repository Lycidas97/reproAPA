"""Preserve separate endpoint and paired-state estimands; audit shared-background context."""
from pathlib import Path
import hashlib,json,itertools
import numpy as np,pandas as pd
from scipy.stats import spearmanr,fisher_exact

P=Path(__file__).resolve().parents[1];V=P/'data/figure4-program-v9'
D=P/'data/figure4-two-contrasts-v12';A=P/'analysis/figure4-two-contrasts-v12';O=P/'outputs/figure4-two-contrasts-v12'
for p in [D,A,O]:p.mkdir(parents=True,exist_ok=True)
inputs=['ortholog_meta_effects.tsv.gz','ortholog_effects.tsv.gz','ortholog_universe.tsv','all_mapped_nodes.tsv.gz','overlap_all_tests.tsv','raw_pm_core_supported50_0.8.tsv','main_native_PAS_source_data.tsv','raw_pm_gene_effects.tsv']
hashes={str(V/x):hashlib.sha256((V/x).read_bytes()).hexdigest() for x in inputs}
meta=pd.read_csv(V/inputs[0],sep='\t');fx=pd.read_csv(V/inputs[1],sep='\t');ort=pd.read_csv(V/inputs[2],sep='\t');nodes=pd.read_csv(V/inputs[3],sep='\t');tests=pd.read_csv(V/inputs[4],sep='\t')
symbols=dict(zip(ort.human_gene_id,ort.human_gene_symbol));hc=['H1','H2','H3','H4'];codes=['W','C']+hc
ps=meta.query("layer=='supported50' and contrast=='SP'").set_index('human_gene_id')[codes]
pm=pd.read_csv(V/inputs[5],sep='\t').set_index('human_gene_id').reindex(columns=codes)
def classify(z,tau):
    mn=z[['W','C']].le(-tau).all(axis=1);mp=z[['W','C']].ge(tau).all(axis=1)
    hn=z[hc].le(-tau).sum(axis=1);hp=z[hc].ge(tau).sum(axis=1)
    return {'both_proximal':mn&(hn>=2)&(hp==0),'both_distal':mp&(hp>=2)&(hn==0),
            'mouse_prox_human_dist':mn&(hp>=2)&(hn==0),'mouse_dist_human_prox':mp&(hn>=2)&(hp==0)}
scatters=[];summaries=[];membership=[];overlaps=[];grids=[]
for con,z in [('SP',ps),('MP',pm)]:
    zz=z.copy();zz['human_symbol']=zz.index.map(symbols);zz['contrast']=con
    zz['mouse_mean']=z[['W','C']].mean(axis=1);zz['human_median']=z[hc].median(axis=1)
    zz['relation']='other';r=classify(z,.02)
    for name,mask in r.items():zz.loc[mask,'relation']=name
    scatters.append(zz.reset_index());corr=spearmanr(zz.mouse_mean,zz.human_median)
    summaries.append(dict(contrast=con,n=len(z),rho=corr.statistic,p=corr.pvalue,**{k:int(v.sum()) for k,v in r.items()}))
    for tau in [.02,.05,.1]:
        for name,mask in classify(z,tau).items():
            for gid in z.index[mask]:membership.append(dict(contrast=con,threshold=tau,relation=name,human_gene_id=gid,human_symbol=symbols[gid]))
    for code in codes:
        for gid,value in z[code].items():grids.append(dict(contrast=con,code=code,human_gene_id=gid,human_symbol=symbols[gid],delta=value))
    if con=='SP':t=tests.query("estimator=='node_means' and layer=='supported50' and contrast=='SP' and threshold==0.02").copy()
    else:t=tests.query("estimator=='raw_paired' and layer=='supported50' and purity==0.8 and contrast=='MP' and threshold==0.02").copy()
    assert len(t)==4
    for row in t.itertuples():
        result=fisher_exact([[row.a,row.b],[row.c,row.d]])
        assert np.isclose(result.pvalue,row.p)
        assert row.N==len(z)
    overlaps.append(t)
pd.concat(scatters).to_csv(D/'panel_A_points.tsv',sep='\t',index=False)
pd.DataFrame(summaries).to_csv(D/'contrast_summary.tsv',sep='\t',index=False)
pd.concat(overlaps).to_csv(D/'panel_B_overlap.tsv',sep='\t',index=False)
mem=pd.DataFrame(membership);mem.to_csv(D/'all_memberships.tsv',sep='\t',index=False)
selected=mem[mem.threshold.eq(.05)].copy();selection=selected.groupby(['human_gene_id','human_symbol']).agg(selected_by=('contrast',lambda x:';'.join(sorted(set(x)))),relations=('relation',lambda x:';'.join(sorted(set(x))))).reset_index()
selection['priority']=selection.selected_by.map({'MP;SP':0,'SP':1,'MP':2});selection=selection.sort_values(['priority','relations','human_symbol']).reset_index(drop=True)
assert len(selection)==18
selection.to_csv(D/'panel_C_selection.tsv',sep='\t',index=False)
# Display all per-study measurable values for selected rows, even if the other contrast
# lacks enough measured studies to enter its primary cross-species background.
context_sp=fx.query("layer=='supported50' and contrast=='SP' and code in @codes and human_gene_id in @selection.human_gene_id")[['human_gene_id','contrast','code','delta']].copy()
context_pm=pd.read_csv(V/'raw_pm_gene_effects.tsv',sep='\t').query("layer=='supported50' and purity==0.8 and quality_pass and code in @codes and human_gene_id in @selection.human_gene_id")[['human_gene_id','code','delta']].copy();context_pm['contrast']='MP'
contexts=pd.concat([context_sp,context_pm],ignore_index=True);contexts['human_symbol']=contexts.human_gene_id.map(symbols)
contexts['in_primary_background']=[gid in (ps.index if con=='SP' else pm.index) for gid,con in zip(contexts.human_gene_id,contexts.contrast)]
assert not contexts.duplicated(['human_gene_id','contrast','code']).any()
contexts.to_csv(D/'panel_C_effects.tsv',sep='\t',index=False)
# Keep the exact same gene and human-study coverage in both contrasts for this diagnostic.
f=fx.query("layer=='supported50' and contrast in ['SP','MP']").pivot(index='human_gene_id',columns=['code','contrast'],values='delta')
meas={h:f[[(h,'SP'),(h,'MP')]].notna().all(axis=1) for h in hc}
ok=f[[('W','SP'),('W','MP'),('C','SP'),('C','MP')]].notna().all(axis=1)&pd.DataFrame(meas).sum(axis=1).ge(2)
fixed=[];fixed_summary=[]
for con in ['SP','MP']:
    z=f.loc[ok,[(c,con) for c in codes]].copy();z.columns=codes
    for h in hc:z.loc[~meas[h][ok],h]=np.nan
    z['contrast']=con;z['mouse_mean']=z[['W','C']].mean(axis=1);z['human_median']=z[hc].median(axis=1)
    fixed.append(z.reset_index());r=spearmanr(z.mouse_mean,z.human_median)
    fixed_summary.append(dict(contrast=con,n=len(z),rho=r.statistic,estimator='same library-balanced node means for both contrasts'))
assert all(x['n']==856 for x in fixed_summary)
pd.concat(fixed).to_csv(D/'fixed_background_diagnostic.tsv',sep='\t',index=False)
pd.DataFrame(fixed_summary).to_csv(D/'fixed_background_summary.tsv',sep='\t',index=False)
# Case extraction uses the whole native site distribution, including observed zeros.
case=ort[ort.human_gene_symbol.isin(['PDLIM7','OTUB1'])]
case.to_csv(D/'example_pairs.tsv',sep='\t',index=False)
nodes[nodes.human_gene_id.isin(case.human_gene_id)].to_csv(D/'case_native_nodes.tsv.gz',sep='\t',index=False)
pd.read_csv(V/'main_native_PAS_source_data.tsv',sep='\t').to_csv(D/'paired_PM_PAS.tsv',sep='\t',index=False)
(A/'input_manifest.json').write_text(json.dumps(dict(source_hashes=hashes,primary_contrasts=['SP','MP'],comparison_specific_backgrounds={'SP':1174,'MP':386},selection='Union of replicated members at absolute effect >=0.05 in either comparison; no P-value gate.',count_tests='Previously defined 288-test BH family, unchanged.',case_role='Main: OTUB1 P-S and CTDSP2 paired P-M. PDLIM7/OTUB1 across-estimator curves are diagnostic; differences do not isolate the reference-state effect.'),indent=2)+'\n')
print(pd.DataFrame(summaries).to_string(index=False));print(selection.to_string(index=False));print(pd.DataFrame(fixed_summary).to_string(index=False))
