"""Prepare Figure 4 tables with unique-ortholog uncertainty for prediction."""
from pathlib import Path
import json, hashlib
import pandas as pd
import numpy as np

P=Path(__file__).resolve().parents[1]
Q=P/'analysis/figure4-supplements-v54';O=P/'outputs/figure4-supplements-v54'
Q.mkdir(exist_ok=True);O.mkdir(exist_ok=True)
sources={};blocks=[]
names={'W':'Schiebinger2019','C':'Guo2019','H1':'Xing2020','H2':'Liu2020','H3':'Guan2022','H4':'Panariello2023'}

def read(folder,name):
    filename=name if name.endswith(('.tsv','.tsv.gz')) else name+'.tsv'
    p=P/'data'/folder/filename
    sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return pd.read_csv(p,sep='\t')

def clean(df):
    df=df.copy()
    for c in ['study','code','omitted']:
        if c in df: df[c]=df[c].replace(names)
    if 'pair' in df:
        df['pair']=df['pair'].map(lambda x:' / '.join(names.get(k,k) for k in str(x).split('-')))
    for c in ['discovery','direction']:
        if c in df:
            df[c]=df[c].replace({'human':'Human to Mouse','mouse':'Mouse to Human','human_to_mouse':'Human to Mouse','mouse_to_human':'Mouse to Human'})
    return df.astype(object).where(pd.notnull(df),None)

def add(sheet,title,df,note,source):
    df=clean(df)
    blocks.append(dict(sheet=sheet,title=title,columns=list(df.columns),rows=df.values.tolist(),note=note,source=source))

# Table S4a: exact values behind every Figure 4A study-pair circle and its state-level test.
a=read('figure4-state-recurrence-v25','state_specific_pair_metrics')
a=a[a.state.isin(['R','E','P'])][['layer','state','pair','n','spearman','pearson','n_both_strong','direction_concordance','direction_score_all']]
add('Table S4','S4a-1. State-specific study-pair correspondence',a,
    'Exact descriptive values behind Figure 4A. n is the jointly measurable one-to-one ortholog count; n_both_strong requires |ΔWPP| ≥ 0.02 in both species. Pair-level P/q values were not computed and are not imputed.',
    'v25 state_specific_pair_metrics')
a=read('figure4-state-recurrence-v25','state_specific_permutation')
a=a[a.state.isin(['R','E','P'])][['layer','state','metric','observed','null_mean','null_lo','null_hi','p_greater','q_within_layer_12']]
a=a.rename(columns={'p_greater':'P_one_sided','q_within_layer_12':'q_12'})
add('Table S4','S4a-2. State-specific aggregate permutation tests',a,
    'Equal weights are assigned to available study pairs. Null bounds are permutation percentiles. q_12 adjusts three states by four metrics within each PAS layer.',
    'v25 state_specific_permutation')

# Table S4b: matched comparisons directly supporting the state-convergence boundary.
a=read('figure4-state-recurrence-v25','matched_state_permutation')
b=read('figure4-m-state-v38','matched_P_vs_M_permutation').rename(columns={'M_value':'state_value','observed_difference_P_minus_M':'observed_difference'})
cols=['layer','comparison','metric','state_value','P_value','observed_difference','null_lo','null_hi','p_P_stronger']
a=pd.concat([a[cols],b[cols]],ignore_index=True)
c=read('figure4-m-state-v38','matched_P_vs_R_E_M_combined_correction')
a=a.merge(c[['layer','comparison','metric','q_within_layer_12']],on=['layer','comparison','metric'],validate='one_to_one')
a=a.rename(columns={'state_value':'reference_state_metric','P_value':'P_state_metric','observed_difference':'P_minus_reference_difference','p_P_stronger':'P_one_sided','q_within_layer_12':'q_12'})
a['comparison']=a.comparison.replace({'P_vs_R':'PS_vs_RS','P_vs_E':'PS_vs_ES','P_vs_M':'PS_vs_MS'})
add('Table S4','S4b. Matched state comparisons',a,
    'P−S and the reference state−S use the same genes and study pairs. PS_vs_MS compares P−S correspondence with M−S correspondence and is distinct from the paired P−M endpoint. Null bounds are permutation percentiles; q_12 adjusts three comparisons by four metrics within each PAS layer.',
    'v25 matched_state_permutation; v38 matched_P_vs_M_permutation and combined correction')

# Table S4c: every plotted shared-event sensitivity setting.
a=read('figure4-supplements-v49','overlap_threshold_layer')
a=a[['layer','tau','contrast','observed','expected','OE','null_lo','null_hi','p','q_two','q_16']]
add('Table S4','S4c. Shared-event sensitivity',a,
    'q_16 adjusts all 16 contrast-by-layer-by-threshold tests and is used for filled symbols in Figure S4A. q_two adjusts P−S and P−M within each setting.',
    'v49 overlap_threshold_layer')

# Table S4d-1: same point estimates, corrected unique-gene uncertainty.
audit_path=P/'analysis/figure4c-reused-controls-v53/thirty_setting_unique_gene_audit.tsv'
sources[str(audit_path)]=hashlib.sha256(audit_path.read_bytes()).hexdigest()
a=pd.read_csv(audit_path,sep='\t')
keep=['layer','contrast','discovery','tau','n_sets','n_unique_controls',
      'n_reused_controls','control_effective_n_by_weights',
      'observed_binary_or','stratified_gene_bootstrap_binary_or_lo','stratified_gene_bootstrap_binary_or_hi',
      'observed_strict_or','stratified_gene_bootstrap_strict_or_lo','stratified_gene_bootstrap_strict_or_hi',
      'observed_continuous_difference','stratified_gene_bootstrap_continuous_lo',
      'stratified_gene_bootstrap_continuous_hi','exact_permutation_p_two_sided','exact_q_four',
      'exact_q_thirty','effect_bin_permutation_p_two_sided','effect_bin_q_thirty']
a=a[keep].sort_values(['layer','contrast','discovery','tau']).reset_index(drop=True)
add('Table S4','S4d-1. Reciprocal prediction sensitivity',a,
    'Fixed-match 95% intervals resample unique orthologs within exact discovery-direction and study-count strata. Reused controls retain their total weight. Continuous P values permute each unique ortholog outcome once; exact_q_four applies to four primary supported50/0.02 settings, exact_q_thirty to all 30 settings. Effect-bin P/q are a stricter exchangeability sensitivity. OR intervals have no P/q. Primary settings use 20,000 draws; other settings use 5,000.',
    'v53 unique-gene audit of v26 frozen matched_comparators')

# Table S4d-2: concise balance audit at the primary threshold.
a=read('figure4-recurrence-validation-v26','matching_balance')
a=a[a.tau.eq(.02)].copy()
rows=[]
for keys,z in a.groupby(['layer','contrast','discovery'],sort=False):
    z=z.assign(abs_smd=z.standardized_difference.abs()).sort_values('abs_smd',ascending=False)
    r=z.iloc[0]
    rows.append(dict(layer=keys[0],contrast=keys[1],discovery=keys[2],tau=.02,
                     max_abs_SMD=r.abs_smd,covariate_at_max=r.covariate,
                     all_abs_SMD_below_0_10=bool((z.abs_smd<.10).all())))
a=pd.DataFrame(rows)
add('Table S4','S4d-2. Primary matching-balance summary',a,
    'Maximum absolute standardized difference over effect magnitude, PAS number, support, libraries and span. Full covariate-level balance remains in the analysis archive.',
    'v26 matching_balance')

# Table S4e: R/E-only technical-confound audit requested for the intermediate-state claim.
obs=read('figure4-state-recurrence-v25','technical_observation_audit.tsv.gz')
obs=obs[obs.layer.eq('supported50')&obs.state.isin(['R','E'])&obs.relation.isin(['concordant','discordant'])].copy()
summary=read('figure4-supplements-v49','technical_RE_only')
rows=[]
for cov in ['min_libraries','min_metacells','min_pas','min_span','max_library_sd','all_nodes_multi_library','relation_stable']:
    z=obs.dropna(subset=[cov]).copy();disc=pd.to_numeric(z.loc[z.relation.eq('discordant'),cov]);conc=pd.to_numeric(z.loc[z.relation.eq('concordant'),cov])
    r=summary[summary.covariate.eq(cov)].iloc[0]
    rows.append(dict(technical_variable=cov,n=int(r.n),n_discordant=int(r.n_discordant),n_concordant=int(r.n_concordant),
                     discordant_mean=disc.mean(),discordant_median=disc.median(),concordant_mean=conc.mean(),concordant_median=conc.median(),
                     raw_mean_difference=disc.mean()-conc.mean(),test_scale=r.kind,test_scale_difference=r.effect,
                     null_lo=r.null_lo,null_hi=r.null_hi,P_two_sided=r.p,q=r.q))
a=pd.DataFrame(rows)
add('Table S4','S4e. Technical audit of intermediate-state concordance',a,
    'Supported50 R/E records with strong cross-species effects only (755 concordant; 777 discordant). Raw group summaries are shown. Permutation tests compare mean log1p values for continuous variables and fractions for binary variables, shuffling concordance labels within state/study-pair strata. q is BH-adjusted separately for five continuous and two binary variables. Repeated genes/studies make this a diagnostic association, not a causal exclusion.',
    'v25 technical_observation_audit; v49 technical_RE_only')

# Table S5: 26 rows rather than 650 long-format records.
cells=read('figure4-integrated-redesign-v29','D_all_plotted_cells')
members=read('figure4-integrated-redesign-v29','D_membership')
assert len(cells)==650 and cells[['gene','study','state','estimator']].duplicated().sum()==0
order={'W':0,'C':1,'H1':2,'H2':3,'H3':4,'H4':5}
cells['study_order']=cells.study.map(order)
cells['key']=cells.apply(lambda r:f"{names[r.study]}_{r.state.replace('−','-')}_{'paired' if r.estimator=='same_library_paired' else 'node'}",axis=1)
col_order=(cells[['key','study_order','state','estimator']].drop_duplicates()
           .assign(state_order=lambda d:d.state.map({'R':0,'E':1,'P':2,'M':3,'P−M':4}))
           .sort_values(['study_order','state_order','estimator']).key.tolist())
wide=cells.pivot(index='gene',columns='key',values='effect').reindex(columns=col_order)
support=(cells[cells.endpoint_support].groupby('gene').apply(lambda z:'; '.join(z.sort_values(['study_order','state']).key),include_groups=False).rename('endpoint_support_records'))
mismatch=(cells[cells.estimator_sign_mismatch].groupby('gene').apply(lambda z:'; '.join(z.sort_values(['study_order','state']).key),include_groups=False).rename('estimator_sign_mismatch_records'))
meta=members.rename(columns={'human_symbol':'gene'}).set_index('gene')[['human_gene_id','in_PS_shared_prox','in_PM_shared_prox','set_relation']]
a=meta.join(wide).join(support).join(mismatch).reset_index()
a['endpoint_support_records']=a.endpoint_support_records.fillna('')
a['estimator_sign_mismatch_records']=a.estimator_sign_mismatch_records.fillna('')
assert len(a)==26 and a[col_order].notna().sum().sum()==cells.effect.notna().sum()
add('Table S5','S5. Gene-level effect matrix',a,
    'Each row is one selected one-to-one ortholog. X_node columns are node-minus-S ΔWPP; X_P-M_paired columns are same-library paired effects. Blank means unavailable, not zero. Record-list columns identify exact supporting or estimator-mismatch cells.',
    'v29 D_membership and D_all_plotted_cells')

(Q/'tables.json').write_text(json.dumps(blocks,ensure_ascii=False,allow_nan=False)+'\n')
(Q/'source_hashes.json').write_text(json.dumps(sources,indent=2)+'\n')
(Q/'removed_from_submission.txt').write_text(
    'Retained in analysis archive, not in manuscript-facing tables:\n'
    '- raw ortholog-permutation draws\n'
    '- recurrence-score groups, trend and interaction tests\n'
    '- leave-one-study detail\n'
    '- covariate-level matching balance beyond the primary max-|SMD| summary\n'
    '- pooled R/E/P technical-association audit beyond the R/E-only manuscript table\n'
    '- library-level, threshold-classification and reference-decomposition audit rows\n')
print(json.dumps({'tables':2,'sections':len(blocks),'rows':sum(len(b['rows']) for b in blocks),'gene_rows':len(a)}))
