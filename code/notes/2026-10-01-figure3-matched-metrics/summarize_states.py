"""Describe matched state relations without treating dependent contrasts as tests."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
KEY=['species','study_id','gene_id']
REQUIRED={'human':['GSE118258'],'mouse':['GSE103221','GSE115943']}
STUDIES={'human':['GSE118258','GSE150311','GSE221739','GSE178325'],'mouse':['GSE103221','GSE115943']}
PATTERNS=['dual','PM_only','PS_only','neither']


def save(frame,name):
    frame.to_csv(HERE/name,sep='\t',index=False,na_rep='NA',float_format='%.17g',
                 compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)


def pattern(ps,pm):
    return 'dual' if ps and pm else 'PM_only' if pm else 'PS_only' if ps else 'neither'


def main():
    m=pd.read_csv(HERE/'recurrence_member_inventory.tsv',sep='\t')
    vectors=pd.read_csv(HERE/'matched_state_vectors.tsv.gz',sep='\t')
    scalar=vectors.drop(columns=[c for c in vectors if 'usage' in c or c in ['PAS_ids','position_weights']])
    configurations=[]
    for gate,wt,mt in [('primary',.02,.1),('stronger',.05,.2)]:
        d=scalar.copy();d['gate']=gate;d['DWUI_threshold']=wt;d['MPRO_threshold']=mt
        for co in ['PS','PM','MS']:
            d['DWUI_sign_'+co]=np.where(d['DWUI_'+co]>=wt,1,np.where(d['DWUI_'+co]<=-wt,-1,0))
            d['MPRO_active_'+co]=d['MPRO_'+co]>=mt
        ps,pm,ms=[d['DWUI_sign_'+c] for c in ['PS','PM','MS']]
        d['DWUI_configuration']=np.select([
            ps*pm==1,ps*pm==-1,(ps==0)&(pm!=0),(ps!=0)&(ms==ps)&(pm==0),
            (ps!=0)&(pm==0)&(ms!=ps),(ps==0)&(pm==0)],
            ['dual_reference_aligned','reference_opposed','PM_event_PS_subthreshold',
             'shared_source_shift_PM_subthreshold','PS_event_PM_subthreshold_otherMS','both_P_contrasts_subthreshold'],default='invalid')
        assert not d.DWUI_configuration.eq('invalid').any()
        d['MPRO_configuration']=[pattern(a,b) for a,b in zip(d.MPRO_active_PS,d.MPRO_active_PM)]
        d['either_summary_PS']=d.DWUI_sign_PS.ne(0)|d.MPRO_active_PS
        d['either_summary_PM']=d.DWUI_sign_PM.ne(0)|d.MPRO_active_PM
        d['either_summary_configuration']=[pattern(a,b) for a,b in zip(d.either_summary_PS,d.either_summary_PM)]
        d['MPRO_PS_with_DWUI_below_gate']=d.MPRO_active_PS&d.DWUI_sign_PS.eq(0)
        d['MPRO_PM_with_DWUI_below_gate']=d.MPRO_active_PM&d.DWUI_sign_PM.eq(0)
        configurations.append(d)
    configs=pd.concat(configurations,ignore_index=True)
    save(configs,'all_matched_state_configurations.tsv.gz')
    # Keep missing studies/genes in the inventory, instead of converting them to zero.
    expected=[]
    for r in m.to_dict('records'):
        for setting in scalar.setting.unique():
            for study in STUDIES[r['species']]:
                expected.append(r|{'setting':setting,'study_id':study})
    grid=pd.DataFrame(expected)
    lookup=configs.drop(columns=['gene_symbol','DWUI_threshold','MPRO_threshold'])
    grid=grid.merge(lookup,on=KEY+['gate','setting'],how='left',validate='many_to_one')
    grid['matched_evaluable']=grid.DWUI_PM.notna()
    grid['source_eligible_study']=[st in REQUIRED[sp] for sp,st in zip(grid.species,grid.study_id)]
    grid['coverage_reason']=np.where(grid.matched_evaluable,'evaluable',
        np.where(grid.study_id.isin(['GSE150311','GSE221739']),'no_qualified_true_source',
        np.where(grid.study_id.eq('GSE178325'),'no_qualified_same_library_PM',
        np.where(grid.setting.eq('supported50_p08_mc5'),'no_complete_triplet_with_three_end_MC5','no_complete_qualified_triplet'))))
    save(grid,'member_study_state_inventory.tsv.gz')
    per_gene=[]
    groupkeys=['species','contrast','gate','setting','group','gene_id']
    for keys,g in grid.groupby(groupkeys,sort=False):
        rec=dict(zip(groupkeys,keys));sp=rec['species']
        sub=g[g.study_id.isin(REQUIRED[sp])]
        assert len(sub)==len(REQUIRED[sp])
        complete=bool(sub.matched_evaluable.all())
        rec.update(gene_symbol=g.gene_symbol.iloc[0],required_source_studies=len(sub),
                   matched_source_studies=int(sub.matched_evaluable.sum()),complete=complete,
                   frozen_DWUI_direction=g.frozen_DWUI_direction.iloc[0],
                   diagnostic_failure=g.diagnostic_failure.iloc[0])
        bools=['DWUI_PM_active_all','DWUI_PS_active_all','DWUI_PM_same_direction_all',
               'DWUI_dual_aligned_within_all','DWUI_dual_same_direction_all','DWUI_frozen_PM_retained_all',
               'DWUI_frozen_dual_retained_all','MPRO_PM_active_all','MPRO_PS_active_all','MPRO_dual_all',
               'either_summary_PM_active_all','either_summary_PS_active_all','either_summary_dual_all',
               'MPRO_PS_with_DWUI_below_gate_any','MPRO_PS_with_DWUI_below_gate_all',
               'MPRO_PM_with_DWUI_below_gate_any','MPRO_PM_with_DWUI_below_gate_all']
        rec.update({b:np.nan for b in bools})
        rec['MPRO_shared_configuration']='unavailable'
        rec['either_summary_shared_configuration']='unavailable'
        if complete:
            ps,pm=sub.DWUI_sign_PS,sub.DWUI_sign_PM
            rec['DWUI_PM_active_all']=bool(pm.ne(0).all())
            rec['DWUI_PS_active_all']=bool(ps.ne(0).all())
            rec['DWUI_PM_same_direction_all']=bool(pm.ne(0).all() and pm.nunique()==1)
            rec['DWUI_dual_aligned_within_all']=bool((ps*pm).eq(1).all())
            rec['DWUI_dual_same_direction_all']=bool((ps*pm).eq(1).all() and pm.nunique()==1)
            if pd.notna(rec['frozen_DWUI_direction']):
                rec['DWUI_frozen_PM_retained_all']=bool(pm.eq(rec['frozen_DWUI_direction']).all())
                rec['DWUI_frozen_dual_retained_all']=bool(pm.eq(rec['frozen_DWUI_direction']).all() and ps.eq(rec['frozen_DWUI_direction']).all())
            for metric,prefix in [('MPRO','MPRO'),('either_summary','either_summary')]:
                a=sub[metric+'_active_PM'] if metric=='MPRO' else sub.either_summary_PM
                b=sub[metric+'_active_PS'] if metric=='MPRO' else sub.either_summary_PS
                rec[prefix+'_PM_active_all']=bool(a.eq(True).all())
                rec[prefix+'_PS_active_all']=bool(b.eq(True).all())
                rec[prefix+'_dual_all']=bool(a.eq(True).all() and b.eq(True).all())
                cats=sub[metric+'_configuration']
                rec[prefix+'_shared_configuration']=cats.iloc[0] if cats.nunique()==1 else 'different_between_studies'
            for co in ['PS','PM']:
                s=sub['MPRO_'+co+'_with_DWUI_below_gate']
                rec['MPRO_'+co+'_with_DWUI_below_gate_any']=bool(s.eq(True).any())
                rec['MPRO_'+co+'_with_DWUI_below_gate_all']=bool(s.eq(True).all())
            rec['min_target_libraries']=int(sub.n_target_libraries.min())
            rec['min_source_libraries']=int(sub.n_source_libraries.min())
            rec['min_state_metacells']=int(sub[['P_min_mc','M_min_mc','S_min_mc']].min().min())
        per_gene.append(rec)
    per_gene=pd.DataFrame(per_gene)
    save(per_gene,'member_cross_study_states.tsv')
    sums=[]
    for keys,g in per_gene.groupby(['species','contrast','gate','setting','group']):
        r=dict(zip(['species','contrast','gate','setting','group'],keys))
        r.update(frozen_candidates=len(g),complete=int(g.complete.sum()),unavailable=int((~g.complete).sum()))
        for b in bools:
            r[b]=np.nan if b.startswith('DWUI_frozen_') and not g[b].notna().any() else int(g[b].eq(True).sum())
        for prefix in ['MPRO','either_summary']:
            for cat in PATTERNS+['different_between_studies']:
                r[prefix+'_'+cat+'_genes']=int(g[prefix+'_shared_configuration'].eq(cat).sum())
        sums.append(r)
    sums=pd.DataFrame(sums)
    save(sums,'group_state_summary.tsv')
    paired_sensitivity=[]
    sensitivity_keys=['species','contrast','gate','group','gene_id']
    comparison=per_gene[per_gene.setting.eq('supported50_p08')].merge(
        per_gene[per_gene.setting.eq('supported50_p08_mc5')],on=sensitivity_keys,
        validate='one_to_one',suffixes=('_primary','_MC5'))
    for keys,g in comparison.groupby(sensitivity_keys[:-1]):
        r=dict(zip(sensitivity_keys[:-1],keys));s=g[g.complete_primary&g.complete_MC5]
        r.update(frozen_candidates=len(g),primary_complete=int(g.complete_primary.sum()),
                 MC5_complete=int(g.complete_MC5.sum()),complete_in_both=len(s))
        for b in ['DWUI_dual_same_direction_all','MPRO_PS_active_all','MPRO_PM_active_all','MPRO_dual_all','either_summary_dual_all']:
            r[b+'_primary_on_common']=int(s[b+'_primary'].eq(True).sum())
            r[b+'_MC5']=int(s[b+'_MC5'].eq(True).sum())
            r[b+'_changed']=int(s[b+'_primary'].ne(s[b+'_MC5']).sum())
        paired_sensitivity.append(r)
    save(pd.DataFrame(paired_sensitivity),'coverage_sensitivity_comparison.tsv')
    bystudy=[]
    for keys,g in grid[grid.source_eligible_study].groupby(['species','contrast','gate','setting','group','study_id']):
        r=dict(zip(['species','contrast','gate','setting','group','study_id'],keys))
        for metric in ['DWUI','MPRO','either_summary']:
            counts=g.loc[g.matched_evaluable,metric+'_configuration'].value_counts()
            for label,n in counts.items():bystudy.append(r|{'metric':metric,'configuration':label,'genes':int(n),'evaluable':int(g.matched_evaluable.sum()),'frozen_candidates':len(g)})
        bystudy.append(r|{'metric':'coverage','configuration':'unavailable','genes':int((~g.matched_evaluable).sum()),'evaluable':int(g.matched_evaluable.sum()),'frozen_candidates':len(g)})
    save(pd.DataFrame(bystudy),'per_study_configuration_counts.tsv')
    # Original C membership and original Figure 3 fractions are locked acceptance targets.
    roster=pd.read_csv(ROOT/'notes/2026-09-27-figure3-single-cell-draft/formal-base-v4/panel_C_native_support.tsv',sep='\t')[['species','gene_id','gene_symbol']].drop_duplicates()
    assert len(roster)==18,len(roster)
    original=configs.merge(roster.drop(columns='gene_symbol'),on=['species','gene_id'],validate='many_to_one')
    original=original.merge(m[m.contrast.eq('P_M')][['species','gene_id','gate','group']],on=['species','gene_id','gate'],how='left',validate='many_to_one')
    save(original,'original_C_matched_metrics.tsv')
    legacy=[]
    legacy_summary=pd.read_csv(ROOT/'notes/2026-09-27-source-state-connection/frozen_cohort_summary.tsv',sep='\t')
    for sp,expected_n,expected_dual in [('mouse',200,97),('human',146,78)]:
        g=per_gene[per_gene.species.eq(sp)&per_gene.contrast.eq('P_M')&per_gene.gate.eq('primary')&per_gene.setting.eq('supported50_p08')&per_gene.group.isin(['both','DWUI_only'])]
        assert g.complete.sum()==expected_n
        reference=legacy_summary[legacy_summary.species.eq(sp)&legacy_summary.setting.eq('supported50_p08')&
            legacy_summary.threshold.eq(.02)&legacy_summary.cohort.eq('P_M')&
            legacy_summary.required_studies.eq(len(REQUIRED[sp]))].iloc[0]
        retained=int(g.DWUI_frozen_PM_retained_all.eq(True).sum())
        assert retained==int(reference.PM_frozen_aligned_genes)
        assert g.DWUI_frozen_dual_retained_all.eq(True).sum()==expected_dual
        legacy.append(dict(species=sp,complete=expected_n,PM_original_direction_and_gate_retained=retained,
                           dual=expected_dual,frozen=len(g)))
    # Sorted solely by coverage; no effect sizes enter the candidate ranking.
    candidate=per_gene[per_gene.contrast.eq('P_M')&per_gene.gate.eq('primary')&per_gene.setting.eq('supported50_p08')].copy()
    selected=set(zip(roster.species,roster.gene_id))
    candidate['original_C']= [(sp,gene) in selected for sp,gene in zip(candidate.species,candidate.gene_id)]
    candidate=candidate.sort_values(['species','group','complete','min_state_metacells','min_target_libraries','gene_id'],ascending=[True,True,False,False,False,True])
    candidate['coverage_rank_within_group']=candidate.groupby(['species','group']).cumcount()+1
    save(candidate,'candidate_coverage_ranking.tsv')
    highlight=configs[configs.gene_symbol.isin(['Slc38a2','Nudt21','Rad23b','Aopep','Tomm20','Fam3c','Mtpn'])]
    save(highlight,'discussed_genes_matched_metrics.tsv')
    # Coverage exclusions and all four unsigned configurations partition each frozen group.
    for _,r in sums.iterrows():
        assert r.complete+r.unavailable==r.frozen_candidates
        assert sum(r['MPRO_'+x+'_genes'] for x in PATTERNS+['different_between_studies'])==r.complete
    (HERE/'state_summary_audit.json').write_text(json.dumps(dict(status='PASS',
        matched_metric_records=len(configs),member_study_records=len(grid),member_gene_records=len(per_gene),
        original_C_genes=len(roster),recovered_original_counts=legacy,missing_never_imputed_zero=True,
        candidate_ranking_uses_coverage_only=True,new_P_values_generated=False),indent=2)+'\n')
    print(sums[sums.gate.eq('primary')&sums.setting.eq('supported50_p08')&sums.contrast.eq('P_M')][[
        'species','group','frozen_candidates','complete','DWUI_dual_same_direction_all','MPRO_dual_all',
        'MPRO_PM_only_genes','MPRO_PS_only_genes','MPRO_different_between_studies_genes']].to_string(index=False))


if __name__=='__main__':main()
