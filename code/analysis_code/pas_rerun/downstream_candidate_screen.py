"""Complete current-cohort screens; descriptive example ranks, never significance ranks."""
import argparse
import numpy as np
import pandas as pd
from downstream_runtime import Context
from downstream_examples import shared_regions


def low_shift(native,regions,well):
    f=native[native.species.eq('mouse')&native.route.isin(['CIP','WOT'])&native.contrast.eq('P_M')&native.design.eq('same_library')]
    rows=[]
    for gid,q in f.groupby('gene_id'):
        assert not q.route.duplicated().any()
        if set(q.route)!={'CIP','WOT'} or q.delta_DWUI.abs().max()>=.02 or q.MPRO.min()<.1:continue
        reg=regions[regions.gene_id.eq(gid)&regions.contrast.eq('MP')]
        w=well[well.gene_id.eq(gid)&well.metric.eq('MPRO')]
        rows.append(dict(gene_id=gid,gene_symbol=q.gene_symbol.iloc[0],min_MPRO=q.MPRO.min(),max_abs_DWUI=q.delta_DWUI.abs().max(),
            min_target_libraries=q.n_target_libraries.min(),strict_region_pairs=len(reg),
            unique_interval_mask_pairs=len(reg.drop_duplicates(['mask_a','mask_b'])),well_MPRO_q=w.q_BH_joint.iloc[0] if len(w) else np.nan,
            **{k+'_'+r.route:getattr(r,k) for r in q.itertuples() for k in ['delta_DWUI','MPRO','n_target_libraries']}))
    cols=['gene_id','gene_symbol','min_MPRO','max_abs_DWUI','min_target_libraries','strict_region_pairs','unique_interval_mask_pairs','well_MPRO_q']+[k+'_'+r for r in ['CIP','WOT'] for k in ['delta_DWUI','MPRO','n_target_libraries']]
    out=pd.DataFrame(rows,columns=cols).sort_values(['min_MPRO','gene_id'],ascending=[False,True]);out.insert(0,'MPRO_only_rank',range(1,len(out)+1))
    return out


def run(ctx):
    matched=ctx.read(ctx.analysis/'matched_pms/matched_state_vectors.tsv.gz');primary=matched[matched.setting.eq('supported50_p08')]
    rec=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz')
    rec=rec[rec.metric.eq('WPP')&rec.threshold.eq(.02)&rec.series_matching&rec.contrast.eq('P_M')]
    ladder=[];dual=[]
    for species,routes in [('human',{'GSE118258'}),('mouse',{'CIP','WOT'})]:
        r=rec[rec.species.eq(species)];assert r.gene_id.is_unique
        measurable=0;directions=[]
        for event in r.itertuples():
            q=primary[primary.species.eq(species)&primary.gene_id.eq(event.gene_id)&primary.route.isin(routes)]
            assert not q.route.duplicated().any()
            if set(q.route)!=routes:continue
            measurable+=1;sign=1 if event.direction=='positive' else -1
            if not ((q[['DWUI_PS','DWUI_PM']]*sign)>=.02).all().all():continue
            directions.append(sign)
            dual.append(dict(species=species,gene_id=event.gene_id,gene_symbol=q.gene_symbol.iloc[0],direction=event.direction,
                minimum_abs_effect=q[['DWUI_PS','DWUI_PM']].abs().min().min(),routes=';'.join(sorted(routes)),
                min_target_libraries=q.n_target_libraries.min(),min_target_dates=q.n_target_dates.min()))
        ladder.append(dict(species=species,route_scope=';'.join(sorted(routes)),recurrent_candidates=len(r),matched_measurable=measurable,
                           dual_reference=len(directions),distal=sum(x==1 for x in directions),proximal=sum(x==-1 for x in directions)))
    dual=pd.DataFrame(dual,columns=['species','gene_id','gene_symbol','direction','minimum_abs_effect','routes','min_target_libraries','min_target_dates']).sort_values(['species','minimum_abs_effect','gene_id'],ascending=[True,False,True])
    ctx.table(dual,'dual_reference_candidates.tsv');ctx.table(ladder,'dual_reference_ladder.tsv')
    sites=ctx.read(ctx.analysis/'site_statistics/recurrent_sites.tsv.gz');regions=shared_regions(sites)
    native=ctx.read(ctx.analysis/'state_vectors/route_effects.primary.tsv.gz');well=ctx.read(ctx.analysis/'well_inference/gene_metric_tests.tsv.gz')
    low=low_shift(native,regions,well);ctx.table(low,'low_shift_candidates.tsv')
    ctx.table(regions,'recurrent_region_pairs.tsv.gz')
    strict=low[low.strict_region_pairs.gt(0)];final=strict[strict.max_abs_DWUI.lt(.01)]
    ctx.table([dict(stage=stage,genes=len(q),interval_pairs=int(q.strict_region_pairs.sum()),gene_ids=';'.join(q.gene_id)) for stage,q in
               [('low_abs_DWUI_lt002_MPRO_ge01',low),('strict_recurrent_interval',strict),('both_abs_DWUI_lt001',final)]],'low_shift_attrition.tsv')
    partial=ctx.read(ctx.analysis/'site_statistics/cross_study_partial_retention.tsv');partial=partial[partial.contrast.isin(['SP','MP'])]
    ctx.table(partial,'partial_retention_relations.tsv.gz')
    rows=[]
    for key,q in partial[partial.exact_partial_recurrence].groupby(['species','contrast','gene_symbol_a']):
        changed=[set(str(v).split(';'))-{'','nan'} for v in q.shared_recurrent_changed_sites]
        retained=[set(str(v).split(';'))-{'','nan'} for v in q.shared_retained_sites]
        rows.append(dict(zip(['species','contrast','gene_symbol'],key))|dict(gene_factor_relationships=len(q),
            changed_site_factor_pairs=sum(map(len,changed)),retained_site_factor_pairs=sum(map(len,retained)),
            unique_changed_coordinates=len(set.union(*changed)),unique_retained_coordinates=len(set.union(*retained))))
    part=pd.DataFrame(rows,columns=['species','contrast','gene_symbol','gene_factor_relationships','changed_site_factor_pairs','retained_site_factor_pairs','unique_changed_coordinates','unique_retained_coordinates'])
    ctx.table(part,'partial_retention_gene_summary.tsv')
    # Selection is explicit and saved; drawing never selects a convenient gene.
    selected=[];limit=ctx.config.get('example_candidates_per_role',3)
    for sp,q in dual.groupby('species'):
        selected.extend(dict(role='dual_reference',species=sp,gene_id=r.gene_id,gene_symbol=r.gene_symbol,
            rank=i,selection='descending minimum absolute effect; descriptive') for i,r in enumerate(q.head(limit).itertuples(),1))
    selected.extend(dict(role='low_shift_interval',species='mouse',gene_id=r.gene_id,gene_symbol=r.gene_symbol,
        rank=i,selection='strict interval then descending minimum MPRO; descriptive') for i,r in enumerate(strict.head(limit).itertuples(),1))
    for sp,q in part.groupby('species'):
        q=q.sort_values(['gene_factor_relationships','gene_symbol'],ascending=[False,True]).drop_duplicates('gene_symbol').head(limit)
        selected.extend(dict(role='partial_retention',species=sp,gene_id='',gene_symbol=r.gene_symbol,rank=i,
            selection='descending exact partial relationships; descriptive') for i,r in enumerate(q.itertuples(),1))
    ctx.table(selected,'example_selection.tsv',columns=['role','species','gene_id','gene_symbol','rank','selection'])
    ctx.finish(dict(full_current_cohorts=True,no_old_gene_or_count_requirement=True,route_scopes_explicit=True,
                    ranks_descriptive_only=True,strict_interval_screen_not_absence_claim=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'candidate_screen',a.out))
