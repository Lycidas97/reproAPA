"""Batch primary selection with the pinned structural and neighbor rules.

The eight tier passes below operate on whole tables, not individual PAS frames.
Only a legacy single-key quicksort tie ordering is reproduced per tied group.
"""
import numpy as np
import pandas as pd


def flags(values):
    return np.where(values, 'true', 'false')


def choose(annotation, pas_map, candidates):
    a=annotation;c=candidates.copy();c['_order']=np.arange(len(c))
    same=c.strand_relation.eq('same_strand')
    informative=same & ~c.candidate_region_class.isin(['unassigned_pas','intergenic_distal_pas'])
    canonical=a._valid_canonical_terminal_candidates(c).index
    core=a._valid_same_strand_core_tes_candidates(c).index
    qualified=a._eligible_same_strand_tes_candidates(c)
    exon=a._same_strand_feature_family_candidates(c,a.EXONIC_CONFLICT_REGION_CLASSES).index
    intron=a._same_strand_feature_family_candidates(c,a.INTRON_REGION_CLASSES).index
    c['_tier']=np.select([c.index.isin(canonical),c.index.isin(core),c.index.isin(exon),
        c.index.isin(intron)&c.candidate_region_class.eq('terminal_intron_pas'),c.index.isin(intron),
        c.index.isin(qualified.index),informative,c.strand_relation.isin(['antisense','opposite_strand'])],
        np.arange(8),default=8)
    # A non-informative same-strand row never outranks any informative row.
    c=c.copy();pool=c[c._tier.eq(c.groupby('pas_id')._tier.transform('min'))].copy()
    pieces=[]
    for tier,q in pool.groupby('_tier',sort=True):
        q=q.copy();q['_ambiguous']=False
        if tier in [0,6]:
            q['_priority']=pd.to_numeric(q.region_priority,errors='coerce')
            q=q[q._priority.eq(q.groupby('pas_id')._priority.transform('min'))].copy()
            q['_ambiguous']=q.groupby('pas_id').candidate_gene_id.transform('nunique').gt(1)
            amb=q[q._ambiguous].sort_values(['pas_id','candidate_id']).drop_duplicates('pas_id')
            unique=q[~q._ambiguous].sort_values(['pas_id','distance_to_transcript_end','overlap_bases','_order'],ascending=[True,True,False,True]).drop_duplicates('pas_id')
            q=pd.concat([amb,unique])
        elif tier in [1,5]:
            q['_distance']=pd.to_numeric(q.downstream_distance_to_transcript_end,errors='raise')
            q=q[q._distance.eq(q.groupby('pas_id')._distance.transform('min'))].copy()
            q['_ambiguous']=q.groupby('pas_id').candidate_gene_id.transform('nunique').gt(1)
            q['_sort_overlap']=np.where(q._ambiguous,0,-q.overlap_bases)
            q=q.sort_values(['pas_id','_sort_overlap','candidate_id']).drop_duplicates('pas_id')
        elif tier in [2,3,4]:
            q['_biotype']=q.candidate_gene_biotype.map(a._gene_biotype_rank)
            q['_distance']=a._model_3p_distance(q)
            q=q[q._biotype.eq(q.groupby('pas_id')._biotype.transform('min'))].copy()
            q=q[q._distance.eq(q.groupby('pas_id')._distance.transform('min'))].copy()
            q['_ambiguous']=q.groupby('pas_id').candidate_gene_id.transform('nunique').gt(1)
            q['_sort_overlap']=np.where(q._ambiguous,0,-q.overlap_bases)
            q=q.sort_values(['pas_id','_sort_overlap','candidate_id']).drop_duplicates('pas_id')
        elif tier==7:
            # pandas' original one-key quicksort is unstable for ties. Reproduce
            # its selected index with NumPy, without constructing per-PAS frames.
            indices=[]
            priority=q.region_priority.to_numpy();positions=q.groupby('pas_id',sort=False).indices
            for ix in positions.values():
                indices.append(ix[np.argsort(priority[ix],kind='quicksort')[0]])
            q=q.iloc[indices]
        else:
            q=q.sort_values('_order').drop_duplicates('pas_id')
        pieces.append(q)
    p=pd.concat(pieces).set_index('pas_id').reindex(pas_map.pas_id).copy()
    assert p._tier.notna().all()
    tier=p._tier;amb=p._ambiguous;priority=pd.to_numeric(p.region_priority)
    region=p.candidate_region_class.copy()
    status=pd.Series(np.select([priority.le(4),priority.isin([5,6,7,8,9,10,11])],
        ['assigned_primary','assigned_sensitivity'],default='unassigned'),index=p.index)
    confidence=pd.Series(np.select([priority.le(4),priority.isin([5,6,7,8,9,10,11])],
        ['medium','low'],default='none'),index=p.index)
    reason=pd.Series(np.select([priority.le(4),priority.isin([5,6,7,8,9,11]),priority.eq(10)],
        ['same_strand_'+region,'sensitivity_'+region,'coding_exon_candidate_retained_for_sensitivity'],default='no_primary_candidate'),index=p.index)
    for mask,st,conf,why in [
        (tier.eq(7),'excluded_antisense','none','antisense_only'),
        (tier.eq(8),'unassigned','none','no_reliable_gene_candidate'),
        (amb & tier.isin([0,6]),'ambiguous_multi_gene','none','multiple_same_priority_gene_candidates'),
        (amb & tier.eq(1),'ambiguous_multi_gene_tes','none','exact_nearest_core_tes_distance_tie_across_genes'),
        (amb & tier.eq(5),'ambiguous_multi_gene_tes','none','exact_nearest_qualified_rescue_tes_distance_tie_across_genes'),
        (amb & tier.eq(2),'ambiguous_multi_gene','none','exact_other_exon_biotype_and_3p_distance_tie_across_genes'),
        (amb & tier.isin([3,4]),'ambiguous_multi_gene','none','exact_intron_biotype_and_3p_distance_tie_across_genes')]:
        status.loc[mask]=st;confidence.loc[mask]=conf;reason.loc[mask]=why
    region.loc[tier.eq(7)]='antisense_pas'
    region.loc[amb & tier.isin([0,6])]='multi_gene_ambiguous_pas'
    # Terminal-exon/CDS check uses the selected transcript (or gene if absent).
    joined=c.merge(p[['candidate_transcript_id','candidate_gene_id']].rename(columns={
        'candidate_transcript_id':'_selected_transcript','candidate_gene_id':'_selected_gene'}),left_on='pas_id',right_index=True,validate='many_to_one').reset_index(drop=True)
    model=np.where(joined._selected_transcript.astype(str).str.len().gt(0),
        joined.candidate_transcript_id.eq(joined._selected_transcript),joined.candidate_gene_id.eq(joined._selected_gene))
    model_rows=joined[model]
    cds=set(model_rows.loc[model_rows.candidate_region_class.eq('coding_exon_pas'),'pas_id'])
    utr=set(model_rows.loc[model_rows.candidate_region_class.eq('annotated_3utr_pas'),'pas_id'])
    overlap=region.eq('terminal_exon_pas') & p.index.isin(cds-utr) & tier.isin([0,6])
    region.loc[overlap]='terminal_exon_cds_overlap_pas';status.loc[overlap]='assigned_sensitivity'
    confidence.loc[overlap]='low';reason.loc[overlap]='terminal_exon_cds_overlap_retained_for_sensitivity'
    has_q=p.index.isin(qualified.pas_id)
    reason.loc[tier.eq(2)&has_q&~amb]='other_exon_over_qualified_tes2kb'
    reason.loc[tier.isin([3,4])&has_q&~amb]='intron_over_qualified_tes2kb'
    # Selected TES policy, including explicit missing-provenance reasons.
    selected_tes=region.eq(a.TES_REGION_CLASS)
    policy=pd.Series('not_selected_tes',index=p.index);policy_reason=pd.Series('primary_region_is_not_tes',index=p.index)
    required={'candidate':p.candidate_id,'gene':p.candidate_gene_id,'transcript':p.candidate_transcript_id,
        'transcript_end':p.candidate_transcript_end,'transcript_strand':p.candidate_strand,
        'strand_relation':p.strand_relation,'pas_strand':pas_map.set_index('pas_id').strand.reindex(p.index),
        'downstream_distance':p.downstream_distance_to_transcript_end}
    missing=pd.Series('',index=p.index)
    for name in sorted(required):
        empty=required[name].map(lambda x:'' if x is None else str(x).strip()).eq('')
        missing.loc[empty]+='_'+name
    distance=pd.to_numeric(p.downstream_distance_to_transcript_end,errors='coerce')
    invalid=~p.candidate_strand.isin(['+','-'])|~required['pas_strand'].isin(['+','-'])|p.strand_relation.ne('same_strand')|p.candidate_strand.ne(required['pas_strand'])
    policy_reason.loc[selected_tes]=np.select([missing.ne(''),invalid,~np.isfinite(distance),distance.lt(0),distance.le(250),distance.le(2000)],
        ['missing'+missing,'invalid_same_strand_tes_provenance','invalid_downstream_distance','negative_downstream_distance','core_tes_at_or_below_250bp','qualified_tes2kb_rescue'],default='distance_above_2000bp')[selected_tes]
    policy.loc[selected_tes]=np.select([policy_reason.eq('core_tes_at_or_below_250bp'),policy_reason.eq('qualified_tes2kb_rescue')],
        ['core_tes','qualified_tes2kb_rescue'],default='policy_ineligible_tes')[selected_tes]
    core_selected=selected_tes & policy.eq('core_tes') & status.ne('ambiguous_multi_gene_tes')
    status.loc[core_selected]='assigned_primary';confidence.loc[core_selected]='medium';reason.loc[core_selected]='core_tes_at_or_below_250bp'
    invalid_selected=selected_tes & policy.eq('policy_ineligible_tes')
    status.loc[invalid_selected]='assigned_sensitivity';confidence.loc[invalid_selected]='low';reason.loc[invalid_selected]='policy_ineligible_tes:'+policy_reason[invalid_selected]
    # Neighbor summaries are joins and grouped aggregations over the full batch.
    same_other=joined[joined.strand_relation.eq('same_strand')&joined.candidate_gene_id.astype(str).str.len().gt(0)&joined.candidate_gene_id.astype(str).ne(joined._selected_gene.astype(str))]
    opposite=joined[joined.strand_relation.isin(['antisense','opposite_strand'])&joined.candidate_gene_id.astype(str).str.len().gt(0)]
    def present(frame): return pd.Series(p.index.isin(frame.pas_id),index=p.index)
    def group_join(frame,column):
        return frame.groupby('pas_id',sort=False)[column].agg(a._sorted_join).reindex(p.index).fillna('')
    intronic_other=present(same_other[same_other.candidate_region_class.isin(a.INTRON_REGION_CLASSES)])
    exonic_other=present(same_other[same_other.candidate_region_class.isin(a.EXONIC_CONFLICT_REGION_CLASSES)])
    opp=present(opposite)
    structural=region.isin(a.EXONIC_CONFLICT_REGION_CLASSES|a.INTRON_REGION_CLASSES)
    displaced=qualified.merge(p[['candidate_id']].rename(columns={'candidate_id':'_selected_id'}),left_on='pas_id',right_index=True).reset_index(drop=True)
    displaced=displaced[displaced.pas_id.isin(p.index[structural])&displaced.candidate_id.astype(str).ne(displaced._selected_id.astype(str))]
    # Structural counterfactual selection: distance first versus biotype first.
    family=pd.Series(np.select([region.isin(a.EXONIC_CONFLICT_REGION_CLASSES),region.isin(a.INTRON_REGION_CLASSES)],['other_exon','intron'],default='not_applicable'),index=p.index)
    struct=c[c.pas_id.isin(p.index[structural])].copy()
    struct['_family']=struct.pas_id.map(family)
    struct=struct[struct.strand_relation.eq('same_strand')&struct.candidate_gene_id.astype(str).str.len().gt(0)&
        ((struct._family.eq('other_exon')&struct.candidate_region_class.isin(a.EXONIC_CONFLICT_REGION_CLASSES))|
         (struct._family.eq('intron')&struct.candidate_region_class.isin(a.INTRON_REGION_CLASSES)))].copy()
    terminal_ids=set(struct.loc[struct.candidate_region_class.eq('terminal_intron_pas'),'pas_id'])
    struct=struct[~struct.pas_id.isin(terminal_ids)|struct.candidate_region_class.eq('terminal_intron_pas')].copy()
    struct['_distance']=a._model_3p_distance(struct)
    nearest=struct[struct._distance.eq(struct.groupby('pas_id')._distance.transform('min'))].copy()
    nearest=nearest[nearest.groupby('pas_id').candidate_gene_id.transform('nunique').eq(1)].sort_values(['pas_id','overlap_bases','candidate_id'],ascending=[True,False,True]).drop_duplicates('pas_id').set_index('pas_id')
    feature_id=p.candidate_id.astype(str).copy();feature_gene=p.candidate_gene_id.astype(str).copy()
    feature_id.loc[structural]=nearest.candidate_id.reindex(p.index[structural]).fillna('')
    feature_gene.loc[structural]=nearest.candidate_gene_id.reindex(p.index[structural]).fillna('')
    biotype_id=p.candidate_id.astype(str).copy();biotype_gene=p.candidate_gene_id.astype(str).copy()
    biotype_conflict=feature_gene.ne('')&biotype_gene.ne('')&feature_gene.ne(biotype_gene)
    ranks=struct.assign(_rank=struct.candidate_gene_biotype.map(a._gene_biotype_rank)).groupby('pas_id')._rank.nunique()
    tie=pd.Series('guarded_canonical_core_or_nonstructural',index=p.index)
    tie.loc[structural]=np.where(ranks.reindex(p.index[structural]).gt(1),'biotype_selected','distance_selected')
    ambiguity=structural&status.eq('ambiguous_multi_gene');biotype_id.loc[ambiguity]='';biotype_gene.loc[ambiguity]='';tie.loc[ambiguity]='exact_cross_gene_tie'
    competing=qualified.groupby('pas_id').candidate_gene_id.nunique().reindex(p.index,fill_value=0).gt(1)
    qother=qualified.merge(p[['candidate_gene_id']].rename(columns={'candidate_gene_id':'_selected_gene'}),left_on='pas_id',right_index=True)
    canonical_over=region.isin(a.CANONICAL_TERMINAL_REGION_CLASSES)&present(qother[qother.candidate_gene_id.astype(str).ne(qother._selected_gene.astype(str))])
    conditions=[ambiguity,status.eq('ambiguous_multi_gene_tes'),canonical_over,
        policy.eq('core_tes')&exonic_other,policy.eq('core_tes'),
        policy.eq('policy_ineligible_tes')&policy_reason.eq('distance_above_2000bp'),policy.eq('policy_ineligible_tes'),
        region.isin(a.EXONIC_CONFLICT_REGION_CLASSES)&has_q,region.isin(a.INTRON_REGION_CLASSES)&has_q,
        policy.eq('qualified_tes2kb_rescue')&exonic_other,policy.eq('qualified_tes2kb_rescue')&competing,
        policy.eq('qualified_tes2kb_rescue')&opp,policy.eq('qualified_tes2kb_rescue')]
    rules=['select_structural_feature_biotype_distance_exact_tie','select_nearest_same_strand_tes_exact_multi_gene_tie','canonical_terminal_other_gene_over_tes2kb',
        'core_tes_exonic_overlap_audit_only','core_tes_assignment_preserved','tes_candidate_outside_rescue_distance','tes_candidate_incomplete_provenance',
        'same_strand_other_exon_over_qualified_tes2kb','same_strand_intron_over_qualified_tes2kb','same_strand_exonic_gene_conflict_to_sensitivity',
        'select_unique_nearest_same_strand_tes','tes2kb_by_pas_transcript_strand_over_opposite_gene','tes2kb_no_disqualifying_neighbor_conflict']
    terminals=['ambiguous_multi_gene_structural','ambiguous_multi_gene_tes','canonical_terminal_other_gene_over_tes2kb',a.CORE_TES_EXONIC_AUDIT_STATE,
        a.CORE_TES_STATE,a.OUT_OF_RANGE_TES_STATE,a.INCOMPLETE_TES_STATE,'other_exon_over_tes2kb','intron_over_tes2kb',
        'sensitivity_same_strand_exonic_gene_conflict','tes2kb_rescue_nearest_competing_tes','tes2kb_rescue_opposite_strand_gene_overlap','tes2kb_rescue_no_neighbor_conflict']
    exclusions=['exact_feature_biotype_3p_distance_tie_across_genes','exact_minimum_tes_distance_tied_across_genes','','','',
        'tes_candidate_outside_rescue_distance','tes_candidate_incomplete_provenance','other_exon_over_tes2kb','intron_over_tes2kb','same_strand_exonic_gene_conflict','','','']
    out=pas_map.set_index('pas_id',drop=False).copy()
    mapping={'primary_candidate_id':'candidate_id','primary_gene_id':'candidate_gene_id','primary_gene_symbol':'candidate_gene_symbol',
        'primary_gene_biotype':'candidate_gene_biotype','primary_transcript_id':'candidate_transcript_id','primary_transcript_biotype':'candidate_transcript_biotype',
        'primary_model_transcript_end':'candidate_model_transcript_end','primary_model_3p_distance':'candidate_model_3p_distance','primary_transcript_strand':'candidate_strand',
        'primary_strand_relation':'strand_relation','primary_transcript_end':'candidate_transcript_end','distance_to_primary_transcript_end':'distance_to_transcript_end',
        'signed_distance_to_primary_transcript_end':'signed_distance_to_transcript_end','downstream_distance_to_primary_transcript_end':'downstream_distance_to_transcript_end'}
    for target,source in mapping.items():out[target]=p[source]
    out['primary_region_class']=region;out['primary_assignment_status']=status;out['primary_assignment_confidence']=confidence;out['primary_assignment_reason']=reason
    out['multi_gene_ambiguous']=flags(status.isin(['ambiguous_multi_gene','ambiguous_multi_gene_tes']))
    for col in c:
        if col.startswith(('same_strand_candidate_','opposite_strand_candidate_','original_tes_','valid_same_strand_core_tes_','qualified_same_strand_tes_')) or col in ['candidate_gene_ids','candidate_region_classes']:
            out[col]=p[col]
    out['same_strand_other_gene_intron_overlap']=flags(intronic_other);out['same_strand_exonic_gene_conflict']=flags(exonic_other);out['opposite_strand_gene_overlap']=flags(opp)
    for prefix,frame in [('neighbor_same_strand',same_other),('neighbor_opposite_strand',opposite)]:
        out[prefix+'_gene_ids']=group_join(frame,'candidate_gene_id');out[prefix+'_region_classes']=group_join(frame,'candidate_region_class')
    out['selected_tes_candidate_id']=np.where(selected_tes,p.candidate_id.astype(str),'')
    out['selected_tes_policy_class']=policy;out['selected_tes_policy_reason']=policy_reason
    out['neighbor_gene_decision_rule']=np.select(conditions,rules,default='region_priority_assignment')
    out['neighbor_gene_terminal_state']=np.select(conditions,terminals,default=status)
    out['neighbor_gene_exclusion_reason']=np.select(conditions,exclusions,default=np.where(status.isin(['assigned_primary','assigned_sensitivity']),'',status))
    for suffix,column in [('candidate_ids','candidate_id'),('gene_ids','candidate_gene_id'),('transcript_ids','candidate_transcript_id')]:
        out['displaced_tes_'+suffix]=group_join(displaced,column)
    for col,value in dict(structural_priority_feature_family=family,feature_first_candidate_id=feature_id,feature_first_gene_id=feature_gene,
        biotype_priority_candidate_id=biotype_id,biotype_priority_gene_id=biotype_gene,biotype_assignment_conflict=flags(biotype_conflict),biotype_tie_break_status=tie).items():out[col]=value
    # Transcript-family summaries shared across all selected rows.
    valid=c[c.strand_relation.eq('same_strand')&c.candidate_gene_id.astype(str).str.len().gt(0)&c.candidate_transcript_id.astype(str).str.len().gt(0)].copy()
    valid['_priority']=pd.to_numeric(valid.region_priority,errors='coerce').fillna(99)
    valid=valid[valid._priority.eq(valid.groupby(['pas_id','candidate_gene_id','candidate_transcript_id'])._priority.transform('min'))].copy()
    valid['_family']=valid.candidate_region_class.map(a._transcript_region_family)
    fam=valid[valid._family.ne('')][['pas_id','candidate_gene_id','_family']].drop_duplicates()
    conflict=fam.groupby(['pas_id','candidate_gene_id']).size().gt(1).groupby('pas_id').any().reindex(p.index,fill_value=False)
    terminal=fam[fam._family.eq('terminal_3end')];intr=fam[fam._family.eq('intronic')]
    within=terminal.merge(intr,on=['pas_id','candidate_gene_id']).pas_id
    across=terminal.merge(intr,on='pas_id',suffixes=('_a','_b'));across=across[across.candidate_gene_id_a.ne(across.candidate_gene_id_b)].pas_id
    out['transcript_model_conflict']=flags(conflict)
    out['transcript_context_class']=np.select([p.index.isin(within),p.index.isin(across),conflict],
        ['alternative_terminal_intron','cross_gene_terminal_intron','alternative_transcript_context'],default='none')
    for col,value in dict(internal_priming_flag='false',external_support_status='unavailable',total_pas_umi='',n_cells_detected='',quantification_status='not_joined').items():out[col]=value
    # One constant-size reference call determines the legacy output column order.
    first=pas_map.iloc[:1];schema=a.choose_primary_assignments(first,candidates[candidates.pas_id.eq(first.pas_id.iloc[0])]).columns
    return out.reset_index(drop=True).reindex(columns=schema).infer_objects()
