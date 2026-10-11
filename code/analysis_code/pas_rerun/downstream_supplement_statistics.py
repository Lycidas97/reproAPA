"""Persist the numerical joins and selections required by the production supplement.

This node consumes the current run only. Legacy numerical functions are replayed
solely to recover unpersisted LOO intermediates, and must reproduce their receipts.
"""
import argparse
import ast
import json
import sys

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

from downstream_runtime import Context
from downstream_continuous import effects_frame
from downstream_profiles import FIT

LOO_SOURCE='github/reproapa-analysis/modules/07_apa_dynamics/code/figure2/analyze_figure2_dynamic_architecture_v1.py'


def loo_tables(ctx):
    members=ctx.read(ctx.analysis/'profiles/gene_membership_k3.tsv.gz')
    nodes=ctx.read(ctx.analysis/'profiles/all_route_node_positions.tsv.gz')
    nodes=nodes[nodes.depth.eq(20)].copy()
    nodes['layer']=nodes.setting.map({'supported50_p08':'supported50','expanded_primary_p08':'expanded_primary'})
    fit=ctx.register(ctx.analysis/'profiles/joint_fit.npz')
    functions=ctx.functions(LOO_SOURCE,['unit_centered','loo_route_validation'],dict(
        FIT=fit,STATES=list('SREP'),CLUSTERS=['C1','C2','C3'],ROUTES=FIT,SEED=91327,
        KMeans=KMeans,linear_sum_assignment=linear_sum_assignment,adjusted_rand_score=adjusted_rand_score))
    fn=functions.loo_route_validation
    tree=ast.parse(open(fn.__code__.co_filename).read())
    block=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==fn.__name__)
    lines=[n.lineno for n in ast.walk(block) if isinstance(n,ast.Expr) and isinstance(n.value,ast.Call)
           and isinstance(n.value.func,ast.Attribute) and isinstance(n.value.func.value,ast.Name)
           and n.value.func.value.id=='rows' and n.value.func.attr=='append']
    assert len(lines)==1,'Reviewed LOO capture point changed'
    records=[];centers=[];seen=set()
    def trace(frame,event,arg):
        if frame.f_code is fn.__code__ and event=='line' and frame.f_lineno==lines[0]:
            v=frame.f_locals;sp,route=v['held_sp'],v['held_route']
            if (sp,route) in seen:return trace
            seen.add((sp,route));ids=v['ids'][v['species']==sp]
            for i,gene in enumerate(ids):
                records.append(dict(species=sp,heldout_route=route,gene_id=gene,
                    full_cluster=f"C{v['full_labels'][i]+1}",projected_cluster=f"C{v['projected'][i]+1}",
                    remaining_route_cluster=f"C{v['training_for_held'][i]+1}" if len(v['same_species_train']) else '',
                    **{f'heldout_{s}':float(v['held'][i,j]) for j,s in enumerate('SREP')}))
            centers.extend(dict(species=sp,heldout_route=route,cluster=f'C{j+1}',
                **dict(zip('SREP',center))) for j,center in enumerate(v['centers']))
        return trace
    previous=sys.gettrace()
    try:
        sys.settrace(trace);summary,metrics,support=fn(members,nodes)
    finally:sys.settrace(previous)
    for name,frame,keys in [('leave_route_projection.tsv',summary,['heldout_species','heldout_route']),
        ('leave_route_centroids.tsv',metrics,['heldout_species','heldout_route','cluster']),
        ('route_shape_support.tsv',support,['species','route','cluster'])]:
        expected=ctx.read(ctx.analysis/'profiles'/name)
        pd.testing.assert_frame_equal(frame.sort_values(keys).reset_index(drop=True),
            expected.sort_values(keys).reset_index(drop=True),check_dtype=False,rtol=1e-9,atol=1e-12)
    assignments=pd.DataFrame(records);display=[]
    assert not assignments.duplicated(['species','heldout_route','gene_id']).any()
    for r in summary.itertuples():
        q=assignments[assignments.species.eq(r.heldout_species)&assignments.heldout_route.eq(r.heldout_route)]
        assert len(q)==r.n_heldout
        np.testing.assert_allclose((q.full_cluster==q.projected_cluster).mean(),r.assignment_agreement_vs_full)
        np.testing.assert_allclose(adjusted_rand_score(q.full_cluster,q.projected_cluster),r.ARI_vs_full)
    for r in centers:
        q=assignments[assignments.species.eq(r['species'])&assignments.heldout_route.eq(r['heldout_route'])&assignments.full_cluster.eq(r['cluster'])]
        obs=q[['heldout_'+s for s in 'SREP']].mean().to_numpy();pred=np.array([r[s] for s in 'SREP'])
        for i,state in enumerate('SREP'):
            display.append(dict(species=r['species'],heldout_route=r['heldout_route'],cluster=r['cluster'],state=state,
                n=len(q),heldout_mean=obs[i],aligned_center=pred[i],difference=obs[i]-pred[i]))
    ctx.table(assignments,'LOO-gene-assignments.tsv.gz')
    ctx.table(centers,'LOO-aligned-centers.tsv')
    ctx.table(display,'LOO-center-display.tsv')
    flow=assignments.groupby(['species','heldout_route','full_cluster','projected_cluster']).size().rename('genes').reset_index()
    flow['route_total']=flow.groupby(['species','heldout_route']).genes.transform('sum')
    ctx.table(flow,'LOO-assignment-flows.tsv')


def dual_points(primary,recurrence):
    primary=primary[primary.setting.eq('supported50_p08')]
    recurrence=recurrence[recurrence.metric.eq('WPP')&recurrence.threshold.eq(.02)&recurrence.series_matching&recurrence.contrast.eq('P_M')]
    parts=[];audit=[]
    for species,routes in [('human',{'GSE118258'}),('mouse',{'CIP','WOT'})]:
        candidates=recurrence[recurrence.species.eq(species)];assert candidates.gene_id.is_unique
        for event in candidates.itertuples():
            q=primary[primary.species.eq(species)&primary.gene_id.eq(event.gene_id)&primary.route.isin(routes)].copy()
            assert not q.route.duplicated().any()
            complete=set(q.route)==routes
            sign=1 if event.direction=='positive' else -1
            qualified=complete and bool(q[['DWUI_PS','DWUI_PM']].mul(sign).ge(.02).all().all())
            audit.append(dict(species=species,gene_id=event.gene_id,matched_measurable=complete,
                dual_qualified_gene=qualified,required_routes=';'.join(sorted(routes)),observed_routes=';'.join(sorted(q.route))))
            if complete:parts.append(q.assign(recurrent_direction=event.direction,dual_qualified_gene=qualified))
    columns=[*primary.columns,'recurrent_direction','dual_qualified_gene']
    return (pd.concat(parts,ignore_index=True) if parts else pd.DataFrame(columns=columns),
        pd.DataFrame(audit,columns=['species','gene_id','matched_measurable','dual_qualified_gene','required_routes','observed_routes']))


def relationship_tables(partial):
    q=partial[partial.exact_partial_recurrence&partial.contrast.isin(['SP','MP'])]
    keys=['species','category','gene_id','factor']
    members=q.groupby(keys,dropna=False).contrast.agg(lambda s:';'.join(sorted(set(s)))).rename('contrasts').reset_index()
    rows=[]
    for species in ['human','mouse']:
        for category in ['miRNA','RBP']:
            g=members[members.species.eq(species)&members.category.eq(category)]
            rows.append(dict(species=species,category=category,relationships=len(g),genes=g.gene_id.nunique(),contrast_scope='SP;MP'))
    return members,pd.DataFrame(rows)


def display_tables(ctx):
    archive=np.load(ctx.register(ctx.analysis/'shape_ortholog/null_draws.npz'))
    tests=ctx.read(ctx.analysis/'shape_ortholog/shape_tests.tsv');rows=[]
    for model in archive.files:
        values=np.trace(archive[model],axis1=1,axis2=2)
        test=tests[tests.model.eq(model)&tests.test.eq('diagonal')].iloc[0]
        np.testing.assert_allclose(values.mean(),test.null_mean)
        counts=pd.Series(values).value_counts().sort_index()
        rows.extend(dict(model=model,same_class_genes=int(n),draws=int(count),total_draws=len(values),
            probability=count/len(values),observed=test.observed,expected=test.expected) for n,count in counts.items())
    ctx.table(rows,'shape-null-PMF.tsv')
    terms=[];term_counts=[]
    for database in ['go','kegg']:
        f=ctx.read(ctx.analysis/f'functional/{database}_ladder.tsv.gz')
        tested=f[f.model.eq('L3_full')&f.comparison.eq('G1_vs_G0')]
        selected=tested[tested.q_database_model.lt(.05)].copy()
        terms.append(selected)
        for key,q in tested.groupby(['database','analysis_id']):
            term_counts.append(dict(database=key[0],analysis_id=key[1],tested_terms=len(q),display_terms=int(q.q_database_model.lt(.05).sum()),
                selection='L3_full;G1_vs_G0;q_database_model<0.05;no_fixed_term_count'))
    terms=pd.concat(terms,ignore_index=True).sort_values(['database','analysis_id','q_database_model','term_id'])
    terms['display_order']=range(len(terms));ctx.table(terms,'endpoint-GO-KEGG-display.tsv')
    ctx.table(term_counts,'endpoint-term-denominators.tsv')
    nodes=ctx.read(ctx.analysis/'cross_overlap/all_mapped_nodes.tsv.gz')
    pairs=ctx.read(ctx.analysis/'cross_continuous/continuous_pair_summary.tsv')
    parts=[];checks=[]
    for r in pairs[pairs.layer.eq('supported50')&pairs.state.isin(list('REP'))].itertuples():
        q=effects_frame(nodes,r.layer,r.mouse,r.human,[r.state])
        rho=q['m_'+r.state].corr(q['h_'+r.state],method='spearman')
        assert len(q)==r.n
        np.testing.assert_allclose(rho,r.rho,rtol=1e-9,atol=1e-12,equal_nan=True)
        q=q.rename(columns={'m_'+r.state:'mouse_effect','h_'+r.state:'human_effect'}).reset_index()
        parts.append(q.assign(layer=r.layer,state=r.state,mouse=r.mouse,human=r.human))
        checks.append(dict(layer=r.layer,state=r.state,mouse=r.mouse,human=r.human,n=len(q),rho=rho))
    ctx.table(pd.concat(parts,ignore_index=True) if parts else [],'continuous-gene-pairs.tsv.gz',
        columns=['human_gene_id','mouse_effect','human_effect','layer','state','mouse','human'])
    ctx.table(checks,'continuous-display-checks.tsv',columns=['layer','state','mouse','human','n','rho'])


def select_partial_sites(partial,exact,features,selection,allow_new_features=False):
    """Choose one reproducible site pair per species, without preferring a seed."""
    chosen=[];statuses=[]
    for species in ['mouse','human']:
        genes=selection[selection.role.eq('partial_retention')&selection.species.eq(species)].sort_values(['rank','gene_symbol'])
        selected=[]
        for gene in genes.itertuples():
            rel=partial[partial.species.eq(species)&partial.gene_symbol_a.str.upper().eq(gene.gene_symbol.upper())&
                partial.contrast.isin(['MP','SP'])&partial.exact_partial_recurrence].sort_values(['contrast','category','factor','gene_id'])
            for r in rel.itertuples():
                pair=[]
                for role,coordinates in [('changed',r.shared_recurrent_changed_sites),('retained',r.shared_retained_sites)]:
                    for coord in sorted(set(str(coordinates).split(';'))-{'','nan'}):
                        chrom,position=coord.rsplit(':',1);start,end=map(int,position.split('-'))
                        q=exact[exact.species.eq(species)&exact.gene_id.eq(r.gene_id)&exact.contrast.eq(r.contrast)&
                            exact.category.eq(r.category)&exact.factor.eq(r.factor)&exact.chrom.eq(chrom)&exact.start0.eq(start)&exact.end0.eq(end)]
                        if q.empty:raise AssertionError('Partial recurrence coordinate absent from exact-site table')
                        assert len(q)==1
                        e=q.iloc[0];assert bool(e.recurrent_strict) if role=='changed' else bool(e.retained_a and e.retained_b)
                        site=[]
                        for side in ['a','b']:
                            route=e['route_'+side]
                            f=features[features.species.eq(species)&features.gene_id.eq(r.gene_id)&features.route.eq(route)&
                                features.kind.eq('site')&features.label.eq(r.factor)&features.start0.eq(start)&features.end0.eq(end)]
                            mask=np.asarray(json.loads(e['mask_'+side]))
                            assert mask.ndim==2 and np.array_equal(mask[:,0],mask[:,1])
                            if len(f)!=1:
                                if not allow_new_features:break
                                assert f.empty
                                feature_id=f'{route}|{r.gene_id}|site|{start}-{end}|{r.factor}'
                                mask_text=json.dumps(mask[:,0].tolist())
                            else:
                                np.testing.assert_array_equal(mask[:,0],json.loads(f.iloc[0]['mask']))
                                feature_id=f.iloc[0].feature_id;mask_text=f.iloc[0]['mask']
                            site.append(dict(species=species,gene_id=r.gene_id,gene_symbol=gene.gene_symbol,contrast=r.contrast,
                                category=r.category,factor=r.factor,site_role=role,route=route,chrom=chrom,strand=e.strand,
                                start0=start,end0=end,feature_id=feature_id,mask=mask_text))
                        if len(site)==2:pair.extend(site);break
                if len(pair)==4:selected=pair;break
            if selected:break
        chosen.extend(selected)
        statuses.append(dict(species=species,status='selected' if selected else 'no_current_exported_exact_site_pair',
            selection='current ranked partial-retention candidates; lexical contrast/factor/coordinate; no favoured seed'))
    cols=['species','gene_id','gene_symbol','contrast','category','factor','site_role','route','chrom','strand','start0','end0','feature_id','mask']
    return pd.DataFrame(chosen,columns=cols),pd.DataFrame(statuses)


def selected_inclusion(ctx,selected,prefix='partial-site'):
    """Late-selected examples use the already quantified terminal models, never BAMs."""
    from downstream_example_inclusion import qualified
    models={};axes=[];masks=[];libraries=[];summaries=[];paired=[]
    for r in selected.itertuples():
        key=(r.species,r.route)
        if key not in models:models[key]=json.loads(ctx.register(ctx.analysis/f'intervals/annotation/{r.species}/{r.route}/models.json').read_text())
        model=models[key][r.gene_id];mask=np.array(json.loads(r.mask));assert len(mask)==len(model['pas'])
        base=dict(feature_id=r.feature_id,species=r.species,gene_id=r.gene_id,gene=r.gene_symbol,route=r.route,site_role=r.site_role,contrast=r.contrast)
        for order,(pas,value) in enumerate(zip(model['pas'],mask)):
            a=dict(species=r.species,gene_id=r.gene_id,route=r.route,pas_order=order,atlas_pas_id=pas['atlas_pas_id'],cut_interbase=pas['cut_interbase'])
            axes.append(a);masks.append(dict(a,feature_id=r.feature_id,site_role=r.site_role,retains_site=int(value)))
        for cutoff in [20,50]:
            for state in 'SREPM':
                values=[]
                for lib,a in qualified(model,state,cutoff).items():
                    counts=np.asarray(a['counts']);assert counts.sum()==a['total'] and len(counts)==len(mask)
                    value=float(counts@mask/counts.sum());values.append(value)
                    libraries.append(dict(base,cutoff=cutoff,state=state,library_key=lib,inclusion=value,family_UMI=a['total'],coverage=a['coverage']))
                qs=np.quantile(values,[.25,.5,.75]) if values else [np.nan]*3
                summaries.append(dict(base,cutoff=cutoff,state=state,n=len(values),mean=np.mean(values) if values else np.nan,
                    q25=qs[0],q50=qs[1],q75=qs[2],status='measured' if values else 'no_qualified_library',uncertainty='empirical_IQR_not_CI'))
            p=qualified(model,'P',cutoff);m=qualified(model,'M',cutoff)
            for lib in sorted(set(p)&set(m)):
                pv=float(np.asarray(p[lib]['counts'])@mask/p[lib]['total']);mv=float(np.asarray(m[lib]['counts'])@mask/m[lib]['total'])
                paired.append(dict(base,cutoff=cutoff,library_key=lib,reference_inclusion=mv,target_inclusion=pv,delta=pv-mv))
    columns=['feature_id','species','gene_id','gene','route','site_role','contrast']
    ctx.table(pd.DataFrame(axes,columns=['species','gene_id','route','pas_order','atlas_pas_id','cut_interbase']).drop_duplicates(),prefix+'-PAS-axis.tsv')
    ctx.table(masks,prefix+'-PAS-masks.tsv',columns=['feature_id','species','gene_id','route','site_role','pas_order','atlas_pas_id','cut_interbase','retains_site'])
    libframe=pd.DataFrame(libraries,columns=columns+['cutoff','state','library_key','inclusion','family_UMI','coverage'])
    summary=pd.DataFrame(summaries,columns=columns+['cutoff','state','n','mean','q25','q50','q75','status','uncertainty'])
    # Where the example producer already exported a feature, require exact agreement.
    exported=ctx.read(ctx.analysis/'examples/inclusion_library_states.tsv.gz')
    common=libframe.merge(exported,on=['feature_id','cutoff','state','library_key'],suffixes=('_new','_old'),validate='one_to_one')
    np.testing.assert_allclose(common.inclusion_new,common.inclusion_old,rtol=0,atol=1e-12)
    ctx.table(libframe,prefix+'-library-inclusion.tsv.gz')
    ctx.table(summary,prefix+'-inclusion-quantiles.tsv')
    ctx.table(summary.rename(columns={'n':'n_libraries','mean':'equal_library'}),prefix+'-state-summary.tsv')
    paired=pd.DataFrame(paired,columns=columns+['cutoff','library_key','reference_inclusion','target_inclusion','delta'])
    ctx.table(paired,prefix+'-paired-inclusion.tsv.gz')
    effects=[]
    for key,q in paired.groupby(['feature_id','cutoff']):
        values=q.delta.to_numpy();n=len(values);influence=(values.sum()-values)/(n-1) if n>1 else np.array([])
        effects.append(dict(feature_id=key[0],cutoff=key[1],n_pairs=n,reference_mean=q.reference_inclusion.mean(),
            target_mean=q.target_inclusion.mean(),delta_mean=values.mean(),leave_one_min=influence.min() if len(influence) else np.nan,
            leave_one_max=influence.max() if len(influence) else np.nan,uncertainty='leave_one_library_influence_not_CI'))
    ctx.table(effects,prefix+'-paired-summary.tsv',columns=['feature_id','cutoff','n_pairs','reference_mean','target_mean','delta_mean','leave_one_min','leave_one_max','uncertainty'])


def examples(ctx,partial):
    selection=ctx.read(ctx.analysis/'candidate_screen/example_selection.tsv')
    features=ctx.read(ctx.analysis/'examples/inclusion_features.tsv')
    exact=ctx.read(ctx.analysis/'site_statistics/exact_site_cross_study.tsv.gz')
    selected,status=select_partial_sites(partial,exact,features,selection,allow_new_features=True)
    ctx.table(selected,'partial-selected-site-pairs.tsv');ctx.table(status,'partial-selection-status.tsv')
    selected_inclusion(ctx,selected)
    low=ctx.read(ctx.analysis/'candidate_screen/low_shift_candidates.tsv')
    strict=low[low.strict_region_pairs.gt(0)].sort_values(['min_MPRO','gene_id'],ascending=[False,True])
    final=strict[strict.max_abs_DWUI.lt(.01)]
    regions=ctx.read(ctx.analysis/'candidate_screen/recurrent_region_pairs.tsv.gz')
    bindings=[]
    for panel in ['D21B','D26A','D26B','D26C','D26D','D26E','D26F','D26G','D27A','D27B','D27C']:
        cohort=panel.startswith('D27') or panel=='D21B'
        genes=(strict if panel=='D27C' else low) if cohort else final.head(1)
        if genes.empty:bindings.append(dict(panel=panel,status='no_eligible_current_gene',role='low_shift_cohort' if cohort else 'strict_low_shift_example',gene_id='',gene_symbol=''))
        for r in genes.itertuples():bindings.append(dict(panel=panel,status='selected',role='low_shift_cohort' if cohort else 'strict_low_shift_example',gene_id=r.gene_id,gene_symbol=r.gene_symbol))
    gids=set(final.head(1).gene_id)
    rr=regions[regions.gene_id.isin(gids)&regions.contrast.eq('MP')].copy()
    rr=rr.sort_values(['gene_id','shared_start0','shared_end0','region_id_a','region_id_b']).drop_duplicates('gene_id')
    ctx.table(rr,'low-shift-example-interval.tsv')
    low_features=[]
    for r in rr.itertuples():
        for side in ['a','b']:
            route=getattr(r,'route_'+side);mask=np.asarray(json.loads(getattr(r,'mask_'+side)))
            assert mask.ndim==2 and np.array_equal(mask[:,0],mask[:,1])
            # All positions in a quantified elementary interval share this mask.
            low_features.append(dict(species='mouse',gene_id=r.gene_id,gene_symbol=r.gene_symbol_a,route=route,
                site_role='strict_interval',contrast='MP',start0=r.shared_start0,end0=r.shared_end0,
                feature_id=f'{route}|{r.gene_id}|interval|{r.shared_start0}-{r.shared_end0}|interval',mask=json.dumps(mask[:,0].tolist())))
    low_features=pd.DataFrame(low_features,columns=['species','gene_id','gene_symbol','route','site_role','contrast','start0','end0','feature_id','mask'])
    ctx.table(low_features,'low-shift-selected-features.tsv');selected_inclusion(ctx,low_features,'low-shift-interval')
    same=ctx.read(ctx.analysis/'axis_extensions/homolog_gene_axes.tsv.gz')
    same=same[same.tolerance.eq(25)&same.design.eq('qualified_P_M')&same.pair.eq('H1_C')&same.mass_min.ge(.8)&same.n_matched.ge(2)]
    structure=ctx.read(ctx.analysis/'geometry/exon-CDS.tsv.gz')
    annotated=set(structure.gene_id)
    same=same[same.human_gene_id.isin(annotated)&same.mouse_gene_id.isin(annotated)]
    # Prefer a gene already chosen for the dual-reference narrative, then stable gene ID order.
    ids=set(selection[selection.role.eq('dual_reference')].gene_id)
    same=same.assign(in_narrative=same.human_gene_id.isin(ids)|same.mouse_gene_id.isin(ids))
    chosen=same.sort_values(['in_narrative','human_gene_id'],ascending=[False,True]).head(1)
    for panel in ['D18G','D18H']:
        if chosen.empty:bindings.append(dict(panel=panel,status='no_eligible_current_gene',role='same_PAS_structure',gene_id='',gene_symbol=''))
        for r in chosen.itertuples():
            for side in ['human','mouse']:bindings.append(dict(panel=panel,status='selected',role='same_PAS_structure',gene_id=getattr(r,side+'_gene_id'),gene_symbol='',species=side))
    sites=ctx.read(ctx.analysis/'axis_extensions/homolog_site_pairs.tsv.gz')
    sites=sites[sites.human_gene_id.isin(chosen.human_gene_id)&sites.tolerance.eq(25)&sites.design.eq('qualified_P_M')&sites.pair.eq('H1_C')].copy()
    ctx.table(sites.sort_values(['human_gene_id','human_signed_cut']),'same-PAS-example-sites.tsv.gz')
    selected_genes=set(chosen.human_gene_id)|set(chosen.mouse_gene_id)
    ctx.table(structure[structure.gene_id.isin(selected_genes)],'same-PAS-example-exon-CDS.tsv.gz')
    checks=ctx.read(ctx.analysis/'geometry/PAS-transcript-coordinate-check.tsv.gz')
    ctx.table(checks[checks.gene_id.isin(selected_genes)],'same-PAS-example-coordinate-check.tsv.gz')
    for panel in ['D22C','D22D','D22E']:
        for r in status.itertuples():
            q=selected[selected.species.eq(r.species)].drop_duplicates('gene_id')
            if q.empty:bindings.append(dict(panel=panel,status=r.status,role='partial_retention',species=r.species,gene_id='',gene_symbol=''))
            for g in q.itertuples():bindings.append(dict(panel=panel,status='selected',role='partial_retention',species=r.species,gene_id=g.gene_id,gene_symbol=g.gene_symbol))
    ctx.table(bindings,'panel-example-bindings.tsv')


def run(ctx):
    from downstream_replicate_diagnostics import run as replicate_diagnostics
    replicate_diagnostics(ctx)
    loo_tables(ctx);display_tables(ctx)
    points,availability=dual_points(ctx.read(ctx.analysis/'matched_pms/matched_state_vectors.tsv.gz'),ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz'))
    ladder=ctx.read(ctx.analysis/'candidate_screen/dual_reference_ladder.tsv')
    for r in ladder.itertuples():
        q=availability[availability.species.eq(r.species)]
        assert len(q)==r.recurrent_candidates and q.matched_measurable.sum()==r.matched_measurable and q.dual_qualified_gene.sum()==r.dual_reference
    ctx.table(points,'dual-reference-all-measurable.tsv.gz');ctx.table(availability,'dual-reference-availability.tsv')
    partial=ctx.read(ctx.analysis/'site_statistics/cross_study_partial_retention.tsv')
    members,counts=relationship_tables(partial)
    ctx.table(members,'partial-relationship-members.tsv');ctx.table(counts,'partial-relationship-counts.tsv')
    examples(ctx,partial)
    ctx.finish(dict(current_run_only=True,LOO_replay_matches_original_statistics=True,all_measurable_dual_candidates_retained=True,
        exact_partial_relationships_deduplicated=True,inclusion_uses_current_terminal_models_and_checks_existing_exports=True,
        no_fixed_gene_or_cohort_count=True,no_v3=True,no_renderer_statistics=True,
        six_well_resolution_recomputed=True,zero_variance_CI_explicitly_unestimable=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'supplement_statistics',a.out))
