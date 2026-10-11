"""Compare new native axes to the frozen old pipeline, with matched contributors."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from downstream_runtime import Context,Study,STUDIES,KEY

SOURCE='notes/2026-10-02-scape-polyapipe-comparison/scripts/'
PAIR=['library_key_a','library_key_b','time_a','time_b','sample_unit_a','sample_unit_b']
SITE=['site_id','chrom','strand','cut_site','primary_axis_position']


def nearest(query,reference):
    """Gene-constrained distances; ties remain measurable distances but never count mappings."""
    rows=[];groups={k:f for k,f in reference.groupby(['primary_gene_id','chrom','strand'])}
    for key,f in query.groupby(['primary_gene_id','chrom','strand'],dropna=False):
        ref=groups.get(key)
        if ref is None:
            rows.extend(dict(site_id=r.atlas_pas_id,nearest_id=None,distance=np.nan,tie=False,status='no_same_gene_reference') for r in f.itertuples());continue
        ref=ref.sort_values(['cut_site','atlas_pas_id']);cuts=ref.cut_site.to_numpy();ids=ref.atlas_pas_id.to_numpy()
        for r in f.itertuples():
            j=np.searchsorted(cuts,r.cut_site);choices=np.unique(np.clip([j-1,j],0,len(cuts)-1));distance=np.abs(cuts[choices]-r.cut_site)
            d=distance.min();ii=choices[distance==d];tied=(np.abs(cuts-r.cut_site)==d).sum()>1
            rows.append(dict(site_id=r.atlas_pas_id,nearest_id=ids[ii[0]] if not tied else None,distance=d,tie=bool(tied),status='nearest_tie' if tied else 'unique_nearest'))
    return pd.DataFrame(rows,columns=['site_id','nearest_id','distance','tie','status'])


def project_old(old,counts,new):
    match=nearest(old,new);ok=match.distance.le(50)&~match.tie
    chosen=match[ok];ids=set(chosen.nearest_id);shared=new[new.atlas_pas_id.isin(ids)].reset_index(drop=True).copy()
    source=pd.Index(old.atlas_pas_id).get_indexer(chosen.site_id);target=pd.Index(shared.atlas_pas_id).get_indexer(chosen.nearest_id)
    mat=sparse.csr_matrix((np.ones(len(chosen),np.int64),(source,target)),shape=(len(old),len(shared)))
    projected=(counts@mat).tocsr()
    assert int(projected.sum())==int(counts[:,source].sum())
    return shared,projected,match


def numerics(ctx):
    funcs=ctx.functions(SOURCE+'atlas_effect_math.py',['hierarchical_weights','state_counts','pair_effect','fixed_original_pairs','collapse_pair_effects'],dict(GROUP=KEY,PAIR=PAIR))
    return ctx.functions(SOURCE+'atlas_effect_engine.py',['effect_class','serializable','gene_effects','agreement_summary'],
        dict(vars(funcs),GROUP=KEY,PAIR=PAIR,VIEWS=['native_terminal_sites','matched_sites_independent_metacells','matched_sites_common_metacells'],
            METRICS=['delta_primary_axis','delta_own_axis','delta_nt','MPRO']))


def compare(ctx,studies=None):
    root=Path(ctx.config['caller_comparison_root']);baseline=root/'comparison_inputs'
    manifest=json.loads(ctx.resource('caller_original_manifest').read_text());verified=json.loads(ctx.resource('caller_source_hashes').read_text())
    frozen=ctx.read(ctx.resource('caller_frozen_pairs'));engine=numerics(ctx)
    distance_summary=[];all_routes=[];all_background=[];availability=[];resource_audit=[]
    for sid in STUDIES if studies is None else studies:
        study=Study(ctx,sid);atlas=study.root.name;folder=baseline/atlas
        entry=next(r for r in manifest['studies'] if r['scientific_study_id']==sid)
        path=ctx.register(entry['source_counts_path']);assert ctx.inputs[str(path)]==verified[str(path)]['sha256']
        meta=ctx.read(folder/'all_metacells.tsv.gz');ix=pd.Index(meta.metacell_id).get_indexer(study.meta.metacell_id)
        assert (ix>=0).all() and np.array_equal(meta.n_cells.to_numpy()[ix],study.meta.n_cells)
        old=ctx.read(folder/'complete_primary_PAS_axes.tsv.gz',low_memory=False)
        assert old.atlas_pas_id.is_unique and (old.peak_end-old.peak_start).eq(1).all() and old.cut_site.eq(old.peak_end).all()
        old['recorded_retained_base_1based']=old.cut_site
        old['cut_site']=old.cut_site-np.where(old.strand.eq('-'),1,0)
        original=sparse.load_npz(path).tocsr()[ix].tocsc()[:,old.feature_index_1based.to_numpy()-1]
        native=study.primary[study.primary.matched_source_count_50nt.ge(1)].copy().reset_index(drop=True)
        native['site_id']=native.atlas_pas_id;native['primary_axis_position']=native.position
        old_supported=old[old.matched_source_count_50nt.ge(1)].copy();old_counts=original[:,np.flatnonzero(old.matched_source_count_50nt.ge(1))]
        shared,legacy_shared,mapping=project_old(old_supported,old_counts,native)
        ctx.table(mapping,f'{atlas}/legacy_to_native_count_mapping.tsv.gz')
        ctx.table(shared,f'{atlas}/common_PAS_axes.tsv.gz')
        poly_shared=study.counts[:,shared.feature_index_1based.to_numpy()-1].tocsr()
        for query,ref,label in [(study.f,old,'native_all_to_legacy_primary'),(old,study.f,'legacy_primary_to_native_all'),
            (native,old_supported,'native_supported_to_legacy_supported'),(old_supported,native,'legacy_supported_to_native_supported')]:
            match=nearest(query,ref);ctx.table(match,f'{atlas}/{label}_nearest.tsv.gz')
            for threshold in [25,50,100]:distance_summary.append(dict(study_id=sid,direction=label,tolerance_nt=threshold,
                full_catalogue_denominator=len(query),gene_assigned=int(query.primary_gene_id.notna().sum()),matched=int(match.distance.le(threshold).sum()),
                unmatched=int(match.distance.gt(threshold).sum()+match.distance.isna().sum()),ties=int((match.distance.le(threshold)&match.tie).sum())))
        native_groups=native.groupby('primary_gene_id').groups;common_groups=shared.groupby('primary_gene_id').groups
        old_pairs=frozen[frozen.study_id.eq(sid)];frozen_groups=dict(tuple(old_pairs.groupby('gene_id')))
        symbols=dict(zip(old.primary_gene_id,old.primary_gene_symbol));genes=sorted(set(old.primary_gene_id)|set(study.primary.primary_gene_id))
        resource_audit.append(dict(study_id=sid,metacells=len(study.meta),old_selected_primary_molecules=int(original.sum()),
            old_common_molecules=int(legacy_shared.sum()),native_common_molecules=int(poly_shared.sum()),common_PAS=len(shared)))
        out=ctx.out/atlas;out.mkdir(exist_ok=True)
        with gzip.open(out/'pair_effects_and_missing.jsonl.gz','wt',compresslevel=3) as handle:
            for i,gene in enumerate(genes):
                ni=np.asarray(native_groups.get(gene,[]),int);ci=np.asarray(common_groups.get(gene,[]),int);old_gene=frozen_groups.get(gene,old_pairs.iloc[:0])
                all_background.append(dict(study_id=sid,gene_id=gene,in_legacy_primary=gene in symbols,in_new_primary=gene in study.genes,
                    new_supported_PAS=len(ni),common_PAS=len(ci),legacy_endpoint_pairs=len(old_gene)))
                if study.pairs.empty:continue
                pp,rr=engine.gene_effects(gene,symbols.get(gene,gene),native.loc[ni,SITE],
                    study.counts[:,native.loc[ni,'feature_index_1based'].to_numpy()-1].tocsr(),shared.loc[ci,SITE],
                    poly_shared[:,ci],legacy_shared[:,ci],study.meta,study.pairs,old_gene,gene in symbols)
                # Preserve contributor evidence without duplicating full sparse matrices.
                for row in pp:handle.write(json.dumps(engine.serializable(row),allow_nan=False,separators=(',',':'))+'\n')
                all_routes.extend(rr)
                if i%500==0:print('caller_compare',sid,i,len(genes),flush=True)
        availability.append(dict(study_id=sid,caller='polyApipe_native',status='recomputed_current_annotation_and_counts',native_result=True))
        scape=root/'scape/study_comparison'/atlas/'union.complete.json'
        availability.append(dict(study_id=sid,caller='SCAPE',status='available_not_in_native_rerun' if scape.exists() else 'not_run_no_complete_study_union',native_result=False))
    routes=pd.DataFrame(all_routes) if all_routes else pd.DataFrame(columns=KEY+['view','contributors','delta_primary_axis'])
    for col in ['site_ids','delta_usage','repro_delta_usage','contributing_library_pairs','pair_weights']:
        if col in routes:routes[col]=routes[col].map(lambda x:json.dumps(engine.serializable(x),allow_nan=False))
    ctx.table(routes,'route_effects_with_full_missing.tsv.gz')
    ctx.table(engine.agreement_summary(routes),'effect_agreement_summary.tsv',
        columns=['species','study_id','route','contrast','design','view','contributors','comparison'])
    ctx.table(distance_summary,'bidirectional_PAS_distances.tsv');ctx.table(all_background,'complete_gene_background.tsv.gz')
    ctx.table(availability,'caller_availability.tsv');ctx.table(resource_audit,'integer_count_and_membership_audit.tsv')
    if studies is None:migration(ctx,routes)
    ctx.finish(dict(new_native_annotation_and_counts=True,old_counts_reused_no_BAM_recount=True,
        full_gene_background_retained=True,common_PAS_MC_pairs_requalified=True,unchanged_GEX_membership=True,
        unavailable_SCAPE_explicit=True,legacy_negative_strand_interbase_conversion=True),
        comparison_scope='polyApipe native rerun versus historical pipeline; SCAPE remains explicitly reserved until separate study unions complete',
        old_coordinate='legacy one-base [peak_start,peak_end) retained base; cut_site=end plus/start minus',
        common_axis='new full-primary normalized positions; legacy supported counts uniquely nearest within 50nt, ties excluded, many-to-one summed')


def migration(ctx,routes):
    old=ctx.read(ctx.resource('caller_legacy_recurrence'))
    if 'level' in old:old=old[old.level.eq('study')]
    ctx.table(old,'legacy_recurrence_input.tsv.gz')
    # Full old recurrent rosters are reported independently of newly selected cohorts.
    recurrent=ctx.read(ctx.resource('caller_legacy_members'));recurrent=recurrent[recurrent.metric.eq('WPP')].copy()
    ctx.table(recurrent,'legacy_members_input.tsv.gz')
    keys=['species','study_id','route','gene_id','contrast','design']
    old_effect=ctx.read(ctx.resource('caller_frozen_routes'))
    current=ctx.read(ctx.analysis/'state_vectors/route_effects.primary.tsv.gz')
    joined=old_effect.merge(current[keys+['delta_DWUI']],on=keys,how='left',validate='one_to_one',suffixes=('_legacy','_new'))
    for threshold in [.02,.05]:
        x=joined.delta_DWUI.to_numpy();y=joined.delta_WPP.to_numpy();direction=np.sign(y)
        joined['migration_'+str(threshold)]=np.select([~np.isfinite(x),x==0,np.sign(x)!=direction,abs(x)>=threshold],
            ['not_measurable','zero','reversed','same_direction_above_threshold'],default='same_direction_below_threshold')
    ctx.table(joined,'legacy_all_route_effect_migration.tsv.gz')
    members=recurrent[['species','contrast','gene_id','threshold','direction']].rename(columns={'threshold':'legacy_threshold','direction':'legacy_recurrent_direction'})
    migrated=members.merge(joined,on=['species','contrast','gene_id'],how='left',validate='many_to_many')
    ctx.table(migrated,'complete_legacy_recurrent_gene_migration.tsv.gz')
    from downstream_examples import EXAMPLES
    examples=set(EXAMPLES)|{'MCM6','BDP1','TACC2'}
    ctx.table(routes[routes.gene_symbol.str.upper().isin(examples)],'main_example_caller_support.tsv.gz')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    p.add_argument('--study',choices=STUDIES)
    a=p.parse_args();compare(Context(a.config,'caller_comparison' if a.study is None else 'caller_comparison_'+a.study,a.out),
        None if a.study is None else [a.study])
