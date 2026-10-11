"""Re-evaluate named examples and expose fresh panels and external query features."""
import argparse
import json

import numpy as np
import pandas as pd

from downstream_runtime import Context,Study,STUDIES,SETTINGS

EXAMPLES=['SLC38A2','NUDT21','SZRD1','RABIF','RTF2','AP1S2','HMGN2','EIF4EBP2','RAD23B']
EXTERNAL=['HMGN2','EIF4EBP2','RAD23B']


def shared_regions(sites):
    keys=['species','contrast','gene_id','chrom','strand','region_id_a','region_id_b']
    x=sites[sites.recurrent_strict.astype(bool)].copy()
    x['shared_start0']=np.maximum(x.region_start0_a,x.region_start0_b)
    x['shared_end0']=np.minimum(x.region_end0_a,x.region_end0_b)
    assert (x.shared_start0.le(x.start0)&x.shared_end0.ge(x.end0)).all()
    x['categories']=x.groupby(keys).category.transform(lambda v:';'.join(sorted(set(v))))
    x=x.drop_duplicates(keys).copy();x['shared_length']=x.shared_end0-x.shared_start0
    assert x.shared_length.gt(0).all()
    return x


def query_features(regions,sites,external=None):
    external=set(EXTERNAL if external is None else external)
    rows=[]
    for r in regions.itertuples():
        if r.species!='mouse' or str(r.gene_symbol_a).upper() not in external:continue
        rows.append(dict(gene_id=r.gene_id,gene=r.gene_symbol_a,kind='shared_region',category=r.categories,
            contrast=r.contrast,chrom=r.chrom,strand=r.strand,start0=int(r.shared_start0),end0=int(r.shared_end0),factor=''))
    for r in sites.itertuples():
        if r.species!='mouse' or str(r.gene_symbol_a).upper() not in external or not r.recurrent_strict:continue
        rows.append(dict(gene_id=r.gene_id,gene=r.gene_symbol_a,kind='site',category=r.category,contrast=r.contrast,
            chrom=r.chrom,strand=r.strand,start0=int(r.start0),end0=int(r.end0),factor=r.factor))
    cols=['gene_id','gene','kind','category','contrast','chrom','strand','start0','end0','factor']
    f=pd.DataFrame(rows,columns=cols).drop_duplicates().sort_values(cols).reset_index(drop=True)
    f['feature_id']=['|'.join(str(r[c]) for c in cols) for r in f.to_dict('records')]
    return f


def run(ctx):
    claims=ctx.read(ctx.analysis/'union_pool/figure5_union_claims.tsv')
    selected=ctx.read(ctx.analysis/'candidate_screen/example_selection.tsv') if 'candidate_screen' in ctx.config['jobs'] else pd.DataFrame(columns=['gene_symbol'])
    examples=sorted(set(EXAMPLES)|set(selected.gene_symbol.str.upper())|{'PDLIM7'})
    external=sorted(set(EXTERNAL)|set(selected.gene_symbol.str.upper()))
    members=ctx.read(ctx.analysis/'cross_trace/member_migration.tsv')
    sites=ctx.read(ctx.analysis/'site_statistics/recurrent_sites.tsv.gz');regions=shared_regions(sites)
    ctx.table(regions,'recurrent_region_pairs.tsv.gz');features=query_features(regions,sites,external)
    ctx.table(features,'external_query_features.tsv')
    candidates=[];axes=[];metacells=[];libraries=[];states=[];effectrows=[]
    for sid in STUDIES:
        study=Study(ctx,sid);study.load_gex();symbols=study.f.primary_gene_symbol.fillna('').str.upper()
        gene_symbols=study.f.loc[symbols.isin(examples),['primary_gene_id','primary_gene_symbol']].drop_duplicates()
        assert not gene_symbols.primary_gene_id.duplicated().any()
        for name in examples:
            genes=gene_symbols[gene_symbols.primary_gene_symbol.str.upper().eq(name)]
            if genes.empty:candidates.append(dict(study_id=sid,example=name,status='absent_from_new_catalog'))
        for gene,symbol in gene_symbols.itertuples(index=False,name=None):
            sp=study.meta.species.iloc[0];cohort=claims[claims.species.eq(sp)&claims.gene_id.eq(gene)]
            current=members[members[sp+'_gene_id'].eq(gene)&members.new_member]
            candidates.append(dict(study_id=sid,species=sp,example=symbol.upper(),gene_id=gene,gene_symbol=symbol,
                status='primary_multiPAS' if gene in study.genes else 'not_primary_multiPAS',
                current_union_contrasts=';'.join(sorted(cohort.contrast.unique())),current_cross_contrasts=';'.join(sorted(current.contrast.unique()))))
            if gene not in study.genes:continue
            for setting in SETTINGS:
                axis=study.axis(gene,setting)
                if len(axis)<2:continue
                axes.append(axis.assign(setting=setting));data=study.gene_counts(axis);depth=data.sum(1)
                usage=np.divide(data,depth[:,None],out=np.zeros(data.shape,float),where=depth[:,None]>0)
                score=usage@axis.position.to_numpy();gex=study.expression(gene,np.arange(len(study.meta)))
                for i,m in study.meta.iterrows():
                    metacells.append(dict(species=sp,study_id=sid,route=m.route,gene_id=gene,gene_symbol=symbol,setting=setting,
                        metacell_id=m.metacell_id,library_key=m.library_key,node=m.node,time=m.experimental_time,n_cells=m.n_cells,purity=m.purity,
                        PAS_UMI=int(depth[i]),qualified=bool(depth[i]>=SETTINGS[setting][3] and m.purity>=SETTINGS[setting][1]),
                        DWUI=score[i] if depth[i] else np.nan,GEX_CPM=gex[i]))
                for (route,lib,node),state in study.profile_states(data,setting).items():
                    for i,p in axis.iterrows():
                        libraries.append(dict(species=sp,study_id=sid,route=route,library_key=lib,node=node,gene_id=gene,gene_symbol=symbol,
                            setting=setting,atlas_pas_id=p.atlas_pas_id,cut_interbase=p.cut_site,position=p.position,usage=state['usage'][i],
                            UMI=int(state['site_counts'][i]),detected_MC=int(state['site_mc'][i]),qualified_MC=state['n_mc'],qualified_cells=state['n_cells']))
            for setting in SETTINGS:
                for kind,target in [('node_positions',states),('route_effects',effectrows)]:
                    frame=study.core(kind,setting);target.append(frame[frame.gene_id.eq(gene)].assign(setting=setting))
    ctx.table(candidates,'example_eligibility.tsv')
    ctx.table(pd.concat(axes,ignore_index=True) if axes else [],'example_PAS_axes.tsv.gz',columns=['study_id','primary_gene_id','setting','atlas_pas_id'])
    ctx.table(metacells,'example_metacell_scores.tsv.gz',columns=['study_id','gene_id','metacell_id','setting','DWUI','qualified'])
    ctx.table(libraries,'example_library_PAS_profiles.tsv.gz',columns=['study_id','gene_id','setting','route','library_key','node','atlas_pas_id','usage'])
    ctx.table(pd.concat(states,ignore_index=True) if states else [],'example_state_profiles.tsv.gz',columns=['study_id','gene_id','setting','node','position_index'])
    ctx.table(pd.concat(effectrows,ignore_index=True) if effectrows else [],'example_endpoint_effects.tsv.gz',columns=['study_id','gene_id','setting','contrast'])
    ctx.table([dict(example=g,status='new_features_available' if features.gene.str.upper().eq(g).any() else 'no_new_recurrent_feature',
                    queries=int(features.gene.str.upper().eq(g).sum())) for g in external],'external_feature_eligibility.tsv')
    from downstream_example_inclusion import run as inclusion_audit
    inclusion_audit(ctx)
    ctx.finish(dict(all_nine_existing_examples_audited=True,new_PAS_and_metacell_profiles=True,
        no_old_intervals_reused=True,missing_examples_explicit=True,new_shared_regions_and_sites=True,
        both_current_selection_and_historical_example_identity_preserved=True,
        native_library_inclusion_and_influence_rebuilt=True,all_ten_historical_loci_reaudited=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'examples',a.out))
