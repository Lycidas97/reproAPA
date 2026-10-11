"""Re-query public GSE99922 counts against the native polyApipe PAS catalogue."""
import argparse
from collections import defaultdict
from itertools import combinations,product
import json

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import spearmanr

from downstream_runtime import Context,Study
from downstream_liftover import lift_points

SOURCE='notes/2026-10-02-gse99922-pilot/'
ENDPOINT=['MEF_R1','MEF_R2','iPSC_R1','iPSC_R2']
METRICS=['tvd','mpro','aligned_shared_axis','aligned_bp','direction_concordance']


def endpoint_samples(manifest,columns):
    # The public table also contains perturbation/sorted-cell columns outside
    # the verified unperturbed sample manifest. Only these four endpoints enter.
    assert manifest.count_column.is_unique, 'Duplicate sample identity'
    assert set(manifest.count_column)<=set(columns), 'Manifest sample absent from counts'
    assert set(ENDPOINT)<=set(manifest.count_column), 'Endpoint sample absent from manifest'
    endpoints=manifest.set_index('count_column').loc[ENDPOINT]
    assert endpoints.group.tolist()==['MEF','MEF','iPSC','iPSC']
    assert endpoints.gsm.is_unique and endpoints.sra_run.is_unique
    return endpoints


def usable_matches(external,reference):
    rows=[]
    for r in external.itertuples():
        out=dict(COORDS=r.COORDS,atlas_pas_id='',distance_nt=np.nan,compatible50=False)
        f=reference[reference.primary_gene_id.eq(r.gene_id)]
        f=f[f.chrom.eq(r.mapped_chrom)&f.strand.eq(r.orientation)]
        if len(f) and not r.coordinate_base_sensitive and r.source_strand==r.orientation:
            d=f.cut_site.sub(r.mapped_cut).abs();a=f.cut_site.sub(r.alternate_cut).abs()
            out['distance_nt']=d.min()
            if d.min()<=50 and a.min()<=50 and d.eq(d.min()).sum()==a.eq(a.min()).sum()==1 and d.idxmin()==a.idxmin():
                out.update(atlas_pas_id=f.loc[d.idxmin(),'atlas_pas_id'],compatible50=True)
        rows.append(out)
    result=external.merge(pd.DataFrame(rows,columns=['COORDS','atlas_pas_id','distance_nt','compatible50']),on='COORDS',validate='one_to_one')
    result['many_to_one_collision']=result.compatible50&result.atlas_pas_id.duplicated(False)
    result['usable']=result.compatible50&~result.many_to_one_collision
    return result


def measurements(ctx,math):
    parts=[]
    for sid in ['GSE103221','GSE115943','GSE171127']:
        folder='WOT_GSE115943' if sid=='GSE115943' else sid
        f=ctx.read(ctx.native/f'studies/{folder}/features.tsv.gz',low_memory=False)
        parts.append(f[f.keep_for_apa_primary.astype(str).str.lower().eq('true')][['atlas_pas_id','primary_gene_id','primary_gene_symbol','chrom','strand','cut_site']])
    ref=pd.concat(parts).drop_duplicates();assert ref.atlas_pas_id.is_unique
    # The same genomic cut detected in two studies is one external matching target.
    keys=['primary_gene_id','chrom','strand','cut_site']
    ctx.table(ref,'native_reference_members.tsv.gz')
    ref=ref.sort_values('atlas_pas_id').drop_duplicates(keys).reset_index(drop=True)
    geometry=ref[['primary_gene_id','primary_gene_symbol','chrom','strand']].drop_duplicates()
    geometry=geometry[geometry.groupby('primary_gene_id').primary_gene_id.transform('size').eq(1)]
    symbols=geometry[geometry.groupby('primary_gene_symbol').primary_gene_id.transform('size').eq(1)].set_index('primary_gene_symbol').primary_gene_id.to_dict()
    geometry=geometry.set_index('primary_gene_id').to_dict('index')
    info=ctx.read(ctx.resource('external_mouse_gene_info'),dtype=str,keep_default_na=False)
    xref,reverse=defaultdict(set),defaultdict(set)
    for row in info[info['#tax_id'].eq('10090')].itertuples():
        for token in row.dbXrefs.split('|'):
            if token.startswith('Ensembl:ENSMUSG'):
                gene=token.split(':',1)[1];xref[row.GeneID].add(gene);reverse[gene].add(row.GeneID)
    raw=ctx.read(ctx.resource('external_pas_counts'),dtype={'IDs':str});columns=list(raw.columns[5:])
    manifest=ctx.read(ctx.resource('external_sample_manifest'),dtype=str)
    endpoint_manifest=endpoint_samples(manifest,columns)
    ctx.table(manifest,'public_sample_manifest.tsv');ctx.table(endpoint_manifest.reset_index(),'endpoint_sample_manifest.tsv')
    ctx.table(pd.DataFrame(dict(count_column=columns,
        in_verified_manifest=[c in set(manifest.count_column) for c in columns],
        used_in_endpoint_comparison=[c in ENDPOINT for c in columns])),'public_count_column_audit.tsv')
    assert raw[columns].ge(0).all().all() and np.equal(raw[columns],np.floor(raw[columns])).all().all()
    assert raw[raw.COORDS.duplicated(False)].groupby('COORDS')[columns].nunique().le(1).all().all()
    ids=[]
    for r in raw[['IDs','GENES']].drop_duplicates().itertuples():
        gene,method,status=math.resolve_gene(r.IDs,r.GENES,xref,reverse,geometry,symbols)
        ids.append(dict(IDs=r.IDs,GENES=r.GENES,gene_id=gene,identity_method=method,identity_status=status))
    ids=pd.DataFrame(ids);ctx.table(ids,'gene_identity_audit.tsv.gz')
    conflict=raw.groupby('COORDS').agg(genes=('GENES','nunique'),types=('TYPES','nunique'))
    keep=conflict.index[conflict.genes.eq(1)&conflict.types.eq(1)]
    sites=raw[raw.COORDS.isin(keep)&raw.TYPES.eq('LastExon')].drop_duplicates('COORDS').merge(ids,on=['IDs','GENES'],validate='many_to_one')
    sites=sites[sites.identity_status.eq('mapped')].copy()
    sites[['chrom','source_strand','pos']]=sites.COORDS.str.extract(r'^(chr.+)([+-])(\d+)$');assert sites.pos.notna().all()
    sites['pos']=sites.pos.astype(int);sites['gene_symbol']=sites.gene_id.map(lambda g:geometry[g]['primary_gene_symbol'])
    for offset in [1,0]:
        q=pd.DataFrame(dict(id=sites.COORDS,chrom=sites.chrom,strand=sites.source_strand,
            cut_site=sites.pos-offset+sites.source_strand.eq('+').astype(int)))
        for a,b in [('mm9','mm10'),('mm10','mm39')]:
            m=lift_points(ctx,q,'chain_'+a+'_'+b,f'offset{offset}_{a}_{b}')
            q=m[m.mapping_status.eq('unique')][['id','mapped_chrom','mapped_strand','mapped_cut']].rename(columns=dict(mapped_chrom='chrom',mapped_strand='strand',mapped_cut='cut_site'))
        prefix='mapped' if offset else 'alternate'
        q=q.rename(columns=dict(id='COORDS',chrom=prefix+'_chrom',strand='orientation' if offset else 'alternate_orientation',cut_site=prefix+'_cut'))
        sites=sites.merge(q,on='COORDS',how='left',validate='one_to_one')
    sites['coordinate_base_sensitive']=(sites.mapped_chrom.ne(sites.alternate_chrom)|sites.orientation.ne(sites.alternate_orientation)|sites.mapped_cut.sub(sites.alternate_cut).abs().gt(1)|sites.mapped_cut.isna()|sites.alternate_cut.isna())
    mapped=usable_matches(sites,ref);ctx.table(mapped,'expanded_site_mapping.tsv.gz')
    depths=raw.drop_duplicates('COORDS')[ENDPOINT].sum().to_numpy();records=[];vectors=[]
    for gene,all_sites in mapped.groupby('gene_id'):
        shared=all_sites[all_sites.usable].copy();x=shared.mapped_cut.to_numpy()*np.where(shared.orientation.eq('+'),1,-1)
        x=x-x.min() if len(x) else x;counts=shared[ENDPOINT].to_numpy().T;totals=counts.sum(1);public=all_sites[ENDPOINT].sum().to_numpy()
        fractions=np.divide(totals,public,out=np.zeros(4),where=public>0)
        p=np.divide(counts,totals[:,None],out=np.full(counts.shape,np.nan),where=totals[:,None]>0)
        span=np.ptp(x) if len(x) else 0;support=int((counts.sum(0)>=10).sum());overlap=len(shared)>1 and shared.pos.sort_values().diff().min()<=80
        values=math.profile_metrics(p,x) if len(x) else np.full(4,np.nan);pooled=p.mean(0)
        row=dict(gene_id=gene,gene_symbol=all_sites.gene_symbol.iloc[0],n_external_last_exon_sites=len(all_sites),n_shared_sites=len(shared),shared_span_bp=span,
            min_shared_count=totals.min(),min_shared_fraction=fractions.min(),overlapping_shared_windows=overlap,n_collision_sites=int(all_sites.many_to_one_collision.sum()),
            mean_log2_last_exon_CPM=float(np.log2(public/depths*1e6+.5).mean()),log2_min_shared_count=np.log2(max(totals.min(),1)),
            log2_external_last_exon_sites=np.log2(len(all_sites)),log2_shared_span_bp=np.log2(max(span,1)),
            pooled_shared_axis_usage=pooled@x/span if span else np.nan,pooled_proportion_impurity=1-np.square(pooled).sum(),
            pas_stratum='2' if len(shared)==2 else '3-4' if len(shared)<=4 else '5+',**dict(zip(['tvd','mpro','delta_shared_axis','delta_bp'],values)))
        for threshold in [20,50,100]:row['measurable'+str(threshold)]=len(shared)>=2 and support>=2 and span>=100 and not overlap and totals.min()>=threshold and fractions.min()>=.8
        for i,sample in enumerate(ENDPOINT):shared['usage__'+sample]=p[i]
        shared['local_cut_bp']=x;records.append(row);vectors.append(shared)
    genes=pd.DataFrame(records);vectors=pd.concat(vectors,ignore_index=True)
    ctx.table(genes,'all_gene_measurements.tsv.gz');ctx.table(vectors,'shared_pas_vectors.tsv.gz')
    return genes.set_index('gene_id'),vectors,mapped


def matched(ctx,math,genes,vectors,roster):
    protocol=json.loads(ctx.resource('external_vector_protocol').read_text());covariates=protocol['matching_covariates']
    matching=ctx.functions(SOURCE+'matched_analysis.py',['match_cases'],dict(linear_sum_assignment=linear_sum_assignment))
    statistic=ctx.functions(SOURCE+'vector_analysis.py',['statistics'])
    profiles={g:(v[['usage__'+s for s in ENDPOINT]].to_numpy().T,v.local_cut_bp.to_numpy()) for g,v in vectors[vectors.gene_id.isin(genes.index[genes.measurable20])].groupby('gene_id')}
    assignments=list(combinations(range(4),2));effects={s:pd.DataFrame({g:math.profile_metrics(p,x,s) for g,(p,x) in profiles.items()}).T for s in assignments}
    rows=[];members=[];balances=[];statuses=[];nulls=[];sens=[]
    for name,threshold,k in [('primary',20,3),('count50',50,3),('count100',100,3),('controls2',20,2)]:
        eligible=genes[genes['measurable'+str(threshold)]];cases=eligible[eligible.index.isin(roster.gene_id)];controls=eligible[~eligible.index.isin(roster.gene_id)]
        groups,balance,status=matching.match_cases(cases[covariates+['pas_stratum']],controls[covariates+['pas_stratum']],k,protocol,covariates)
        statuses.append(roster.merge(status,on='gene_id',how='left',validate='one_to_one').assign(configuration=name))
        if len(balance):balances.append(balance.assign(configuration=name))
        if not groups:
            rows.extend(dict(configuration=name,statistic=m,status='no_complete_matched_sets',candidate_n=len(roster),measurable_n=len(cases),matched_n=0,estimate=np.nan,whole_library_upper_tail_diagnostic=np.nan) for m in METRICS)
            continue
        signs=roster.set_index('gene_id').loc[[g[0] for g in groups],'direction_sign'].to_numpy()
        observed,individual,_=statistic.statistics(effects[(0,1)],groups,signs)
        null=[];weights=[]
        for source in assignments:
            value,*_=statistic.statistics(effects[source],groups,signs);null.append(value)
            nulls.extend(dict(configuration=name,source_libraries=';'.join(ENDPOINT[i] for i in source),statistic=m,value=v) for m,v in zip(METRICS,value))
        for sw,tw in product([0.,.5,1.],repeat=2):
            values={}
            for gene in {g for group in groups for g in group}:
                p,x=profiles[gene];d=tw*p[2]+(1-tw)*p[3]-sw*p[0]-(1-sw)*p[1]
                values[gene]=[abs(d).sum()/2,np.ptp(d),d@x/np.ptp(x),d@x]
            value,*_=statistic.statistics(pd.DataFrame(values).T,groups,signs);weights.append(value)
            sens.extend(dict(configuration=name,source_weight=sw,target_weight=tw,statistic=m,value=v) for m,v in zip(METRICS,value))
        null=np.asarray(null);weights=np.asarray(weights);retention=len(groups)/len(cases);smd=balance.smd_after.abs().max()
        for i,m in enumerate(METRICS):rows.append(dict(configuration=name,statistic=m,status='descriptive_whole_library_diagnostic',candidate_n=len(roster),measurable_n=len(cases),matched_n=len(groups),controls_n=len(groups)*k,
            retention=retention,max_abs_smd=smd,design_gate_pass=retention>=.8 and smd<=.2,estimate=observed[i],
            whole_library_upper_tail_diagnostic=float((null[:,i]>=observed[i]-1e-12).mean()),label_assignments=6,sensitivity_min=weights[:,i].min(),sensitivity_max=weights[:,i].max()))
        for group,direction in zip(groups,signs):members.extend(dict(configuration=name,case_gene=group[0],gene_id=g,role='candidate' if j==0 else 'control',direction_sign=direction) for j,g in enumerate(group))
    for data,path,cols in [(rows,'vector_summary.tsv',['configuration','statistic','status']),
        (members,'vector_matched_members.tsv',['configuration','case_gene','gene_id','role']),
        (pd.concat(balances) if balances else [],'vector_matching_balance.tsv',['configuration','covariate','smd_after']),
        (pd.concat(statuses),'vector_candidate_status.tsv',['gene_id','configuration']),
        (nulls,'whole_library_label_assignments.tsv',['configuration','statistic','value']),
        (sens,'whole_library_sensitivity.tsv',['configuration','statistic','value'])]:ctx.table(data,path,columns=cols)


def direct(ctx,mapped,roster):
    rows=[];site_rows=[]
    for sid,route in [('GSE103221','CIP'),('GSE115943','WOT')]:
        study=Study(ctx,sid);design=study.design.set_index(['library_key','state'])
        for candidate in roster.itertuples():
            gid=candidate.gene_id;r=dict(study_id=sid,route=route,gene_id=gid,gene_symbol=candidate.gene_symbol,status='no_native_or_external_axis')
            ext=mapped[mapped.gene_id.eq(gid)]
            if gid not in study.genes or ext.empty:rows.append(r);continue
            axis=study.axis(gid);counts=study.gene_counts(axis);states=study.profile_states(counts);records=[]
            for (rr,lib,node),s in states.items():
                if rr!=route or node not in ['S','P']:continue
                meta=design.loc[lib,node]
                if node=='S' and meta.time!='D0':continue
                records.append(dict(node=node,time=meta.time,usage=s['usage']))
            if not {'S','P'}<={x['node'] for x in records}:r['status']='no_qualified_source_or_target';rows.append(r);continue
            m=usable_matches(ext.drop(columns=['atlas_pas_id','distance_nt','compatible50','many_to_one_collision','usable']),axis)
            shared=m[m.usable];r['n_common_pas']=len(shared)
            if len(shared)<2:r['status']='fewer_than_two_common_pas';rows.append(r);continue
            ix=pd.Index(axis.atlas_pas_id).get_indexer(shared.atlas_pas_id);assert (ix>=0).all()
            x=axis.cut_site.to_numpy()[ix]*np.where(axis.strand.to_numpy()[ix]=='+',1,-1);x=x-x.min();span=np.ptp(x)
            c=shared[ENDPOINT].to_numpy().T;total=c.sum(1);public=ext[ENDPOINT].sum().to_numpy();fraction=np.divide(total,public,out=np.zeros(4),where=public>0)
            native_fraction=min(p['usage'][ix].sum() for p in records)
            r.update(shared_span_bp=span,min_external_common_fraction=fraction.min(),min_native_common_fraction=native_fraction,min_external_common_count=total.min())
            if span<100 or shared.pos.sort_values().diff().min()<=80 or (c.sum(0)>=10).sum()<2 or total.min()<20:r['status']='external_coverage_or_windows_fail';rows.append(r);continue
            if fraction.min()<.8 or native_fraction<.8:r['status']='common_mass_below_80percent';rows.append(r);continue
            def mean(node):
                p=[p for p in records if p['node']==node]
                u=np.mean([np.mean([r['usage'] for r in p if r['time']==t],axis=0) for t in sorted({r['time'] for r in p})],axis=0)[ix]
                return u/u.sum()
            source,target=mean('S'),mean('P');p=c/total[:,None];delta=target-source;external=p[2:].mean(0)-p[:2].mean(0);norm=np.linalg.norm(delta)*np.linalg.norm(external)
            r.update(status='comparable',native_delta_axis=delta@x/span,external_delta_axis=external@x/span,native_tvd=abs(delta).sum()/2,external_tvd=abs(external).sum()/2,pas_change_cosine=delta@external/norm if norm else np.nan,same_position_direction=bool((delta@x)*(external@x)>0));rows.append(r)
            site_rows.extend(dict(study_id=sid,gene_id=gid,COORDS=coord,atlas_pas_id=pas,native_usage_change=a,external_usage_change=b) for coord,pas,a,b in zip(shared.COORDS,shared.atlas_pas_id,delta,external))
    ctx.table(rows,'native_comparison_candidates.tsv',columns=['study_id','gene_id','status']);ctx.table(site_rows,'native_comparison_pas.tsv',columns=['study_id','gene_id','COORDS','atlas_pas_id','native_usage_change','external_usage_change'])
    summary=[]
    for sid,g in (pd.DataFrame(rows) if rows else pd.DataFrame(columns=['study_id','status'])).groupby('study_id'):
        q=g[g.status.eq('comparable')];summary.append(dict(study_id=sid,candidate_n=len(g),comparable_n=len(q),same_direction_n=int(q.same_position_direction.sum()) if len(q) else 0,
            descriptive_spearman=spearmanr(q.native_delta_axis,q.external_delta_axis).statistic if len(q)>=3 else np.nan))
    ctx.table(summary,'native_comparison_summary.tsv',columns=['study_id','candidate_n','comparable_n','same_direction_n','descriptive_spearman'])


def run(ctx):
    math=ctx.functions(SOURCE+'expanded_measurements.py',['resolve_gene','profile_metrics'])
    genes,vectors,mapped=measurements(ctx,math)
    rec=ctx.read(ctx.analysis/'recurrence/recurrence_members.tsv.gz')
    roster=rec[rec.species.eq('mouse')&rec.contrast.eq('P_S')&rec.metric.eq('WPP')&rec.threshold.eq(.02)&rec.series_matching][['gene_id','gene_symbol','direction']].drop_duplicates()
    assert roster.gene_id.is_unique
    roster['direction_sign']=roster.direction.map({'negative':-1,'positive':1});assert roster.direction_sign.notna().all()
    ctx.table(roster,'native_candidate_roster.tsv');matched(ctx,math,genes,vectors,roster);direct(ctx,mapped,roster)
    old=ctx.read(ctx.resource('external_legacy_candidates'));migration=old.merge(roster,on='gene_id',how='outer',suffixes=('_legacy','_native'),indicator=True)
    migration['scope']='historical_roster_tracking_only; current matching uses native DWUI PS roster'
    ctx.table(migration,'historical_roster_migration.tsv')
    ctx.finish(dict(new_native_catalogue_matched=True,public_integer_counts_reused=True,both_coordinate_origins_checked=True,
        unique_gene_site_mapping=True,new_recurrence_roster=True,outcome_blind_disjoint_controls=True,six_whole_library_assignments=True,
        no_confirmatory_gene_as_replicate_p_values=True,new_common_PAS_comparisons=True,public_sample_roles_verified_against_pinned_manifest=True),
        reference_independence='Sample-level reference provenance unresolved; external effect measurements only',legacy_six_rule='historical audit only, not native roster')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'external_pas',a.out))
