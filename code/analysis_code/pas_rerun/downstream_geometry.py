"""Current low-shift cohort geometry and structure, without a fixed gene count."""
import argparse
import gzip
import json
import re
import hashlib
from pathlib import Path
import numpy as np
import pandas as pd
from downstream_runtime import Context


def decompose(ctx):
    inputs={};read=ctx.read
    def write(name,frame):ctx.table(frame,name)
    def js(path):return json.loads(ctx.register(path).read_text())
    cohort=read(ctx.analysis/'candidate_screen/low_shift_candidates.tsv')
    if cohort.empty:
        for name,columns in [('axis-decomposition.tsv',['gene_id','route','family_status']),('native-sites.tsv.gz',['gene_id','route','pas_id']),
          ('native-contributors.tsv.gz',['gene_id','route']),('interval-eligibility.tsv',['gene_id','route','status']),('model-depth.tsv',['gene_id','route'])]:
            ctx.table([],name,columns=columns)
        ctx.table(cohort,'low-shift-cohort.tsv');return
    effects=read(ctx.analysis/'state_vectors/route_effects.primary.tsv.gz')
    q=effects[effects.setting.eq('supported50_p08')&effects.route.isin(['CIP','WOT'])&effects.contrast.eq('P_M')&effects.design.eq('same_library')&effects.gene_id.isin(cohort.gene_id)]
    assert len(q)==2*len(cohort), 'Every selected low-shift gene must have both routes'
    axes=read(ctx.analysis/'state_vectors/PAS_axes.primary.tsv.gz',usecols=['study_id','setting','primary_gene_id','atlas_pas_id','cut_site','strand','primary_region_class','position','native_span'])
    eligibility=read(ctx.analysis/'intervals/eligibility.tsv'); mfx=read(ctx.analysis/'intervals/model_effects.tsv')
    models={r:js(ctx.analysis/f'intervals/annotation/mouse/{r}/models.json') for r in ['CIP','WOT']}
    results=[];sites=[]
    for a in q.itertuples():
        c=cohort.set_index('gene_id').loc[a.gene_id];ids=json.loads(a.PAS_ids);P=np.array(json.loads(a.target_usage));M=np.array(json.loads(a.reference_usage))
        ax=axes[axes.study_id.eq(a.study_id)&axes.setting.eq(a.setting)&axes.primary_gene_id.eq(a.gene_id)].set_index('atlas_pas_id').loc[ids]
        cuts=ax.cut_site.to_numpy();strand=ax.strand.iloc[0];signed=cuts*(1 if strand=='+' else -1);span=float(ax.native_span.iloc[0]);assert ax.native_span.nunique()==1
        np.testing.assert_allclose([P.sum(),M.sum()],[1,1],atol=1e-12)
        np.testing.assert_allclose((P-M)@((signed-signed.min())/span),a.delta_DWUI,atol=1e-12)
        np.testing.assert_allclose(a.delta_DWUI,c['delta_DWUI_'+a.route],atol=1e-12)
        row=dict(gene_id=a.gene_id,gene_symbol=a.gene_symbol,route=a.route,study_id=a.study_id,native_pas=len(ids),axis_span_nt=span,supported_pas_span_nt=float(np.ptp(signed)),n_native_pairs=a.n_pairs,native_delta_DWUI=a.delta_DWUI,native_MPRO=a.MPRO,strict_region_pairs=c.strict_region_pairs)
        e=eligibility[eligibility.route.eq(a.route)&eligibility.gene_id.eq(a.gene_id)&eligibility.contrast.eq('MP')]
        row['interval_status']=';'.join(e.status.unique()) if len(e) else 'no_interval_eligibility_record'
        m=models[a.route].get(a.gene_id);mask=np.zeros(len(ids),bool)
        if m:
            fam=m['family'];lo=fam['common_start0'];hi=fam['common_end0']
            # Family bounds are genomic half-open UTR intervals. A cut retains the
            # base at cut on + and cut+1 on -, in 1-based annotation coordinates.
            retained=cuts+(strand=='-');mask=(retained>lo)&(retained<=hi)
            origin=lo if strand=='+' else hi;x=(cuts-origin)*(1 if strand=='+' else -1)
            mp=P[mask].sum();mm=M[mask].sum()
            row.update(family_id=fam['family_id'],family_start0=lo,family_end0=hi,terminal_native_pas=int(mask.sum()),terminal_pas_span_nt=float(np.ptp(signed[mask])) if mask.any() else np.nan,terminal_mass_P=mp,terminal_mass_M=mm)
            row['family_status']='measurable_native_family' if mask.sum()>=2 and min(mp,mm)>0 else 'fewer_than_two_native_family_sites_or_zero_mass'
            if mask.sum()>=2 and min(mp,mm)>0:
                up=P[mask]@x[mask]/mp;um=M[mask]@x[mask]/mm
                within=(mp+mm)/2*(up-um)/span;mass=(mp-mm)*(up+um)/2/span;outside=(P[~mask]-M[~mask])@x[~mask]/span
                np.testing.assert_allclose(within+mass+outside,a.delta_DWUI,atol=1e-12)
                row.update(terminal_position_P_nt=up,terminal_position_M_nt=um,terminal_shift_nt=up-um,within_family_DWUI=within,family_mass_DWUI=mass,nonterminal_DWUI=outside,nonterminal_mass_delta=mm-mp,decomposition_residual=within+mass+outside-a.delta_DWUI)
            for threshold in [20,50]:
                pp={v['library'] for v in m['profiles'].get('P',{}).get('libraries',[]) if v['total']>=threshold and v['coverage']>=.5}
                rr={v['library'] for v in m['profiles'].get('M',{}).get('libraries',[]) if v['total']>=threshold and v['coverage']>=.5}
                row[f'interval_paired_libraries_{threshold}']=len(pp&rr)
        else:
            row['family_status']='no_frozen_terminal_family';x=signed-signed.min()
        for i,pas in enumerate(ids):
            sites.append(dict(gene_id=a.gene_id,gene_symbol=a.gene_symbol,route=a.route,pas_id=pas,cut_interbase=int(cuts[i]),strand=strand,region_class=ax.primary_region_class.iloc[i],terminal_family=bool(mask[i]),relative_to_family_origin_nt=float(x[i]),position=float(ax.position.iloc[i]),P_usage=P[i],M_usage=M[i],delta_usage=P[i]-M[i],native_contribution=(P[i]-M[i])*ax.position.iloc[i]))
        results.append(row)
    res=pd.DataFrame(results);write('axis-decomposition.tsv',res);write('native-sites.tsv.gz',sites)
    write('low-shift-cohort.tsv',cohort);write('interval-eligibility.tsv',eligibility[eligibility.gene_id.isin(cohort.gene_id)])
    write('model-depth.tsv',mfx[mfx.gene_id.isin(cohort.gene_id)])
    pairs=[]
    for chunk in pd.read_csv(ctx.register(ctx.analysis/'state_vectors/pair_effects.primary.tsv.gz'),sep='\t',chunksize=20000):
        z=chunk[chunk.gene_id.isin(cohort.gene_id)&chunk.route.isin(['CIP','WOT'])&chunk.setting.eq('supported50_p08')&chunk.contrast.eq('P_M')&chunk.design.eq('same_library')]
        pairs.append(z)
    pairpath=ctx.analysis/'state_vectors/pair_effects.primary.tsv.gz';inputs[str(pairpath)]=hashlib.sha256(pairpath.read_bytes()).hexdigest();write('native-contributors.tsv.gz',pd.concat(pairs))


def state_references(ctx):
    f=ctx.read(ctx.analysis/'profiles/all_route_node_positions.tsv.gz')
    f=f[f.setting.eq('supported50_p08')&f.depth.eq(20)]
    rows=[];members=[]
    for route,q in f.groupby('route'):
        p=q.pivot(index='gene_id',columns='node',values='position_index')
        observed={s:(p[s]-p.S).dropna() for s in 'REPM' if s in p and 'S' in p}
        common=set.intersection(*(set(z.index) for z in observed.values())) if observed else set()
        for state in 'REPM':
            z=observed.get(state,pd.Series(dtype=float))
            for background,ids in [('available',set(z.index)),('fixed_common_observed_states_this_route',common&set(z.index))]:
                values=z.loc[sorted(ids)]
                rows.append(dict(route=route,species=q.species.iloc[0],study_id=q.study_id.iloc[0],state=state,background=background,
                                 n=len(ids),mean=values.mean(),status='measured' if len(ids) else 'unavailable'))
                members.extend(dict(route=route,state=state,background=background,gene_id=g,delta=float(z[g])) for g in sorted(ids))
    ctx.table(rows,'state-reference-summary.tsv');ctx.table(members,'state-reference-members.tsv.gz',columns=['route','state','background','gene_id','delta'])
    fit=np.load(ctx.register(ctx.analysis/'profiles/joint_fit.npz'),allow_pickle=True)
    centers=ctx.read(ctx.analysis/'profiles/shape_centers.tsv');centers=centers[centers.k.eq(3)].sort_values('cluster')[list('SREP')].to_numpy()
    distances=np.linalg.norm(fit['features'][:,None,:]-centers[None,:,:],axis=2)
    assigned=distances.argmin(1);distances=np.sort(distances,axis=1)
    assert np.array_equal(assigned,fit['labels_k3'])
    ctx.table(pd.DataFrame(dict(species=fit['species'],gene_id=fit['gene_ids'],cluster=['C'+str(i+1) for i in assigned],
        assigned_distance=distances[:,0],next_distance=distances[:,1],margin=distances[:,1]-distances[:,0])),'k3-assignment-margins.tsv')


def structures(ctx):
    selected=ctx.read(ctx.analysis/'candidate_screen/example_selection.tsv')
    cross=ctx.read(ctx.analysis/'cross_trace/member_migration.tsv');cross=cross[cross.new_member]
    axes=ctx.read(ctx.analysis/'state_vectors/PAS_axes.primary.tsv.gz')
    # Include current manuscript loci for audit, without requiring eligibility.
    symbols=set(selected.gene_symbol.str.upper())|{'PDLIM7','EIF4EBP2','RAD23B'}
    wanted=set(axes.loc[axes.primary_gene_symbol.fillna('').str.upper().isin(symbols),'primary_gene_id'])
    for side in ['human','mouse']:wanted|=set(cross[side+'_gene_id'])
    txs={};segments=[]
    for species in ['human','mouse']:
        with gzip.open(ctx.resource(species+'_ensembl_gtf'),'rt') as handle:
            for line in handle:
                if line.startswith('#'):continue
                v=line.rstrip().split('\t')
                if len(v)!=9 or v[2] not in ['transcript','exon','CDS','stop_codon']:continue
                a=dict(re.findall(r'(\w+) "([^"]*)"',v[8]));gene=a.get('gene_id','').split('.')[0]
                if gene not in wanted:continue
                tx=a.get('transcript_id','');assert tx
                row=dict(species=species,gene_id=gene,transcript_id=tx,chrom=v[0],strand=v[6],start0=int(v[3])-1,end0=int(v[4]),kind=v[2],annotation_source='Ensembl release 116; static structural annotation')
                if v[2]=='transcript':txs[species,tx]=row
                else:segments.append(row)
    cols=['species','gene_id','transcript_id','chrom','strand','start0','end0','kind','annotation_source']
    seg=pd.DataFrame(segments,columns=cols);tx=pd.DataFrame(txs.values(),columns=cols)
    ctx.table(tx,'transcript-models.tsv.gz');ctx.table(seg,'exon-CDS.tsv.gz')
    checks=[]
    cutcols=['primary_gene_id','atlas_pas_id','chrom','strand','cut_site']
    cuts=axes[axes.setting.eq('supported50_p08')&axes.primary_gene_id.isin(wanted)][cutcols].drop_duplicates()
    groups={g:q for g,q in seg.groupby('gene_id')}
    for r in cuts.itertuples():
        retained=int(r.cut_site)+(r.strand=='-')
        for tid,q in groups.get(r.primary_gene_id,seg.iloc[:0]).groupby('transcript_id'):
            overlap=q.chrom.astype(str).str.removeprefix('chr').eq(str(r.chrom).removeprefix('chr'))&q.strand.eq(r.strand)&q.start0.lt(retained)&q.end0.ge(retained)
            checks.append(dict(gene_id=r.primary_gene_id,pas_id=r.atlas_pas_id,transcript_id=tid,cut_interbase=r.cut_site,
                retained_base1=retained,overlaps_exon=bool((overlap&q.kind.eq('exon')).any()),overlaps_CDS=bool((overlap&q.kind.eq('CDS')).any())))
    ctx.table(checks,'PAS-transcript-coordinate-check.tsv.gz',columns=['gene_id','pas_id','transcript_id','cut_interbase','retained_base1','overlaps_exon','overlaps_CDS'])
    ctx.table(cuts,'selected-native-PAS.tsv.gz')


def screen_failures(ctx):
    cohort=ctx.read(ctx.analysis/'candidate_screen/low_shift_candidates.tsv')
    r=ctx.read(ctx.analysis/'intervals/all_coordinate_regions.tsv.gz')
    r=r[r.gene_id.isin(cohort.gene_id)&r.contrast.eq('MP')&r.route.isin(['CIP','WOT'])].copy()
    r['fails_magnitude']=r.effect20.abs().lt(.10)|r.effect50.abs().lt(.10)
    r['fails_same_sign']=r.effect20*r.effect50<=0
    r['fails_paired_included_excluded_support']=~r.paired_observed_both5
    r['missing_estimate']=r.effect20.isna()|r.effect50.isna()
    ctx.table(r,'all-regions-and-failure-flags.tsv.gz')
    merged=ctx.read(ctx.out/'axis-decomposition.tsv');rows=[]
    for x in cohort.itertuples():
        z=merged[merged.gene_id.eq(x.gene_id)];a=r[r.gene_id.eq(x.gene_id)]
        rows.append(dict(gene_id=x.gene_id,gene_symbol=x.gene_symbol,strict_region_pairs=x.strict_region_pairs,
            outcome='strict_recurrent_interval' if x.strict_region_pairs else 'not_admitted_to_strict_recurrent_interval',
            route_eligibility=';'.join(t.route+':'+t.interval_status for t in z.itertuples()),recorded_regions=len(a),
            single_route_strict_regions=int(a.strict_change.sum()),note='Eligibility/failure flags do not establish absence of biological interval change'))
    ctx.table(rows,'screen-outcomes.tsv',columns=['gene_id','gene_symbol','strict_region_pairs','outcome','route_eligibility','recorded_regions','single_route_strict_regions','note'])


def run(ctx):
    decompose(ctx);state_references(ctx);structures(ctx);screen_failures(ctx)
    ctx.finish(dict(full_new_low_shift_cohort=True,exact_DWUI_decomposition=True,unavailable_not_zero=True,
        within_route_state_background=True,new_PAS_annotation_intersections=True),
        structure_scope='Exon/CDS compatibility; no measured protein length, abundance or activity inferred')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'geometry',a.out))
