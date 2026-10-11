"""Fresh annotation PAS, reciprocal architecture and native-count driver features."""
import argparse
from collections import defaultdict
import json

import numpy as np
import pandas as pd

from downstream_runtime import Context,Study
from downstream_cross import CODES
from downstream_regulatory import SOURCES,reference_lookup,F5
from downstream_liftover import lift_points
from downstream_cis_models import models,CIS,KEYS

RANK={'AATAAA':3,'ATTAAA':2,'AGTAAA':2,'AATACA':2,'AATATA':2,'CATAAA':2,'TATAAA':2,
      'AAGAAA':1,'AATAAT':1,'AACAAA':1,'AATAGA':1,'AATGAA':1,'GATAAA':1}


def select_pairs(candidates):
    if candidates.empty:return candidates.copy()
    order=candidates.sort_values(['reciprocal25','reciprocal50','forward_distance','backward_distance','human_id','mouse_id'],
        ascending=[False,False,True,True,True,True])
    used_h=set();used_m=set();rows=[]
    for r in order.to_dict('records'):
        if r['human_id'] in used_h or r['mouse_id'] in used_m:continue
        used_h.add(r['human_id']);used_m.add(r['mouse_id']);rows.append(r)
    return pd.DataFrame(rows,columns=candidates.columns)


def catalog(ctx,orth):
    chooser=ctx.functions(F5+'site-recurrence-v1/figure5-publication-v14r4/prepare_annotated_tracks.py',['choose'])
    mappings=[];usage=[]
    for sp in ['human','mouse']:
        lookup=reference_lookup(ctx,sp);genes=set(orth[sp+'_gene_id'])
        for code,(sid,route) in CODES.items():
            if (code in ['W','C'])!=(sp=='mouse'):continue
            study=Study(ctx,sid);f=study.primary[study.primary.primary_gene_id.isin(genes)&study.primary.matched_source_count_50nt.ge(1)]
            by_pas={}
            for r in f.itertuples():
                candidates=[]
                for source in SOURCES:
                    distance=getattr(r,source+'_distance_nt')
                    if pd.isna(distance) or not 0<=distance<=50:continue
                    cut=int(getattr(r,source+'_nearest_cut'));identity=lookup[source,r.chrom,r.strand,cut]
                    candidates.append(dict(annotation_id=identity,annotation_cut=int(identity.split(':')[5]),source=source,distance_nt=int(distance)))
                selected=chooser.choose(candidates);assert selected is not None
                best=selected[0];aid='|'.join(map(str,[sp,r.primary_gene_id,r.chrom,r.strand,best['annotation_cut']]))
                by_pas[r.atlas_pas_id]=aid
                mappings.append(dict(species=sp,code=code,gene_id=r.primary_gene_id,chrom=r.chrom,strand=r.strand,
                    native_pas_id=r.atlas_pas_id,native_cut=r.cut_site,id=aid,cut_site=best['annotation_cut'],
                    annotation_id=best['annotation_id'],source_distance=best['distance_nt'],representative_distance=abs(r.cut_site-best['annotation_cut']),
                    n_sources=len({c['source'] for c in selected}),sources=';'.join(sorted({c['source'] for c in selected}))))
            for gene in sorted(set(f.primary_gene_id)):
                axis=study.axis(gene);data=study.gene_counts(axis);states=study.profile_states(data)
                grouped=defaultdict(list)
                for (rr,lib,node),s in states.items():
                    if rr==route:grouped[node].append(s['usage'])
                for node,values in grouped.items():
                    mean=np.mean(values,axis=0);merged=defaultdict(float)
                    for pas,value in zip(axis.atlas_pas_id,mean):merged[by_pas[pas]]+=float(value)
                    assert abs(sum(merged.values())-1)<1e-12
                    usage.extend(dict(species=sp,code=code,gene_id=gene,node=node,id=aid,usage=u,n_libraries=len(values)) for aid,u in merged.items())
            print('cis native usage',sp,code,len(by_pas),flush=True)
    mapping=pd.DataFrame(mappings)
    cat=mapping.groupby(['species','gene_id','id','chrom','strand','cut_site'],as_index=False).agg(n_studies=('code','nunique'),n_detected=('native_pas_id','size'),max_sources=('n_sources','max'))
    assert cat.id.is_unique
    usage=pd.DataFrame(usage,columns=['species','code','gene_id','node','id','usage','n_libraries'])
    ctx.table(mapping,'detected_to_annotation.tsv.gz');ctx.table(cat,'annotation_catalog.tsv.gz');ctx.table(usage,'annotated_state_usage.tsv.gz')
    return cat,usage


def architecture(ctx,cat,orth):
    lifts={};frames={}
    for sp,a,b in [('human','hg38','mm39'),('mouse','mm39','hg38')]:
        f=cat[cat.species.eq(sp)]
        q=lift_points(ctx,f[['id','chrom','strand','cut_site']],'chain_'+a+'_'+b,sp+'_reciprocal')
        lifts[sp]=q.set_index('id');frames[sp]=f.set_index('id')
    mgenes={g:v for g,v in frames['mouse'].groupby('gene_id')};partners=orth.set_index('human_gene_id').mouse_gene_id.to_dict();rows=[]
    for r in frames['human'].itertuples():
        a=lifts['human'].loc[r.Index]
        if a.mapping_status!='unique':continue
        m=mgenes.get(partners[r.gene_id])
        if m is None:continue
        m=m[m.chrom.eq(a.mapped_chrom)&m.strand.eq(a.mapped_strand)&m.cut_site.sub(a.mapped_cut).abs().le(50)]
        for v in m.itertuples():
            b=lifts['mouse'].loc[v.Index];forward=abs(v.cut_site-a.mapped_cut)
            backward=abs(b.mapped_cut-r.cut_site) if b.mapping_status=='unique' and b.mapped_chrom==r.chrom and b.mapped_strand==r.strand else np.inf
            rows.append(dict(human_gene_id=r.gene_id,mouse_gene_id=v.gene_id,human_id=r.Index,mouse_id=v.Index,
                forward_distance=forward,backward_distance=backward,reciprocal25=bool(max(forward,backward)<=25),
                reciprocal50=bool(max(forward,backward)<=50),legacy_asymmetric25=bool(forward<=50 and backward<=25)))
    candidates=pd.DataFrame(rows,columns=['human_gene_id','mouse_gene_id','human_id','mouse_id','forward_distance','backward_distance','reciprocal25','reciprocal50','legacy_asymmetric25'])
    one=select_pairs(candidates);ctx.table(candidates,'reciprocal_PAS_candidates.tsv.gz');ctx.table(one,'reciprocal_one_to_one.tsv.gz')
    rows=[]
    for tolerance in [25,50]:
        selected=one[one['reciprocal'+str(tolerance)]];conserved=set(selected.human_id)|set(selected.mouse_id)
        lookup={sp:{gid:g for gid,g in cat[cat.species.eq(sp)].groupby('gene_id')} for sp in ['human','mouse']}
        for r in orth.itertuples():
            row=dict(human_gene_id=r.human_gene_id,mouse_gene_id=r.mouse_gene_id,tolerance=tolerance)
            nh=nm=0
            for sp,gid,stem in [('human',r.human_gene_id,'h'),('mouse',r.mouse_gene_id,'m')]:
                f=lookup[sp].get(gid,cat.iloc[:0]);n=len(f)
                if sp=='human':nh=n
                else:nm=n
                row['n_'+stem+'_pas']=n;row[stem+'_span']=float(f.cut_site.max()-f.cut_site.min()) if n else np.nan
                if n:
                    signed=f.cut_site*np.where(f.strand.eq('+'),1,-1);order=f.loc[signed.sort_values().index]
                    row[stem+'_proximal_conserved']=int(order.id.iloc[0] in conserved);row[stem+'_distal_conserved']=int(order.id.iloc[-1] in conserved)
                else:row[stem+'_proximal_conserved']=row[stem+'_distal_conserved']=np.nan
            matched=int(selected.human_gene_id.eq(r.human_gene_id).sum());row.update(n_matched_pas=matched,A_g=2*matched/(nh+nm) if nh+nm else np.nan)
            rows.append(row)
    result=pd.DataFrame(rows);ctx.table(result,'gene_architecture.tsv.gz')
    return one,result


def sequence_features(ctx,cat):
    import pysam
    import pyBigWig
    fa={sp:pysam.FastaFile(str(ctx.resource(key))) for sp,key in [('human','cis_hg38_fasta'),('mouse','cis_mm10_fasta')]}
    bw=pyBigWig.open(str(ctx.resource('cis_phastcons')));sizes=bw.chroms();rows=[]
    def motif(fasta,chrom,cut,strand):
        if chrom not in fasta.references:return np.nan
        start,end=(cut-35,cut-9) if strand=='+' else (cut+9,cut+35)
        if start<0 or end>fasta.get_reference_length(chrom):return np.nan
        seq=fasta.fetch(chrom,start,end).upper()
        if strand=='-':seq=seq.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
        if len(seq)!=26 or set(seq)-set('ACGT'):return np.nan
        return max([v for word,v in RANK.items() if word in seq]+[0])
    for sp in ['human','mouse']:
        f=cat[cat.species.eq(sp)];q=f[['id','chrom','strand','cut_site']]
        chains=[('hg38','mm10')] if sp=='human' else [('mm39','hg38'),('hg38','mm10')]
        for a,b in chains:
            m=lift_points(ctx,q,'chain_'+a+'_'+b,sp+'_sequence_'+a+'_'+b)
            q=m[m.mapping_status.eq('unique')][['id','mapped_chrom','mapped_strand','mapped_cut']].rename(columns=dict(mapped_chrom='chrom',mapped_strand='strand',mapped_cut='cut_site'))
        lookup=q.set_index('id')
        for r in f.itertuples():
            mm=lookup.loc[r.id] if r.id in lookup.index else None;phast=np.nan
            if mm is not None and mm.chrom in sizes:
                a=max(0,int(mm.cut_site)-100);b=min(sizes[mm.chrom],int(mm.cut_site)+50)
                values=np.asarray(bw.values(mm.chrom,a,b),float)
                if np.isfinite(values).any():phast=float(np.nanmean(values))
            rank=motif(fa[sp],r.chrom,int(r.cut_site),r.strand) if sp=='human' else motif(fa[sp],mm.chrom,int(mm.cut_site),mm.strand) if mm is not None else np.nan
            rows.append(dict(id=r.id,species=sp,gene_id=r.gene_id,phastcons=phast,motif_rank=rank,
                sequence_mapping_status='unique' if mm is not None else 'unavailable_or_ambiguous'))
    for f in fa.values():f.close()
    bw.close();result=pd.DataFrame(rows);ctx.table(result,'site_sequence_features.tsv.gz');return result


def drivers(ctx,usage,cat,one,seq,orth):
    lookup=cat.set_index('id');sequence=seq.set_index('id');rows=[]
    m2h=orth.set_index('mouse_gene_id').human_gene_id.to_dict()
    for tolerance in [25,50]:
        conserved=set(one.loc[one['reciprocal'+str(tolerance)],'human_id'])|set(one.loc[one['reciprocal'+str(tolerance)],'mouse_id'])
        for contrast,reference in [('P-S','S'),('P-M','M')]:
            for sp in ['human','mouse']:
                f=usage[usage.species.eq(sp)&usage.node.isin(['P',reference])]
                wide=f.pivot(index=['gene_id','id','code'],columns='node',values='usage').reindex(columns=['P',reference]).dropna()
                wide['delta']=wide.P-wide[reference]
                d=wide.reset_index().groupby(['gene_id','id'],as_index=False).agg(delta=('delta','mean'),n_studies=('code','nunique'))
                d['conserved']=d.id.isin(conserved)
                for gid,g in d.groupby('gene_id'):
                    g=g.sort_values('id');gain=g.loc[g.delta.idxmax()];loss=g.loc[g.delta.idxmin()]
                    s=sequence.loc[[gain.id,loss.id]]
                    rows.append(dict(human_gene_id=gid if sp=='human' else m2h[gid],gene_id=gid,discovery=sp,contrast=contrast,tolerance=tolerance,
                        n_discovery_pas=len(g),frac_conserved_discovery_pas=float(g.conserved.mean()),
                        gain_id=gain.id,loss_id=loss.id,gain_delta=gain.delta,loss_delta=loss.delta,driver_conserved=int(gain.conserved)+int(loss.conserved),
                        driver_same_site=gain.id==loss.id,driver_min_studies=min(gain.n_studies,loss.n_studies),
                        phastcons_driver_mean=s.phastcons.mean(),motif_rank_driver_mean=s.motif_rank.mean(),
                        phastcons_available_drivers=int(s.phastcons.notna().sum()),motif_available_drivers=int(s.motif_rank.notna().sum())))
    out=pd.DataFrame(rows);ctx.table(out,'driver_sites.tsv.gz');return out


def descriptive_controls(ctx,data):
    """Retain the prespecified matched S/R/B populations without new matching."""
    primary=data[data.layer.eq('supported50')&data.tau.eq(.02)].copy()
    matches=ctx.read(ctx.analysis/'cross_continuous/matched_comparators.tsv.gz')
    if len(matches):
        matches=matches[matches.layer.eq('supported50')&matches.tau.eq(.02)].copy()
        matches=matches.merge(primary[KEYS+['human_gene_id']+CIS],on=KEYS+['human_gene_id'],validate='many_to_one')
        matches['class']=np.where(matches.role.eq('control'),'Matched background',
            np.where(matches.validation_same_direction,'Shared','Recurrent-only'))
    else:matches=pd.DataFrame(columns=KEYS+['human_gene_id','match_id','role','match_weight','class']+CIS)
    ctx.table(matches,'descriptive_matched_features.tsv.gz')
    counts=[];distributions=[];correlations=[]
    classes=['Shared','Recurrent-only','Matched background']
    for contrast,discovery in [('P-S','human'),('P-S','mouse'),('P-M','human'),('P-M','mouse')]:
        identity=dict(contrast=contrast,discovery=discovery)
        p=primary[primary.contrast.eq(contrast)&primary.discovery.eq(discovery)]
        m=matches[matches.contrast.eq(contrast)&matches.discovery.eq(discovery)]
        counts.append(dict(identity,model_n=len(p),model_Y1=int(p.Y.sum()),
            **{name:int(m.loc[m['class'].eq(name),'human_gene_id'].nunique()) for name in classes}))
        for feature in CIS:
            for name in classes:
                z=m[m['class'].eq(name)].dropna(subset=[feature])
                z=z.groupby(['human_gene_id',feature],as_index=False).match_weight.sum()
                z=z[z.match_weight.gt(0)].sort_values(feature)
                if len(z):
                    positions=(z.match_weight.cumsum()-.5*z.match_weight)/z.match_weight.sum()
                    quantiles=np.interp([.25,.5,.75],positions,z[feature])
                else:quantiles=[np.nan]*3
                distributions.append(dict(identity,feature=feature,**{'class':name},n_observed=len(z),
                    weight_sum=z.match_weight.sum(),q25=quantiles[0],median=quantiles[1],q75=quantiles[2]))
            for other in CIS:
                z=p[[feature,other]].dropna() if feature!=other else p[[feature]].dropna()
                rho=z[feature].corr(z[other],method='spearman') if len(z)>1 else np.nan
                correlations.append(dict(identity,feature1=feature,feature2=other,n_pairwise=len(z),spearman=rho))
    ctx.table(counts,'descriptive_class_counts.tsv');ctx.table(distributions,'feature_distribution_summary.tsv')
    ctx.table(correlations,'feature_correlations.tsv')


def run(ctx):
    orth=ctx.read(ctx.resource('ortholog'));cat,usage=catalog(ctx,orth);one,arch=architecture(ctx,cat,orth)
    seq=sequence_features(ctx,cat);driver=drivers(ctx,usage,cat,one,seq,orth)
    universe=ctx.read(ctx.analysis/'cross_continuous/reciprocal_gene_universe.tsv.gz')
    primary=driver[driver.tolerance.eq(25)].drop(columns='tolerance')
    data=universe.merge(primary,on=['human_gene_id','contrast','discovery'],how='left',validate='many_to_one')
    data=data.merge(arch[arch.tolerance.eq(25)].drop(columns=['tolerance','mouse_gene_id']),on='human_gene_id',how='left',validate='many_to_one')
    data['Y']=data.validation_same_direction.astype(int);data['abs_effect']=data.disc_meta.abs()
    ctx.table(data,'model_input.tsv.gz');descriptive_controls(ctx,data);models(ctx,data)
    ctx.finish(dict(new_catalogue_from_native_features=True,source_to_annotation_mapping_verified=True,
        interbase_strand_aware_liftover=True,unique_reciprocal_one_to_one=True,new_driver_usage_from_native_counts=True,
        sequence_requeried_at_new_PAS=True,new_model_universe=True,Firth_fits_and_conditional_bootstrap=True,paired_out_of_fold_comparisons=True,
        original_matched_control_protocol_and_weighted_diagnostics=True),
        method_corrections=['reciprocal25 requires both forward and backward <=25; old asymmetric flag preserved in audit',
            'minus PAS lifted as retained [cut,cut+1); strand changes propagated into motif sequence',
            'proximal and distal architecture anchors follow transcript strand',
            'phastCons uses exact base mean rather than bigWig zoom approximation'])


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);a=p.parse_args();run(Context(a.config,'cis',a.out))
