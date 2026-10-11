"""Frozen, strand-aware UGUA susceptibility tests on independent strict-UTR backgrounds."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
TABLES=ROOT/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables'
LABELS=['proximal','distal','nonrecurrent']
WINDOWS={'upstream40_100':(-100,-40),'upstream0_100':(-100,0),'upstream40_80':(-80,-40)}
MOTIFS=['TGTA','TGAT','CGTA']
MODELS=['joint_covariate_median','assay_architecture','residualized_assay_architecture']
RANK={'AATAAA':3,'ATTAAA':2,'AGTAAA':2,'AATACA':2,'AATATA':2,'CATAAA':2,'TATAAA':2,'AAGAAA':1,'AATAAT':1,'AACAAA':1,'AATAGA':1,'AATGAA':1,'GATAAA':1}
NPERM=20000


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def save(frame,name):
    frame.to_csv(HERE/name,sep='\t',index=False,na_rep='NA',float_format='%.17g',compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)


def bh(values):
    p=np.asarray(values,float);order=np.argsort(p,kind='stable');ranked=p[order]*len(p)/(np.arange(len(p))+1)
    q=np.minimum.accumulate(ranked[::-1])[::-1].clip(0,1);out=np.empty(len(p));out[order]=q;return out


def seed(*values):
    return int(hashlib.sha256('|'.join(map(str,values)).encode()).hexdigest()[:16],16)%2**32


def motif_count(sequence,motif):
    return sum(sequence[i:i+len(motif)]==motif for i in range(len(sequence)-len(motif)+1))


def window_sequence(upstream100,bounds):
    lo,hi=bounds;assert -100<=lo<hi<=0
    return upstream100[100+lo:100+hi]


def feature_table():
    inventory=pd.read_csv(HERE/'background_inventory.tsv.gz',sep='\t',float_precision='round_trip')
    s=pd.read_csv(HERE/'upstream100_sequences.tsv',sep='\t')
    assert s.sequence.str.len().eq(100).all()
    seqs=s.set_index(['species','gene_id','anchor']).sequence.to_dict()
    eligible=inventory[inventory.status.eq('eligible')].copy()
    rows=[];windows=[]
    for r in eligible.itertuples():
        proximal=seqs[r.species,r.gene_id,'proximal'];distal=seqs[r.species,r.gene_id,'distal']
        assert not(set(proximal+distal)-set('ACGT'))
        core=lambda seq:max([v for k,v in RANK.items() if k in window_sequence(seq,(-35,-9))]+[0])
        row=dict(species=r.species,gene_id=r.gene_id,GC_mean=(sum(c in 'GC' for c in proximal+distal))/200,proximal_core_rank=core(proximal),distal_core_rank=core(distal))
        row['core_rank_difference']=row['distal_core_rank']-row['proximal_core_rank']
        for wn,bounds in WINDOWS.items():
            lo,hi=bounds
            for anchor,sequence,cut in [('proximal',proximal,int(r.proximal_cut)),('distal',distal,int(r.distal_cut))]:
                a0,b0=(cut+lo,cut+hi) if r.strand=='+' else (cut-hi,cut-lo)
                windows.append(dict(species=r.species,gene_id=r.gene_id,anchor=anchor,window=wn,chrom=r.chrom,strand=r.strand,cut_site=cut,start0=a0,end0=b0,sequence=window_sequence(sequence,bounds)))
            for motif in MOTIFS:
                pc=motif_count(window_sequence(proximal,bounds),motif);dc=motif_count(window_sequence(distal,bounds),motif)
                row[f'{wn}_{motif}_proximal_count']=pc;row[f'{wn}_{motif}_distal_count']=dc;row[f'{wn}_{motif}_difference']=dc-pc
        rows.append(row)
    save(pd.DataFrame(windows),'motif_windows.tsv.gz')
    features=eligible.merge(pd.DataFrame(rows),on=['species','gene_id'],validate='one_to_one')
    # Join response labels only after all structural and complete-sequence filters.
    f=pd.read_csv(TABLES/'recurrence_controls/recurrence_frozen_covariates.tsv.gz',sep='\t',float_precision='round_trip')
    f=f[f.contrast.eq('P_M')&f.threshold.eq(.02)]
    roster=pd.read_csv(ROOT/'notes/2026-10-10-directional-recurrence/observed_members.tsv.gz',sep='\t')
    roster=roster[roster.contrast.eq('P_M')&roster.threshold.eq(.02)]
    labels=[]
    for (sp,gene),part in f.groupby(['species','gene_id'],sort=True):
        nprox=int(part.event.eq(-1).sum());ndist=int(part.event.eq(1).sum())
        assert not(nprox>=2 and ndist>=2)
        labels.append(dict(species=sp,gene_id=gene,direction='proximal' if nprox>=2 else 'distal' if ndist>=2 else 'nonrecurrent',no_contributor_event=bool(part.event.eq(0).all()),proximal_support=nprox,distal_support=ndist))
    labs=pd.DataFrame(labels)
    recurrent=labs[labs.direction.ne('nonrecurrent')]
    assert set(zip(recurrent.species,recurrent.gene_id,recurrent.direction))==set(zip(roster.species,roster.gene_id,roster.direction))
    features=features.merge(labs,on=['species','gene_id'],validate='one_to_one')
    assert features.notna()[['GC_mean','PAS_number','baseline','anchor_span','mean_GEX','APA_coverage','core_rank_difference','direction']].all().all()
    save(features,'gene_features_and_labels.tsv.gz')
    return features


def quantiles(values,n):
    assert np.isfinite(values).all()
    return pd.qcut(values,n,labels=False,duplicates='drop').fillna(0).astype(int)


def partition(frame,model):
    z=frame.copy()
    for name,n in [('PAS_number',3),('baseline',2),('anchor_span',2),('GC_mean',2),('mean_GEX',2),('APA_coverage',2)]:z[name+'_bin']=quantiles(z[name],n)
    z['core_rank_difference_sign']=np.sign(z.core_rank_difference).astype(int)
    cols=['measured_pattern','proximal_capacity','distal_capacity','PAS_number_bin','baseline_bin','anchor_span_bin']
    if model=='joint_covariate_median':cols+=['GC_mean_bin','mean_GEX_bin','APA_coverage_bin','core_rank_difference_sign']
    cells=[np.asarray(i,int) for i in z.groupby(cols,sort=True).indices.values()]
    z['stratum']=-1
    for j,ii in enumerate(cells):z.loc[ii,'stratum']=j
    assert z.stratum.ge(0).all()
    return z,cells,cols


def residualize(values,frame):
    data=pd.DataFrame(dict(baseline=frame.baseline,log_PAS=np.log1p(frame.PAS_number),log_span=np.log(frame.anchor_span),GC=frame.GC_mean,log_GEX=np.log1p(frame.mean_GEX),log_coverage=np.log1p(frame.APA_coverage),core_difference=frame.core_rank_difference))
    continuous=data.to_numpy(float);sd=continuous.std(0);continuous[:,sd>0]=(continuous[:,sd>0]-continuous[:,sd>0].mean(0))/sd[sd>0]
    categorical=pd.get_dummies(frame[['measured_pattern','proximal_capacity','distal_capacity']].astype(str),drop_first=True,dtype=float).to_numpy()
    x=np.column_stack([np.ones(len(frame)),continuous,categorical]);coef=np.linalg.lstsq(x,values,rcond=None)[0]
    resid=values-x@coef
    assert np.allclose(x.T@resid,0,atol=1e-9,rtol=0)
    return resid,int(np.linalg.matrix_rank(x))


def label_draws(labels,cells,random):
    out=np.broadcast_to(labels,(NPERM,len(labels))).copy()
    for ii in cells:
        if len(np.unique(labels[ii]))>1:out[:,ii]=random.permuted(out[:,ii],axis=1)
        for label in range(3):assert (out[:,ii]==label).sum(1).tolist()==[int((labels[ii]==label).sum())]*NPERM
    return out


def stat(values,labels,other):
    n0=int((labels==0).sum());n1=int((labels==other).sum())
    if not n0 or not n1:return np.nan
    return float(values[labels==0].mean()-values[labels==other].mean())


def statistic_draws(values,draws,other):
    n0=int((draws[0]==0).sum());n1=int((draws[0]==other).sum())
    if not n0 or not n1:return np.full(NPERM,np.nan)
    return (draws==0)@values/n0-(draws==other)@values/n1


def expectation(values,labels,cells,other):
    n0=int((labels==0).sum());n1=int((labels==other).sum())
    if not n0 or not n1:return np.nan
    score=0.
    for ii in cells:
        score+=float(values[ii].sum())*(float((labels[ii]==0).mean())/n0-float((labels[ii]==other).mean())/n1)
    return score


def null_summary(observed,expected,values):
    if not np.isfinite(observed):return dict(observed=np.nan,analytic_expected=np.nan,null_mean=np.nan,null_sd=np.nan,null_p025=np.nan,null_p975=np.nan,excess=np.nan,p_enrichment=1.,p_two_sided=1.,tail_upper=-1,tail_two_sided=-1,informative=False,status='empty_comparison_group')
    sd=float(values.std(ddof=1));assert abs(values.mean()-expected)<=6*sd/np.sqrt(len(values))+1e-10
    upper=int((values>=observed-1e-12).sum());both=int((abs(values-expected)>=abs(observed-expected)-1e-12).sum())
    informative=bool(sd>1e-12)
    return dict(observed=observed,analytic_expected=expected,null_mean=float(values.mean()),null_sd=sd,null_p025=float(np.quantile(values,.025)),null_p975=float(np.quantile(values,.975)),excess=observed-expected,p_enrichment=(1+upper)/(len(values)+1) if informative else 1.,p_two_sided=(1+both)/(len(values)+1) if informative else 1.,tail_upper=upper,tail_two_sided=both,informative=informative,status='estimable' if informative else 'fixed_under_matched_null')


def main():
    locked=json.loads((HERE/'protocol.json').read_text());assert locked['protocol_sha256']==digest(HERE/'PROTOCOL.zh.md')
    features=feature_table();rows=[];drawstore={};strata=[];stratum_members=[];summaries=[];balance=[];paired={}
    cov=['PAS_number','baseline','anchor_span','GC_mean','mean_GEX','APA_coverage','core_rank_difference']
    for background in ['paired_strict','species_strict']:
        for species in ['human','mouse']:
            part=features[features.species.eq(species)&(features.paired_eligible if background=='paired_strict' else True)].copy().sort_values('gene_id').reset_index(drop=True)
            context=dict(background=background,species=species,n_genes=len(part))
            label=part.direction.map({v:i for i,v in enumerate(LABELS)}).to_numpy(np.int8)
            for direction in LABELS+['all_contributors_no_event']:
                subset=part[part.direction.eq(direction)] if direction in LABELS else part[part.no_contributor_event]
                base=context|dict(direction=direction,n_group=len(subset),n_no_event=int(subset.no_contributor_event.sum()))
                for wn in WINDOWS:
                    for motif in MOTIFS:
                        pc=subset[f'{wn}_{motif}_proximal_count'];dc=subset[f'{wn}_{motif}_distal_count'];diff=dc-pc
                        summaries.append(base|dict(window=wn,motif=motif,proximal_mean=pc.mean(),distal_mean=dc.mean(),difference_mean=diff.mean(),difference_median=diff.median(),distal_gt_proximal_fraction=diff.gt(0).mean(),distal_has_motif_fraction=dc.gt(0).mean(),proximal_has_motif_fraction=pc.gt(0).mean(),n_overlapping_upstream100=int(subset.upstream100_overlap_nt.gt(0).sum())))
                for c in cov:balance.append(base|dict(covariate=c,mean=subset[c].mean(),median=subset[c].median(),sd=subset[c].std(ddof=1)))
            for model in MODELS:
                z,cells,cols=partition(part,model);stratum_members.append(z.assign(background=background,model=model));random=np.random.default_rng(seed('CFIm_sequence',background,species,model));draws=label_draws(label,cells,random)
                probability=np.zeros((len(z),3));mobile=np.zeros(len(z),bool)
                for j,ii in enumerate(cells):
                    counts=np.bincount(label[ii],minlength=3);probability[ii]=counts/len(ii);mobile[ii]=np.count_nonzero(counts)>1
                    strata.append(context|dict(model=model,stratum=j,n=len(ii),proximal=counts[0],distal=counts[1],nonrecurrent=counts[2],mobile=bool(mobile[ii].any()),key='|'.join(map(str,z.loc[ii[0],cols].to_list()))))
                for wn in WINDOWS:
                    for motif in MOTIFS:
                        raw=z[f'{wn}_{motif}_difference'].to_numpy(float);values,design_rank=(residualize(raw,z) if model=='residualized_assay_architecture' else (raw,0))
                        for other,comparison in [(2,'proximal_vs_nonrecurrent'),(1,'proximal_vs_distal')]:
                            obs=stat(values,label,other);nul=statistic_draws(values,draws,other);exp=expectation(values,label,cells,other)
                            key='|'.join([background,species,model,wn,motif,comparison]);drawstore[key]=nul
                            rows.append(context|dict(model=model,window=wn,motif=motif,comparison=comparison,n_proximal=int((label==0).sum()),n_control=int((label==other).sum()),n_strata=len(cells),n_mobile_genes=int(mobile.sum()),n_mobile_proximal=int(((label==0)&mobile).sum()),n_mobile_control=int(((label==other)&mobile).sum()),proximal_assignment_informative=bool(((label==0)&mobile).any()),control_assignment_informative=bool(((label==other)&mobile).any()),residual_design_rank=design_rank,permutations=NPERM,seed=seed('CFIm_sequence',background,species,model),primary=(model=='joint_covariate_median' and wn=='upstream40_100' and motif=='TGTA'),**null_summary(obs,exp,nul)))
                            if background=='paired_strict' and wn=='upstream40_100' and motif=='TGTA':paired[species,model,comparison]=(obs,exp,nul,rows[-1])
            print(background,species,'n',len(part),'labels',np.bincount(label,minlength=3).tolist(),flush=True)
    results=pd.DataFrame(rows);results['q_family']=np.nan
    results.loc[results.primary,'q_family']=bh(results.loc[results.primary,'p_enrichment'])
    results.loc[~results.primary,'q_family']=bh(results.loc[~results.primary,'p_enrichment'])
    results['family']=np.where(results.primary,'primary_8','exploratory_208')
    assert results.primary.sum()==8 and (~results.primary).sum()==208
    save(results,'sequence_signature_tests.tsv');save(pd.DataFrame(summaries),'motif_group_summary.tsv');save(pd.DataFrame(balance),'covariate_group_summary.tsv');save(pd.DataFrame(strata),'stratum_margins.tsv.gz');save(pd.concat(stratum_members,ignore_index=True),'stratum_gene_members.tsv.gz')
    interactions=[]
    for model in MODELS:
        for comparison in ['proximal_vs_nonrecurrent','proximal_vs_distal']:
            ho,he,hn,hr=paired['human',model,comparison];mo,me,mn,mr=paired['mouse',model,comparison]
            obs=mo-ho;exp=me-he;nul=mn-hn;key='interaction|'+model+'|'+comparison;drawstore[key]=nul
            row=dict(model=model,comparison=comparison,n_paired_orthologs=hr['n_genes'],human_n_proximal=hr['n_proximal'],human_n_control=hr['n_control'],mouse_n_proximal=mr['n_proximal'],mouse_n_control=mr['n_control'],primary=model=='joint_covariate_median',species_interaction_proximal_assignment_informative=bool(hr['proximal_assignment_informative'] and mr['proximal_assignment_informative']),human_n_mobile_proximal=hr['n_mobile_proximal'],mouse_n_mobile_proximal=mr['n_mobile_proximal'],permutations=NPERM,**null_summary(obs,exp,nul))
            if not hr['informative'] or not mr['informative']:
                row.update(informative=False,status='at_least_one_species_uninformative',p_enrichment=1.,p_two_sided=1.)
            interactions.append(row)
    interactions=pd.DataFrame(interactions);interactions['q_two_sided']=np.nan
    for primary,g in interactions.groupby('primary'):interactions.loc[g.index,'q_two_sided']=bh(g.p_two_sided)
    save(interactions,'paired_structure_interaction.tsv')
    np.savez_compressed(HERE/'permutation_draws.npz',**drawstore)
    inputs=json.loads((HERE/'input_hashes.json').read_text());inputs.update({str(p.relative_to(ROOT)):digest(p) for p in [HERE/'prepare_inputs.py',HERE/'collect_sequences.py',HERE/'extract_remote.py',HERE/'analyze.py',HERE/'background_inventory.tsv.gz',HERE/'upstream100_sequences.tsv',ROOT/'notes/2026-10-10-directional-recurrence/observed_members.tsv.gz']})
    (HERE/'analysis_audit.json').write_text(json.dumps(dict(protocol_verified=True,eligibility_before_labels=True,roster_verified=True,main_tests=8,exploratory_tests=208,paired_interaction_tests=6,inputs=inputs),indent=2)+'\n')
    print(results[results.primary][['background','species','comparison','n_genes','n_proximal','n_control','n_mobile_genes','observed','analytic_expected','null_sd','p_enrichment','q_family','status']].to_string(index=False),flush=True)


if __name__=='__main__':main()
