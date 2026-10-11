"""Fixed target windows, explicit background ladder, and conditional power."""
import sys, json
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from common import *
sys.path.insert(0,str(ROOT/'notes/2026-10-02-reviewer-followup'))
from motif_background import fasta, runs, U, MODELS, ROUTES

def build_context():
    regions=read(U/'all_coordinate_regions.tsv.gz');contexts={};changes=[];hitgroups={}
    for species,routes in ROUTES.items():
        models={}
        for route in routes:
            p=MODELS/species/route/'models.json';INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();models[route]=json.loads(p.read_text())
        p=U/species/'query.fa';INPUTS[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest();seqs=fasta(p)
        hits=read(U/f'full-cohort-v2/data/{species}/site_effects.tsv.gz',usecols=['gene_id','category','chrom','strand','start0','end0'])
        hits=hits[hits.category.isin(['miRNA','RBP'])].drop_duplicates()
        hitgroups[species]=dict(tuple(hits.groupby(['gene_id','category'])))
        for (contrast,gene),group in regions[regions.species.eq(species)].groupby(['contrast','gene_id']):
            if set(group.route)!=set(routes):continue
            families=[models[r][gene]['family'] for r in routes]
            if len({(f['chrom'],f['strand']) for f in families})!=1:continue
            lo=max(f['common_start0'] for f in families);hi=min(f['common_end0'] for f in families)
            if hi<=lo:continue
            family=families[0]
            available=[m for m in family['models'] if m['region_id'] in seqs and m['start0']<=lo and m['end0']>=hi]
            if not available:continue
            m=available[0];seq=seqs[m['region_id']]
            if family['strand']=='-':seq=seq[::-1]
            seq=seq[lo-m['start0']:hi-m['start0']].upper();assert len(seq)==hi-lo
            gc=np.r_[0,np.cumsum(np.array(list(seq))=='G')+np.cumsum(np.array(list(seq))=='C')]
            cuts=sorted({lo,hi}|{max(lo,min(hi,int(x))) for x in np.r_[group.start0,group.end0]});parts=[]
            for a,b in zip(cuts[:-1],cuts[1:]):
                changed=True;stable=True;effects20=[];effects50=[]
                for route in routes:
                    block=group[group.route.eq(route)&group.start0.le(a)&group.end0.ge(b)]
                    if block.empty:
                        pos=np.array([p['cut_site']-1 for p in models[route][gene]['pas']])
                        constitutive=b<=pos.min() if family['strand']=='+' else a>=pos.max()
                        changed=False;stable &= constitutive
                        effects20.append(0. if constitutive else np.nan);effects50.append(0. if constitutive else np.nan)
                    else:
                        assert len(block)==1;row=block.iloc[0]
                        effects20.append(row.effect20);effects50.append(row.effect50)
                        changed &= bool(row.observed_both5 and row.strict_change and row.eligible)
                        stable &= bool(row.eligible and np.isfinite(row.effect20+row.effect50) and max(abs(row.effect20),abs(row.effect50))<.05)
                effect=np.r_[effects20,effects50]
                changed &= np.isfinite(effect).all() and (min(effect)>=.1 or max(effect)<=-.1)
                parts.append((a,b,'changed' if changed else 'stable' if stable else 'other'))
                if changed:
                    mid=((a+b)/2-lo)/(hi-lo)
                    changes.append(dict(species=species,contrast=contrast,gene_id=gene,start0=a,end0=b,length=b-a,
                        GC=(gc[b-lo]-gc[a-lo])/(b-a),relative_position=mid if family['strand']=='+' else 1-mid,
                        effect20=np.mean(effects20),effect50=np.mean(effects50),min_effect=min(effect),max_effect=max(effect)))
            contexts[(species,contrast,gene)]=dict(lo=lo,hi=hi,gc=gc,strand=family['strand'],stable=runs(parts,'stable'))
    return contexts,pd.DataFrame(changes),hitgroups

def candidates(context,length,target_gc,target_position):
    lo=context['lo'];hi=context['hi'];gc=context['gc']
    blocks=[np.arange(a,b-length+1,dtype=int) for a,b in context['stable'] if b-a>=length]
    if not blocks:return np.empty((0,4))
    starts=np.concatenate(blocks);fraction=(gc[starts+length-lo]-gc[starts-lo])/length
    positions=(starts+length/2-lo)/(hi-lo)
    if context['strand']=='-':positions=1-positions
    return np.column_stack([starts,fraction,abs(fraction-target_gc),abs(positions-target_position)])

def count_hits(hits,species,gene,category,a,b):
    h=hits[species].get((gene,category))
    return 0 if h is None else int((h.start0.ge(a)&h.end0.le(b)).sum())

def paired(x,key):
    n=len(x);nz=x[abs(x)>1e-12];row=dict(n_genes=n,n_nonzero=len(nz),mean_difference=x.mean() if n else np.nan)
    if n<10:return row|dict(p=1.,lo=np.nan,hi=np.nan,status='fewer_than_10_genes')
    random=rng('motif_summary',*key);boot=x[random.integers(n,size=(5000,n))].mean(1)
    if len(nz)<=20:
        sums=np.array([0.])
        for value in nz:sums=np.r_[sums+value,sums-value]
        p=float((abs(sums)>=abs(x.sum())-1e-10).mean())
    else:
        null=random.choice([-1,1],size=(20000,n))@x
        p=(1+int((abs(null)>=abs(x.sum())-1e-10).sum()))/20001
    return row|dict(p=p,lo=np.quantile(boot,.025),hi=np.quantile(boot,.975),status='paired_gene_signs')

def power_analysis():
    old=read(ROOT/'notes/2026-10-02-reviewer-followup/motif-background/gene_density.tsv');out=[]
    for key,f in old.groupby(['species','contrast','category']):
        x=f.sort_values('gene_id').density_difference.to_numpy();n=len(x)
        random=rng('motif_power',*key);noise=random.choice([-1,1],size=(1000,n))*abs(x)
        flips=random.choice([-1,1],size=(4095,n)).astype(float)
        sparse=np.zeros(n);sparse[random.choice(n,size=max(1,int(np.ceil(.25*n))),replace=False)]=1
        for prevalence,mask in [('all',np.ones(n)),('quarter',sparse)]:
            for delta in [0,2,5,10,20]:
                sim=noise+delta*mask
                if n>=10:
                    randomized=flips@sim.T;observed=sim.sum(1)
                    p=(1+(abs(randomized)>=abs(observed)[None,:]-1e-10).sum(0))/4096
                else:p=np.ones(1000)
                for alpha in [.05,.05/8]:
                    out.append(dict(species=key[0],contrast=key[1],category=key[2],n_genes=n,n_nonzero_noise=int((abs(x)>1e-12).sum()),
                        prevalence=prevalence,injected_genes=int(mask.sum()),injected_sites_per_kb=delta,
                        population_mean_shift=float(delta*mask.mean()),alpha=alpha,**rate(p<=alpha)))
    save(pd.DataFrame(out),'motif_power.tsv')

def main():
    context,changed,hits=build_context()
    old=read(ROOT/'notes/2026-10-02-reviewer-followup/motif-background/matched_windows.tsv.gz')
    old=old[old.category.eq('miRNA')]
    keys=['species','contrast','gene_id','start0','end0']
    target=old[keys].merge(changed,on=keys,validate='one_to_one');assert len(target)==len(old)
    attrition=read(ROOT/'notes/2026-10-02-null-result-diagnostics/motif_attrition.tsv')
    assert set(map(tuple,changed[keys].values))==set(map(tuple,attrition[keys].values))
    windows=[];assignment=[];edges=[]
    for (species,contrast),t in target.groupby(['species','contrast']):
        all_targets=set(changed[changed.species.eq(species)&changed.contrast.eq(contrast)].gene_id)
        donor_ids=sorted(g for sp,c,g in context if sp==species and c==contrast and g not in all_targets and context[(sp,c,g)]['stable'])
        # Randomized order is fixed without looking at counts or effects.
        random=rng('donor_assignment',species,contrast);donor_ids=list(random.permutation(donor_ids))
        query_ids=list(random.permutation(sorted(t.gene_id.unique())))
        adjacency=np.zeros((len(query_ids),len(donor_ids)),dtype=np.int8)
        for i,gene in enumerate(query_ids):
            for j,donor in enumerate(donor_ids):
                valid=True
                for q in t[t.gene_id.eq(gene)].itertuples():
                    c=candidates(context[(species,contrast,donor)],q.length,q.GC,q.relative_position)
                    if not len(c) or not ((c[:,2]<=.05)&(c[:,3]<=.2)).any():valid=False;break
                adjacency[i,j]=valid
                if valid:edges.append(dict(species=species,contrast=contrast,gene_id=gene,donor_gene=donor))
        match=maximum_bipartite_matching(csr_matrix(adjacency),perm_type='column')
        assert len(match)==len(query_ids)
        assigned=[donor_ids[j] for j in match if j>=0];assert len(assigned)==len(set(assigned)) and not set(assigned)&all_targets
        for i,gene in enumerate(query_ids):
            donor=donor_ids[match[i]] if match[i]>=0 else None
            assignment.append(dict(species=species,contrast=contrast,gene_id=gene,donor_gene=donor,
                compatible_donors=int(adjacency[i].sum()),matched=donor is not None))
            for q in t[t.gene_id.eq(gene)].itertuples():
                for source,control_gene in [('same_gene',gene),('cross_gene',donor)]:
                    if control_gene is None:continue
                    c=candidates(context[(species,contrast,control_gene)],q.length,q.GC,q.relative_position)
                    priority=rng('window_priority',species,contrast,gene,control_gene,q.start0,q.end0).random(len(c))
                    for level,mask in [('length',np.ones(len(c),bool)),('GC',c[:,2]<=.05),('GC_position',(c[:,2]<=.05)&(c[:,3]<=.2))]:
                        assert mask.any();selected=np.flatnonzero(mask)[np.argmin(priority[mask])];a=int(c[selected,0]);b=a+q.length
                        for category in ['miRNA','RBP']:
                            nc=count_hits(hits,species,gene,category,q.start0,q.end0);ns=count_hits(hits,species,control_gene,category,a,b)
                            windows.append(dict(species=species,contrast=contrast,gene_id=gene,category=category,source=source,level=level,
                                common_cohort=donor is not None,start0=q.start0,end0=q.end0,length=q.length,control_gene=control_gene,
                                control_start0=a,control_end0=b,GC_difference=c[selected,2],position_difference=c[selected,3],
                                changed_sites=nc,control_sites=ns,density_difference=1000*(nc-ns)/q.length))
        print('MOTIF',species,contrast,'targets',len(query_ids),'donors',len(donor_ids),'matched',len(assigned),flush=True)
    w=pd.DataFrame(windows);summaries=[];gene_parts=[]
    for cohort in ['same_gene_original_fixed','six_background_common']:
        subset=w[w.source.eq('same_gene')] if cohort=='same_gene_original_fixed' else w[w.common_cohort]
        gene=subset.groupby(['species','contrast','category','source','level','gene_id']).density_difference.mean().reset_index()
        gene_parts.append(gene.assign(cohort=cohort))
        for key,g in gene.groupby(['species','contrast','category','source','level']):
            summaries.append(dict(zip(['species','contrast','category','source','level'],key))|dict(cohort=cohort)|paired(g.density_difference.to_numpy(),(cohort,*key)))
    s=pd.DataFrame(summaries);s['q_all_backgrounds']=s.groupby('cohort').p.transform(bh)
    s['q_eight']=s.groupby(['cohort','source','level']).p.transform(bh)
    save(w,'motif_windows.tsv.gz');save(pd.concat(gene_parts),'motif_genes.tsv.gz');save(s,'motif_ladder.tsv')
    save(pd.DataFrame(assignment),'motif_donor_assignment.tsv');save(pd.DataFrame(edges),'motif_donor_candidates.tsv.gz')
    # Quantity of predicted sites affected by usage; no density-enrichment test.
    retention=[]
    for q in changed.itertuples():
        for category in ['miRNA','RBP']:
            n=count_hits(hits,q.species,q.gene_id,category,q.start0,q.end0)
            h=hits[q.species].get((q.gene_id,category))
            overlap=0 if h is None else int((h.start0.lt(q.end0)&h.end0.gt(q.start0)).sum())
            retention.append(dict(species=q.species,contrast=q.contrast,gene_id=q.gene_id,category=category,start0=q.start0,end0=q.end0,
                length=q.length,sites=n,overlapping_not_contained=overlap-n,effect20=q.effect20,effect50=q.effect50,
                weighted_sites20=n*q.effect20,weighted_sites50=n*q.effect50))
    save(changed,'motif_all_changed_regions.tsv.gz');save(pd.DataFrame(retention),'motif_site_retention.tsv.gz')
    power_analysis();manifest('motif_manifest.json',status='PASS',target_windows=len(target),all_changed_regions=len(changed),
        inference='paired gene/unique donor blocks; conditional exchangeability not proven',power='statistic-level symmetric-noise injection')

if __name__=='__main__':main()
