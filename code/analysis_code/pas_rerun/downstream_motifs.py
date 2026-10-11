"""New atomic windows, fixed motif-blind controls and gene-level inference."""
import argparse
import json
import pickle
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import maximum_bipartite_matching
from scipy.stats import beta

from downstream_runtime import Context,bh,seed
from downstream_regulatory import ROUTES,observed_regions

GRID=[(g,p) for g in [.05,.10,.20,np.inf] for p in [.20,.40,np.inf]]
BGS='notes/2026-10-02-reviewer-followup/motif_background.py'
CTL='notes/2026-10-02-null-model-controls/motif.py'


def label(g,p):return f'GC{g:g}_POS{p:g}'
def rng(*keys):return np.random.default_rng(seed(*keys))


def helpers(ctx):
    num=ctx.functions(CTL,['candidates','count_hits','paired'],{'rng':rng})
    read=ctx.functions(BGS,['fasta','runs'])
    choose=ctx.functions('notes/2026-10-02-null-result-diagnostics/motif_sensitivity.py',['choose'])
    return num,read,choose


def build_context(ctx):
    num,read,choose=helpers(ctx);regions=observed_regions(ctx.read(ctx.analysis/'intervals/all_coordinate_regions.tsv.gz'))
    regions=regions[regions.contrast.isin(['SP','MP'])];contexts={};changed=[];hitgroups={}
    for sp,routes in ROUTES.items():
        names=[r for sid,r in routes]
        models={r:json.loads(ctx.register(ctx.analysis/f'intervals/annotation/{sp}/{r}/models.json').read_text()) for r in names}
        seqs=read.fasta(ctx.register(ctx.analysis/f'site_scan/{sp}/query.fa'))
        hits=ctx.read(ctx.analysis/f'site_scan/full-cohort-v2/data/{sp}/site_effects.tsv.gz',usecols=['gene_id','category','chrom','strand','start0','end0']).drop_duplicates()
        hitgroups[sp]=dict(tuple(hits.groupby(['gene_id','category'])))
        for (co,gene),group in regions[regions.species.eq(sp)].groupby(['contrast','gene_id']):
            if set(group.route)!=set(names):continue
            families=[models[r][gene]['family'] for r in names]
            if len({(f['chrom'],f['strand']) for f in families})!=1:continue
            lo=max(f['common_start0'] for f in families);hi=min(f['common_end0'] for f in families)
            if hi<=lo:continue
            family=families[0];available=[m for m in family['models'] if m['region_id'] in seqs and m['start0']<=lo and m['end0']>=hi]
            if not available:continue
            m=available[0];seq=seqs[m['region_id']]
            # Reversal preserves GC, and makes array indexing genomic ascending.
            if family['strand']=='-':seq=seq[::-1]
            seq=seq[lo-m['start0']:hi-m['start0']].upper();assert len(seq)==hi-lo
            letters=np.array(list(seq));gc=np.r_[0,np.cumsum((letters=='G')|(letters=='C'))]
            cuts=sorted({lo,hi}|{max(lo,min(hi,int(v))) for v in np.r_[group.start0,group.end0]});parts=[]
            for a,b in zip(cuts[:-1],cuts[1:]):
                change=True;stable=True;effects20=[];effects50=[]
                for route in names:
                    block=group[group.route.eq(route)&group.start0.le(a)&group.end0.ge(b)]
                    if block.empty:
                        pos=np.array([p['cut_interbase'] for p in models[route][gene]['pas']])
                        constitutive=b<=pos.min() if family['strand']=='+' else a>=pos.max()
                        change=False;stable &= constitutive
                        effects20.append(0. if constitutive else np.nan);effects50.append(0. if constitutive else np.nan)
                    else:
                        assert len(block)==1;row=block.iloc[0];effects20.append(row.effect20);effects50.append(row.effect50)
                        change &= bool(row.observed_both5 and row.strict_change)
                        stable &= bool(np.isfinite(row.effect20+row.effect50) and max(abs(row.effect20),abs(row.effect50))<.05)
                effect=np.r_[effects20,effects50];change &= np.isfinite(effect).all() and (min(effect)>=.1 or max(effect)<=-.1)
                parts.append((a,b,'changed' if change else 'stable' if stable else 'other'))
                if change:
                    mid=((a+b)/2-lo)/(hi-lo)
                    changed.append(dict(species=sp,contrast=co,gene_id=gene,start0=a,end0=b,length=b-a,GC=(gc[b-lo]-gc[a-lo])/(b-a),
                        relative_position=mid if family['strand']=='+' else 1-mid,effect20=np.mean(effects20),effect50=np.mean(effects50),min_effect=min(effect),max_effect=max(effect)))
            contexts[sp,co,gene]=dict(lo=lo,hi=hi,gc=gc,strand=family['strand'],stable=read.runs(parts,'stable'))
    return contexts,pd.DataFrame(changed,columns=['species','contrast','gene_id','start0','end0','length','GC','relative_position','effect20','effect50','min_effect','max_effect']),hitgroups


def density_row(q,c,cat,num,hits,control_gene=None):
    control_gene=q.gene_id if control_gene is None else control_gene;start=int(c[0]);end=start+q.length
    nc=num.count_hits(hits,q.species,q.gene_id,cat,q.start0,q.end0);ns=num.count_hits(hits,q.species,control_gene,cat,start,end)
    return dict(species=q.species,contrast=q.contrast,gene_id=q.gene_id,category=cat,start0=q.start0,end0=q.end0,length=q.length,
        control_gene=control_gene,control_start0=start,control_end0=end,GC_change=q.GC,GC_control=c[1],GC_difference=c[2],relative_position_distance=c[3],
        changed_sites=nc,control_sites=ns,changed_density=1000*nc/q.length,control_density=1000*ns/q.length,density_difference=1000*(nc-ns)/q.length)


WINDOW=['species','contrast','gene_id','category','start0','end0','length','control_gene','control_start0','control_end0','GC_change','GC_control','GC_difference','relative_position_distance','changed_sites','control_sites','changed_density','control_density','density_difference']


def primary(ctx):
    num,read,choose=helpers(ctx);contexts,changed,hits=build_context(ctx);rows=[];attrition=[]
    for q in changed.itertuples():
        c=num.candidates(contexts[q.species,q.contrast,q.gene_id],q.length,q.GC,q.relative_position);picked=choose.choose(c,.05,.2)
        context=contexts[q.species,q.contrast,q.gene_id]
        status=('no_stable_sequence' if not context['stable'] else 'no_equal_length_window' if not len(c) else
                'no_GC_match' if not (c[:,2]<=.05).any() else 'no_joint_position_match' if picked is None else 'matched')
        attrition.append(dict(species=q.species,contrast=q.contrast,gene_id=q.gene_id,start0=q.start0,end0=q.end0,length=q.length,status=status,
            stable_bp=sum(b-a for a,b in context['stable']),longest_stable=max([b-a for a,b in context['stable']],default=0),
            equal_length_windows=len(c),gc_alone_windows=int((c[:,2]<=.05).sum()),position_alone_windows=int((c[:,3]<=.2).sum()),joint_windows=int(((c[:,2]<=.05)&(c[:,3]<=.2)).sum())))
        if picked is not None:
            for cat in ['miRNA','RBP']:rows.append(density_row(q,picked,cat,num,hits))
    frame=pd.DataFrame(rows,columns=WINDOW);genes=frame.groupby(['species','contrast','category','gene_id'])[['density_difference','changed_density','control_density']].mean().reset_index()
    results=[]
    for sp in ROUTES:
        for co in ['SP','MP']:
            for cat in ['miRNA','RBP']:
                g=genes[genes.species.eq(sp)&genes.contrast.eq(co)&genes.category.eq(cat)];x=g.density_difference.to_numpy();random=rng('motif',sp,co,cat)
                row=dict(species=sp,contrast=co,category=cat,n_genes=len(x),mean_difference=x.mean() if len(x) else np.nan,
                    mean_changed_density=g.changed_density.mean(),mean_control_density=g.control_density.mean(),lo=np.nan,hi=np.nan,p=1.,status='fewer_than_10_genes_descriptive_only')
                if len(x)>=10:
                    boot=x[random.integers(len(x),size=(5000,len(x)))].mean(1);null=(random.choice([-1,1],size=(10000,len(x)))*x).mean(1)
                    row.update(lo=np.quantile(boot,.025),hi=np.quantile(boot,.975),p=(1+(abs(null)>=abs(x.mean())-1e-12).sum())/10001,status='matched_gene_permutation')
                results.append(row)
    out=pd.DataFrame(results);out['q_8_tests']=bh(out.p)
    ctx.table(frame,'matched_windows.tsv.gz');ctx.table(genes,'gene_density.tsv');ctx.table(out,'density_summary.tsv')
    ctx.table(attrition,'motif_attrition.tsv',columns=['species','contrast','gene_id','start0','end0','length','status']);ctx.table(changed,'all_changed_regions.tsv.gz')
    with (ctx.out/'window_context.pickle').open('wb') as handle:pickle.dump((contexts,changed,hits),handle)
    ctx.finish(dict(new_atomic_intervals=True,interbase_boundaries=True,matching_blind_to_motifs=True,all_unmatched_reported=True,planned_eight_tests=True))


def controls(ctx):
    num,read,choose=helpers(ctx)
    with ctx.register(ctx.analysis/'motif_background/window_context.pickle').open('rb') as handle:context,changed,hits=pickle.load(handle)
    windows=[];transfer=[]
    for q in changed.itertuples():
        candidates=num.candidates(context[q.species,q.contrast,q.gene_id],q.length,q.GC,q.relative_position)
        baseline=choose.choose(candidates,.05,.2) is not None
        for g,p in GRID:
            c=choose.choose(candidates,g,p)
            if c is None:continue
            for cat in ['miRNA','RBP']+(['RBP_transfer'] if q.species=='mouse' else []):
                row=density_row(q,c,cat,num,hits)|dict(setting=label(g,p),baseline_window=baseline)
                (transfer if cat=='RBP_transfer' else windows).append(row)
    w=pd.DataFrame(windows+transfer,columns=WINDOW+['setting','baseline_window']);ctx.table(w,'grid_windows.tsv.gz')
    results=[];gene_frames=[]
    for cohort in ['dynamic','fixed_baseline_windows']:
        f=w if cohort=='dynamic' else w.loc[w.baseline_window.astype(bool)]
        genes=f.groupby(['setting','species','contrast','category','gene_id'])[['density_difference','changed_density','control_density']].mean().reset_index();gene_frames.append(genes.assign(cohort=cohort))
        for g,p in GRID:
            for sp in ROUTES:
                for co in ['SP','MP']:
                    for cat in ['miRNA','RBP']+(['RBP_transfer'] if sp=='mouse' else []):
                        z=genes[genes.setting.eq(label(g,p))&genes.species.eq(sp)&genes.contrast.eq(co)&genes.category.eq(cat)];x=z.density_difference.to_numpy()
                        result=dict(cohort=cohort,setting=label(g,p),species=sp,contrast=co,category=cat,**num.paired(x,(cohort,g,p,sp,co,cat)),
                            mean_changed_density=z.changed_density.mean(),mean_control_density=z.control_density.mean(),n_positive=int((x>1e-12).sum()),n_negative=int((x< -1e-12).sum()),
                            loo_min=(x.sum()-x).min()/(len(x)-1) if len(x)>1 else np.nan,loo_max=(x.sum()-x).max()/(len(x)-1) if len(x)>1 else np.nan)
                        results.append(result)
    s=pd.DataFrame(results);s['resource_family']=np.where(s.category.eq('RBP_transfer'),'transfer24','direct96')
    s['q_all_settings']=s.groupby(['cohort','resource_family']).p.transform(bh);s['q_within_setting']=s.groupby(['cohort','resource_family','setting']).p.transform(bh)
    ctx.table(s,'grid_summary.tsv');ctx.table(pd.concat(gene_frames),'grid_gene_density.tsv.gz')
    targets=ctx.read(ctx.analysis/'motif_background/matched_windows.tsv.gz');keys=['species','contrast','gene_id','start0','end0']
    target=targets[targets.category.eq('miRNA')][keys].merge(changed,on=keys,validate='one_to_one');assert len(target)==len(targets[targets.category.eq('miRNA')])
    backgrounds=[];assignment=[];edges=[]
    for (sp,co),t in target.groupby(['species','contrast']):
        excluded=set(changed[changed.species.eq(sp)&changed.contrast.eq(co)].gene_id);random=rng('donor_assignment',sp,co)
        donors=list(random.permutation(sorted(g for species,contrast,g in context if species==sp and contrast==co and g not in excluded and context[species,contrast,g]['stable'])))
        genes=list(random.permutation(sorted(t.gene_id.unique())));adj=np.zeros((len(genes),len(donors)),dtype=np.int8)
        for i,gene in enumerate(genes):
            for j,donor in enumerate(donors):
                okay=True
                for q in t[t.gene_id.eq(gene)].itertuples():
                    c=num.candidates(context[sp,co,donor],q.length,q.GC,q.relative_position)
                    if not len(c) or not ((c[:,2]<=.05)&(c[:,3]<=.2)).any():okay=False;break
                adj[i,j]=okay
                if okay:edges.append(dict(species=sp,contrast=co,gene_id=gene,donor_gene=donor))
        match=maximum_bipartite_matching(csr_matrix(adj),perm_type='column') if len(donors) else np.full(len(genes),-1)
        assigned=[donors[j] for j in match if j>=0];assert len(assigned)==len(set(assigned)) and not set(assigned)&excluded
        for i,gene in enumerate(genes):
            donor=donors[match[i]] if match[i]>=0 else None
            assignment.append(dict(species=sp,contrast=co,gene_id=gene,donor_gene=donor,compatible_donors=int(adj[i].sum()),matched=donor is not None))
            for q in t[t.gene_id.eq(gene)].itertuples():
                for source,control_gene in [('same_gene',gene),('cross_gene',donor)]:
                    if control_gene is None:continue
                    c=num.candidates(context[sp,co,control_gene],q.length,q.GC,q.relative_position)
                    priority=rng('window_priority',sp,co,gene,control_gene,q.start0,q.end0).random(len(c))
                    for level,mask in [('length',np.ones(len(c),bool)),('GC',c[:,2]<=.05),('GC_position',(c[:,2]<=.05)&(c[:,3]<=.2))]:
                        assert mask.any();pick=np.flatnonzero(mask)[np.argmin(priority[mask])]
                        for cat in ['miRNA','RBP']:backgrounds.append(density_row(q,c[pick],cat,num,hits,control_gene)|dict(source=source,level=level,common_cohort=donor is not None))
    back=pd.DataFrame(backgrounds,columns=WINDOW+['source','level','common_cohort']);summaries=[];gene_parts=[]
    for cohort in ['same_gene_original_fixed','six_background_common']:
        subset=back[back.source.eq('same_gene')] if cohort=='same_gene_original_fixed' else back.loc[back.common_cohort.astype(bool)]
        g=subset.groupby(['species','contrast','category','source','level','gene_id']).density_difference.mean().reset_index();gene_parts.append(g.assign(cohort=cohort))
        for sp in ROUTES:
            for co in ['SP','MP']:
                for cat in ['miRNA','RBP']:
                    for source in ['same_gene'] if cohort=='same_gene_original_fixed' else ['same_gene','cross_gene']:
                        for level in ['length','GC','GC_position']:
                            z=g[g.species.eq(sp)&g.contrast.eq(co)&g.category.eq(cat)&g.source.eq(source)&g.level.eq(level)]
                            key=dict(cohort=cohort,species=sp,contrast=co,category=cat,source=source,level=level)
                            summaries.append(key|num.paired(z.density_difference.to_numpy(),tuple(key.values())))
    summary=pd.DataFrame(summaries);summary['q_all_backgrounds']=summary.groupby('cohort').p.transform(bh);summary['q_eight']=summary.groupby(['cohort','source','level']).p.transform(bh)
    ctx.table(back,'background_windows.tsv.gz');ctx.table(pd.concat(gene_parts),'background_genes.tsv.gz');ctx.table(summary,'background_ladder.tsv')
    ctx.table(assignment,'donor_assignment.tsv',columns=['species','contrast','gene_id','donor_gene','compatible_donors','matched'])
    ctx.table(edges,'donor_candidates.tsv.gz',columns=['species','contrast','gene_id','donor_gene'])
    power(ctx)
    ctx.finish(dict(new_full_grid=True,planned_direct96_and_transfer24=True,unique_external_donors=True,matching_blind_to_motifs=True,
        same_targets_six_backgrounds=True,planned_48_background_tests=True,full_and_quarter_signal_injection=True))


def power(ctx):
    gene=ctx.read(ctx.analysis/'motif_background/gene_density.tsv');out=[]
    rates=ctx.functions('notes/2026-10-02-null-model-controls/common.py',['rate'],{'beta':beta})
    for sp in ROUTES:
        for co in ['SP','MP']:
            for cat in ['miRNA','RBP']:
                x=gene[gene.species.eq(sp)&gene.contrast.eq(co)&gene.category.eq(cat)].sort_values('gene_id').density_difference.to_numpy();n=len(x)
                random=rng('motif_power',sp,co,cat);noise=random.choice([-1,1],size=(1000,n))*abs(x);flips=random.choice([-1,1],size=(4095,n)).astype(float)
                quarter=np.zeros(n)
                if n:quarter[random.choice(n,size=max(1,int(np.ceil(.25*n))),replace=False)]=1
                for prevalence,mask in [('all',np.ones(n)),('quarter',quarter)]:
                    for delta in [0,2,5,10,20]:
                        sim=noise+delta*mask
                        p=(1+(abs(flips@sim.T)>=abs(sim.sum(1))[None,:]-1e-10).sum(0))/4096 if n>=10 else np.ones(1000)
                        for alpha in [.05,.05/8]:out.append(dict(species=sp,contrast=co,category=cat,n_genes=n,n_nonzero_noise=int((abs(x)>1e-12).sum()),
                            prevalence=prevalence,injected_genes=int(mask.sum()),injected_sites_per_kb=delta,population_mean_shift=float(delta*mask.mean()) if n else np.nan,
                            alpha=alpha,**rates.rate(p<=alpha)))
    ctx.table(out,'motif_power.tsv')


def retention(ctx):
    num,read,choose=helpers(ctx)
    with ctx.register(ctx.analysis/'motif_background/window_context.pickle').open('rb') as handle:context,changed,hits=pickle.load(handle)
    rows=[]
    for q in changed.itertuples():
        for cat in ['miRNA','RBP']:
            n=num.count_hits(hits,q.species,q.gene_id,cat,q.start0,q.end0);h=hits[q.species].get((q.gene_id,cat))
            overlap=0 if h is None else int((h.start0.lt(q.end0)&h.end0.gt(q.start0)).sum())
            rows.append(dict(species=q.species,contrast=q.contrast,gene_id=q.gene_id,category=cat,start0=q.start0,end0=q.end0,length=q.length,
                sites=n,overlapping_not_contained=overlap-n,effect20=q.effect20,effect50=q.effect50,weighted_sites20=n*q.effect20,weighted_sites50=n*q.effect50))
    f=pd.DataFrame(rows,columns=['species','contrast','gene_id','category','start0','end0','length','sites','overlapping_not_contained','effect20','effect50','weighted_sites20','weighted_sites50'])
    gene=f.groupby(['species','contrast','category','gene_id'])[['sites','weighted_sites20','weighted_sites50']].sum().reset_index();summary=[]
    for sp in ROUTES:
        for co in ['SP','MP']:
            for cat in ['miRNA','RBP']:
                z=gene[gene.species.eq(sp)&gene.contrast.eq(co)&gene.category.eq(cat)];eligible=z[z.sites.gt(0)]
                summary.append(dict(species=sp,contrast=co,category=cat,genes=len(z),genes_with_complete_sites=len(eligible),
                    positive=int(eligible.weighted_sites20.gt(0).sum()),negative=int(eligible.weighted_sites20.lt(0).sum()),zero=int(eligible.weighted_sites20.eq(0).sum()),
                    net_weighted20=eligible.weighted_sites20.sum(),net_weighted50=eligible.weighted_sites50.sum()))
    ctx.table(f,'interval_retention.tsv.gz');ctx.table(gene,'gene_retention.tsv');ctx.table(summary,'retention_summary.tsv')
    ctx.finish(dict(all_changed_intervals_retained=True,matched_and_unmatched_included=True,complete_sites_only=True,zero_gene_signs_separate=True))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True);p.add_argument('--node',choices=['motif_background','motif_controls','weighted_retention'],required=True)
    a=p.parse_args();{'motif_background':primary,'motif_controls':controls,'weighted_retention':retention}[a.node](Context(a.config,a.node,a.out))
