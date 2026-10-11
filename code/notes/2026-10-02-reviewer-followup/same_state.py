"""Literal S-vs-S and P-vs-P halves through the frozen library weighting."""
import argparse
import numpy as np
import pandas as pd
from apa_followup import HERE, OLD, BASE, SPECS, read, save, load_study, legacy, rng_for


def extract(study):
    meta,axes,counts=load_study(study)
    pairs=read(OLD / 'pair_effects.tsv.gz');pairs=pairs[pairs.study_id.eq(study)]
    groups={g:np.asarray(ix) for g,ix in axes.groupby('primary_gene_id').groups.items()}
    mrows={k:np.flatnonzero(meta.route.eq(k[0])&meta.library_key.eq(k[1])&meta.node.eq(k[2])&meta.purity.ge(.8))
           for k in {(r.route,r.library_key_a,'P') for r in pairs.itertuples()} | {(r.route,r.library_key_b,r.contrast[-1]) for r in pairs.itertuples()}}
    output=[]
    for gene,gp in pairs.groupby('gene_id'):
        columns=groups[gene];c=counts[:,columns].toarray();w=axes.position.iloc[columns].to_numpy();cache={}
        for k in {(r.route,r.library_key_a,'P') for r in gp.itertuples()} | {(r.route,r.library_key_b,r.contrast[-1]) for r in gp.itertuples()}:
            ix=mrows[k];v=c[ix];cells=meta.n_cells.iloc[ix].to_numpy();keep=v.sum(1)>=20;v=v[keep];cells=cells[keep]
            if len(v)<6 or cells.sum()<60 or (((v.sum(0)>=10)&((v>0).sum(0)>=2)).sum()<2):
                continue
            p=v/v.sum(1)[:,None];rng=rng_for('literal_halves',study,gene,*k)
            order=rng.permuted(np.broadcast_to(np.arange(len(v)),(100,len(v))),axis=1);n=len(v)//2;a,b=order[:,:n],order[:,n:]
            eligible=(cells[a].sum(1)>=30)&(cells[b].sum(1)>=30)
            delta=p[a].mean(1)-p[b].mean(1);delta[~eligible]=np.nan
            cache[k]=delta
        for key,part in gp.groupby(BASE):
            for which in ['reference','target']:
                # Both halves come from the SAME real state; each state draw is reused across comparisons.
                keys=[(r.route,r.library_key_b,r.contrast[-1]) if which=='reference' else (r.route,r.library_key_a,'P') for r in part.itertuples()]
                keep=np.array([k in cache for k in keys]);q=part.loc[keep].copy()
                if q.empty:continue
                q=legacy.pair_weights(q.drop(columns='weight'));pw=q.weight.to_numpy();d=sum(ww*cache[k] for ww,k in zip(pw,[k for k,ok in zip(keys,keep) if ok]))
                for rep in range(100):
                    output.append(dict(zip(BASE,key)) | dict(pseudo_kind=which,replicate=rep,delta_WPP=float(d[rep]@w),MPRO=float(np.ptp(d[rep])),n_pairs=len(q)))
    save(pd.DataFrame(output),HERE / f'same-state/{study}.tsv.gz')
    print(study,len(output),flush=True)


def summarize():
    pseudo=pd.concat([read(p) for p in sorted((HERE / 'same-state').glob('GSE*.tsv.gz'))],ignore_index=True)
    observed=read(OLD / 'route_effects.tsv.gz')
    rows,noise=[] ,[]
    for kind in ['reference','target']:
        for (species,contrast),p in pseudo[pseudo.pseudo_kind.eq(kind)].groupby(['species','contrast']):
            for metric,thresholds in [('WPP',[.02,.05]),('MPRO',[.1,.2])]:
                field='delta_WPP' if metric=='WPP' else 'MPRO'
                for rep,frame in p.groupby('replicate'):
                    frame=frame[np.isfinite(frame[field])]
                    # All original series must have eligible halves; no partial-series recollapse.
                    expected=observed[observed.species.eq(species)&observed.contrast.eq(contrast)].groupby(['study_id','gene_id']).size()
                    available=frame.groupby(['study_id','gene_id']).size();good=available[available.eq(expected.reindex(available.index))].index
                    idx=pd.MultiIndex.from_frame(frame[['study_id','gene_id']]);frame=frame[idx.isin(good)]
                    sub=observed.merge(frame[BASE],on=BASE,validate='one_to_one')
                    for threshold in thresholds:
                        def recurrent(d):
                            f=d.groupby(['study_id','gene_id'])[field].agg(['min','max'])
                            plus=f['min'].ge(threshold);minus=f['max'].le(-threshold) if metric=='WPP' else pd.Series(False,index=f.index)
                            votes=pd.DataFrame({'plus':plus,'minus':minus}).groupby('gene_id').sum();nunits=f.reset_index().groupby('gene_id').study_id.nunique()
                            return int(((votes.plus>=2)|(votes.minus>=2)).sum()),int((nunits>=2).sum())
                        n,eligible=recurrent(frame);obs,oe=recurrent(sub);assert eligible==oe
                        rows.append(dict(pseudo_kind=kind,species=species,contrast=contrast,metric=metric,threshold=threshold,replicate=rep,
                                         eligible_cross_study_genes=eligible,pseudo_recurrent=n,observed_on_same_background=obs))
                values=p[field].dropna().to_numpy();noise.append(dict(pseudo_kind=kind,species=species,contrast=contrast,metric=metric,
                    n_pseudo_routes=len(values),absolute_median=float(np.median(abs(values))),absolute_p90=float(np.quantile(abs(values),.9)),absolute_p95=float(np.quantile(abs(values),.95)),
                    lower_threshold_fraction=float((abs(values)>=(.02 if metric=='WPP' else .1)).mean())))
    result=pd.DataFrame(rows);save(result,HERE / 'same-state/recurrence_draws.tsv.gz');save(pd.DataFrame(noise),HERE / 'same-state/noise_thresholds.tsv')
    summary=[]
    keys=['pseudo_kind','species','contrast','metric','threshold']
    for key,f in result.groupby(keys):
        summary.append(dict(zip(keys,key)) | dict(eligible_median=float(f.eligible_cross_study_genes.median()),pseudo_median=float(f.pseudo_recurrent.median()),pseudo_p95=float(f.pseudo_recurrent.quantile(.95)),
            observed_same_background_median=float(f.observed_on_same_background.median()),observed_exceeds_pseudo_fraction=float((f.observed_on_same_background>f.pseudo_recurrent).mean())))
    save(pd.DataFrame(summary),HERE / 'same-state/recurrence_summary.tsv');print(pd.DataFrame(summary).to_string(index=False))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['extract','summarize']);p.add_argument('--study',choices=list(SPECS));a=p.parse_args()
    extract(a.study) if a.command=='extract' else summarize()
