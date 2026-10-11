"""Separate resource sensitivity on the exact same matched mouse coordinates."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT/'notes/2026-10-02-reviewer-followup'))
from apa_followup import read,save,rng_for,nullmod,INPUTS


def main():
    w=read(HERE/'motif_windows.tsv.gz');w=w[w.species.eq('mouse')&w.category.eq('RBP')].copy()
    h=read(ROOT/'notes/2026-10-01-figure5-union-update/full-cohort-v2/data/mouse/site_effects.tsv.gz',usecols=['gene_id','category','chrom','strand','start0','end0'])
    h=h[h.category.eq('RBP_transfer')].drop_duplicates();groups=dict(tuple(h.groupby('gene_id')))
    rows=[]
    for r in w.itertuples():
        g=groups.get(r.gene_id,h.iloc[:0])
        nc=int((g.start0.ge(r.start0)&g.end0.le(r.end0)).sum())
        ns=int((g.start0.ge(r.control_start0)&g.end0.le(r.control_end0)).sum())
        rows.append(dict(species='mouse',contrast=r.contrast,setting=r.setting,gene_id=r.gene_id,start0=r.start0,end0=r.end0,
            control_start0=r.control_start0,control_end0=r.control_end0,baseline_window=r.baseline_window,
            changed_sites=nc,control_sites=ns,length=r.length,density_difference=1000*(nc-ns)/r.length))
    w=pd.DataFrame(rows);save(w,HERE/'transfer_windows.tsv.gz')
    result=[]
    for cohort in ['dynamic','fixed_baseline_windows']:
        f=w if cohort=='dynamic' else w[w.baseline_window]
        g=f.groupby(['setting','contrast','gene_id']).density_difference.mean().reset_index()
        for (setting,contrast),b in g.groupby(['setting','contrast']):
            x=b.density_difference.to_numpy();nz=x[abs(x)>1e-12]
            row=dict(cohort=cohort,setting=setting,contrast=contrast,category='RBP_transfer',n_genes=len(x),n_nonzero=len(nz),mean_difference=x.mean())
            if len(x)>=10:
                rng=rng_for('RBP_transfer',setting,contrast)
                boot=x[rng.integers(len(x),size=(5000,len(x)))].mean(1)
                if len(nz)<=20:
                    sums=np.array([0.])
                    for v in nz:sums=np.r_[sums+v,sums-v]
                    p=float((abs(sums)>=abs(x.sum())-1e-10).mean())
                else:
                    null=(rng.choice([-1,1],size=(10000,len(x)))*x).mean(1)
                    p=(1+int((abs(null)>=abs(x.mean())-1e-12).sum()))/10001
                row.update(lo=np.quantile(boot,.025),hi=np.quantile(boot,.975),p=p)
            else:row.update(lo=np.nan,hi=np.nan,p=1.)
            result.append(row)
    d=pd.DataFrame(result);assert len(d)==48
    d['q_24']=d.groupby('cohort').p.transform(lambda p:nullmod.bh(p));save(d,HERE/'transfer_summary.tsv')
    (HERE/'transfer_input_sha256.json').write_text(json.dumps(INPUTS,indent=2)+'\n')
    print(d[d.cohort.eq('dynamic')&d.setting.isin(['GC0.05_POS0.2','GCinf_POSinf'])].to_string(index=False))
    print('Minimum p',d.p.min(),'minimum q',d.q_24.min())


if __name__=='__main__':main()
