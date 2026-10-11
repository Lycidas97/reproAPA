"""Unique one-base mappings with explicit interbase and strand semantics."""
import subprocess

import numpy as np
import pandas as pd


def lift_points(ctx,frame,chain,tag):
    assert frame.id.is_unique
    work=ctx.out/'liftover';work.mkdir(exist_ok=True)
    source=work/(tag+'.input.bed');mapped=work/(tag+'.mapped.bed');unmapped=work/(tag+'.unmapped.bed')
    start=frame.cut_site.astype(int)-frame.strand.eq('+').astype(int)
    assert start.ge(0).all()
    pd.DataFrame(dict(chrom=frame.chrom,start=start,end=start+1,name=frame.id,score=0,strand=frame.strand)).to_csv(source,sep='\t',index=False,header=False)
    command=[str(ctx.resource('liftover')), '-multiple','-minMatch=1',str(source),str(ctx.resource(chain)),str(mapped),str(unmapped)]
    with (work/(tag+'.log')).open('w') as log:
        subprocess.run(command,stdout=log,stderr=log,check=True)
    columns=['mapped_chrom','start0','end0','id','score','mapped_strand']
    m=pd.read_csv(mapped,sep='\t',names=columns) if mapped.stat().st_size else pd.DataFrame(columns=columns)
    assert set(m.id)<=set(frame.id)
    duplicate=set(m.loc[m.id.duplicated(False),'id'])
    invalid=set(m.loc[m.end0.sub(m.start0).ne(1),'id'])
    good=m[~m.id.isin(duplicate|invalid)].copy()
    good['mapped_cut']=np.where(good.mapped_strand.eq('+'),good.end0,good.start0)
    out=frame.merge(good[['id','mapped_chrom','mapped_strand','mapped_cut']],on='id',how='left',validate='one_to_one')
    out['mapping_status']=np.select([out.id.isin(duplicate),out.id.isin(invalid),out.mapped_cut.notna()],
        ['ambiguous','not_one_base','unique'],default='unmapped')
    ctx.table(out,'liftover/'+tag+'.audit.tsv.gz')
    return out
