"""Estimability flags and fixed-cohort, six-well sign-resolution diagnostics."""
import json
import numpy as np
import pandas as pd
from downstream_runtime import bh

METRICS=['DWUI','RWUI','PDUI','PPUI','MPRO']


def confidence_intervals(frame):
    """A zero observed variance does not define a finite-df Student interval."""
    out=frame.copy()
    supplied=out.SE.notna()
    valid=supplied&out.SE.gt(0)&np.isfinite(out.df.to_numpy(dtype=float))&out.df.gt(0)
    out['CI_status']=np.where(valid,'estimable',np.where(supplied,
        'degenerate_observed_between_unit_variance','not_provided'))
    out.loc[~valid,['CI95_low','CI95_high']]=np.nan
    return out


def six_well_indices(frame,date_priority,unit_priority):
    groups=frame.groupby('time_a').indices
    dates=sorted(groups,key=date_priority.get)[:4]
    assert len(dates)==4
    ix=[];weights=[]
    for day,count in zip(dates,[2,2,1,1]):
        selected=sorted(groups[day],key=lambda i:unit_priority[frame.sample_unit_a.iloc[i]])[:count]
        assert len(selected)==count
        ix.extend(selected);weights.extend([1/4/count]*count)
    return np.asarray(ix),np.asarray(weights),dates


def exact_vector_test(delta,weights,linear):
    patterns=2*((np.arange(64)[:,None]>>np.arange(6))&1)-1
    v=delta*weights[:,None];mean=v.sum(0);effect=np.r_[mean@linear,np.ptp(mean)]
    transformed=patterns@v
    null=np.column_stack([np.abs(transformed@linear),np.ptp(transformed,axis=1)])
    p=(null>=np.r_[np.abs(effect[:4]),effect[4]]-1e-12).mean(0)
    assert (p>=2/64).all()
    return effect,p


def run(ctx,draws=100):
    baseline=ctx.read(ctx.analysis/'well_inference/gene_metric_tests.tsv.gz')
    audited=confidence_intervals(baseline)
    ctx.table(audited,'well-metric-tests-audited.tsv.gz')
    issues=audited[audited.CI_status.eq('degenerate_observed_between_unit_variance')]
    ctx.table(issues[['gene_id','gene_symbol','metric','effect','SE','df','p','q_BH_joint','CI_status']],
        'well-CI-estimability-audit.tsv')
    wells=ctx.read(ctx.analysis/'well_inference/qualified_well_pairs.tsv.gz')
    axes=ctx.read(ctx.analysis/'state_vectors/PAS_axes.primary.tsv.gz')
    axes=axes[axes.study_id.eq('GSE115943')].set_index(['primary_gene_id','atlas_pas_id'])
    data=[]
    for gene,q in wells.groupby('gene_id',sort=True):
        q=q.reset_index(drop=True)
        assert not q.sample_unit_a.duplicated().any() and q.replicate_status_a.eq('biological_independent').all()
        ids=json.loads(q.PAS_ids.iloc[0]);assert q.PAS_ids.nunique()==1
        mat=axes.loc[[(gene,p) for p in ids],['WPP','RWUI','PDUI','PPUI']].to_numpy()
        delta=np.stack(q.target_usage.map(json.loads))-np.stack(q.reference_usage.map(json.loads))
        assert np.allclose(delta.sum(1),0,atol=1e-12)
        data.append((gene,q,delta,mat))
    assert {g for g,*_ in data}==set(baseline.gene_id)
    rng=np.random.default_rng(20261005);summary=[];members=[];results=[]
    for draw in range(draws if data else 0):
        dp={x:i for i,x in enumerate(rng.permutation(wells.time_a.unique()))}
        up={x:i for i,x in enumerate(rng.permutation(wells.sample_unit_a.unique()))};rows=[]
        for gene,q,delta,mat in data:
            ix,weights,dates=six_well_indices(q,dp,up);effect,p=exact_vector_test(delta[ix],weights,mat)
            rows.extend(dict(draw=draw,gene_id=gene,metric=m,effect=e,p=pv) for m,e,pv in zip(METRICS,effect,p))
            members.append(dict(draw=draw,gene_id=gene,sample_units=';'.join(q.sample_unit_a.iloc[ix]),dates=';'.join(dates),weights=json.dumps(weights.tolist())))
        f=pd.DataFrame(rows);f['q_joint']=bh(f.p);results.append(f)
        for metric,q in f.groupby('metric'):
            old=baseline[baseline.metric.eq(metric)]
            summary.append(dict(draw=draw,metric=metric,genes=len(q),units=6,dates=4,
                full_significant=int(old.joint_q05.sum()),downsample_significant=int(q.q_joint.lt(.05).sum()),
                min_p=q.p.min(),min_joint_q=q.q_joint.min(),median_abs_effect=q.effect.abs().median(),status='estimated'))
    columns=['draw','metric','genes','units','dates','full_significant','downsample_significant','min_p','min_joint_q','median_abs_effect','status']
    if not data:
        summary=[dict(draw=-1,metric=m,genes=0,units=6,dates=4,full_significant=0,downsample_significant=0,
            min_p=np.nan,min_joint_q=np.nan,median_abs_effect=np.nan,status='no_eligible_current_well_genes') for m in METRICS]
    ctx.table(summary,'well_six_unit_downsample_summary.tsv',columns=columns)
    ctx.table(members,'well_six_unit_downsample_members.tsv.gz',columns=['draw','gene_id','sample_units','dates','weights'])
    ctx.table(pd.concat(results) if results else pd.DataFrame(columns=['draw','gene_id','metric','effect','p','q_joint']),
        'well_six_unit_downsample_tests.tsv.gz')
