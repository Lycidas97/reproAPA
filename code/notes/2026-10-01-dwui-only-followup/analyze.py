"""Inspect DWUI-only calls without changing frozen effects or FDR families."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
DATA=ROOT/'notes/2026-09-30-multimetric-statistics'
FOLLOWUP=ROOT/'notes/2026-10-01-mpro-biology'
INPUTS={}

def read(path):
    INPUTS[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
    return pd.read_csv(path,sep='\t')

def save(frame,name):
    frame.to_csv(HERE/name,sep='\t',index=False,float_format='%.17g')

def main():
    tests=read(DATA/'metric_inference/gene_metric_tests.tsv.gz')
    inventory=read(FOLLOWUP/'matched_well_candidate_inventory.tsv')
    axis=read(DATA/'metric_inference/frozen_supported50_PAS.tsv.gz')
    wells=read(DATA/'metric_inference/frozen_well_vectors.tsv.gz')
    routes=read(DATA/'route_effects.tsv.gz')
    recurrence=read(DATA/'nulls/recurrence_members.tsv.gz')
    eligible=read(FOLLOWUP/'native/mouse/eligibility.tsv')
    chosen=inventory[inventory.DWUI_q.lt(.05)&inventory.MPRO_q.ge(.05)].copy()
    assert len(chosen)==15
    summaries,site_rows=[],[]
    for row in chosen.itertuples():
        # Serialized vectors retain source-file PAS order, not transcript rank.
        a=axis[axis.primary_gene_id.eq(row.gene_id)].copy()
        w=wells[wells.gene_id.eq(row.gene_id)]
        t=tests[tests.gene_id.eq(row.gene_id)].set_index('metric')
        d=np.asarray(json.loads(t.loc['DWUI','weighted_delta_usage']))
        matrix=np.stack(w.delta_usage.map(json.loads))
        weights=1/w.time.nunique()/w.date_replicates.to_numpy()
        assert np.isclose(weights.sum(),1)
        assert np.allclose(weights@matrix,d,atol=1e-12)
        assert np.isclose(d@a.WPP.to_numpy(),row.DWUI,atol=1e-12)
        assert np.isclose(np.ptp(d),row.MPRO,atol=1e-12)
        assert np.allclose(matrix@a.WPP.to_numpy(),w.delta_DWUI,atol=1e-12)
        dates=w.groupby('time').delta_DWUI.mean()
        s=dict(gene_id=row.gene_id,gene_symbol=row.gene_symbol,
            DWUI=row.DWUI,DWUI_p=t.loc['DWUI','p'],DWUI_q=row.DWUI_q,
            MPRO=row.MPRO,MPRO_p=t.loc['MPRO','p'],MPRO_q=row.MPRO_q,
            n_PAS=row.n_PAS,n_wells=len(w),n_dates=len(dates),
            wells_in_mean_direction=int((np.sign(w.delta_DWUI)==np.sign(row.DWUI)).sum()),
            dates_in_mean_direction=int((np.sign(dates)==np.sign(row.DWUI)).sum()),
            DWUI_effect_ge_0p02=abs(row.DWUI)>=.02)
        for route in ['CIP','WOT']:
            e=routes[routes.gene_id.eq(row.gene_id)&routes.route.eq(route)&routes.contrast.eq('P_M')]
            if len(e):
                assert len(e)==1
                s[f'{route}_native_PM_DWUI']=e.iloc[0].delta_WPP
                s[f'{route}_native_PM_MPRO']=e.iloc[0].MPRO
            f=eligible[eligible.gene_id.eq(row.gene_id)&eligible.route.eq(route)&eligible.contrast.eq('MP')]
            s[f'{route}_F5_status']=f.iloc[0].status if len(f) else 'not_admitted_for_MP'
        summaries.append(s)
        a['gene_symbol']=row.gene_symbol
        a['weighted_delta_usage']=d
        a['DWUI_contribution']=d*a.WPP.to_numpy()
        a['negative_wells']=(matrix<0).sum(0)
        a['positive_wells']=(matrix>0).sum(0)
        a['delta_usage_sd']=matrix.std(0,ddof=1)
        site_rows.append(a)
    summary=pd.DataFrame(summaries).sort_values('DWUI',key=abs,ascending=False)
    save(summary,'DWUI_only_significant_gene_audit.tsv')
    save(pd.concat(site_rows,ignore_index=True),'DWUI_only_PAS_geometry.tsv')
    save(wells[wells.gene_id.isin(chosen.gene_id)],'DWUI_only_frozen_well_vectors.tsv')
    additions=[]
    for sp in ['mouse','human']:
        for co in ['P_M','P_S']:
            a=recurrence[recurrence.species.eq(sp)&recurrence.contrast.eq(co)&recurrence.metric.eq('WPP')&recurrence.threshold.eq(.02)]
            b=recurrence[recurrence.species.eq(sp)&recurrence.contrast.eq(co)&recurrence.metric.eq('MPRO')&recurrence.threshold.eq(.1)]
            additions.append(a[~a.gene_id.isin(b.gene_id)].assign(DWUI_threshold=.02,MPRO_threshold=.1))
    additions=pd.concat(additions,ignore_index=True)
    save(additions,'DWUI_only_recurrence_members.tsv')
    save(routes.merge(additions[['species','contrast','gene_id']],on=['species','contrast','gene_id'],validate='many_to_one'),'DWUI_only_recurrence_route_effects.tsv')
    assert summary.DWUI_effect_ge_0p02.sum()==3
    tom=summary[summary.gene_symbol.eq('Tomm20')].iloc[0]
    assert tom.wells_in_mean_direction==tom.n_wells==22 and tom.dates_in_mean_direction==tom.n_dates==11
    for path,digest in INPUTS.items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    report=dict(status='PASS',significant_only_genes=len(summary),absolute_DWUI_ge_0p02=3,
        MPRO_q_between_0p05_and_0p10=int(summary.MPRO_q.lt(.1).sum()),
        recurrence_only_counts=additions.groupby(['species','contrast']).size().to_dict(),
        original_tests_changed=False,new_tests_performed=False,input_sha256=INPUTS)
    report['recurrence_only_counts']={f'{sp}/{co}':n for (sp,co),n in report['recurrence_only_counts'].items()}
    (HERE/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(summary[['gene_symbol','DWUI','DWUI_q','MPRO_q','wells_in_mean_direction','n_wells']].to_string(index=False))

if __name__=='__main__':main()
