import numpy as np
import pandas as pd
from common import *


def run():
    d = pd.read_parquet(ROOT/'A06/matched_cell_diagnostics.parquet')
    x = obs()
    assert np.array_equal(d.cell_id, x.cell_id)
    out = RUN/'A06'
    progress('A06','RUNNING',stage='UMI_stratification_and_matched_random_exclusion')
    # Deterministic UMI deciles within each library, with cell ID as the tie-breaker.
    order = d.sort_values(['sample_key','rna_total_umi','cell_id']).index
    rank = d.loc[order].groupby('sample_key').cumcount()
    sizes = d.groupby('sample_key').size()
    strata = np.empty(len(d),int)
    strata[order] = np.minimum(9, (10*rank/d.loc[order,'sample_key'].map(sizes)).astype(int))
    d['umi_decile'] = strata
    d['umi_stratum'] = d.sample_key + '|' + pd.Series(strata,index=d.index).astype(str)
    groups = list(d.groupby('umi_stratum',sort=True).indices.values())
    masks = {'automatic_call': d.automatic_call.to_numpy(bool)}
    for fraction in [.05,.10]:
        mask = np.zeros(len(d),bool)
        for lib, idx in d.groupby('sample_key',sort=True).indices.items():
            ranked = d.iloc[idx].sort_values(['doublet_score','cell_id'],ascending=[False,True]).index
            mask[ranked[:int(np.ceil(len(idx)*fraction))]] = True
        masks[f'top_{int(fraction*100):02d}_percent'] = mask
    summary=[]; retained=[]
    keys=['species_key','true_study','state_key','display_group','sample_key']
    for arm, target in list(masks.items()):
        counts = np.array([target[ix].sum() for ix in groups])
        for seed in range(20):
            rng=np.random.default_rng(SEED+seed)
            control=np.zeros(len(d),bool)
            for ix,number in zip(groups,counts):
                control[rng.choice(ix,int(number),replace=False)]=True
            assert control.sum()==target.sum()
            assert all(control[ix].sum()==n for ix,n in zip(groups,counts))
            masks[f'{arm}_matched_random_{seed:02d}']=control
    for name,mask in masks.items():
        z=d[keys].copy();z['removed']=mask;z['remaining']=~mask
        v=z.groupby(keys).agg(n_cells=('removed','size'),removed=('removed','sum'),remaining=('remaining','sum')).reset_index()
        v['variant']=name;retained.append(v)
        summary.append(dict(variant=name,removed_cells=int(mask.sum()),remaining_cells=int((~mask).sum())))
    np.savez_compressed(out/'exclusion_masks.npz',cell_ids=d.cell_id.to_numpy(str),
                        **{k:np.packbits(v) for k,v in masks.items()},n_cells=np.array(len(d)))
    pd.concat(retained,ignore_index=True).to_csv(out/'state_library_retention.tsv.gz',sep='\t',index=False)
    pd.DataFrame(summary).to_csv(out/'exclusion_summary.tsv',sep='\t',index=False)
    # UMI-matched mean-score contrast, state versus other cells in the same stratum.
    totals=d.groupby('umi_stratum').agg(total=('doublet_score_percentile','sum'),n=('cell_id','size'))
    v=d.groupby(keys+['umi_stratum']).agg(total=('doublet_score_percentile','sum'),n=('cell_id','size')).reset_index()
    v['other_n']=v.umi_stratum.map(totals.n)-v.n
    v['other_total']=v.umi_stratum.map(totals.total)-v.total
    v['delta']=v.total/v.n-v.other_total/v.other_n.replace(0,np.nan)
    v['valid_cells']=np.where(v.other_n.gt(0),v.n,0)
    v['weighted_delta']=v.delta.fillna(0)*v.valid_cells
    a=v.groupby(keys).agg(valid_cells=('valid_cells','sum'),weighted_delta=('weighted_delta','sum')).reset_index()
    a['UMI_stratified_percentile_difference']=a.weighted_delta/a.valid_cells.replace(0,np.nan)
    a.to_csv(out/'UMI_adjusted_state_library_scores.tsv',sep='\t',index=False)
    receipt([ROOT/'A06/matched_cell_diagnostics.parquet',out/'exclusion_masks.npz'],out/'receipt.json')
    finish('A06','exclusions',variants=len(masks),matched_random_controls=60,cells=len(d),
           next_steps=['A04 expression support after exclusion','A05 APA reaggregation'],
           interpretation='Stress-test masks and score contrasts; not confirmed doublet identities.')


if __name__=='__main__': execute('A06',run)
