import numpy as np
import pandas as pd
from common import *


def run():
    out=RUN/'A05'; records=[]; summaries=[]; receipts=[]
    for species,expected in [('human',1194),('mouse',1569)]:
        member=ROOT/'prior_evidence/APA'/species/'aligned_gene_membership.tsv'
        m=pd.read_csv(member,sep='\t');assert len(m)==expected and m.gene_id.is_unique
        root=MOUSE/'results_apa_shared_local_v1'/f'{species}_supported50'
        votes_path=root/'all_study_votes.tsv.gz'
        local_path=root/'all_gene_contrast_evidence.tsv.gz'
        effect_path=root/'all_route_effects.tsv.gz'
        receipts += [member,votes_path,local_path,effect_path]
        v=pd.read_csv(votes_path,sep='\t');v=v[v.gene_id.isin(m.gene_id)]
        baseline=pd.read_csv(local_path,sep='\t').set_index(['gene_id','contrast'])
        effects=pd.read_csv(effect_path,sep='\t');effects=effects[effects.gene_id.isin(m.gene_id)]
        assert not v.duplicated(['gene_id','contrast','study_id']).any()
        progress('A05','RUNNING',stage='fixed_candidate_leave_one_study_out',species=species)
        studies=sorted(v.study_id.unique())
        for contrast,z in v.groupby('contrast',sort=True):
            study_list=sorted(z.study_id.unique())
            matrix=z.pivot(index='gene_id',columns='study_id',values='study_vote').reindex(index=m.gene_id,columns=study_list)
            assert matrix.notna().all().all()
            raw=effects[effects.contrast.eq(contrast)].groupby(['gene_id','study_id']).delta.mean().unstack('study_id').reindex(index=m.gene_id,columns=study_list)
            base=baseline.reindex(pd.MultiIndex.from_arrays([m.gene_id,np.repeat(contrast,len(m))]))
            def direction(a):
                need=a.shape[1]//2+1
                return np.where((a==1).sum(1)>=need,1,np.where((a==-1).sum(1)>=need,-1,0))
            assert np.array_equal(direction(matrix.to_numpy()),base.direction.to_numpy())
            for excluded in ['NONE',*studies]:
                columns=[s for s in study_list if s!=excluded]
                remaining=matrix[columns].to_numpy()
                sign=direction(remaining)
                assessable=len(columns)>=2
                rows=pd.DataFrame(dict(species=species,gene_id=m.gene_id.to_numpy(),contrast=contrast,
                    excluded_study=excluded,n_required_studies=len(columns),cross_study_assessable=assessable,
                    baseline_direction=base.direction.to_numpy(),recomputed_direction=sign,
                    baseline_mean_effect=raw.mean(axis=1).to_numpy(),
                    leave_study_mean_effect=raw[columns].mean(axis=1).to_numpy(),
                    n_evaluable_study_effects=raw[columns].notna().sum(axis=1).to_numpy()))
                rows['status']=np.select([~rows.cross_study_assessable,rows.n_evaluable_study_effects.eq(0),
                    rows.recomputed_direction.eq(0),rows.recomputed_direction.eq(rows.baseline_direction)],
                    ['fewer_than_two_studies','missing_effect','no_majority','same_direction'],default='changed_direction')
                records.append(rows)
                summaries.append(dict(species=species,contrast=contrast,excluded_study=excluded,
                    fixed_candidate_denominator=len(m),cross_study_assessable=assessable,
                    statuses=rows.status.value_counts().to_dict()))
    table=pd.concat(records,ignore_index=True)
    table.to_csv(out/'fixed_candidate_leave_study_out.tsv.gz',sep='\t',index=False)
    write_json(out/'leave_study_out_summary.json',summaries)
    receipt(receipts,out/'input_receipt.json')
    finish('A05','fixed_candidate_LOSO',rows=len(table),
        next_steps=['A04 alternative annotation membership','A07 alternative partitions','A06 cell-level APA reaggregation'],
        boundary='Uses frozen route effects; does not claim cell-level APA reaggregation or biological replicate validation.')


if __name__=='__main__': execute('A05',run)
