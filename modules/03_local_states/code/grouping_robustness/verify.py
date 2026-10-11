"""Independent checks of cell exclusions, direct PAS-weighted means and outcome accounting."""
import os
os.environ.update(OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2')
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parents[1];D=ROOT/'data'
BASE=Path('/root/nfsdata/apadetect/analysis_projects');PRIOR=BASE/'reproapa_pas_reliability_2026_09_13'
def read(p):return pd.read_csv(p,sep='\t',low_memory=False)
spec=json.loads((D/'apa_frozen_variants.json').read_text())
atlas=pd.read_parquet(BASE/'reproapa_integration_validation_20260913/execution_v1/shared/atlas.parquet')
ann=read(ROOT/'inputs/S5-local-state-annotations.tsv').set_index('original_state_id')
atlas['family']=atlas.state_key.map(ann.original_family_code.replace({'T':'P'}))
removed=read(D/'removed_cell_library_counts.tsv');deletion=[]
for (study,states),q in removed.groupby(['study','excluded_states']):
    expected=int(atlas.state_key.isin(states.split(';')).sum())
    assert int(q.removed_cells.sum())==expected
    assert (q.all_cells==q.removed_cells+q.remaining_cells).all()
    deletion.append(dict(study=study,states=states,removed_cells=expected))
values=read(D/'affected_library_positions.tsv.gz')
cases=[(spec['studies'][0],'exclude_WOT_L00','WOT_L00'),(spec['studies'][1],'exclude_CIP_C012','CIP_C012'),
       (spec['studies'][2],'exclude_human_C005','human_C005')]
direct=[]
for study,variant,target in cases:
    state=spec['targets'][target]
    # Sampling is fixed by removed-cell count, not APA effect or agreement.
    libs=removed[(removed.study==study)&(removed.excluded_states==state)].sort_values(['removed_cells','library_key'],ascending=[False,True]).head(3)
    f=read(PRIOR/'revision/data'/study/'features.tsv.gz');meta=read(PRIOR/'data'/study/'selected_metacells.tsv.gz')
    for lib in libs.library_key:
        folder=PRIOR/'revision/data'/study/'recount'/lib
        binding=read(folder/'cells.tsv.gz');raw=sparse.load_npz(folder/'original_cells.npz').tocsr()
        local=atlas[atlas.library_key.eq(lib)|atlas.sample_key.eq(lib)].set_index('barcode').loc[binding.barcode]
        keep=~local.state_key.eq(state).to_numpy()
        mc=pd.Index(meta.metacell_id);indices=mc.get_indexer(binding.metacell_id)
        selected=indices>=0;ix=indices[selected]
        counts=np.bincount(ix,weights=keep[selected],minlength=len(meta))
        correct=keep[selected]&(local.family.to_numpy()[selected]==meta.node.to_numpy()[ix])
        same=np.bincount(ix,weights=correct,minlength=len(meta))
        purity=np.divide(same,counts,out=np.zeros(len(meta)),where=counts>0)
        eligible=(counts>0)&(purity>=.8)&meta.library_key.eq(lib).to_numpy()
        expected=values[(values.study==study)&values.variant.eq(variant)&values.policy.eq('reapply_purity')&values.library_key.eq(lib)]
        # First 20 sorted genes with a reported value; this checks numerical estimation, not route eligibility.
        genes=sorted(expected.loc[expected.position_index.notna(),'gene_id'].unique())[:20]
        for gene in genes:
            cols=np.flatnonzero(f.primary_gene_id.eq(gene)&f.supported50)
            totals_per_cell=np.asarray(raw[:,cols].sum(1)).ravel()
            numerator_per_cell=np.asarray(raw[:,cols]@f.relative.to_numpy()[cols]).ravel()
            # np.add.at is separate from the production sparse aggregation and does not import its estimator.
            total=np.zeros(len(meta));num=np.zeros(len(meta))
            np.add.at(total,ix,totals_per_cell[selected]*keep[selected])
            np.add.at(num,ix,numerator_per_cell[selected]*keep[selected])
            for row in expected[expected.gene_id.eq(gene)&expected.position_index.notna()].itertuples():
                group=eligible&meta.node.eq(row.node).to_numpy()&meta.route.eq(row.route).to_numpy()
                assert group.sum()>=3 and counts[group].sum()>=30
                good=group&(total>=20);assert good.sum()>=3
                actual=float(np.mean(num[good]/total[good]))
                assert np.isclose(actual,row.position_index,rtol=0,atol=1e-12)
                direct.append(dict(study=study,library_key=lib,variant=variant,gene_id=gene,node=row.node,
                    direct_position=actual,reported_position=row.position_index,qualified_metacells=int(good.sum())))
pd.DataFrame(direct).to_csv(D/'independent_direct_PA_checks.tsv',sep='\t',index=False)
assert len(direct)>=100
# The fixed original effects must reproduce an independently existing seven-study result.
major=read(D/'gene_majority_outcomes.tsv.gz');old=read(PRIOR/'revision/data/local_majority_events.tsv.gz')
old=old[old.population.eq('current_figure2')][['species','gene_id','contrast','original_direction']].drop_duplicates()
b=major[major.variant.eq('baseline')&major.policy.eq('reapply_purity')]
check=b.merge(old,on=['species','gene_id','contrast'],validate='one_to_one')
assert len(check)==len(b) and np.array_equal(check.baseline_direction,check.original_direction)
shape=read(D/'gene_shape_assignments.tsv.gz');summary=read(D/'shape_retention.tsv')
for r in summary.itertuples():
    q=shape[(shape.species==r.species)&(shape.variant==r.variant)&(shape.policy==r.policy)]
    assert len(q)==r.fixed_genes and int(q.evaluable.sum())==r.evaluable
    assert int(q.fixed_center_cluster.eq(q.original_cluster).sum())==r.same_fixed_center
    assert int(q.refitted_cluster.eq(q.original_cluster).sum())==r.same_refitted_cluster
    if r.variant=='exclude_human_candidates' and r.species=='mouse':
        baseline=shape[(shape.species=='mouse')&shape.variant.eq('baseline')&shape.policy.eq(r.policy)].set_index('gene_id')
        residual=baseline.reindex(q.gene_id).feature_distance.to_numpy()
        # Both reconstructions can differ from the saved fit by floating-point roundoff.
        assert np.allclose(q.feature_distance,residual,atol=1e-15,rtol=0) and r.same_fixed_center==540
        assert r.same_refit_on_identical_eligible_set==540 and r.refit_species_present=='mouse'
effect=read(D/'local_effect_summary.tsv')
assert (effect.baseline_effects_ge001==effect.same_direction+effect.weakened+effect.opposite_direction+effect.new_unavailable).all()
assert not (major.fixed_studies<major.required_studies).any()
# A stability score computed against broad families is not an independent label transfer.
runs=read(D/'partition_runs.tsv');assert len(runs)==110 and runs[runs.role.eq('baseline_replay')].ARI.eq(1).all()
print('INDEPENDENT_VERIFICATION_PASS',len(direct),len(check),flush=True)
(D/'independent_verification.json').write_text(json.dumps(dict(status='PASS',cell_exclusion_counts=deletion,
    direct_PA_estimates=len(direct),maximum_direct_PA_difference=float(np.max(np.abs(pd.DataFrame(direct).direct_position-pd.DataFrame(direct).reported_position))),
    independently_existing_baseline_majority_rows=len(check),shape_denominators_checked=True,
    effect_status_denominators_checked=True,missing_human_refit_composition_control_checked=True),indent=2)+'\n')
