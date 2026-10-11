"""Independently check cell coverage and temporal provenance of the display layer."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

P=Path(__file__).resolve().parents[1]
D=P/'data/state-flow-v1'


def main():
    cells=pd.read_csv(P/'data/results_v1/prepared/primary_cells.tsv.gz',sep='\t',low_memory=False)
    plan=pd.read_csv(P/'data/results_v1/prepared/temporal_edge_plan.tsv',sep='\t')
    original=pd.read_csv(P/'data/results_temporal_v1/layout/state_condition_flows.tsv.gz',sep='\t')
    mapping=pd.read_csv(D/'local_to_coarse_state.tsv',sep='\t')
    coarse=pd.read_csv(D/'all_coarse_time_flows.tsv',sep='\t')
    conditions=pd.read_csv(D/'condition_preserving_coarse_flows.tsv',sep='\t')
    selected=pd.read_csv(D/'figure_temporal_rows.tsv',sep='\t')
    assertions=[]
    def check(name,value):
        assertions.append(dict(check=name,passed=bool(value)))
        assert value,name
    check('all_primary_cells_once',cells.cell_id.is_unique and len(cells)==215114)
    check('all_frozen_states_once',mapping.state_id.is_unique and set(mapping.state_id)==set(cells.original_state_id))
    expected=cells.original_state_id.value_counts().sort_index()
    check('each_state_cell_count_exact',mapping.set_index('state_id').n_cells.sort_index().equals(expected))
    check('no_primary_cells_lost',mapping.n_cells.sum()==215114)
    check('all_34_time_pairs_retained',set(original.pair)==set(coarse.pair)==set(range(34)))
    for pair,row in plan.iterrows():
        g=original[original.pair.eq(pair)]
        check(f'pair_{pair}_route_and_order',g.route_id.eq(row.route_id).all() and g.time_from.eq(row.time_from).all() and g.time_to.eq(row.time_to).all())
        for side,column in [('source','source_libraries'),('target','target_libraries')]:
            cc=cells[cells.library_key.isin(json.loads(row[column]))].copy()
            observed=cc.groupby(cc.original_state_id+'|'+cc.condition,dropna=False).size().sort_index()
            recorded=g.drop_duplicates(side+'_state_condition').set_index(side+'_state_condition')[side+'_cells'].sort_index()
            check(f'pair_{pair}_{side}_registered_cell_counts',observed.equals(recorded))
        check(f'pair_{pair}_mass_preserved',np.isclose(g.mass_fraction.sum(),coarse[coarse.pair.eq(pair)].mass_fraction.sum(),rtol=0,atol=1e-10))
    check('no_GSE150311_time_edges',not original.source_state_condition.str.startswith('GSE150311').any())
    check('no_GSE221739_fabricated_D0',not coarse[coarse.route_id.eq('GSE221739')].time_from.eq('D0').any())
    check('conditional_source_normalization',np.allclose(coarse.groupby(['pair','coarse_source']).source_fraction.sum(),1))
    check('conditions_preserve_pair_mass',np.allclose(conditions.groupby('pair').mass_fraction.sum(),1))
    check('figure_rows_exist_in_full_table',set(selected.edge_id)<=set(coarse.edge_id))
    check('all_display_edges_cross_coarse',selected.coarse_source.ne(selected.coarse_target).all())
    check('no_naive_time_incoming',not selected.coarse_target.eq('N').any())
    check('chemical_endpoint_edge_is_0618',selected[selected.route_id.str.startswith('gse178325') & selected.coarse_target.eq('P')].route_id.eq('gse178325_0618').all())
    check('minor_keratin_probability_exact',np.isclose(selected.set_index('edge_id').loc['H0045','source_fraction'],0.03061604,atol=1e-8))
    check('source_annotations_preserved',mapping.original_annotation_preserved.all())
    check('cross_study_equivalence_not_invented',not mapping.cross_study_equivalence_validated.any())
    for path,digest in json.loads((D/'data_review.json').read_text())['source_hashes'].items():
        check('unchanged_'+path,hashlib.sha256((P/path).read_bytes()).hexdigest()==digest)
    expression_checks=pd.read_csv(D/'input_checks.tsv',sep='\t')
    check('all_expression_cells_validated',expression_checks.n_cells.sum()==215114 and expression_checks.exact_cell_set.all() and expression_checks.exact_frozen_cluster_labels.all())
    report=dict(status='PASS',checks=len(assertions),all_passed=True,items=assertions,
                interpretation='Coverage and provenance validation, not biological lineage truth')
    (P/'analysis/human-state-flow-validation.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print('PASS',len(assertions))


if __name__=='__main__':main()
