"""Check fixed denominators and ensure matching failures never become flips."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent

def read(name):return pd.read_csv(HERE/name,sep='\t')

checks={}
for caller in ['polyapipe','apadetect_p50']:
    species=read(caller+'/dual_species_summary.tsv')
    route=read(caller+'/dual_route_summary.tsv')
    for frame in [species,route]:
        assert np.array_equal(frame.candidates,frame.retained+frame.flipped+frame.lost_alignment+frame.unmeasurable)
        assert np.array_equal(frame.measurable,frame.retained+frame.flipped+frame.lost_alignment)
    f=species[species.cohort.eq('caller_own')&species.subset.eq('all_candidates')]
    assert f.set_index('species').candidates.to_dict()==({'mouse':395,'human':250} if caller=='polyapipe' else {'mouse':203,'human':174})
    b=species[species.cohort.eq('caller_own')&species.subset.eq('baseline_dual')]
    assert b.set_index('species').candidates.to_dict()==({'mouse':222,'human':128} if caller=='polyapipe' else {'mouse':112,'human':85})
    triplets=read(caller+'/dual_triplets.tsv.gz');good=triplets.status.eq('measurable')
    assert triplets.loc[~good,['matched_PS','matched_PM','matched_MS']].isna().all().all()
    assert np.allclose(triplets.loc[good,'matched_PS'],triplets.loc[good,'matched_PM']+triplets.loc[good,'matched_MS'],rtol=0,atol=1e-10)
    for state in 'PMS':
        assert triplets.loc[good,state+'_n_mc'].ge(3).all()
        assert triplets.loc[good,state+'_n_cells'].ge(30).all()
    means=read(caller+'/pm_library_means.tsv')
    common=means[means.version.eq('raw_common')].set_index(['species','study_id','route','library_key'])
    matched=means[means.version.eq('cycle_matched')].set_index(['species','study_id','route','library_key'])
    assert np.array_equal(common.n_genes,matched.n_genes)
    for file in ['dual_alignment_tests.tsv','pm_route_mean_tests.tsv']:
        tests=read(caller+'/'+file);assert tests[tests.n_libraries.lt(2)].p.isna().all()
        assert tests[tests.n_libraries.lt(2)][['ci_low','ci_high']].isna().all().all()
        assert tests[tests.p.notna()].p.between(0,1).all()
    checks[caller]=dict(fixed_candidate_denominators=True,failures_separate_from_flips=True,
        shared_vector_additivity=True,minimum_coverage_enforced=True,pm_common_gene_denominators_equal=True,
        single_library_inference_withheld=True)
audit=json.loads((HERE/'audit.json').read_text())
assert all(x['max_raw_dual_usage_error']<1e-10 and x['max_raw_pm_effect_error']<1e-10 for x in audit['checks'])
checks['raw_reconstruction_agrees_with_frozen_inputs']=True
coverage=json.loads((HERE/'score_coverage.json').read_text())
assert coverage['all_seven_studies_covered'] and not coverage['excluded_studies']
score_states=read('cycle_score_library_states.tsv');score_tests=read('cycle_score_tests.tsv')
assert score_states.study_id.nunique()==7
assert score_tests[score_tests.n_libraries.lt(2)].p.isna().all()
assert np.isfinite(score_states[score_states.n_metacells.gt(0)]['mean']).all()
checks['all_seven_study_scores_covered_without_cell_pseudoreplication']=True
(HERE/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
