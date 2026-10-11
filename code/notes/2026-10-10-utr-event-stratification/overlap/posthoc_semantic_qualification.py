"""Post-hoc ALE endpoint UTR-annotation qualification; no inferential tests."""
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
spec = importlib.util.spec_from_file_location('event_overlap_runtime', HERE / 'analyze.py')
P = importlib.util.module_from_spec(spec)
spec.loader.exec_module(P)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    qualification = HERE.parent / 'classification/both_ALE_terminal_UTR_coverage_pairs.tsv'
    assert sha(qualification) == 'b76ececf9964ada1d25448c2bdaf89667bfa6de46995bdd9cf50847cddce95e2'
    source_members = P.SCAN / 'measured_ortholog_members.tsv.gz'
    source_features = P.SCAN / 'ortholog_assay_features.tsv.gz'
    classification = HERE.parent / 'classification/event_class_pairs.tsv'
    assert sha(classification) == '5ba39cd3a5439557c7c3edf604ec6c9ca03c88acb0249996f4748572c17842c3'
    paths=[qualification,source_members,source_features,classification,HERE/'analyze.py',
           P.SCAN/'analyze.py',P.OLD/'analyze.py']
    hashes={str(p.relative_to(ROOT)):sha(p) for p in paths}
    q=pd.read_csv(qualification,sep='\t')
    gate='both_species_all_routes_both_ends_UTR_covered'
    assert len(q)==52 and q[gate].notna().all()
    assert not q.duplicated(P.KEY).any()
    eligible=q[q[gate]][P.KEY]
    assert len(eligible)==15
    members=pd.read_csv(source_members,sep='\t',float_precision='round_trip')
    members=members[members.contrast.eq('P_M')]
    classes=pd.read_csv(classification,sep='\t')
    joined=members.merge(classes[P.KEY+['pair_class']],on=P.KEY,validate='many_to_one')
    selected=joined.merge(eligible,on=P.KEY,validate='many_to_one')
    assert selected.pair_class.eq('both_ALE').all()
    features=pd.read_csv(source_features,sep='\t',float_precision='round_trip')
    features=features[features.contrast.eq('P_M')]
    rows,mobility=[],[]
    checks=0
    for threshold in P.THRESHOLDS:
        pair=selected[selected.threshold.eq(threshold)].reset_index(drop=True)
        assert len(pair)==15
        for model in P.MODELS:
            fs={sp:P.select_features(features,pair,sp,threshold) for sp in ['human','mouse']}
            ls={sp:pair[sp+'_direction'].map(dict(zip(P.LABELS,range(3)))).to_numpy(np.int8)
                for sp in ['human','mouse']}
            cs={sp:P.A.partitions(fs[sp],model) for sp in ['human','mouse']}
            probabilities={sp:P.A.class_probabilities(ls[sp],cs[sp]) for sp in ['human','mouse']}
            expected=probabilities['human'].T@probabilities['mouse']
            # Independent expectation construction from explicit equality keys.
            columns=['measured_pattern','proximal_capacity','distal_capacity','PAS_number_bin']
            if model=='full_assay':
                columns+=['mean_GEX_bin','APA_coverage_bin','baseline_bin']
            pp={}
            for sp in ['human','mouse']:
                pp[sp]=np.zeros((15,3))
                keys=list(map(tuple,fs[sp][columns].to_numpy()))
                for i,key in enumerate(keys):
                    idx=np.array([j for j,other in enumerate(keys) if other==key])
                    for direction in range(3):
                        pp[sp][i,direction]=float((ls[sp][idx]==direction).sum()/len(idx))
                context=dict(contrast='P_M',pair_class='both_ALE_UTR_covered_posthoc',
                             threshold=threshold,model=model,species=sp)
                _,summary=P.S.feature_diagnostics(fs[sp],ls[sp],cs[sp],context)
                mobility.append(summary)
            independent=pp['human'].T@pp['mouse']
            np.testing.assert_allclose(expected,independent,atol=0,rtol=0)
            checks+=9
            for a,hd in enumerate(P.LABELS):
                for b,md in enumerate(P.LABELS):
                    observed=int(((ls['human']==a)&(ls['mouse']==b)).sum())
                    rows.append(dict(contrast='P_M',pair_class='both_ALE_UTR_covered_posthoc',
                        threshold=threshold,model=model,measured_orthologs=15,
                        human_direction=hd,mouse_direction=md,observed=observed,
                        analytic_expected=float(expected[a,b]),excess=float(observed-expected[a,b]),
                        fold=float(observed/expected[a,b]) if expected[a,b] else np.nan))
    result=pd.DataFrame(rows)
    assert not any(c.startswith('p_') or c.startswith('q_') for c in result)
    P.save(result,'strict_ale_utr_descriptive.tsv')
    P.save(selected,'strict_ale_utr_members.tsv.gz')
    P.save(pd.DataFrame(mobility),'strict_ale_utr_mobility.tsv')
    assert all(sha(ROOT/name)==value for name,value in hashes.items())
    outputs=['strict_ale_utr_descriptive.tsv','strict_ale_utr_members.tsv.gz','strict_ale_utr_mobility.tsv']
    validation=dict(status='PASS',qualification='POSTHOC_SEMANTIC_DESCRIPTION_ONLY',
        parent_authorized_addendum='Both endpoints covered by terminal UTR annotations in all routes of both species.',
        background=15,contexts=18,expectation_checks=checks,skipped=0,
        no_new_p_values_or_families=True,original_main_analysis_unmodified=True,
        qualification_does_not_require_same_transcript_or_CDS=True,
        input_sha256=hashes,output_sha256={name:sha(HERE/name) for name in outputs},
        script_sha256=sha(Path(__file__)))
    (HERE/'posthoc_semantic_qualification_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(result[result.human_direction.eq('proximal')&result.mouse_direction.eq('proximal')].to_string(index=False))


if __name__=='__main__':
    main()
