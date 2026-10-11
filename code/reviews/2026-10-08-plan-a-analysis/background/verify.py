#!/usr/bin/env python3
"""Independent counts, source hashes and null moments for A1/A7."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASES = {
    'polyapipe': ROOT / 'notes/2026-10-06-supplement-polyapipe/contract-r2/bundle/tables',
    'apadetect-p50': ROOT / 'notes/2026-10-06-dual-main-plot-data/apadetect-p50/supplement/tables',
}


def main():
    output = pd.read_csv(HERE / 'gene_background.tsv.gz', sep='\t')
    strata = pd.read_csv(HERE / 'strata.tsv.gz', sep='\t')
    tests = pd.read_csv(HERE / 'tests.tsv', sep='\t')
    report = []
    for pas, base in BASES.items():
        frozen = pd.read_csv(base / 'matched_pms/matched_state_vectors.tsv.gz', sep='\t')
        frozen = frozen[frozen.setting.eq('supported50_p08')]
        rec = pd.read_csv(base / 'recurrence/recurrence_members.tsv.gz', sep='\t')
        rec = rec[rec.series_matching & rec.metric.eq('WPP') & rec.threshold.eq(.02) & rec.contrast.eq('P_M')]
        dual = pd.read_csv(base / 'candidate_screen/dual_reference_candidates.tsv', sep='\t')
        ladder = pd.read_csv(base / 'candidate_screen/dual_reference_ladder.tsv', sep='\t')
        for sp, routes in [('human', ['GSE118258']), ('mouse', ['CIP', 'WOT'])]:
            z = frozen[frozen.species.eq(sp) & frozen.route.isin(routes)]
            ps = z.pivot(index='gene_id', columns='route', values='DWUI_PS')[routes].dropna()
            pm = z.pivot(index='gene_id', columns='route', values='DWUI_PM')[routes].dropna()
            measurable = ps.index.intersection(pm.index)
            ps = ps.loc[measurable]; pm = pm.loc[measurable]
            any_positive = (ps.ge(.02) & pm.ge(.02)).all(axis=1)
            any_negative = (ps.le(-.02) & pm.le(-.02)).all(axis=1)
            r = rec[rec.species.eq(sp)].set_index('gene_id')
            qualifying = measurable.intersection(r.index)
            original_positive = r.loc[qualifying].direction.eq('positive')
            fixed_support = (any_positive.loc[qualifying] & original_positive) | (any_negative.loc[qualifying] & ~original_positive)
            o = output[output.pas_source.eq(pas) & output.species.eq(sp)]
            assert set(o.gene_id) == set(measurable)
            assert int(o.is_pm_recurrent.sum()) == len(qualifying)
            assert set(o.loc[o.is_dual_reference, 'gene_id']) == set(fixed_support[fixed_support].index)
            assert set(o.loc[o.is_dual_reference, 'gene_id']) == set(dual.loc[dual.species.eq(sp), 'gene_id'])
            assert int(o.aligned_any_direction.sum()) == int((any_positive | any_negative).sum())
            known = ladder[ladder.species.eq(sp)].iloc[0]
            assert (len(qualifying), int(fixed_support.sum())) == (known.matched_measurable, known.dual_reference)
            for model in ['reference_capacity', 'unmatched_signed_labels', 'conditional_PS_reference_capacity']:
                s = strata[strata.pas_source.eq(pas) & strata.species.eq(sp) & strata.model.eq(model)]
                direct_expected = ((s.positive_candidates * s.positive_aligned + s.negative_candidates * s.negative_aligned) / s.genes).sum()
                row = tests[tests.test_id.eq(f'{pas}|{sp}|{model}')].iloc[0]
                assert abs(direct_expected - row.null_mean) < 1e-8
                with np.load(HERE / f'null-{pas}-{sp}-{model}.npz') as archive:
                    prob = archive['exact_probability']
                    assert abs(prob.sum() - 1) < 1e-10
                    assert abs(np.arange(len(prob)) @ prob - direct_expected) < 1e-8
                    assert abs(prob[int(fixed_support.sum()):].sum() - row.p_value) < 1e-10
            report.append(dict(pas_source=pas, species=sp, measurable=len(measurable), pm_candidates=len(qualifying),
                               fixed_direction_dual=int(fixed_support.sum()), unsigned_aligned=int((any_positive | any_negative).sum()),
                               opposite_candidates=int(o.opposite_recurrent_alignment.sum()), status='PASS'))
    provenance = json.loads((HERE / 'provenance.json').read_text())
    for path, expected in provenance['input_files'].items():
        with Path(path).open('rb') as stream:
            assert hashlib.file_digest(stream, 'sha256').hexdigest() == expected['sha256']
    eroute = pd.read_csv(HERE / 'E_route_gene_effects.tsv.gz', sep='\t')
    estrict = pd.read_csv(HERE / 'E_same_library_gene_effects.tsv.gz', sep='\t')
    assert np.max(np.abs(eroute.DWUI_PE - (eroute.P_DWUI - eroute.E_DWUI))) < 1e-10
    assert not eroute.duplicated(['pas_source', 'route', 'gene_id']).any()
    assert not estrict.duplicated(['pas_source', 'route', 'gene_id']).any()
    assert set(estrict.route) == {'CIP'}
    # Independent frozen physical-library support audit: WOT cannot become a same-library E control.
    for pas, base in BASES.items():
        triplets = pd.read_csv(base / 'matched_pms/weighted_triplets.tsv.gz', sep='\t',
                              usecols=['species', 'setting', 'route', 'gene_id', 'library_key_a'])
        pe = pd.read_csv(base / 'intermediate/paired_library_effects.tsv.gz', sep='\t')
        pe = pe[pe.species.eq('mouse') & pe.setting.eq('supported50_p08') & pe.contrast.eq('P_E')]
        triplets = triplets[triplets.species.eq('mouse') & triplets.setting.eq('supported50_p08')]
        support = triplets.merge(pe.rename(columns={'library_key': 'library_key_a'}),
                                 on=['species', 'setting', 'route', 'gene_id', 'library_key_a'])
        assert set(support.route) == {'CIP'}
        counts = pd.read_csv(HERE / 'E_summary.tsv', sep='\t')
        for rt in ['CIP', 'WOT']:
            row = counts[counts.pas_source.eq(pas) & counts.design.eq('same_library_P_S_M_E') & counts.route.eq(rt)].iloc[0]
            assert row.measurable_background == len(estrict[estrict.pas_source.eq(pas) & estrict.route.eq(rt)])
    result = dict(status='PASS', checks=report, source_hashes_checked=len(provenance['input_files']),
                  independent_null_expectation='signed candidate margins vs alignment margins',
                  E_same_library_WOT='verified unavailable', skipped_checks=[])
    (HERE / 'verification.json').write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
