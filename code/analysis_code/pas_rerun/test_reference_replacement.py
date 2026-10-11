"""Reference changes must not move count axes or admit ambiguous liftOver hits."""
import unittest
import json
from pathlib import Path
import tempfile

import numpy as np
import pandas as pd

from replace_pas_reference import BED, SOURCES, accept_mappings, cuts, support_fields, source_records,carry_membership_audit
from polyapipe_prepare import digest


class IntegratedReference(unittest.TestCase):
    def test_reference_rebind_requires_identical_counts_before_inheriting_membership_evidence(self):
        with tempfile.TemporaryDirectory() as tmp:
            old=Path(tmp)/'old';new=Path(tmp)/'new'
            for root in [old,new]:(root/'studies/A').mkdir(parents=True)
            studies=[dict(study_id='A',cells=2,metacells=1,libraries=1)]
            (old/'intake.native.json').write_text(json.dumps(dict(studies=studies)))
            (new/'intake.native.json').write_text(json.dumps(dict(studies=studies,previous_native_intake_sha256=digest(old/'intake.native.json'))))
            audit=dict(status='PASS',cells=2,metacells=1,libraries=1,checks=dict(exact_mapping=True))
            (old/'frozen_native_membership_audit.json').write_text(json.dumps(audit))
            (old/'frozen_native_membership_audit.tsv').write_text('verified mapping\n')
            hashes={}
            for name in ['metacells.tsv.gz','metacell_by_pas_counts.npz']:
                for root in [old,new]:(root/'studies/A'/name).write_bytes(b'frozen')
                hashes[name]=digest(old/'studies/A'/name)
            (new/'annotation.complete.json').write_text(json.dumps(dict(studies=[dict(study_id='A',matrix_and_metadata_reused=hashes)])))
            counts=new/'studies/A/metacell_by_pas_counts.npz';counts.write_bytes(b'changed')
            with self.assertRaises(AssertionError):carry_membership_audit(old,new)
            self.assertFalse((new/'frozen_native_membership_audit.json').exists())
            counts.write_bytes(b'frozen');carry_membership_audit(old,new)
            self.assertEqual(digest(old/'frozen_native_membership_audit.json'),digest(new/'frozen_native_membership_audit.json'))

    def test_retained_base_has_strand_specific_interbase_boundary(self):
        np.testing.assert_array_equal(cuts([99, 99], [100, 100], ['+', '-']), [100, 99])
        with self.assertRaises(AssertionError):
            cuts([99], [101], ['+'])

    def test_ambiguous_mapping_is_rejected_before_contig_filter(self):
        source = pd.DataFrame([['chr1', 9, 10, 'a', 0, '+'], ['chr1', 20, 21, 'b', 0, '-'],
                               ['chr1', 30, 31, 'c', 0, '+']], columns=BED)
        mapped = pd.DataFrame([['chr1', 12, 13, 'a', 0, '+'], ['missing', 12, 13, 'a', 0, '+'],
                               ['chr1', 22, 23, 'b', 0, '+']], columns=BED)
        result = accept_mappings(source, mapped, {'1': 100}).set_index('pas_id')
        self.assertEqual(result.loc['a', 'mapping_status'], 'ambiguous')
        self.assertEqual(result.loc['b', 'target_cut_site'], 23)
        self.assertEqual(result.loc['b', 'target_strand'], '+')
        self.assertEqual(result.loc['b', 'target_chrom'], '1')
        self.assertEqual(result.loc['c', 'mapping_status'], 'unmapped')

    def test_support_uses_integrated_representative_and_inclusive_window(self):
        frame = pd.DataFrame(dict(chrom=['chr1']*4, strand=['+', '+', '-', '+'],
                                  cut_site=[150, 151, 100, 100], feature_index_1based=[1, 2, 3, 4],
                                  primary_gene_id=['g']*4, keep_for_apa_primary=[True]*4))
        reference = pd.DataFrame(dict(target_chrom=['chr1'], target_strand=['+'],
                                      target_cut_site=[100], source_family=['polyasite_broad']))
        result = support_fields(frame, reference)
        self.assertEqual(result.matched_source_count_50nt.tolist(), [1, 0, 0, 1])
        pd.testing.assert_frame_equal(result[frame.columns], frame)
        self.assertEqual(result.loc[2, 'polyasite_broad_distance_nt'], -1)

    def test_target_coordinate_collision_does_not_duplicate_source_support(self):
        frame = pd.DataFrame(dict(chrom=['chr1'], strand=['-'], cut_site=[100]))
        ref = pd.DataFrame(dict(target_chrom=['chr1']*3, target_strand=['-']*3,
                                target_cut_site=[100]*3, source_family=['polyasite_broad']*3))
        result = support_fields(frame, ref)
        self.assertEqual(result.matched_source_count_50nt.iloc[0], 1)

    def test_composite_metadata_preserves_each_database_source(self):
        ref = pd.DataFrame(dict(target_chrom=['chr1'], target_strand=['+'], target_cut_site=[100],
            supporting_sources=['GENCODE_polyA_site,PolyASite,PolyA_DB4.1Main']))
        records = source_records(ref)
        self.assertEqual(set(records.source_family), set(SOURCES.values()))
        query = pd.DataFrame(dict(chrom=['chr1'], strand=['+'], cut_site=[150]))
        self.assertEqual(support_fields(query, records).matched_source_count_50nt.iloc[0], 3)


if __name__ == '__main__':
    unittest.main()
