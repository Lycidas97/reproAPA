"""Protect scientific identity across the new caller boundary."""
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from apadetect_prepare import aggregate_library, convert_features, library_metacells, receipt, validate_handoff, write_json


class IntakeTests(unittest.TestCase):
    def test_frozen_key_join_preserves_missing_or_different_GEX_library_ids(self):
        meta = pd.DataFrame(dict(metacell_id=['m1', 'm2'], library_key=['study::A', 'study::B']))
        self.assertEqual(library_metacells(meta, 'study::A').metacell_id.tolist(), ['m1'])
        meta['library_id'] = ['original_GEX_A', 'original_GEX_B']
        selected = library_metacells(meta, 'study::B')
        self.assertEqual(selected.metacell_id.tolist(), ['m2'])
        self.assertEqual(selected.library_id.tolist(), ['original_GEX_B'])

    def test_frozen_key_join_never_guesses_prefix_or_suffix(self):
        meta = pd.DataFrame(dict(metacell_id=['m1'], library_key=['study::A']))
        with self.assertRaisesRegex(ValueError, 'exact library key'):
            library_metacells(meta, 'A')

    def feature_fixture(self):
        f = pd.DataFrame(dict(name=['p', 'm', 'empty', 'boundary'], strand=['+', '-', '+', '+'],
            chrom=['chr1'] * 4, score=[100, 200, 300, 1000], cut_site=[100, 200, 300, np.nan],
            coordinate_supported=[True, True, True, False], quantifiable=[True, True, False, True],
            terminal_cut_boundary0=[101, 200, 301, np.nan], pas_id=['p', 'm', 'empty', 'boundary'],
            min_cut=[100, 200, 300, 1000], max_cut=[100, 200, 300, 1000],
            member_count=[1] * 4, library_count=[1] * 4,
            reject_original_point=[False, True, np.nan, np.nan]))
        tiers = pd.DataFrame(dict(pas_id=f.pas_id, keep_for_apa_primary=['true'] * 4,
                                  keep_for_sensitivity=['true'] * 4, keep_for_pas_quant=['true'] * 4))
        ref = pd.DataFrame(dict(target_chrom=['1', '1'], target_strand=['+', '-'],
                                target_cut_site=[151, 200], source_family=['polyasite_broad'] * 2))
        return f, ref, tiers

    def test_base_to_boundary_is_strand_specific_and_support_is_rebuilt(self):
        f, ref, tiers = self.feature_fixture()
        f['matched_source_count_50nt'] = 99
        out = convert_features(f, ref, tiers)
        self.assertEqual(out.cut_site.tolist(), [101, 200, 301, 1001])
        self.assertEqual(out.matched_source_count_50nt.tolist(), [1, 1, 0, 0])
        self.assertEqual(out.matched_source_count_25nt.tolist(), [0, 1, 0, 0])

    def test_empty_and_unsupported_sites_stay_on_axis_but_never_define_primary(self):
        out = convert_features(*self.feature_fixture())
        self.assertEqual(len(out), 4)
        self.assertEqual(out.feature_index_1based.tolist(), [1, 2, 3, 4])
        self.assertEqual(out.keep_for_apa_primary.tolist(), [True, True, False, False])
        self.assertEqual(out.measurement_exclusion.iloc[2], 'unavailable_count_window')
        self.assertEqual(out.measurement_exclusion.iloc[3], 'unavailable_coordinate')

    def test_sensitivity_is_separate_and_unknown_is_not_retained(self):
        args = self.feature_fixture()
        a = convert_features(*args)
        b = convert_features(*args, view='reject_original_point')
        self.assertTrue(a.keep_for_apa_primary.iloc[1])
        self.assertFalse(b.keep_for_apa_primary.iloc[1])
        self.assertEqual(a.atlas_pas_id.tolist(), b.atlas_pas_id.tolist())

    def test_independent_coordinate_check_rejects_off_by_one(self):
        f, ref, tiers = self.feature_fixture()
        f.loc[0, 'terminal_cut_boundary0'] = 100
        with self.assertRaisesRegex(ValueError, 'coordinate conversion'):
            convert_features(f, ref, tiers)

    def matrix_fixture(self):
        x = sparse.csr_matrix([[3, 0, 99], [0, 2, 88]], dtype=np.int64)
        members = pd.DataFrame(dict(cell_id=['c1', 'c2'], barcode=['AA-1', 'BB-1'], metacell_id=['M', 'M']))
        mc = pd.DataFrame(dict(metacell_id=['M'], n_cells=[2]))
        return x, ['p', 'm'], ['p', 'm'], ['AA-1', 'BB-1', 'other-1'], members, mc

    def test_selects_exact_frozen_cells_and_preserves_per_pas_umi(self):
        out = aggregate_library(*self.matrix_fixture())
        np.testing.assert_array_equal(out.toarray(), [[3, 2]])

    def test_missing_barcode_does_not_use_suffix_guessing(self):
        args = list(self.matrix_fixture()); args[3] = ['AA-2', 'BB-1', 'other-1']
        with self.assertRaisesRegex(ValueError, 'suffix guessing'):
            aggregate_library(*args)

    def test_permuted_pas_axis_is_rejected_even_when_dimensions_match(self):
        args = list(self.matrix_fixture()); args[1] = ['m', 'p']
        with self.assertRaisesRegex(ValueError, 'row order'):
            aggregate_library(*args)

    def test_fractional_counts_are_not_treated_as_umi(self):
        args = list(self.matrix_fixture()); args[0] = args[0].astype(float)
        with self.assertRaisesRegex(ValueError, 'integer UMI'):
            aggregate_library(*args)

    def test_receipt_cannot_silently_reuse_modified_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp); (p / 'a.tsv').write_text('changed')
            write_json(p / 'complete.json', dict(status='COMPLETE', artifacts={'a.tsv': '0' * 64}))
            with self.assertRaisesRegex(ValueError, 'Changed input'):
                receipt(p, {})

    def test_running_release_cannot_publish_statistics_intake(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp)
            write_json(p / 'plan.json', {})
            write_json(p / 'active_run.json', {'status_file': 'status_v2.json'})
            write_json(p / 'status_v2.json', {'status': 'RUNNING', 'failures': []})
            with self.assertRaisesRegex(ValueError, 'still running'):
                validate_handoff(p, {})


if __name__ == '__main__':
    unittest.main()
