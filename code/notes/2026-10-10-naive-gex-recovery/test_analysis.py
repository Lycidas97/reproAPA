"""Intent checks for genuine GEX, weighting, and absent naive-library P cells."""
import unittest

import numpy as np
import pandas as pd

from analyze import csr_dense_rows, expression_rows, contrast_rows


class RecoveryIntentTests(unittest.TestCase):
    def test_csr_row_selection_conserves_duplicate_column_counts(self):
        archive = dict(format=np.array(b'csr'), shape=np.array([2, 3]),
                       data=np.array([2, 3, 7]), indices=np.array([1, 1, 2]), indptr=np.array([0, 2, 3]))
        counts = csr_dense_rows(archive, np.array([1, 0]))
        np.testing.assert_array_equal(counts, [[0, 0, 7], [0, 5, 0]])

    def test_GEX_uses_all_gene_denominator_and_MC_equal_mean(self):
        meta = pd.DataFrame(dict(library_key=['N6'] * 3, sample_accession=['N6'] * 3,
                                 node=['N'] * 3, condition=['naive'] * 3,
                                 experimental_time=['endpoint'] * 3, n_cells=[10, 10, 10]))
        counts = np.array([[1, 1], [1, 19], [0, 10]])
        rows = expression_rows(meta, counts, counts.sum(1), pd.DataFrame(dict(gene_id=['A', 'B'], gene_symbol=['a', 'b'])))
        a = rows.set_index('gene_id').loc['A']
        self.assertTrue(a.qualified)
        self.assertAlmostEqual(a.mean_log2CPM1, np.log2(np.array([500000, 50000, 0]) + 1).mean())
        self.assertAlmostEqual(a.pooled_CPM, 2 / 32 * 1e6)
        self.assertNotAlmostEqual(a.mean_log2CPM1, np.log2(a.pooled_CPM + 1))
        # A requested factor with zero observed GEX is observed zero, not unmapped.
        self.assertTrue(rows.mapped.all())
        meta['n_cells'] = 9
        self.assertFalse(expression_rows(meta, counts, counts.sum(1), pd.DataFrame(dict(gene_id=['A', 'B']))).qualified.any())

    def test_absent_N6_P_cannot_be_imputed_or_called_complete(self):
        lib = pd.DataFrame(dict(gene_id=['A'] * 3, sample_accession=['N6', 'N8', 'N8'], state=['N', 'N', 'P'],
                                qualified=[True] * 3, mean_log2CPM1=[4, 7, 5], mean_CPM=[16, 128, 32], pooled_CPM=[16, 128, 32]))
        f = contrast_rows(lib, {'same_library': [('N6', 'N', 'N6', 'P'), ('N8', 'N', 'N8', 'P')]})
        self.assertEqual(f.n_pairs.iloc[0], 1)
        self.assertEqual(f.planned_pairs.iloc[0], 2)
        self.assertFalse(f.complete_contributors.iloc[0])
        self.assertEqual(f.target_libraries.iloc[0], 'N8')
        self.assertEqual(f.delta_log2CPM1.iloc[0], 2)


if __name__ == '__main__':
    unittest.main()
