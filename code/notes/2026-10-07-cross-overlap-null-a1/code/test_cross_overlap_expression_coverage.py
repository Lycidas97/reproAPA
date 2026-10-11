"""The sensitivity must retain the original estimands, labels and backgrounds."""
import ast
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from cross_overlap_expression_coverage import BASES, paired_weights, summarize_contributors, terciles
from downstream_cross import feature_frame


class CrossOverlapSensitivityTests(unittest.TestCase):
    def test_paired_measurement_count_is_not_the_available_node_reference_count(self):
        orth = pd.DataFrame(dict(human_gene_id=['h0', 'h1'], mouse_gene_id=['m0', 'm1']))
        nodes = pd.DataFrame([
            dict(layer='supported50', measurement_pass=True, gene_id=gene, code=code,
                 node=node, position_index=.5, supported_pas_in_route=2)
            for gene, code in [('h0', 'H1'), ('h1', 'H1'), ('m0', 'W'), ('m1', 'W')]
            for node in ['P', 'S']
        ])
        paired = pd.DataFrame([
            dict(layer='supported50', gene_id=gene, code=code, mean_delta=-.03)
            for gene, codes in [('h0', ['H1', 'H2']), ('h1', ['H1']),
                                ('m0', ['W', 'C']), ('m1', ['W'])]
            for code in codes
        ])
        frame = feature_frame(orth, nodes, paired, 'supported50')
        for side in ['h', 'm']:
            self.assertEqual(frame[side + '_PM_effect_n'].tolist(), [0, 0])
            self.assertEqual(frame[side + '_PM_actual_effect_n'].tolist(), [2, 1])
            self.assertEqual(frame[side + '_PM'].tolist(), [True, False])

    def test_extended_nulls_cannot_place_events_on_one_study_genes(self):
        path = Path(__file__).resolve().parent / 'numerics/common.py'
        if not path.exists():
            path = Path(__file__).resolve().parents[2] / 'notes/2026-10-02-null-model-controls/common.py'
        functions = [item for item in ast.parse(path.read_text()).body
                     if isinstance(item, ast.FunctionDef) and item.name in ['groups', 'binary_draw', 'mobility']]
        ns = dict(np=np, pd=pd)
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), ns)
        for model, fields in BASES.items():
            columns = fields + ['mean_GEX_q3', 'APA_coverage_q3']
            frame = pd.DataFrame({column: [0, 0, 0, 0] for column in columns})
            frame['actual_effect_n'] = [2, 2, 1, 0]
            labels = np.array([True, False, False, False])
            cells = ns['groups'](frame, columns)
            probabilities, _ = ns['mobility'](labels, cells)
            draws = ns['binary_draw'](labels, cells, 101, np.random.default_rng(7))
            with self.subTest(model=model):
                self.assertTrue(np.array_equal(probabilities, [.5, .5, 0., 0.]))
                self.assertFalse(draws[:, 2:].any())
                self.assertTrue((draws.sum(1) == 1).all())
                self.assertTrue(draws[:, 0].any() and draws[:, 1].any())

    def test_capture_multiplicity_cannot_change_sampling_unit_or_date_weight(self):
        frame = pd.DataFrame(dict(time=[0, 0, 0, 0, 1], sample_unit=['a', 'a', 'a', 'b', 'c']))
        weights = paired_weights(frame)
        self.assertTrue(np.allclose(weights, [1/12, 1/12, 1/12, 1/4, 1/2]))
        self.assertAlmostEqual(weights[:3].sum(), weights[3])
        self.assertAlmostEqual(weights[:4].sum(), weights[4])

    def test_observed_source_nodes_keep_independent_equal_library_averages(self):
        frame = pd.DataFrame([
            dict(library_key='a', node='P', position=.8, GEX_CPM=10, APA_coverage=20, n_mc=3),
            dict(library_key='b', node='P', position=.4, GEX_CPM=30, APA_coverage=60, n_mc=30),
            dict(library_key='c', node='S', position=.2, GEX_CPM=40, APA_coverage=80, n_mc=100),
            dict(library_key='a', node='M', position=.5, GEX_CPM=99, APA_coverage=99, n_mc=9),
        ])
        row, selected = summarize_contributors(frame, 'PS')
        self.assertAlmostEqual(row['effect'], .4)
        self.assertAlmostEqual(row['mean_GEX'], 30)
        self.assertAlmostEqual(row['APA_coverage'], 60)
        self.assertEqual(set(selected.node), {'P', 'S'})
        self.assertTrue(np.allclose(selected.weight, [.25, .25, .5]))

    def test_unpaired_state_cannot_enter_paired_PM_covariates(self):
        frame = pd.DataFrame([
            dict(library_key='a', node='P', position=.2, GEX_CPM=20, APA_coverage=40, n_mc=3),
            dict(library_key='a', node='M', position=.4, GEX_CPM=40, APA_coverage=80, n_mc=3),
            dict(library_key='b', node='M', position=.99, GEX_CPM=999, APA_coverage=999, n_mc=90),
        ])
        design = pd.DataFrame(dict(library_key=['a', 'b'], time=[0, 1], sample_unit=['a', 'b']))
        row, selected = summarize_contributors(frame, 'PM', design)
        self.assertAlmostEqual(row['effect'], -.2)
        self.assertAlmostEqual(row['mean_GEX'], 30)
        self.assertAlmostEqual(row['APA_coverage'], 60)
        self.assertEqual(set(selected.library_key), {'a'})

    def test_missing_expression_stays_missing_without_dropping_contributors(self):
        frame = pd.DataFrame([
            dict(library_key='a', node='P', position=.2, GEX_CPM=np.nan, APA_coverage=40, n_mc=3),
            dict(library_key='a', node='M', position=.4, GEX_CPM=0, APA_coverage=80, n_mc=3),
        ])
        row, selected = summarize_contributors(
            frame, 'PM', pd.DataFrame(dict(library_key=['a'], time=[0], sample_unit=['a'])))
        self.assertTrue(np.isnan(row['mean_GEX']))
        self.assertAlmostEqual(row['APA_coverage'], 60)
        self.assertEqual(len(selected), 2)
        codes, _ = terciles(pd.Series([np.nan, 0., 0., 1., 2., 3.]))
        self.assertEqual(codes.iloc[0], -1)
        self.assertNotEqual(codes.iloc[1], -1)
        self.assertEqual(codes.iloc[1], codes.iloc[2])


if __name__ == '__main__':
    unittest.main()
