"""The original three nulls must never assign events to unmeasurable genes."""
import ast
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np
import pandas as pd

from downstream_cross import feature_frame, overlap


class OriginalNullEligibilityTests(unittest.TestCase):
    def test_both_species_measurement_counts_constrain_all_original_nulls(self):
        orth = pd.DataFrame(dict(human_gene_id=['h0', 'h1', 'h2'], mouse_gene_id=['m0', 'm1', 'm2']))
        nodes = pd.DataFrame([
            dict(layer=layer, measurement_pass=True, gene_id=gene, code=code,
                 node=node, position_index=.5, supported_pas_in_route=2)
            for layer in ['supported50', 'expanded_primary']
            for gene, code in [('h0', 'H1'), ('h1', 'H1'), ('h2', 'H1'),
                               ('m0', 'W'), ('m1', 'W'), ('m2', 'W')]
            for node in ['P', 'S']
        ])
        paired = pd.DataFrame([
            dict(layer=layer, gene_id=gene, code=code, mean_delta=effect)
            for layer in ['supported50', 'expanded_primary']
            for gene, codes, effect in [('h0', ['H1', 'H2'], -.03), ('h1', ['H1'], -.03),
                                       ('h2', ['H1', 'H2'], 0.), ('m0', ['W', 'C'], -.03),
                                       ('m1', ['W', 'C'], 0.), ('m2', ['W'], -.03)]
            for code in codes
        ])
        frozen = feature_frame(orth, nodes, paired, 'supported50')
        frozen = frozen.drop(columns=[column for column in frozen if column.endswith('_actual_effect_n')])
        import downstream_cross
        path = Path(downstream_cross.__file__).parent / 'numerics/figure4_matched_overlap.py'
        if not path.exists():
            path = Path(downstream_cross.__file__).resolve().parents[2] / 'notes/2026-09-29-core-sensitivity/figure4_matched_overlap.py'

        class Context:
            analysis = Path('/frozen')

            def resource(self, key):
                return Path('/orthologs')

            def read(self, source):
                return {'orthologs': orth, 'ortholog_matched_features.tsv.gz': frozen,
                        'all_mapped_nodes.tsv.gz': nodes, 'paired_pm_study_effects.tsv.gz': paired}[source.name].copy()

            def functions(self, relative, names, namespace):
                functions = [item for item in ast.parse(path.read_text()).body
                             if isinstance(item, ast.FunctionDef) and item.name in names]
                ns = dict(np=np, pd=pd, **namespace)
                ns['B'] = 101
                exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), ns)
                return SimpleNamespace(**{name: ns[name] for name in names})

            def table(self, frame, name):
                if name == 'overlap_null_results.tsv':
                    self.results = frame

            def finish(self, checks, **details):
                assert all(checks.values())

        with tempfile.TemporaryDirectory() as temporary:
            context = Context()
            context.out = Path(temporary)
            overlap(context, cached=True)
            rows = context.results[context.results.contrast.eq('PM')].set_index('model')
            self.assertEqual(rows.loc['joint_measurement_PAS', 'analytic_null_mean'], 1.)
            self.assertEqual(rows.loc['independent_q5_reference_capacity', 'analytic_null_mean'], .25)
            self.assertEqual(rows.loc['independent_range10_reference_capacity', 'analytic_null_mean'], .25)
            # These genes share all legacy strata; omitting either real count changes the expectation.
            pd.testing.assert_series_equal(frozen.h_PM, pd.Series([True, False, False], name='h_PM'))
            pd.testing.assert_series_equal(frozen.m_PM, pd.Series([True, False, False], name='m_PM'))


if __name__ == '__main__':
    unittest.main()
