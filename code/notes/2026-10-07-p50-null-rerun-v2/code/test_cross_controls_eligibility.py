"""Calibration must use the same measurable event support as observed sharing."""
import ast
from pathlib import Path
import unittest

import numpy as np
import pandas as pd

from downstream_cross import MODELS


class CrossControlsEligibilityTests(unittest.TestCase):
    def test_every_ladder_step_and_injection_keeps_both_species_eligible(self):
        path = Path(__file__).resolve().parents[2] / 'notes/2026-10-02-null-model-controls/common.py'
        if not path.exists():
            path = Path(__file__).resolve().parent / 'numerics/common.py'
        functions = [item for item in ast.parse(path.read_text()).body
                     if isinstance(item, ast.FunctionDef) and item.name in ['groups', 'binary_draw', 'mobility']]
        ns = dict(np=np, pd=pd)
        exec(compile(ast.Module(body=functions, type_ignores=[]), str(path), 'exec'), ns)
        for model, fields in MODELS.items():
            frame = pd.DataFrame({field: [0, 0, 0, 0] for field in fields})
            # Reference nodes miss a valid paired effect, so effect_n alone cannot constrain eligibility.
            frame['actual_effect_n'] = [2, 2, 1, 0]
            labels = np.array([True, False, False, False])
            cells = ns['groups'](frame, fields)
            probability, _ = ns['mobility'](labels, cells)
            draws = ns['binary_draw'](labels, cells, 101, np.random.default_rng(7))
            with self.subTest(model=model):
                self.assertTrue(np.array_equal(probability, [.5, .5, 0., 0.]))
                self.assertFalse(draws[:, 2:].any())
                self.assertTrue((draws.sum(1) == 1).all())
                self.assertTrue(draws[:, 0].any() and draws[:, 1].any())


if __name__ == '__main__':
    unittest.main()
