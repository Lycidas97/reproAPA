"""Guard biological time identity and the scope of display-only short-edge removal."""
import unittest
from types import SimpleNamespace

import numpy as np
import pandas as pd

from render_compact_state_trajectory_panel import prepare_display_edges


def prepare(rows):
    arrows = pd.DataFrame([
        dict(species='mouse', route=route, source=a, target=b,
             source_x=x, source_y=0., target_x=y, target_y=0.)
        for route, a, b, t1, t2, x, y in rows
    ])
    times = pd.DataFrame([
        dict(route=route, source=a, target=b, day_from=t1, day_to=t2,
             source_fraction=.5, target_fraction=.5)
        for route, a, b, t1, t2, x, y in rows
    ])
    axis = SimpleNamespace(transData=SimpleNamespace(transform=lambda xy: np.asarray(xy)*100))
    return prepare_display_edges(arrows, times, axis, 100)


class TimedDisplayTests(unittest.TestCase):
    def test_same_state_at_different_times_does_not_create_a_fork(self):
        display, audit = prepare([
            ('r', 'A', 'B', 0, 1, 0., 1.),
            ('r', 'B', 'C', 2, 3, 2., 3.),
            ('r', 'B', 'D', 4, 5, 4., 5.),
        ])
        self.assertNotEqual(display.iloc[0].target, display.iloc[1].source)
        self.assertNotEqual(display.iloc[1].source, display.iloc[2].source)
        self.assertEqual(set(audit.loc[audit.state.eq('B'), 'day']), {1, 2, 4})
        self.assertEqual(list(display.display_source_x), [0., 2., 4.])

    def test_observed_same_time_fork_is_preserved(self):
        display, _ = prepare([
            ('r', 'A', 'B', 1, 3, 0., 1.),
            ('r', 'A', 'C', 1, 3, 0., 2.),
        ])
        self.assertEqual(display.iloc[0].source, display.iloc[1].source)
        self.assertTrue(display.displayed.all())

    def test_short_cleanup_does_not_break_an_existing_chain(self):
        display, _ = prepare([
            ('isolated', 'A', 'B', 0, 1, 0., .01),
            ('connected', 'A', 'B', 0, 1, 0., .01),
            ('connected', 'B', 'C', 1, 2, .01, 1.),
        ])
        self.assertEqual(list(display.displayed), [False, True, True])

    def test_nonforward_time_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'forward time'):
            prepare([('r', 'A', 'B', 3, 1, 0., 1.)])


if __name__ == '__main__':
    unittest.main()
