"""Mouse example panels must never borrow a human ortholog's measurements."""
import unittest
from types import SimpleNamespace

import matplotlib.pyplot as plt
import pandas as pd

from render_regulatory import inclusion_example
from render_external import inclusion_sensitivity


class ExampleSpeciesTests(unittest.TestCase):
    def tearDown(self):plt.close('all')

    def test_mouse_site_panel_excludes_same_named_human_sites(self):
        rows=[]
        for species,route,starts,value in [('mouse','CIP',[10,20,30],.2),('human','GSE118258',[40,50],.9)]:
            for start in starts:
                for state in 'SREPM':
                    rows.append(dict(species=species,route=route,gene='EIF4EBP2',cutoff=20,
                        kind='site',label='PCBP1',start0=start,end0=start+7,state=state,
                        equal_library=value,minimum=value,maximum=value))
        reader=SimpleNamespace(read=lambda _:pd.DataFrame(rows));fig,ax=plt.subplots()
        inclusion_example(reader,ax,'EIF4EBP2')
        self.assertEqual([t.get_text() for t in ax.get_legend().get_texts()],['10–17','20–27','30–37'])
        for line in ax.lines:
            self.assertFalse(any(v==90 for v in line.get_ydata()))

    def test_unavailable_mouse_locus_does_not_use_human_ortholog(self):
        frame=pd.DataFrame([dict(species='human',gene='HMGN2',cutoff=20,kind='interval')])
        reader=SimpleNamespace(read=lambda _:frame);fig,ax=plt.subplots()
        inclusion_example(reader,ax,'HMGN2')
        self.assertEqual(len(ax.lines),0)
        self.assertTrue(any('No newly qualified' in t.get_text() for t in ax.texts))

    def test_mouse_influence_panel_excludes_human_rows(self):
        rows=[dict(species=sp,route=route,gene='EIF4EBP2',contrast='SP',kind='site',label='PCBP1',
            start0=start,end0=start+7,cutoff=20,equal_library=.2,time_equal=.2,leave_one_min=.1,leave_one_max=.3)
            for sp,route,start in [('mouse','CIP',10),('human','GSE118258',40)]]
        reader=SimpleNamespace(read=lambda _:pd.DataFrame(rows));fig,ax=plt.subplots()
        inclusion_sensitivity(reader,ax,'EIF4EBP2')
        self.assertEqual([t.get_text() for t in ax.get_yticklabels()],['CIP 10–17'])


if __name__=='__main__':unittest.main()
