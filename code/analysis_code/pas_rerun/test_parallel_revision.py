import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd
import openpyxl
from downstream_runtime import Context
from downstream_apa_gex import permutations,permutation_task
from downstream_parallel import ordered_map
from downstream_external_structure import published_longread_rows


class ParallelRevisionTests(unittest.TestCase):
    def test_apa_gex_process_results_preserve_nulls_and_input_precision(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);config=root/'config.json'
            config.write_text(json.dumps(dict(project_root=tmp,native_run=tmp,frozen_gex=tmp,
                analysis_root=tmp,run_id='test',node_workers={'apa_gex':2})))
            ctx=Context(config,'apa_gex',root/'out')
            rng=np.random.default_rng(991);n=18
            f=pd.DataFrame(dict(species=['mouse']*n,guan_series_pattern=['not_applicable_mouse']*n))
            for suffix in ['1','2']:
                v=rng.choice([-.08,0,.08],n)
                for col in ['alo','ahi']:f[col+'_'+suffix]=v
                for col in ['glo','ghi','GEX']:f[col+'_'+suffix]=rng.choice([-.5,0,.5],n)
                f['meas_'+suffix]=True;f['depth_'+suffix]=rng.uniform(1,100,n)
                f['PAS_'+suffix]=3
                for col in ['reference_WPP','reference_WPP_min','reference_WPP_max']:f[col+'_'+suffix]=.5
            ctx.table(f,'permutation_input.tsv.gz')
            pd.testing.assert_frame_equal(ctx.read(ctx.out/'permutation_input.tsv.gz'),f,check_exact=True)
            tasks=[('mouse',.02,m,3,40) for m in ['response_identity','direction_specific']]
            expected=[permutations(f,*t) for t in tasks]
            actual=ordered_map(ctx,permutation_task,tasks)
            for a,b in zip(expected,actual):
                for x,y in zip(a,b):pd.testing.assert_frame_equal(pd.DataFrame(x),pd.DataFrame(y),check_exact=True)

    def workbook(self,root,rows):
        p=root/'published.xlsx';w=openpyxl.Workbook();s=w.active;s.title='Classification-file'
        s.append(['Transcript name ID','Supported either by junction or TSS'])
        for row in rows:s.append(row)
        w.save(p);return p

    def test_new_candidate_is_selected_from_complete_classification(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=self.workbook(Path(tmp),[('Rad23b.old',True),('Coa5.G485.18',True),('unrelated',False)])
            f=published_longread_rows(p,['Coa5.G485.18','Rad23b.old'])
            self.assertEqual(set(f.index),{'Coa5.G485.18','Rad23b.old'})
            self.assertTrue(f['Supported either by junction or TSS'].all())
            with self.assertRaisesRegex(AssertionError,'Missing published'):published_longread_rows(p,['missing'])

    def test_duplicate_classification_cannot_be_silently_selected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=self.workbook(Path(tmp),[('new',True),('new',True)])
            with self.assertRaisesRegex(AssertionError,'Duplicate published'):published_longread_rows(p,['new'])


if __name__=='__main__':unittest.main()
