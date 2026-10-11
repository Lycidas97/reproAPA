"""Unknown QC outside the primary axis must not invalidate eligible intervals."""
import ast
from pathlib import Path
import unittest
import numpy as np
import pandas as pd
from downstream_regulatory import strict_interval_pas


class IntervalQCTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path=Path(__file__).resolve().parents[2]/'notes/2026-10-01-mpro-biology/legacy_prepare_pilot.py'
        tree=ast.parse(path.read_text())
        nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='strict_pas']
        ns={};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(path),'exec'),ns)
        cls.strict_pas=staticmethod(ns['strict_pas'])

    def frame(self):
        return pd.DataFrame(dict(keep_for_apa_primary=[True,True,True,False],
            transcript_model_conflict=[False,True,False,False],multi_gene_ambiguous=[False]*4,
            internal_priming_flag=[False,False,True,False],sequence_qc_status=['computed']*4,
            primary_gene_id=['a','b','c','d']),index=[10,20,30,40])

    def test_complete_annotation_preserves_exact_legacy_selection(self):
        f=self.frame();pd.testing.assert_series_equal(strict_interval_pas(f,self.strict_pas),self.strict_pas(f))

    def test_excluded_unannotated_rows_remain_excluded_without_filling_qc(self):
        f=self.frame()
        for col in ['transcript_model_conflict','multi_gene_ambiguous','internal_priming_flag','sequence_qc_status','primary_gene_id']:
            f[col]=f[col].astype(object);f.loc[40,col]=np.nan
        before=f.copy(deep=True)
        with self.assertRaisesRegex(ValueError,'transcript_model_conflict'):self.strict_pas(f)
        self.assertEqual(strict_interval_pas(f,self.strict_pas).tolist(),[True,False,False,False])
        pd.testing.assert_frame_equal(f,before)

    def test_unknown_qc_on_any_primary_row_still_fails(self):
        for col in ['transcript_model_conflict','multi_gene_ambiguous','internal_priming_flag']:
            for value in [np.nan,'unknown','false']:
                with self.subTest(col=col,value=value):
                    f=self.frame();f[col]=f[col].astype(object);f.loc[10,col]=value
                    with self.assertRaisesRegex(ValueError,col):strict_interval_pas(f,self.strict_pas)

    def test_primary_membership_cannot_be_unknown_or_string(self):
        for value in [np.nan,'false']:
            f=self.frame();f['keep_for_apa_primary']=f.keep_for_apa_primary.astype(object);f.loc[40,'keep_for_apa_primary']=value
            with self.assertRaisesRegex(ValueError,'keep_for_apa_primary'):strict_interval_pas(f,self.strict_pas)

    def test_empty_primary_axis_is_a_false_mask_with_original_index(self):
        f=self.frame();f['keep_for_apa_primary']=False
        f['transcript_model_conflict']=np.nan
        result=strict_interval_pas(f,self.strict_pas)
        self.assertFalse(result.any());self.assertEqual(result.index.tolist(),f.index.tolist())


if __name__=='__main__':unittest.main()
