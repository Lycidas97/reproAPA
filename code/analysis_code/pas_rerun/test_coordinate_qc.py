"""Coordinate denominators must use the exact primary QC axis used by intervals."""
import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from downstream_coordinate_states import prepare_family


class CoordinateQCTests(unittest.TestCase):
    def fixture(self,unknown_primary=False):
        f=pd.DataFrame(dict(atlas_pas_id=['p1','p2','excluded'],feature_index_1based=[1,2,3],
            primary_gene_id=['g','g',None],matched_source_count_50nt=[1,1,0],
            keep_for_apa_primary=[True,True,False],transcript_model_conflict=[np.nan if unknown_primary else False,False,np.nan],
            multi_gene_ambiguous=[False,False,np.nan],internal_priming_flag=[False,False,np.nan],sequence_qc_status=['computed','computed',None]))
        study=SimpleNamespace(f=f,counts=csr_matrix([[2,3,900],[5,7,800]]))
        raw=dict(gene_id='g',pas=[{'atlas_pas_id':'p1'},{'atlas_pas_id':'p2'}])
        annotation=dict(pas=[{'atlas_pas_id':'a'}]);mapping=pd.DataFrame(dict(native_pas_id=['p1','p2'],annotation_id=['a','a']))
        def strict(x):
            assert x.transcript_model_conflict.notna().all()
            return x.keep_for_apa_primary&~x.transcript_model_conflict&~x.multi_gene_ambiguous&~x.internal_priming_flag
        return study,raw,annotation,mapping,strict

    def test_excluded_unknown_qc_does_not_enter_family_denominator(self):
        family,total,projection=prepare_family(*self.fixture())
        np.testing.assert_array_equal(total,[5,12]);np.testing.assert_array_equal(family,[[2,3],[5,7]])
        np.testing.assert_array_equal(family@projection,[[5],[12]])

    def test_unknown_primary_qc_remains_a_hard_failure(self):
        with self.assertRaisesRegex(ValueError,'transcript_model_conflict'):prepare_family(*self.fixture(True))

if __name__=='__main__':unittest.main()
