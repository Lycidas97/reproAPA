import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
from scipy import sparse

from polyapipe_prepare import clusters, save
from polyapipe_effects import weights, metric_weights, state, reason, study


class PolyApipeIntentTests(unittest.TestCase):
    def test_adjacent_merge_cannot_chain_beyond_fifty(self):
        self.assertEqual(clusters(np.array([0,25,50,75,100,126])).tolist(), [0,0,0,1,1,2])

    def test_capture_replication_does_not_change_source_date_weight(self):
        f=pd.DataFrame(dict(time_a=[1,1,1,2],sample_unit_a=['A','A','B','C'],
                            time_b=[0]*4,sample_unit_b=['S']*4))
        np.testing.assert_allclose(weights(f), [.125,.125,.25,.5])

    def test_minus_strand_endpoint_weights_follow_transcription(self):
        f=pd.DataFrame(dict(cut_site=[100,200,300],strand=['-']*3,position=[1,.5,0]))
        w=metric_weights(f)
        np.testing.assert_allclose(w[:,1],[1,2/3,1/3])
        np.testing.assert_array_equal(w[:,2],[1,0,0])
        np.testing.assert_array_equal(w[:,3],[0,0,1])

    def test_gene_qualified_cells_define_thirty_cell_gate(self):
        data=np.array([[20,0],[20,0],[0,20],[9,9]])
        meta=pd.DataFrame(dict(n_cells=[9,9,9,100]))
        s=state(data,np.arange(4),meta)
        self.assertEqual(s['n_cells'],27)
        self.assertEqual(reason(s,s,3),'insufficient_gene_qualified_cells')

    def test_support_filter_does_not_rescale_full_primary_axis(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            f=pd.DataFrame(dict(atlas_pas_id=['p0','p1','p2'],feature_index_1based=[1,2,3],
                primary_gene_id=['new_gene']*3,primary_gene_symbol=['New']*3,
                keep_for_apa_primary=[True]*3,chrom=['1']*3,strand=['+']*3,
                cut_site=[100,200,300],matched_source_count_50nt=[0,1,1],species=['human']*3))
            save(f,p/'features.tsv.gz')
            meta=pd.DataFrame(dict(species=['human']*6,study_id=['GSE_TEST']*6,
                metacell_id=[f'm{i}' for i in range(6)],library_key=['target']*3+['source']*3,
                route=['R']*6,node=['P']*3+['S']*3,purity=[.85]*6,n_cells=[10]*6,matrix_row=np.arange(6)))
            save(meta,p/'metacells.tsv.gz');save(meta,p/'selected_metacells.tsv.gz')
            sparse.save_npz(p/'metacell_by_pas_counts.npz',sparse.csr_matrix(np.array([[0,6,54]]*3+[[0,54,6]]*3)))
            design=pd.DataFrame([dict(species='human',study_id='GSE_TEST',route='R',library_key_a='target',
                library_key_b='source',time_a='D16',time_b='D0',sample_unit_a='T',sample_unit_b='S',contrast='P_S',design='true_source')])
            save(design,p/'endpoint_library_design.tsv');save(design,p/'sampling_design.tsv')
            study(p,[('supported50_p08',True,.8,3)])
            r=pd.read_csv(p/'effects/supported50_p08_route_effects.tsv.gz',sep='\t').iloc[0]
            self.assertEqual(r.gene_id,'new_gene')  # No frozen old gene roster required.
            self.assertAlmostEqual(r.delta_DWUI,.4)
            self.assertAlmostEqual(r.MPRO,1.6)
            self.assertEqual(json.loads(r.PAS_ids),['p1','p2'])


if __name__ == '__main__':
    unittest.main()
