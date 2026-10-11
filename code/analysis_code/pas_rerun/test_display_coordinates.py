import json
from types import SimpleNamespace
import unittest
import numpy as np
import pandas as pd
from downstream_display import pas_coordinate_lookup,coordinates_for_row

class DisplayCoordinateTests(unittest.TestCase):
    def test_opaque_native_ids_use_frozen_coordinates_and_vector_order(self):
        f=pd.DataFrame(dict(study_id=['study']*3,primary_gene_id=['gene']*3,
            atlas_pas_id=['chr4_129537773_r','opaque-second','third'],chrom=['chr4']*3,
            strand=['-']*3,cut_site=[100,500,300]))
        row=SimpleNamespace(study_id='study',gene_id='gene',PAS_ids=json.dumps(['opaque-second','chr4_129537773_r','third']))
        x=coordinates_for_row(pas_coordinate_lookup(pd.concat([f,f])),row)
        np.testing.assert_array_equal(x,[-500,-100,-300])
        usage=np.array([.2,.3,.5]);order=np.argsort(x)
        np.testing.assert_allclose(np.cumsum(usage[order]),[.2,.7,1])

    def test_mismatched_study_or_conflicting_axis_fails(self):
        f=pd.DataFrame(dict(study_id=['s'],primary_gene_id=['g'],atlas_pas_id=['id'],chrom=['chr1'],strand=['+'],cut_site=[42]))
        with self.assertRaisesRegex(AssertionError,'absent'):
            coordinates_for_row(pas_coordinate_lookup(f),SimpleNamespace(study_id='different',gene_id='g',PAS_ids='["id"]'))
        with self.assertRaisesRegex(AssertionError,'Ambiguous'):
            pas_coordinate_lookup(pd.concat([f,f.assign(cut_site=43)]))

if __name__=='__main__':unittest.main()
