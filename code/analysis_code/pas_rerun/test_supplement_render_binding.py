"""A supplement redraw cannot silently reuse an old PAS table or denominator."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
import pandas as pd
from pipeline import digest
from render_theme_bundle import SupplementBundle
from package_theme import resolve_sources


class SupplementBindingTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.adapter=SupplementBundle.__new__(SupplementBundle)
        self.adapter.reader=SimpleNamespace(root=self.root,manifest={},identity={'source_key':'current'},used={})
        self.adapter.renderer=SimpleNamespace(inputs={},pins={},cache={},sources=set(),fid='D04',active='D04A',evidence={'D04A':{'sources':{}}})
        self.adapter.bindings=[]

    def tearDown(self):self.tmp.cleanup()

    def table(self,name,rows):
        p=self.root/'tables'/name;p.parent.mkdir(parents=True,exist_ok=True)
        pd.DataFrame(rows).to_csv(p,sep='\t',index=False)
        self.adapter.reader.manifest['tables/'+name]={'sha256':digest(p)}
        return p

    def test_existing_old_file_is_not_a_fallback(self):
        with self.assertRaisesRegex(AssertionError,'Unbound supplement input'):
            self.adapter.read('profiles/missing.tsv')
        with self.assertRaises(KeyError):
            self.adapter.read('base:old-result/unknown.tsv')

    def test_removed_v3_cannot_reenter_final_display(self):
        self.table('same_pas/v3_ge62_comparison.tsv',[{'effect':1}])
        with self.assertRaisesRegex(AssertionError,'v3 is excluded'):
            self.adapter.read('hyp:C_reference/v3_ge62_comparison.tsv')

    def test_changed_value_is_rejected_before_plotting(self):
        p=self.table('example/current.tsv',[{'effect':1}]);p.write_text('effect\n9\n')
        with self.assertRaisesRegex(AssertionError,'Changed table'):
            self.adapter.read('example/current.tsv')

    def test_support_windows_keep_native_denominator(self):
        self.table('pas_qc/catalog_summary.tsv',[dict(study_id='A',species='mouse',layer='primary',
            sites=100,supported25_sites=20,supported50_sites=40,supported100_sites=70)])
        x=self.adapter.read('reference_comparison/PAS_support_comparison.tsv')
        self.assertEqual(x.denominator.tolist(),[100,100,100])
        self.assertEqual(x.supported_v2.tolist(),[20,40,70])
        self.assertEqual(self.adapter.bindings[0]['table'],'pas_qc/catalog_summary.tsv')

    def test_tandem_display_cannot_mix_other_contrasts_or_settings(self):
        rows=[dict(contrast=c,design='true_source',setting=s,background='available',n=n)
              for c,s,n in [('P_S','supported50_p08',8),('P_M','supported50_p08',99),('P_S','expanded_primary_p08',77)]]
        self.table('panel_statistics/tandem-comparison.tsv',rows)
        x=self.adapter.read('base:review-data/E03b-tandem-comparison.tsv')
        self.assertEqual(x.n.tolist(),[8])

    def test_reused_dataframe_remains_traceable_in_later_panels(self):
        ref='base:review-data/E14d-axis-decomposition-128.tsv'
        bound=dict(table='geometry/axis-decomposition.tsv',path='/verified/current.tsv',
                   sha256='current-hash',design_reference=ref,logical_panel='D27A')
        receipt={'figure':'D27','panels':{'D27A':{'sources':{ref:bound}},'D27B':{'sources':{}},'D27C':{'sources':{}}}}
        index={('D27','D27A',ref):bound}
        resolved=resolve_sources(receipt,index)
        for panel in ['D27B','D27C']:
            self.assertEqual(resolved[panel][ref]['sha256'],'current-hash')
            self.assertEqual(resolved[panel][ref]['reused_from'],'D27A')
        with self.assertRaises(AssertionError):
            resolve_sources({'figure':'D27','panels':{'D27B':{'sources':{}}}}, {})


if __name__=='__main__':unittest.main()
