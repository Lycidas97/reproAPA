"""Intent checks: protect the biological comparison, not only implementation mirroring."""
import unittest

import numpy as np
import pandas as pd

import analyze as a
from prepare_inputs import native_order,upstream_bounds


class CFImIntentTests(unittest.TestCase):
    def test_ugua_is_read_upstream_on_each_transcript_strand(self):
        # The same biologically oriented upstream RNA must give the same score;
        # moving the cleavage interbase boundary by one base loses the final TGTA.
        seq='A'*96+'TGTA';rc=lambda s:s.translate(str.maketrans('ACGT','TGCA'))[::-1]
        plus='C'*20+seq+'G'*120
        cut=120;lo,hi=upstream_bounds(cut,'+');self.assertEqual(plus[lo:hi],seq)
        minus='C'*120+rc(seq)+'G'*20
        lo,hi=upstream_bounds(cut,'-');self.assertEqual(rc(minus[lo:hi]),seq)
        self.assertEqual(a.motif_count(plus[20:120],'TGTA'),1)
        self.assertEqual(a.motif_count(plus[19:119],'TGTA'),0)
        self.assertEqual(a.window_sequence(seq,(-100,-40)),seq[:60])

    def test_motifs_must_be_fully_inside_fixed_window(self):
        s=list('A'*100)
        for i in [0,56,60]:s[i:i+4]=list('TGTA')
        seq=''.join(s)
        self.assertEqual(a.motif_count(a.window_sequence(seq,(-100,-40)),'TGTA'),2)
        self.assertEqual(a.motif_count(seq,'TGTA'),3)
        self.assertEqual(a.motif_count('TGTATGTA','TGTA'),2)

    def test_anchors_follow_strand_not_response_driver(self):
        cuts=[1050,1000,1070,1050]
        self.assertEqual(native_order(cuts,'+'),[1000,1050,1070])
        self.assertEqual(native_order(cuts,'-'),[1070,1050,1000])

    def test_shuffling_cannot_change_direction_margins(self):
        labels=np.array([0,0,1,2,0,2,2],np.int8);cells=[np.array([0,1,2,3]),np.array([4,5,6])]
        draws=a.label_draws(labels,cells,np.random.default_rng(71))
        # Fixed margins prevent a proximal-heavy mouse roster itself from creating
        # an apparent sequence signature during randomization.
        for ii in cells:
            for label in range(3):self.assertTrue(np.all((draws[:,ii]==label).sum(1)==int((labels[ii]==label).sum())))
        x=np.array([2,0,-1,1,0,3,1],float)
        expected=a.expectation(x,labels,cells,2);values=a.statistic_draws(x,draws,2)
        self.assertLess(abs(values.mean()-expected),6*values.std(ddof=1)/np.sqrt(len(values)))

    def test_absent_variation_does_not_mean_no_cfIm_mechanism(self):
        summary=a.null_summary(1.,1.,np.ones(a.NPERM))
        self.assertFalse(summary['informative']);self.assertEqual(summary['p_enrichment'],1)
        self.assertEqual(summary['status'],'fixed_under_matched_null')
        empty=a.null_summary(np.nan,np.nan,np.full(a.NPERM,np.nan))
        self.assertEqual(empty['status'],'empty_comparison_group')

    def test_covariate_residualization_is_label_independent(self):
        n=20;z=pd.DataFrame(dict(baseline=np.linspace(0,1,n),PAS_number=np.arange(n)%4+2,anchor_span=np.arange(n)*50+100,GC_mean=np.linspace(.1,.6,n),mean_GEX=np.arange(n)+1,APA_coverage=np.arange(n)+20,core_rank_difference=np.arange(n)%3-1,measured_pattern=['A;B']*n,proximal_capacity=[2]*n,distal_capacity=[2]*n))
        vals=3*z.baseline.to_numpy()+np.sin(np.arange(n));r,rank=a.residualize(vals,z)
        self.assertGreater(rank,1);self.assertLess(abs(r.mean()),1e-12)
        z['direction']=['proximal']*n;r2,_=a.residualize(vals,z)
        np.testing.assert_array_equal(r,r2)

    def test_complete_window_stays_inside_verified_terminal_utr(self):
        p=a.HERE/'background_inventory.tsv.gz'
        if not p.exists():self.fail('Sequence-qualified inventory required')
        z=pd.read_csv(p,sep='\t');z=z[z.status.eq('eligible')]
        self.assertGreater(len(z),0)
        for r in z.itertuples():
            self.assertTrue(r.full_primary_axis_tandem)
            self.assertIn(r.chosen_transcript,r.compatible_transcripts.split(';'))
            for cut in [r.proximal_cut,r.distal_cut]:
                lo,hi=upstream_bounds(int(cut),r.strand)
                self.assertGreaterEqual(lo,r.utr_start0);self.assertLessEqual(hi,r.utr_end0)
        self.assertFalse(z.gene_symbol.str.upper().eq('PDLIM7').any())


if __name__=='__main__':unittest.main(verbosity=2)
