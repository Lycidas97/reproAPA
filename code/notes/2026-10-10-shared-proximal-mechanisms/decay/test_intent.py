"""Tests protect biological interpretation, coordinate orientation and unknown handling."""
import unittest
import pandas as pd
from prepare_events import classify,retained,HERE
from analyze import occurrence_offsets

class IntentTests(unittest.TestCase):
 def structures(self,strand='+',utr=True):
  u=pd.DataFrame([dict(chrom='1',strand=strand,utr_start=101,utr_end=200,transcript_id='tx')]) if utr else pd.DataFrame(columns=['chrom','strand','utr_start','utr_end','transcript_id'])
  m=pd.DataFrame([dict(chrom='1',strand=strand,start0=0,end0=200,transcript_id='tx',kind='transcript'),dict(chrom='1',strand=strand,start0=0,end0=50,transcript_id='tx',kind='exon'),dict(chrom='1',strand=strand,start0=100,end0=200,transcript_id='tx',kind='exon')])
  return u,m
 def test_retained_nucleotide_tracks_the_long_RNA_boundary(self):
  self.assertEqual(retained(150,'+'),149);self.assertEqual(retained(150,'-'),150)
 def test_plus_final_UTR_nucleotide_is_kept_but_next_cut_is_not_UTR(self):
  u,m=self.structures();self.assertEqual(classify(150,200,'chr1','+',u,m)[0],'strict_tandem_UTR');self.assertNotEqual(classify(150,201,'chr1','+',u,m)[0],'strict_tandem_UTR')
 def test_minus_boundary_uses_cut_plus_one_in_GTF_coordinates(self):
  u,m=self.structures('-');self.assertEqual(classify(150,199,'chr1','-',u,m)[0],'strict_tandem_UTR');self.assertNotEqual(classify(150,200,'chr1','-',u,m)[0],'strict_tandem_UTR')
 def test_same_terminal_exon_without_CDS_boundary_evidence_is_not_UTR(self):
  u,m=self.structures(utr=False);self.assertEqual(classify(150,190,'chr1','+',u,m)[0],'terminal_exon_nonUTR')
 def test_intronic_cleavage_is_not_long_UTR_loss(self):
  u,m=self.structures(utr=False);self.assertEqual(classify(80,190,'chr1','+',u,m)[0],'intronic_compatible')
 def test_different_terminal_exons_are_not_tandem_UTR(self):
  u,m=self.structures(utr=False);m=pd.concat([m,pd.DataFrame([dict(chrom='1',strand='+',start0=210,end0=300,transcript_id='other',kind='exon'),dict(chrom='1',strand='+',start0=0,end0=300,transcript_id='other',kind='transcript')])],ignore_index=True)
  self.assertEqual(classify(190,280,'chr1','+',u,m)[0],'alternative_terminal_exon_compatible')
 def test_RNA_seed_is_not_scanned_in_its_own_orientation(self):
  self.assertEqual(occurrence_offsets('AAGTGCT','AGCACTT'),[]);self.assertEqual(occurrence_offsets('AGCACTT','AGCACTT'),[0])
 def test_8mer_is_nested_in_both_7mer_classes_without_pseudoreplication(self):
  s='AGCACTTA';self.assertEqual(occurrence_offsets(s,'AGCACTT'),[0]);self.assertEqual(occurrence_offsets(s,'GCACTTA'),[1]);self.assertEqual(occurrence_offsets(s,'AGCACTTA'),[0])
 def test_U_rich_runs_use_overlapping_binding_windows(self):
  self.assertEqual(occurrence_offsets('TTTTTT','TTTT'),[0,1,2])
 def test_noneligible_core_is_not_imputed_as_a_seed_negative(self):
  c=pd.read_csv(HERE/'paired_coverage.tsv.gz',sep='\t');b=pd.read_csv(HERE/'paired_motif_background.tsv.gz',sep='\t');self.assertEqual(int(c.foreground_002.sum()),27);self.assertEqual(int(b.foreground_002.sum()),8);self.assertFalse(c.loc[~c.paired_eligible].human_gene_id.isin(b.human_gene_id).any())
 def test_gene_wide_net_shift_does_not_force_each_driver_to_the_same_direction(self):
  # Multiple sites can shift the weighted mean proximally while the largest gain/loss pair is distal.
  positions=[0.,.1,.2,1.];changes=[.12,-.13,.14,-.13]
  self.assertLess(sum(p*d for p,d in zip(positions,changes)),0)
  gain=max(range(4),key=lambda k:changes[k]);loss=min(range(4),key=lambda k:changes[k])
  self.assertGreater(positions[gain],positions[loss])
 def test_negative_results_do_not_claim_that_the_pathway_is_excluded(self):
  r=pd.read_csv(HERE/'enrichment_tests.tsv',sep='\t');self.assertTrue(r.q18.eq(1).all());self.assertFalse(r.positive_evidence_gate.any());self.assertTrue((r.max_abs_SMD.gt(.2)).all())
if __name__=='__main__':unittest.main(verbosity=2)
