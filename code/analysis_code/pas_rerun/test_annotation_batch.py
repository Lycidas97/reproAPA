"""Scientific selection boundaries must survive the batch implementation.

Run with PAS_ANNOTATION_CODE pointing to the pinned postprocessing package.
"""
import os
import sys
import unittest

import numpy as np
import pandas as pd

sys.path.insert(0, os.environ['PAS_ANNOTATION_CODE'])
from postprocessing.core import annotation as ann
from polyapipe_assignment_batch import choose


def hit(region, gene='g1', transcript='t1', distance=100, biotype='protein_coding', **extra):
    priorities={'annotated_3utr_pas':2,'terminal_exon_pas':3,'transcript_end_zone_pas':4,
        'terminal_intron_pas':7,'internal_intron_pas':8,'internal_exon_pas':9,
        'coding_exon_pas':10,'gene_downstream_pas':11}
    return dict(candidate_region_class=region,candidate_gene_id=gene,candidate_transcript_id=transcript,
        candidate_strand='+',strand_relation='same_strand',candidate_transcript_end=10000,
        candidate_gene_biotype=biotype,region_priority=priorities[region],
        distance_to_transcript_end=abs(distance),downstream_distance_to_transcript_end=distance,
        candidate_model_3p_distance=abs(distance),overlap_bases=51,**extra)


def tables(cases):
    pmap=ann.normalize_pas_features(pd.DataFrame([
        dict(chrom='1',start=9999+i,end=10000+i,name=name,score=10000+i,strand='+')
        for i,name in enumerate(cases)]),'human','GRCh38_2024_A')
    records=[]
    for pas,rows in zip(pmap.itertuples(index=False),cases.values()):
        for i,update in enumerate(rows):
            record=ann._unassigned_candidate(pas)
            record.update(update);record['candidate_id']=f'{pas.pas_id}:c{i:03d}'
            records.append(record)
    return pmap,ann.annotate_neighbor_gene_conflicts(pd.DataFrame(records))


class BatchPolicy(unittest.TestCase):
    def compare(self,cases):
        pmap,candidates=tables(cases)
        expected=ann.choose_primary_assignments(pmap,candidates)
        actual=choose(ann,pmap,candidates)
        pd.testing.assert_frame_equal(expected,actual,check_exact=True)
        return actual.set_index('raw_peak_id')

    def test_tes_boundaries_and_incomplete_provenance(self):
        cases={f'd{d}':[hit('transcript_end_zone_pas',distance=d)] for d in [-1,0,250,251,2000,2001]}
        for field in ['candidate_gene_id','candidate_transcript_id','candidate_transcript_end',
                      'candidate_strand','downstream_distance_to_transcript_end']:
            row=hit('transcript_end_zone_pas');row[field]='';cases[field]=[row]
        row=hit('transcript_end_zone_pas');row['downstream_distance_to_transcript_end']=float('inf');cases['infinite']=[row]
        result=self.compare(cases)
        self.assertEqual(result.loc['d250','selected_tes_policy_class'],'core_tes')
        self.assertEqual(result.loc['d251','selected_tes_policy_class'],'qualified_tes2kb_rescue')
        self.assertEqual(result.loc['d2000','selected_tes_policy_class'],'qualified_tes2kb_rescue')
        self.assertEqual(result.loc['d2001','selected_tes_policy_class'],'policy_ineligible_tes')

    def test_guarded_tiers_biotype_and_exact_gene_ties(self):
        tes=lambda d,g='g1':hit('transcript_end_zone_pas',gene=g,distance=d)
        exon=lambda g,d,b='protein_coding':hit('internal_exon_pas',gene=g,distance=d,biotype=b)
        cases={
            'canonical_over_core':[hit('annotated_3utr_pas'),tes(0,'g2')],
            'canonical_tie':[hit('annotated_3utr_pas'),hit('annotated_3utr_pas',gene='g2')],
            'core_over_exon':[tes(250),exon('g2',1)],
            'core_tie':[tes(250),tes(250,'g2')],
            'rescue_tie':[tes(251),tes(251,'g2')],
            'near_rescue':[tes(251),tes(252,'g2')],
            'exon_over_rescue':[tes(251),exon('g2',1000)],
            'intron_over_rescue':[tes(251),hit('internal_intron_pas',gene='g2')],
            'terminal_intron_precedes_internal':[hit('terminal_intron_pas',biotype='pseudogene'),hit('internal_intron_pas',gene='g2')],
            'biotype_over_distance':[exon('g1',1,'lncRNA'),exon('g2',1000)],
            'structural_tie':[exon('g1',1),exon('g2',1)],
            'terminal_cds':[hit('terminal_exon_pas'),hit('coding_exon_pas')],
            'terminal_cds_with_utr':[hit('terminal_exon_pas'),hit('coding_exon_pas'),hit('annotated_3utr_pas')],
            'same_gene_alternative_intron':[hit('terminal_exon_pas'),hit('internal_intron_pas',transcript='t2')],
            'different_gene_intron':[hit('terminal_exon_pas'),hit('internal_intron_pas',gene='g2')],
            'unassigned':[{}]}
        result=self.compare(cases)
        self.assertEqual(result.loc['biotype_over_distance','primary_gene_id'],'g2')
        self.assertEqual(result.loc['biotype_over_distance','biotype_assignment_conflict'],'true')
        self.assertEqual(result.loc['structural_tie','primary_assignment_status'],'ambiguous_multi_gene')
        self.assertEqual(result.loc['same_gene_alternative_intron','transcript_context_class'],'alternative_terminal_intron')
        self.assertEqual(result.loc['different_gene_intron','transcript_context_class'],'cross_gene_terminal_intron')

    def test_opposite_strand_retains_legacy_tie_order(self):
        rows=[hit('internal_exon_pas',gene=f'g{i}') for i in range(50)]
        for r in rows:r.update(strand_relation='opposite_strand',candidate_strand='-')
        result=self.compare({'opposite':rows})
        self.assertEqual(result.loc['opposite','primary_assignment_status'],'excluded_antisense')

    def test_empty_model_families_preserve_unassigned_sites(self):
        result=self.compare({'unassigned':[{}]})
        self.assertEqual(result.loc['unassigned','primary_assignment_status'],'unassigned')
        self.assertEqual(result.loc['unassigned','transcript_model_conflict'],'false')

    def test_random_mixed_models_against_pinned_selector(self):
        rng=np.random.default_rng(20261003)
        classes=['annotated_3utr_pas','terminal_exon_pas','transcript_end_zone_pas',
            'terminal_intron_pas','internal_intron_pas','internal_exon_pas','coding_exon_pas']
        cases={}
        for i in range(200):
            rows=[]
            for _ in range(rng.integers(1,15)):
                row=hit(rng.choice(classes),gene=f'g{rng.integers(4)}',transcript=f't{rng.integers(8)}',
                    distance=int(rng.choice([0,1,250,251,2000,2001])),biotype=rng.choice(['protein_coding','lncRNA','pseudogene']))
                if rng.random()<.25:row.update(strand_relation='opposite_strand',candidate_strand='-')
                rows.append(row)
            cases[f'random{i}']=rows
        self.compare(cases)


if __name__=='__main__':unittest.main()
