import csv
import gzip
import json
from pathlib import Path
import tempfile
import unittest

import pandas as pd

from polyapipe_finalize_annotation import finalize,sha


def save(path,rows):
    with gzip.open(path,'wt',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(rows)


class AnnotationAssemblyIntentTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.ids={s:[f'polyApipe:{s}:chr{c}:+:100' for c in [1,2]] for s in ['A','B']}
        lock={}
        for study,ids in self.ids.items():
            d=self.root/'studies'/study;d.mkdir(parents=True)
            p=d/'coordinate_features.tsv.gz'
            save(p,[dict(atlas_pas_id=pas,feature_index_1based=i+1,chrom=f'chr{i+1}',cut_site=100) for i,pas in enumerate(ids)])
            lock[str(p)]=sha(p)
        (self.root/'annotation_inputs_bedtools.json').write_text(json.dumps(lock))
        (self.root/'intake.native.json').write_text(json.dumps(dict(studies=[dict(study_id=s,atlas_pas=2) for s in self.ids])))
        self.chunks=[]
        for c in range(2):
            d=self.root/'annotation_chunks/human'/f'{c:06d}';d.mkdir(parents=True)
            rows=[dict(atlas_pas_id=self.ids[s][c],species='human',chrom=str(c+1),cut_site=999,
                       primary_assignment_status='assigned_primary',primary_gene_id=f'gene_{s}',keep_for_apa_primary='true') for s in self.ids]
            self.chunks.append((d,rows));self.write_chunk(d,rows)

    def tearDown(self):self.tmp.cleanup()

    def write_chunk(self,d,rows):
        p=d/'features.tsv.gz';save(p,rows)
        for name in ['sequence_qc.tsv.gz','annotation_candidates.tsv.gz']:(d/name).touch()
        (d/'complete.json').write_text(json.dumps(dict(species='human',block=int(d.name),rows=len(rows),features_sha256=sha(p))))

    def test_interleaved_studies_preserve_matrix_order_and_coordinate_source(self):
        annotations=pd.concat([pd.read_csv(d/'features.tsv.gz',sep='\t') for d,_ in self.chunks],ignore_index=True)
        expected={}
        for study in self.ids:
            raw=pd.read_csv(self.root/'studies'/study/'coordinate_features.tsv.gz',sep='\t')
            extra=annotations.drop(columns=[c for c in annotations if c in raw and c!='atlas_pas_id'])
            expected[study]=raw.merge(extra,on='atlas_pas_id',validate='one_to_one',how='left')
        result=finalize(self.root)
        self.assertEqual(result['rows'],4)
        for study in self.ids:
            actual=pd.read_csv(self.root/'studies'/study/'features.tsv.gz',sep='\t')
            pd.testing.assert_frame_equal(actual,expected[study],check_exact=True)

    def test_duplicate_pas_cannot_silently_replace_missing_matrix_column(self):
        d,rows=self.chunks[1];rows[0]['atlas_pas_id']=self.ids['A'][0];self.write_chunk(d,rows)
        with self.assertRaisesRegex(AssertionError,'Missing, duplicate or reordered'):finalize(self.root)
        self.assertFalse((self.root/'annotation.complete.json').exists())
        self.assertFalse((self.root/'studies/A/features.tsv.gz').exists())

    def test_missing_pas_prevents_publishing_study_tables(self):
        d,rows=self.chunks[1];self.write_chunk(d,rows[1:])
        with self.assertRaisesRegex(AssertionError,'Study row counts differ'):finalize(self.root)
        self.assertFalse((self.root/'annotation.complete.json').exists())

    def test_changed_chunk_is_rejected_before_its_rows_are_published(self):
        d,rows=self.chunks[1];rows[0]['primary_gene_id']='changed';save(d/'features.tsv.gz',rows)
        with self.assertRaisesRegex(AssertionError,'Chunk checksum changed'):finalize(self.root)

    def test_changed_matrix_axis_input_is_rejected(self):
        p=self.root/'studies/A/coordinate_features.tsv.gz'
        with p.open('ab') as f:f.write(b'changed')
        with self.assertRaisesRegex(AssertionError,'Annotation input changed'):finalize(self.root)


if __name__=='__main__':unittest.main()
