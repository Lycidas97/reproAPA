"""Protect common-catalogue identities, zero rows and integer conservation."""
import gzip
import json
from pathlib import Path
import subprocess
import tempfile
import unittest

import numpy as np
import pandas as pd
from polyapipe_common_prepare import catalogue, digest, member_axis


class CommonIntake(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp=tempfile.TemporaryDirectory();cls.root=Path(cls.tmp.name);cls.exe=cls.root/'aggregate'
        subprocess.run(['g++','-O2','-std=c++17',str(Path(__file__).with_name('polyapipe_common_aggregate.cpp')),
                        '-lz','-o',str(cls.exe)],check=True)

    @classmethod
    def tearDownClass(cls):cls.tmp.cleanup()

    def aggregate(self, text, passtext='p1\t0\np2\t1\np0\t2\n', cbtext='A\t0\nB\t1\n'):
        with tempfile.TemporaryDirectory(dir=self.root) as d:
            p=Path(d)
            with gzip.open(p/'counts.gz','wt') as f:f.write(text)
            (p/'pas').write_text(passtext);(p/'cb').write_text(cbtext)
            r=subprocess.run([str(self.exe),str(p/'counts.gz'),str(p/'pas'),str(p/'cb'),'3','3',str(p/'matrix'),str(p/'audit')],capture_output=True,text=True)
            if r.returncode:return r,None,None
            data=np.fromfile(p/'matrix',dtype=[('row','<u4'),('col','<u4'),('count','<u8')])
            return r,data,json.loads((p/'audit').read_text())

    def test_integer_sum_zero_rows_and_zero_pas_preserved(self):
        r,d,a=self.aggregate('gene\tcell\tcount\np2\tA\t3\np1\tA\t5\np1\tB\t7\n')
        self.assertEqual(r.returncode,0,r.stderr);self.assertEqual(a,dict(records=3,umi=15,observed_cells=2,nnz=3))
        self.assertEqual(d.tolist(),[(0,0,5),(0,1,3),(1,0,7)])

    def test_unknown_pas_and_barcode_fail_instead_of_dropping_counts(self):
        for line in ['unknown\tA\t1','p1\tother_library_CB\t1']:
            r,_,_=self.aggregate('gene\tcell\tcount\n'+line+'\n');self.assertNotEqual(r.returncode,0)

    def test_nonpositive_noninteger_overflow_fail(self):
        for value in ['0','-1','1.5','9223372036854775808']:
            r,_,_=self.aggregate('gene\tcell\tcount\np1\tA\t'+value+'\n');self.assertNotEqual(r.returncode,0)
        r,_,_=self.aggregate('gene\tcell\tcount\np1\tA\t9223372036854775807\np2\tA\t1\n')
        self.assertNotEqual(r.returncode,0)

    def test_axis_duplicates_and_changed_header_fail(self):
        r,_,_=self.aggregate('gene\tcell\tcount\np1\tA\t1\n',passtext='p1\t0\np1\t1\n');self.assertNotEqual(r.returncode,0)
        r,_,_=self.aggregate('gene\tbarcode\tcount\np1\tA\t1\n');self.assertNotEqual(r.returncode,0)

    def test_common_peaks_not_reclustered_and_base_interbase_distinct(self):
        p=self.root/'peaks.gff';p.write_text('##gff\n1\tpolyAends\tpolyAends\t101\t350\t.\t+\t.\tpeak="p1"; peakdepth="1"; misprime="False";\n1\tpolyAends\tpolyAends\t350\t599\t.\t-\t.\tpeak="p2"; peakdepth="2"; misprime="True";\n')
        f=catalogue(dict(catalogue=dict(path=str(p),sha256=digest(p)),native_pas_count=2,study_id='study',species='human'))
        self.assertEqual(len(f),2);self.assertEqual(f.cut_site.tolist(),[350,349])
        self.assertEqual(f.annotation_cut_site_0based_base.tolist(),[349,349])
        self.assertEqual(f.native_peak_id.tolist(),['p1','p2'])

    def test_library_local_normalization_collision_rejected(self):
        m=pd.DataFrame(dict(cell_id=['s|l|A-1','s|l|A-2'],barcode=['A-1','A-2'],metacell_id=['mc','mc']))
        meta=pd.DataFrame(dict(metacell_id=['mc'],n_cells=[2]))
        with self.assertRaisesRegex(AssertionError,'collision'):
            member_axis(m,meta,dict(cells=2,metacells=1,study_id='s',library_id='l'))

    def test_frozen_metacell_order_and_cell_sizes_authoritative(self):
        m=pd.DataFrame(dict(cell_id=['s|l|A-1','s|l|B-1'],barcode=['A-1','B-1'],metacell_id=['m1','m0']))
        meta=pd.DataFrame(dict(metacell_id=['m0','m1'],n_cells=[1,1]))
        axis=member_axis(m,meta,dict(cells=2,metacells=2,study_id='s',library_id='l'))
        self.assertEqual(axis.row.tolist(),[1,0])
        meta.loc[0,'n_cells']=2
        with self.assertRaises(AssertionError):member_axis(m,meta,dict(cells=2,metacells=2,study_id='s',library_id='l'))

    def test_wot_folder_alias_keeps_original_cell_identity(self):
        m=pd.DataFrame(dict(cell_id=['GSE115943|l|A-1'],barcode=['A-1'],metacell_id=['mc']))
        meta=pd.DataFrame(dict(metacell_id=['mc'],n_cells=[1]))
        axis=member_axis(m,meta,dict(cells=1,metacells=1,study_id='WOT_GSE115943',library_id='l'))
        self.assertEqual(axis.row.tolist(),[0])

    def test_cip_frozen_library_key_resolves_legacy_library_alias(self):
        m=pd.DataFrame(dict(cell_id=['GSE118258|gse118258_D0|A-1'],barcode=['A-1'],metacell_id=['mc']))
        meta=pd.DataFrame(dict(metacell_id=['mc'],n_cells=[1]))
        axis=member_axis(m,meta,dict(cells=1,metacells=1,study_id='GSE118258',
            library_id='gse118258_human_GSM3323825',library_key='GSE118258::gse118258_D0'))
        self.assertEqual(axis.row.tolist(),[0])


if __name__=='__main__':unittest.main()
