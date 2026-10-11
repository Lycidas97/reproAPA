"""Small BAM fixture: molecule rules, cell identity and strand endpoint semantics."""
from pathlib import Path
import tempfile
import unittest

import numpy as np
import pandas as pd
import pysam

from downstream_recount import gff,make_windows,recount_library


class NativeWindowRecountTests(unittest.TestCase):
    def test_one_umi_per_cell_and_atlas_site_without_erasing_sample_overlap(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);path=root/'calls.gff'
            path.write_text('1\ttool\tpeak\t81\t110\t.\t+\t.\tpeak="a";\n'
                            '1\ttool\tpeak\t91\t120\t.\t+\t.\tpeak="b";\n'
                            '1\ttool\tpeak\t101\t130\t.\t-\t.\tpeak="c";\n')
            sample=gff(path).assign(atlas_column_0based=[0,0,1],target_index=[0,0,1])
            self.assertEqual(sample.endpoint.tolist(),[109,119,100])
            feature=pd.DataFrame(dict(atlas_column_0based=[0,1],atlas_pas_id=['plus','minus'],chrom=['chr1']*2,
                                      strand=['+','-'],cut_site=[115,100],keep_for_apa_primary=[True,True]))
            windows=make_windows(feature,sample);windows['target_index']=windows.atlas_column_0based
            bam=root/'reads.bam'
            with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'1','LN':1000}]}) as writer:
                for i,(start,reverse,cb,ub,flag) in enumerate([(95,False,'A-1','X',0),(95,False,'A-1','X',0),
                    (95,False,'B-1','X',0),(95,False,'Z-1','Y',0),(95,False,'A-1','Y',1024),
                    (95,False,'A-1','Y',256),(100,True,'A-1','X',0),(100,True,'A-1','X',0)]):
                    read=pysam.AlignedSegment();read.query_name=str(i);read.query_sequence='A'*10
                    read.flag=flag+(16 if reverse else 0);read.reference_id=0;read.reference_start=start;read.mapping_quality=60
                    read.cigar=[(0,10)];read.set_tag('CB',cb);read.set_tag('UB',ub);writer.write(read)
            pysam.index(str(bam))
            meta=pd.DataFrame(dict(metacell_id=['mc'],n_cells=[2]))
            binding=pd.DataFrame(dict(barcode=['A-1','B-1'],metacell_id=['mc']*2))
            counts,audit=recount_library(str(bam),binding,meta,feature,sample,windows)
            np.testing.assert_array_equal(counts['original'].toarray(),[[4,1]])
            for mode in ['atlas_dedup','transferred_observed','transferred','primary_clipped','all_clipped']:
                np.testing.assert_array_equal(counts[mode].toarray(),[[2,1]])
            self.assertEqual(audit['eligible_reads'],5)

    def test_unselected_neighbor_clips_target_and_cannot_add_molecules(self):
        feature=pd.DataFrame(dict(atlas_column_0based=[0,1,2],atlas_pas_id=['a','b','c'],chrom=['chr1']*3,
            strand=['+']*3,cut_site=[100,140,200],keep_for_apa_primary=[True,False,True]))
        origin=pd.DataFrame(dict(atlas_column_0based=[0,1,2],start=[49,89,149],end=[150,190,250],endpoint=[99,139,199]))
        windows=make_windows(feature,origin)
        full=windows[windows['mode'].eq('all_clipped')].set_index('atlas_pas_id')
        primary=windows[windows['mode'].eq('primary_clipped')].set_index('atlas_pas_id')
        self.assertLess(full.loc['a','end'],primary.loc['a','end'])
        self.assertGreater(full.loc['c','start'],primary.loc['c','start'])


if __name__=='__main__':unittest.main()
