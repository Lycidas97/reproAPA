"""Independently verify native aggregation used the frozen cell-to-MC mapping."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from polyapipe_prepare import digest, save


def verify(source, frozen, run):
    libraries=pd.read_csv(source/'input_libraries.tsv',sep='\t')
    bindings=pd.read_csv(frozen,sep='\t')
    records=[]
    for row in libraries.itertuples(index=False):
        path=source/'inputs'/row.study_id/row.library_id/'members.tsv.gz'
        assert digest(path)==row.members_sha256
        native=pd.read_csv(path,sep='\t')
        binding=bindings[bindings.library_id.eq(row.library_id)&bindings.role.eq('membership')].iloc[0]
        assert digest(binding.path)==binding.sha256
        original=pd.read_csv(binding.path,sep='\t',compression='gzip')
        cols=['cell_id','barcode','metacell_id']
        a=native[cols].sort_values(cols).reset_index(drop=True)
        b=original[cols].sort_values(cols).reset_index(drop=True)
        assert a.equals(b), (row.library_id,'changed cell or metacell assignment')
        bare=native.barcode.str.split('-').str[0]
        assert bare.is_unique and len(native)==row.cells
        folder=source/'polyapipe'/row.study_id/row.library_id
        audit=json.loads((folder/'count_audit.json').read_text())
        counts=sparse.load_npz(folder/'metacell_by_site_counts.npz').tocsr()
        ids=pd.read_csv(folder/'metacell_index.tsv',sep='\t').metacell_id
        assert np.issubdtype(counts.dtype,np.integer) and (counts.data>=0).all()
        assert set(ids)==set(native.metacell_id) and counts.shape[0]==len(ids)
        assert audit['input_cells']==len(native) and audit['input_metacells']==len(ids)
        assert int(counts.sum())==audit['integer_molecules']
        assert audit['unknown_barcode_or_peak_entries']==0 and audit['normalized_barcode_mapping_one_to_one']
        records.append(dict(study_id=row.study_id,library_id=row.library_id,cells=len(native),metacells=len(ids),
            frozen_members_sha256=binding.sha256,native_members_sha256=row.members_sha256,
            native_counts_sha256=audit['native_counts_sha256'],
            mc_matrix_sha256=digest(folder/'metacell_by_site_counts.npz'),
            cell_id_barcode_metacell_mapping_identical=True,barcode_normalization_one_to_one=True))
    save(pd.DataFrame(records),run/'frozen_native_membership_audit.tsv')
    (run/'frozen_native_membership_audit.json').write_text(json.dumps(dict(status='PASS',libraries=len(records),
        cells=sum(x['cells'] for x in records),metacells=sum(x['metacells'] for x in records),
        raw_native_counts_reparsed_locally=False,
        counts_evidence='Completed native count audits plus integer MC matrix validation; raw cell-peak counts remain on analysis host',
        checks=dict(exact_cell_barcode_MC_mapping=True,native_member_hashes=True,integer_conservation=True)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--frozen',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
    a=p.parse_args();verify(a.source,a.frozen,a.run)
