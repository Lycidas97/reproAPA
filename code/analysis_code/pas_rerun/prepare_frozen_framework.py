"""Export PAS-independent majority states from the frozen raw-cell labels."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from polyapipe_prepare import digest,save


def main(a):
    a.out.mkdir(parents=True,exist_ok=False)
    manifest=pd.read_csv(a.atlas/'input_libraries.tsv',sep='\t');inputs={str(a.states):digest(a.states)};parts=[]
    for row in manifest.itertuples():
        p=a.atlas/'inputs'/row.study_id/row.library_id/'members.tsv.gz'
        assert digest(p)==row.members_sha256
        parts.append(pd.read_csv(p,sep='\t',usecols=['cell_id','metacell_id']))
        inputs[str(p)]=row.members_sha256
    membership=pd.concat(parts,ignore_index=True)
    states=pd.read_csv(a.states,sep='\t',usecols=['cell_id','state_group','library_key'])
    cells=membership.merge(states,on='cell_id',validate='one_to_one')
    assert len(cells)==int(manifest.cells.sum()) and cells.cell_id.is_unique
    assert cells.groupby('metacell_id').library_key.nunique().eq(1).all()
    counts=pd.crosstab(cells.metacell_id,cells.state_group)
    result=pd.DataFrame(dict(metacell_id=counts.index,state_group=counts.idxmax(axis=1).to_numpy(),
        state_purity=(counts.max(axis=1)/counts.sum(axis=1)).to_numpy(),
        library_key=cells.groupby('metacell_id').library_key.first().reindex(counts.index).to_numpy()))
    assert len(result)==manifest.metacells.sum()
    save(result,a.out/'metacell_states.tsv.gz')
    assert all(digest(p)==h for p,h in inputs.items())
    (a.out/'complete.json').write_text(json.dumps(dict(status='complete',cells=len(cells),metacells=len(result),
        inputs=inputs,outputs={'metacell_states.tsv.gz':digest(a.out/'metacell_states.tsv.gz')}),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for k in ['atlas','states','out']:p.add_argument('--'+k,type=Path,required=True)
    main(p.parse_args())
