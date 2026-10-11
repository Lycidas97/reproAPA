"""Run the prespecified DWUI source-range sensitivity as a separate receipt."""
import argparse

import numpy as np
import pandas as pd

from downstream_runtime import Context
from downstream_recurrence import numerical,cross_species_current


def run(ctx):
    events=ctx.read(ctx.analysis/'recurrence/study_event_covariates.tsv.gz')
    ortholog=ctx.read(ctx.resource('ortholog'))
    assert ortholog.human_gene_id.is_unique and ortholog.mouse_gene_id.is_unique
    diagnostics,archive=cross_species_current(ctx,numerical(ctx),events,ortholog,10000,20261003,model='Range10')
    ctx.table(pd.concat(diagnostics) if diagnostics else pd.DataFrame(columns=['species','metric','human_study',
        'mouse_study','threshold','stratum','genes']),'null_mobility.tsv.gz')
    np.savez_compressed(ctx.out/'null_draws.npz',**archive)
    ctx.finish(dict(new_ortholog_backgrounds=True,source_min_max_bins_rebuilt=True,
        same_DWUI_primary_event_definition=True,separate_8_16_24_families=True),
        permutations=10000,scope='Prespecified DWUI Range10 source-reference sensitivity; Q5 multimetric results remain separate')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--config',required=True);p.add_argument('--out',required=True)
    a=p.parse_args();run(Context(a.config,'cross_sharing_range10',a.out))
