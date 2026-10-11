"""Collect read-only fixed query windows after effect-independent local preparation."""
import json
import subprocess

import pandas as pd

from prepare_inputs import HERE,CONFIG,save


def main():
    config=json.loads(CONFIG.read_text());query=pd.read_csv(HERE/'sequence_queries.tsv',sep='\t')
    code='CONFIG = '+repr(config)+'\nQUERIES = '+repr(query.to_dict('records'))+'\n'+(HERE/'extract_remote.py').read_text()
    partial=HERE/'upstream100_sequences.partial.tsv'
    command='OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 /root/miniconda3/envs/scanpy/bin/python -'
    with partial.open('wb') as out,(HERE/'extract.log').open('wb') as log:
        subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=20','213docker',command],input=code.encode(),stdout=out,stderr=log,check=True)
    partial.rename(HERE/'upstream100_sequences.tsv');print((HERE/'extract.log').read_text(),flush=True)
    seq=pd.read_csv(HERE/'upstream100_sequences.tsv',sep='\t');assert len(seq)==len(query)
    assert not seq.duplicated(['species','gene_id','anchor']).any()
    seq['complete_ACGT']=seq.sequence.str.len().eq(100)&seq.sequence.str.fullmatch('[ACGT]{100}');save(seq,'upstream100_sequences.tsv')
    inventory=pd.read_csv(HERE/'background_before_sequence.tsv.gz',sep='\t',float_precision='round_trip')
    valid=seq.groupby(['species','gene_id']).complete_ACGT.all()
    for i,r in inventory[inventory.status.eq('sequence_pending')].iterrows():inventory.loc[i,'status']='eligible' if valid.loc[r.species,r.gene_id] else 'ambiguous_or_incomplete_sequence'
    inventory['paired_eligible']=False
    orth=pd.read_csv(HERE.parents[1]/'notes/2026-10-06-dual-main-plot-data/polyapipe/main/tables/cross_overlap/ortholog_matched_features.tsv.gz',sep='\t')[['human_gene_id','mouse_gene_id','human_gene_symbol','mouse_gene_symbol','homology_type','is_high_confidence']]
    eligible={sp:set(inventory[inventory.species.eq(sp)&inventory.status.eq('eligible')].gene_id) for sp in ['human','mouse']}
    pairs=orth[orth.human_gene_id.isin(eligible['human'])&orth.mouse_gene_id.isin(eligible['mouse'])].copy()
    for sp in ['human','mouse']:inventory.loc[inventory.species.eq(sp)&inventory.gene_id.isin(pairs[sp+'_gene_id']),'paired_eligible']=True
    save(inventory,'background_inventory.tsv.gz');save(pairs,'paired_ortholog_background.tsv')
    print('complete strict:',{sp:len(v) for sp,v in eligible.items()},'paired:',len(pairs),flush=True)


if __name__=='__main__':main()
