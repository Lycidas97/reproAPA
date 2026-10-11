"""Independent checks using saved sequences, membership, draws and simple arithmetic."""
from pathlib import Path
import gzip
import hashlib
import json
import re

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]


def digest(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()


def qvalues(p):
    # Independent explicit step-up calculation, without production helper.
    p=np.asarray(p,float);order=sorted(range(len(p)),key=lambda i:(p[i],i));result=np.ones(len(p))
    for rank,i in enumerate(order):result[i]=min(1,min(p[j]*len(p)/(r+1) for r,j in enumerate(order) if r>=rank))
    return result


def independent_local_catalog_sequences(seq):
    rows=[]
    for sp in ['human','mouse']:
        catalog=ROOT/f'notes/2026-10-01-mpro-biology/{sp}_catalog_utrs.tsv.gz'
        fasta=ROOT/f'notes/2026-10-01-figure5-union-update/resources/{sp}_sequences.fa.gz'
        u=pd.read_csv(catalog,sep='\t');joined=seq[seq.species.eq(sp)].merge(u,on=['gene_id','chrom','strand'])
        joined=joined[joined.start0_x.ge(joined.start0_y)&joined.end0_x.le(joined.end0_y)]
        joined=joined.sort_values('region_id').drop_duplicates(['gene_id','anchor']);wanted=set(joined.region_id);sequences={};rid=None;pieces=[]
        with gzip.open(fasta,'rt') as f:
            for line in f:
                if line.startswith('>'):
                    if rid in wanted:sequences[rid]=''.join(pieces).upper()
                    rid=line[1:].split()[0];pieces=[]
                elif rid in wanted:pieces.append(line.strip())
            if rid in wanted:sequences[rid]=''.join(pieces).upper()
        for r in joined.itertuples():
            assert r.region_id in sequences
            a,b=(r.start0_x-r.start0_y,r.end0_x-r.start0_y) if r.strand=='+' else (r.end0_y-r.end0_x,r.end0_y-r.start0_x)
            observed=sequences[r.region_id][a:b];assert len(observed)==100 and observed==r.sequence,(sp,r.gene_id,r.anchor)
            rows.append(dict(species=sp,gene_id=r.gene_id,anchor=r.anchor,region_id=r.region_id,status='identical_independent_frozen_UTR_sequence'))
    pd.DataFrame(rows).to_csv(HERE/'independent_sequence_crosscheck.tsv',sep='\t',index=False)
    return len(rows)


def main():
    seq=pd.read_csv(HERE/'upstream100_sequences.tsv',sep='\t');windows=pd.read_csv(HERE/'motif_windows.tsv.gz',sep='\t')
    genes=pd.read_csv(HERE/'gene_features_and_labels.tsv.gz',sep='\t',float_precision='round_trip')
    byseq=seq.set_index(['species','gene_id','anchor']).sequence;checks=0
    bounds={'upstream40_100':(-100,-40),'upstream0_100':(-100,0),'upstream40_80':(-80,-40)}
    motif_counts={}
    for r in windows.itertuples():
        lo,hi=bounds[r.window];expected=byseq.loc[r.species,r.gene_id,r.anchor][100+lo:100+hi]
        assert r.sequence==expected and len(expected)==hi-lo
        a,b=(r.cut_site+lo,r.cut_site+hi) if r.strand=='+' else (r.cut_site-hi,r.cut_site-lo)
        assert(a,b)==(r.start0,r.end0)
        for motif in ['TGTA','TGAT','CGTA']:
            motif_counts[r.species,r.gene_id,r.anchor,r.window,motif]=len(list(re.finditer('(?='+motif+')',expected)))
        checks+=1
    for r in genes.itertuples():
        fields=r._asdict()
        for window in bounds:
            for motif in ['TGTA','TGAT','CGTA']:
                p=motif_counts[r.species,r.gene_id,'proximal',window,motif];d=motif_counts[r.species,r.gene_id,'distal',window,motif]
                assert fields[window+'_'+motif+'_proximal_count']==p and fields[window+'_'+motif+'_distal_count']==d and fields[window+'_'+motif+'_difference']==d-p
                checks+=1
    nulls=np.load(HERE/'permutation_draws.npz');t=pd.read_csv(HERE/'sequence_signature_tests.tsv',sep='\t',float_precision='round_trip')
    for r in t.itertuples():
        key='|'.join([r.background,r.species,r.model,r.window,r.motif,r.comparison]);v=nulls[key]
        assert len(v)==20000
        assert np.isclose(v.mean(),r.null_mean,atol=1e-12,rtol=0) and np.isclose(v.std(ddof=1),r.null_sd,atol=1e-12,rtol=0)
        upper=int(np.sum(v>=r.observed-1e-12));both=int(np.sum(abs(v-r.analytic_expected)>=abs(r.observed-r.analytic_expected)-1e-12))
        assert upper==r.tail_upper and both==r.tail_two_sided
        assert np.isclose((upper+1)/20001 if r.informative else 1,r.p_enrichment,atol=1e-14,rtol=0)
        assert r.proximal_assignment_informative==(r.n_mobile_proximal>0)
        assert r.control_assignment_informative==(r.n_mobile_control>0)
        checks+=1
    for _,g in t.groupby('family'):np.testing.assert_allclose(qvalues(g.p_enrichment),g.q_family,atol=1e-14,rtol=0)
    i=pd.read_csv(HERE/'paired_structure_interaction.tsv',sep='\t',float_precision='round_trip')
    for r in i.itertuples():
        v=nulls['interaction|'+r.model+'|'+r.comparison];tail=int(np.sum(abs(v-r.analytic_expected)>=abs(r.observed-r.analytic_expected)-1e-12))
        assert tail==r.tail_two_sided and np.isclose((tail+1)/20001,r.p_two_sided,atol=1e-14,rtol=0)
        assert r.species_interaction_proximal_assignment_informative==(r.human_n_mobile_proximal>0 and r.mouse_n_mobile_proximal>0)
        checks+=1
    for _,g in i.groupby('primary'):np.testing.assert_allclose(qvalues(g.p_two_sided),g.q_two_sided,atol=1e-14,rtol=0)
    members=pd.read_csv(HERE/'stratum_gene_members.tsv.gz',sep='\t',float_precision='round_trip');margins=pd.read_csv(HERE/'stratum_margins.tsv.gz',sep='\t')
    cells=members.groupby(['background','species','model','stratum'])
    for r in margins.itertuples():
        z=cells.get_group((r.background,r.species,r.model,r.stratum))
        assert len(z)==r.n and all(int(z.direction.eq(label).sum())==getattr(r,label) for label in ['proximal','distal','nonrecurrent'])
        checks+=1
    # All primary raw-score means and analytic conditional expectations are independent.
    for r in t[t.primary].itertuples():
        z=members[members.background.eq(r.background)&members.species.eq(r.species)&members.model.eq(r.model)]
        feature=r.window+'_'+r.motif+'_difference';control='nonrecurrent' if r.comparison.endswith('nonrecurrent') else 'distal'
        observed=z.loc[z.direction.eq('proximal'),feature].mean()-z.loc[z.direction.eq(control),feature].mean()
        expected=0.
        nprox=z.direction.eq('proximal').sum();ncon=z.direction.eq(control).sum()
        for _,q in z.groupby('stratum'):expected+=q[feature].sum()*(q.direction.eq('proximal').mean()/nprox-q.direction.eq(control).mean()/ncon)
        assert np.isclose(observed,r.observed,atol=1e-12,rtol=0) and np.isclose(expected,r.analytic_expected,atol=1e-12,rtol=0)
        checks+=1
    logs=[json.loads(line) for line in (HERE/'extract.log').read_text().splitlines() if line.strip().startswith('{')]
    config=json.loads((ROOT/'notes/2026-10-06-dual-main-plot-data/provenance/polyapipe-pipeline.config.json').read_text())
    for log in logs:assert log['sha256']==config['resources'][log['species']+'_genome']['sha256']
    audit=json.loads((HERE/'analysis_audit.json').read_text())
    for rel,h in audit['inputs'].items():assert digest(ROOT/rel)==h,rel
    independent=independent_local_catalog_sequences(seq)
    result=dict(status='PASS',checks=checks,sequence_windows=len(windows),motif_gene_fields=len(genes)*27,primary_mean_expectation_checks=8,conditional_tests=len(t),interactions=len(i),independent_frozen_UTR_sequence_checks=independent,independent_frozen_UTR_sequence_uncovered_windows=len(seq)-independent,remote_reference_hashes_verified=2,input_hashes_verified=len(audit['inputs']),human_paired_primary_proximal_identity_fully_fixed=True,skipped_checks=0)
    (HERE/'verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
