"""Independent numerical and eligibility checks for the frozen U1 candidate analysis."""
from pathlib import Path
import gzip
import hashlib
import json
import math

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
checks = []


def check(name, value):
    assert bool(value), name
    checks.append(name)


def sha(p):
    with Path(p).open('rb') as fi:
        return hashlib.file_digest(fi,'sha256').hexdigest()


def check_structure():
    tests = pd.read_csv(HERE/'structure_tests.tsv',sep='\t',float_precision='round_trip')
    strata = pd.read_csv(HERE/'structure_strata.tsv.gz',sep='\t')
    draws = np.load(HERE/'structure_null_draws.npz')
    check('12 frozen structure tests',len(tests)==12)
    for row in tests.itertuples():
        z = pd.read_csv(HERE/(row.species+'_structure_inputs.tsv.gz'),sep='\t',float_precision='round_trip')
        m = strata[(strata.species==row.species)&(strata.model==row.model)]
        d = z.merge(m[['gene_id','stratum','mobile']],on='gene_id',validate='one_to_one')
        n1=int(d.foreground_002.sum());n0=len(d)-n1
        vals=d[row.feature]
        observed=vals[d.foreground_002].sum()/n1-vals[~d.foreground_002].sum()/n0
        expected=0.;variance=0.
        for _,part in d.groupby('stratum'):
            N=len(part);n=int(part.foreground_002.sum());value=part[row.feature].to_numpy(float)
            expected += value.sum()*(n/N/n1-(1-n/N)/n0)
            if N>1:
                variance += (1/n1+1/n0)**2*n*(N-n)/N*float(value.var(ddof=1))
        key='|'.join([row.species,row.model,row.feature]);null=draws[key]
        tail=sum(abs(float(v)-expected)>=abs(observed-expected)-1e-12 for v in null)
        p=(tail+1)/(len(null)+1) if np.std(null,ddof=1)>1e-12 else 1.
        check(key+': observed',abs(observed-row.observed_difference)<1e-12)
        check(key+': analytic expectation',abs(expected-row.analytic_expected_difference)<1e-12)
        check(key+': Monte Carlo expectation',abs(null.mean()-expected)<=6*math.sqrt(variance)/math.sqrt(len(null))+1e-10)
        check(key+': Monte Carlo variance',abs(null.var(ddof=1)-variance)<=.08*max(variance,1e-12))
        check(key+': two-tail p',abs(p-row.p_two_sided)<1e-15)
        check(key+': null length',len(null)==20000)
    p=tests.p_two_sided.to_numpy();q=np.empty(len(p))
    # Independent BH using the minimum adjusted p among all thresholds >= this p.
    order=np.argsort(p);rank=np.empty(len(p),int);rank[order]=np.arange(1,len(p)+1)
    for i in range(len(p)):
        q[i]=min(1.,min(p[j]*len(p)/rank[j] for j in range(len(p)) if rank[j]>=rank[i]))
    check('independent BH twelve family',np.max(abs(q-tests.q_12))<1e-15)
    sequences=pd.read_csv(ROOT/'notes/2026-10-10-cfim-sequence-signature/upstream100_sequences.tsv',sep='\t')
    sequences=sequences.set_index(['species','gene_id','anchor']).sequence.to_dict()
    ranks={'AATAAA':3,'ATTAAA':2,'AGTAAA':2,'AATACA':2,'AATATA':2,'CATAAA':2,'TATAAA':2,
           'AAGAAA':1,'AATAAT':1,'AACAAA':1,'AATAGA':1,'AATGAA':1,'GATAAA':1}
    for species in ['human','mouse']:
        z=pd.read_csv(HERE/(species+'_structure_inputs.tsv.gz'),sep='\t',float_precision='round_trip')
        for row in z.itertuples():
            computed=[]
            for anchor in ['proximal','distal']:
                text=sequences[species,row.gene_id,anchor][65:91]
                hits=[ranks.get(text[i:i+6],0) for i in range(len(text)-5)]
                value=max(hits+[0]);computed.append(value)
                check(species+':'+row.gene_id+':'+anchor+': core rank independently from sequence',value==getattr(row,anchor+'_core_rank'))
            check(species+':'+row.gene_id+': rank difference',computed[1]-computed[0]==row.core_rank_difference)
            check(species+':'+row.gene_id+': native strand distance',abs(row.distal_cut-row.proximal_cut)==row.anchor_span)


def check_probes():
    inv=pd.read_csv(HERE/'external_probe_inventory.tsv.gz',sep='\t',float_precision='round_trip')
    x=pd.read_csv(HERE.parent/'decay/background_transcript_exon_CDS.tsv.gz',sep='\t',dtype={'chrom':str})
    e=pd.read_csv(HERE.parent/'event_inventory.tsv.gz',sep='\t',dtype={'chrom':str})
    eligible=inv[inv.status.eq('static_probe_eligible')]
    check('fixed 751 probe inventory',len(inv)==751 and inv.gene_id.is_unique)
    check('fixed 27 core inventory',int(inv.foreground_002.sum())==27)
    check('222 static probes 11 core',len(eligible)==222 and int(eligible.foreground_002.sum())==11)
    for row in eligible.itertuples():
        part=x[x.species.eq('human')&x.gene_id.eq(row.gene_id)&x.transcript_id.eq(row.transcript_id)]
        exon=part[part.kind.eq('exon')].sort_values('start0').iloc[-1 if row.strand=='+' else 0]
        cds=part[part.kind.eq('CDS')]
        a=max(int(exon.start0),int(cds[cds.end0.gt(exon.start0)&cds.start0.lt(exon.end0)].start0.min()))
        b=min(int(exon.end0),int(cds[cds.end0.gt(exon.start0)&cds.start0.lt(exon.end0)].end0.max()))
        expected=(b,int(exon.end0)) if row.strand=='+' else (int(exon.start0),a)
        check(row.gene_id+': exact terminal CDS/UTR',a==row.cds_start0 and b==row.cds_end0 and expected==(row.utr_start0,row.utr_end0))
        eg=e[e.species.eq('human')&e.gene_id.eq(row.gene_id)]
        check(row.gene_id+': driver all strict',eg.driver_class.eq('strict_tandem_UTR').all())
        check(row.gene_id+': transcript compatible across routes',all(row.transcript_id in str(t).split(';') for t in eg.driver_compatible_transcripts))
    check('plus UTR downstream intent',eligible[eligible.strand.eq('+')].utr_start0.eq(eligible[eligible.strand.eq('+')].cds_end0).all())
    check('minus UTR upstream intent',eligible[eligible.strand.eq('-')].utr_end0.eq(eligible[eligible.strand.eq('-')].cds_start0).all())


def check_coverage():
    tracks=['AMO_control','AMO_12_5','AMO_62_5','OE_control','OE_1','OE_1_5']
    coverage_files=[HERE/('coverage_'+track+'.tsv.gz') for track in tracks]
    if not all(p.exists() for p in coverage_files):
        return 'awaiting_six_tracks'
    for p in coverage_files:
        cov=pd.read_csv(p,sep='\t',float_precision='round_trip')
        items=pd.read_csv(HERE/p.name.replace('coverage_','intervals_'),sep='\t',float_precision='round_trip')
        check(p.name+': 444 queries',len(cov)==444 and not cov.duplicated(['gene_id','region']).any())
        parts={key:g for key,g in items.groupby(['gene_id','region'],sort=False)}
        for row in cov[cov.status.eq('measured')].itertuples():
            g=parts.get((row.gene_id,row.region))
            total=0.
            if g is not None:
                # Weighted numpy interval clipping differs from production scalar integration.
                lengths=np.minimum(g.end0.to_numpy(),row.end0)-np.maximum(g.start0.to_numpy(),row.start0)
                total=float(np.maximum(lengths,0)@g.value.to_numpy())
            check(p.name+':'+row.gene_id+':'+row.region+': exact integral',abs(total-row.integral)<=1e-8*max(1,abs(total)))
            check(p.name+':'+row.gene_id+':'+row.region+': gaps zero length denominator',abs(total/(row.end0-row.start0)-row.mean_coverage)<=1e-10*max(1,abs(row.mean_coverage)))
    if not (HERE/'external_gene_effects.tsv.gz').exists():
        return 'awaiting_summary'
    z=pd.read_csv(HERE/'external_gene_effects.tsv.gz',sep='\t',float_precision='round_trip')
    cols=['AMO_12_5_log2_fold','AMO_62_5_log2_fold','OE_1_log2_fold','OE_1_5_log2_fold']
    for tr,ctrl in [('AMO_12_5','AMO_control'),('AMO_62_5','AMO_control'),('OE_1','OE_control'),('OE_1_5','OE_control')]:
        exp=np.log2(z[tr+'_UTR_coverage']/z[tr+'_CDS_coverage']/z[ctrl+'_UTR_coverage']*z[ctrl+'_CDS_coverage'])
        finite=np.isfinite(exp)&np.isfinite(z[tr+'_log2_fold'])
        check(tr+': independent log-ratio',np.max(abs(exp[finite]-z.loc[finite,tr+'_log2_fold']))<1e-12)
    counts=pd.read_csv(HERE/'external_direction_counts.tsv',sep='\t')
    for row in counts.itertuples():
        q=z[(z.minimum_control_coverage>=row.control_coverage_gate)&np.isfinite(z[cols]).all(axis=1)&z.complete_all_queries]
        q=q[q.foreground_002] if row.group=='core002' else q[q.foreground_005] if row.group=='core005' else q[~q.foreground_002] if row.group=='other' else q
        t=row.log2_threshold
        a=((q[cols[0]]<-t)&(q[cols[1]]<-t)) if t==0 else ((q[cols[0]]<=-t)&(q[cols[1]]<=-t))
        b=((q[cols[2]]>t)&(q[cols[3]]>t)) if t==0 else ((q[cols[2]]>=t)&(q[cols[3]]>=t))
        reverse=((q[cols[0]]>t)&(q[cols[1]]>t)&(q[cols[2]]<-t)&(q[cols[3]]<-t)) if t==0 else ((q[cols[0]]>=t)&(q[cols[1]]>=t)&(q[cols[2]]<=-t)&(q[cols[3]]<=-t))
        check(str((row.control_coverage_gate,row.group,t))+': count orientation',len(q)==row.n and int(a.sum())==row.AMO_both_proximal and int(b.sum())==row.OE_both_distal and int((a&b).sum())==row.reciprocal_expected_direction and int(reverse.sum())==row.reciprocal_reverse_direction)
    allgene=pd.read_csv(HERE/'external_all_751_coverage.tsv.gz',sep='\t')
    for row in allgene.itertuples():
        expected=row.status
        if row.status=='static_probe_eligible':
            if not bool(row.complete_all_queries):
                expected='query_failed_or_missing'
            elif row.minimum_control_coverage < .05:
                expected='control_coverage_below_0.05'
            elif not bool(row.finite_all_effects):
                expected='treatment_or_control_nonfinite_ratio'
            else:
                expected='primary_external_eligible'
        check(row.gene_id+': independent explicit missing status',row.primary_external_status==expected)
    core=pd.read_csv(HERE/'external_all_27_coverage.tsv',sep='\t')
    check('all27 attrition retains all core',len(core)==27 and core.foreground_002.all())
    check('all27 primary coverage agrees main counts',core.primary_external_status.eq('primary_external_eligible').sum()==10)
    return 'complete'


if __name__=='__main__':
    p=json.loads((HERE/'protocol.json').read_text())
    check('protocol unchanged',sha(HERE/'PROTOCOL.zh.md')==p['protocol_sha256'])
    for path,expected in json.loads((HERE/'input_hashes.json').read_text()).items():
        check(path+': source unchanged',sha(ROOT/path)==expected)
    check_structure();check_probes();status=check_coverage()
    (HERE/'verification.json').write_text(json.dumps(dict(status=status,n_checks=len(checks),failed=0,checks=checks),indent=2)+'\n')
    print(json.dumps(dict(status=status,n_checks=len(checks),failed=0)))
