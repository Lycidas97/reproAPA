"""Frozen, conditional PAS architecture tests and fixed-probe U1 coverage summaries."""
from pathlib import Path
import bisect
import concurrent.futures
import gzip
import hashlib
import json
import re
import sys
import time

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PARENT = HERE.parent
ROSTER = ROOT / 'notes/2026-10-10-rna-module-enrichment/ortholog_roster.tsv.gz'
FEATURES = ROOT / 'notes/2026-10-10-cfim-sequence-signature/gene_features_and_labels.tsv.gz'
GTF = ROOT / 'notes/2026-10-03-polyapipe-rerun/statistics_package_v4r1/resources/human_ensembl_gtf__Homo_sapiens.GRCh38.116.gtf.gz'
FILES = {
    'AMO_control': ('GSM4173953', 'GSM4173953_cAMO_for_U1AMO.norm.bw'),
    'AMO_12_5': ('GSM4173954', 'GSM4173954_U1AMO_12.5_pmol.norm.bw'),
    'AMO_62_5': ('GSM4173955', 'GSM4173955_U1AMO_62.5_pmol.norm.bw'),
    'OE_control': ('GSM4173956', 'GSM4173956_hela_overexpression_control.norm.bw'),
    'OE_1': ('GSM4173957', 'GSM4173957_hela_overexpression_1ug.norm.bw'),
    'OE_1_5': ('GSM4173958', 'GSM4173958_hela_overexpression_1.5ug.norm.bw'),
}


def digest(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(d, name):
    d.to_csv(HERE / name, sep='\t', index=False, na_rep='NA', float_format='%.17g',
             compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)


def locked():
    p = json.loads((HERE / 'protocol.json').read_text())
    assert p['protocol_sha256'] == digest(HERE / 'PROTOCOL.zh.md')


def roster():
    z = pd.read_csv(ROSTER, sep='\t')
    z = z[z.contrast.eq('P_M')].copy()
    assert len(z) == 751 and z.human_gene_id.is_unique and int(z.foreground_002.sum()) == 27
    return z


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p, kind='stable')
    q = np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1]
    out = np.empty(len(p)); out[order] = q.clip(0, 1)
    return out


def seed(context):
    return int(hashlib.sha256(context.encode()).hexdigest()[:16], 16) % 2**32


def structure():
    locked(); r = roster()
    f = pd.read_csv(FEATURES, sep='\t', float_precision='round_trip')
    rows = []; memberships = []; draws = {}; coverage = []
    for species, col in [('human', 'human_gene_id'), ('mouse', 'mouse_gene_id')]:
        z = f[f.species.eq(species) & f.gene_id.isin(r[col])].copy()
        z = z.merge(r[[col, 'foreground_002', 'foreground_005']].rename(columns={col: 'gene_id'}), on='gene_id', validate='one_to_one')
        z = z.sort_values('gene_id').reset_index(drop=True)
        z['log_anchor_span'] = np.log(z.anchor_span)
        coverage.append(dict(species=species, n_background=len(z), n_core=int(z.foreground_002.sum()),
                             n_strict_core=int(z.foreground_005.sum()), n_paired=int(z.paired_eligible.sum()),
                             n_paired_core=int((z.paired_eligible & z.foreground_002).sum())))
        save(z, species + '_structure_inputs.tsv.gz')
        for model in ['assay_architecture', 'assay_expression']:
            zz = z.copy()
            columns = ['measured_pattern', 'proximal_capacity', 'distal_capacity']
            for name in ['PAS_number', 'baseline'] + (['mean_GEX', 'APA_coverage'] if model == 'assay_expression' else []):
                zz[name + '_bin'] = pd.qcut(zz[name], 2, labels=False, duplicates='drop').fillna(0).astype(int)
                columns.append(name + '_bin')
            groups = [np.asarray(v, int) for v in zz.groupby(columns, sort=True).indices.values()]
            labels = zz.foreground_002.to_numpy(np.int8)
            shuffled = np.broadcast_to(labels, (20000, len(labels))).copy()
            rng = np.random.default_rng(seed(species + '|' + model))
            mobile = np.zeros(len(z), bool)
            for j, ii in enumerate(groups):
                if len(np.unique(labels[ii])) > 1:
                    shuffled[:, ii] = rng.permuted(shuffled[:, ii], axis=1); mobile[ii] = True
                assert (shuffled[:, ii].sum(1) == labels[ii].sum()).all()
                for i in ii:
                    memberships.append(dict(species=species, model=model, stratum=j, gene_id=zz.loc[i, 'gene_id'],
                                            foreground=bool(labels[i]), mobile=bool(mobile[i])))
            n1 = int(labels.sum()); n0 = len(labels) - n1
            for feature in ['proximal_core_rank', 'core_rank_difference', 'log_anchor_span']:
                values = zz[feature].to_numpy(float)
                obs = float(values[labels == 1].mean() - values[labels == 0].mean())
                expected = sum(float(values[ii].sum()) * (labels[ii].mean() / n1 - (1 - labels[ii].mean()) / n0) for ii in groups)
                null = shuffled @ values / n1 - (1 - shuffled) @ values / n0
                sd = float(null.std(ddof=1)); informative = sd > 1e-12
                assert abs(null.mean() - expected) <= 6 * sd / np.sqrt(len(null)) + 1e-10
                tail = int((abs(null - expected) >= abs(obs - expected) - 1e-12).sum())
                context = species + '|' + model + '|' + feature
                draws[context] = null
                rows.append(dict(species=species, model=model, feature=feature, n_background=len(z), n_core=n1,
                                 core_mean=float(values[labels == 1].mean()), other_mean=float(values[labels == 0].mean()),
                                 observed_difference=obs, analytic_expected_difference=expected, excess=obs-expected,
                                 null_sd=sd, null_mean=float(null.mean()), two_sided_tail=tail,
                                 p_two_sided=(tail + 1) / 20001 if informative else 1., informative=informative,
                                 mobile_background=int(mobile.sum()), mobile_core=int((mobile & (labels == 1)).sum())))
    results = pd.DataFrame(rows); assert len(results) == 12
    results['q_12'] = bh(results.p_two_sided)
    save(results, 'structure_tests.tsv'); save(pd.DataFrame(coverage), 'structure_coverage.tsv')
    save(pd.DataFrame(memberships), 'structure_strata.tsv.gz')
    np.savez_compressed(HERE / 'structure_null_draws.npz', **draws)
    print(results.to_string(index=False), flush=True)


def terminal_probe(part):
    """CDS/UTR half-open intervals for the terminal exon of a fixed transcript."""
    ex = part[part.kind.eq('exon')]; cds = part[part.kind.eq('CDS')]
    if ex.empty or cds.empty:
        return None
    strand = str(part.strand.iloc[0])
    terminal = ex.sort_values('start0').iloc[-1 if strand == '+' else 0]
    lo, hi = int(terminal.start0), int(terminal.end0)
    blocks = sorted(set((max(lo, int(v.start0)), min(hi, int(v.end0))) for v in cds.itertuples()
                        if max(lo, int(v.start0)) < min(hi, int(v.end0))))
    if len(blocks) != 1:
        return None
    a, b = blocks[0]
    u0, u1 = (b, hi) if strand == '+' else (lo, a)
    if b - a < 50 or u1 - u0 < 50:
        return None
    return dict(chrom='chr' + str(terminal.chrom).removeprefix('chr'), strand=strand,
                cds_start0=a, cds_end0=b, utr_start0=u0, utr_end0=u1,
                CDS_length=int(sum(int(v.end0)-int(v.start0) for v in cds.itertuples())),
                last_CDS_length=b-a, last_UTR_length=u1-u0)


def prepare_probes():
    locked(); r = roster()
    e = pd.read_csv(PARENT / 'event_inventory.tsv.gz', sep='\t', dtype={'chrom': str})
    x = pd.read_csv(PARENT / 'decay/background_transcript_exon_CDS.tsv.gz', sep='\t', dtype={'chrom': str})
    x = x[x.species.eq('human')]
    transcript_groups = {key: part for key, part in x.groupby(['gene_id', 'transcript_id'], sort=False)}
    rows = []
    for row in r.itertuples():
        gene = row.human_gene_id; eg = e[e.species.eq('human') & e.gene_id.eq(gene)]
        base = dict(gene_id=gene, gene_symbol=row.human_gene_symbol, foreground_002=row.foreground_002,
                    foreground_005=row.foreground_005, n_routes=len(eg), status='pending')
        if eg.empty or not eg.driver_class.eq('strict_tandem_UTR').all():
            rows.append(base | {'status': 'driver_not_all_strict_tandem'}); continue
        common = set.intersection(*(set(str(v).split(';')) for v in eg.driver_compatible_transcripts))
        candidates = []
        for transcript in sorted(common):
            part = transcript_groups.get((gene, transcript))
            if part is None:
                continue
            probe = terminal_probe(part)
            if probe is not None:
                candidates.append(probe | {'transcript_id': transcript})
        if not candidates:
            rows.append(base | {'status': 'no_common_terminal_CDS_UTR_probe'}); continue
        pick = sorted(candidates, key=lambda v: (-v['CDS_length'], -v['last_UTR_length'], v['transcript_id']))[0]
        rows.append(base | pick | {'status': 'static_probe_pending_overlap', 'n_candidate_transcripts': len(candidates)})
    z = pd.DataFrame(rows)
    # Reject any selected unstranded probe overlapping an exon of another gene, genome-wide.
    query_by_chr = {}
    for idx, q in z[z.status.eq('static_probe_pending_overlap')].iterrows():
        for region in ['cds', 'utr']:
            query_by_chr.setdefault(str(q.chrom).removeprefix('chr'), []).append((int(q[region+'_start0']), int(q[region+'_end0']), idx, q.gene_id, region))
    indexes = {}
    for chrom, values in query_by_chr.items():
        values.sort(); indexes[chrom] = (values, [v[0] for v in values])
    bad = {}; count = 0
    pattern = re.compile(r'gene_id "([^"]+)"')
    with gzip.open(GTF, 'rt') as fi:
        for line in fi:
            if line.startswith('#'):
                continue
            fields = line.rstrip('\n').split('\t')
            if fields[2] != 'exon' or fields[0].removeprefix('chr') not in indexes:
                continue
            g = pattern.search(fields[8]).group(1); lo, hi = int(fields[3])-1, int(fields[4])
            values, starts = indexes[fields[0].removeprefix('chr')]
            for qlo, qhi, idx, own, region in values[:bisect.bisect_left(starts, hi)]:
                if qhi > lo and own != g:
                    bad.setdefault(idx, set()).add((g, region))
            count += 1
    z['other_gene_exon_overlap'] = ''
    for idx, values in bad.items():
        z.loc[idx, 'status'] = 'other_gene_exon_overlap'
        z.loc[idx, 'other_gene_exon_overlap'] = ';'.join(sorted(g+':'+region for g, region in values))
    z.loc[z.status.eq('static_probe_pending_overlap'), 'status'] = 'static_probe_eligible'
    save(z, 'external_probe_inventory.tsv.gz')
    inputs = [ROSTER, FEATURES, GTF, PARENT/'event_inventory.tsv.gz', PARENT/'event_contract.json',
              PARENT/'decay/background_transcript_exon_CDS.tsv.gz', HERE/'PROTOCOL.zh.md', HERE/'protocol.json']
    (HERE/'input_hashes.json').write_text(json.dumps({str(p.relative_to(ROOT)): digest(p) for p in inputs}, indent=2)+'\n')
    print(z.groupby('status').agg(n=('gene_id','size'), core=('foreground_002','sum')).to_string(), flush=True)


def integrate(intervals, lo, hi):
    return sum((min(b, hi) - max(a, lo)) * v for a, b, v in intervals if min(b, hi) > max(a, lo))


def coverage_track(track, info, shard=0, n_shards=1):
    import pyBigWig
    import urllib.request
    accession, filename = info
    url = 'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4173nnn/' + accession + '/suppl/' + filename
    p = pd.read_csv(HERE/'external_probe_inventory.tsv.gz', sep='\t')
    p = p[p.status.eq('static_probe_eligible')].sort_values(['chrom','cds_start0'])
    p = p.iloc[shard::n_shards]
    bw = pyBigWig.open(url); chroms = bw.chroms(); header = bw.header()
    remote = {'head_metadata_unavailable':'NCBI optional HEAD request returned HTTP503 in initial collection; file sizes archived in GEO filelist'}
    rows = []; intervals_out = []
    for q in p.itertuples():
        for region in ['cds', 'utr']:
            lo, hi = int(getattr(q, region+'_start0')), int(getattr(q, region+'_end0'))
            base = dict(track=track, gene_id=q.gene_id, transcript_id=q.transcript_id, region=region,
                        chrom=q.chrom, start0=lo, end0=hi, length=hi-lo)
            if q.chrom not in chroms or not 0 <= lo < hi <= chroms[q.chrom]:
                rows.append(base | dict(status='chrom_or_bounds_missing', integral=np.nan, mean_coverage=np.nan, n_intervals=0)); continue
            try:
                items = bw.intervals(q.chrom, lo, hi) or ()
                value = integrate(items, lo, hi)
                assert np.isfinite(value) and value >= 0
                rows.append(base | dict(status='measured', integral=value, mean_coverage=value/(hi-lo), n_intervals=len(items)))
                intervals_out.extend(dict(track=track, gene_id=q.gene_id, region=region, start0=a, end0=b, value=v) for a,b,v in items)
            except Exception as exc:
                rows.append(base | dict(status='query_failed:'+str(exc), integral=np.nan, mean_coverage=np.nan, n_intervals=0))
        if len(rows) % 40 == 0:
            print(track, shard, len(rows), 'queries checkpoint', flush=True)
    bw.close()
    suffix = track + ('_shard' + str(shard) if n_shards > 1 else '')
    save(pd.DataFrame(rows), 'coverage_' + suffix + '.tsv.gz')
    save(pd.DataFrame(intervals_out), 'intervals_' + suffix + '.tsv.gz')
    (HERE/('source_'+suffix+'.json')).write_text(json.dumps(dict(url=url, accession=accession, header=header, chroms=chroms,
        http_headers=remote, complete_file_sha256=None, note='HTTPS range queries; source interval data archived and hashed, full 2.56GB source not downloaded'),indent=2)+'\n')
    print(track, shard, len(rows), len(intervals_out), flush=True)


def external_collect():
    locked()
    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as pool:
        futures = [pool.submit(coverage_track, track, info, shard, 3) for track, info in FILES.items() for shard in range(3)]
        for future in concurrent.futures.as_completed(futures):
            future.result()
    for track in FILES:
        for prefix in ['coverage_', 'intervals_']:
            parts = [pd.read_csv(HERE/(prefix+track+'_shard'+str(j)+'.tsv.gz'),sep='\t',float_precision='round_trip') for j in range(3)]
            joined = pd.concat(parts,ignore_index=True)
            if prefix == 'coverage_':
                assert len(joined)==444 and not joined.duplicated(['gene_id','region']).any()
            save(joined,prefix+track+'.tsv.gz')
        meta=json.loads((HERE/('source_'+track+'_shard0.json')).read_text())
        meta['collection_shards']=3
        (HERE/('source_'+track+'.json')).write_text(json.dumps(meta,indent=2)+'\n')


def external_summary():
    locked(); r = roster()
    inv = pd.read_csv(HERE/'external_probe_inventory.tsv.gz', sep='\t')
    coverage = pd.concat([pd.read_csv(HERE/('coverage_'+track+'.tsv.gz'), sep='\t', float_precision='round_trip') for track in FILES], ignore_index=True)
    assert not coverage.duplicated(['track','gene_id','region']).any()
    cov = coverage.pivot(index='gene_id', columns=['track','region'], values='mean_coverage')
    z = inv[inv.status.eq('static_probe_eligible')].copy().set_index('gene_id')
    query_states = coverage.assign(ok=coverage.status.eq('measured')).groupby('gene_id').ok.all()
    z['complete_all_queries'] = query_states
    for track in FILES:
        z[track+'_CDS_coverage'] = cov[track, 'cds']
        z[track+'_UTR_coverage'] = cov[track, 'utr']
        z[track+'_LECDS'] = cov[track,'utr']/cov[track,'cds']
    for track, control in [('AMO_12_5','AMO_control'),('AMO_62_5','AMO_control'),('OE_1','OE_control'),('OE_1_5','OE_control')]:
        z[track+'_log2_fold'] = np.log2(z[track+'_LECDS']/z[control+'_LECDS'])
    effects = ['AMO_12_5_log2_fold','AMO_62_5_log2_fold','OE_1_log2_fold','OE_1_5_log2_fold']
    z['finite_all_effects'] = np.isfinite(z[effects]).all(axis=1)
    z['AMO_mean_log2_fold'] = z[effects[:2]].mean(axis=1)
    z['OE_mean_log2_fold'] = z[effects[2:]].mean(axis=1)
    z['reciprocal_score'] = np.minimum(-z.AMO_mean_log2_fold, z.OE_mean_log2_fold)
    controlcols = [s+'_'+region+'_coverage' for s in ['AMO_control','OE_control'] for region in ['CDS','UTR']]
    z['minimum_control_coverage'] = z[controlcols].min(axis=1)
    save(z.reset_index(), 'external_gene_effects.tsv.gz')
    rows=[]; stats=[]; attrition=[]
    for gate in [.01,.05,.1]:
        eligible = z[z.minimum_control_coverage.ge(gate) & z.finite_all_effects & z.complete_all_queries]
        attrition.append(dict(control_coverage_gate=gate, static_probes=len(z), static_core=int(z.foreground_002.sum()),
             n_background=len(eligible), n_core=int(eligible.foreground_002.sum()), n_strict_core=int(eligible.foreground_005.sum()),
             coverage_failed=int(z.minimum_control_coverage.lt(gate).sum()), nonfinite_effect=int((~z.finite_all_effects).sum()),query_failed=int((~z.complete_all_queries).sum())))
        for group in ['core002','core005','other','all_background']:
            part = eligible[eligible.foreground_002] if group=='core002' else eligible[eligible.foreground_005] if group=='core005' else eligible[~eligible.foreground_002] if group=='other' else eligible
            for feat in effects+['AMO_mean_log2_fold','OE_mean_log2_fold','reciprocal_score']:
                stats.append(dict(control_coverage_gate=gate,group=group,feature=feat,n=len(part),mean=part[feat].mean(),median=part[feat].median()))
            for threshold in [0.,.1,.2,.5]:
                amo = part[effects[:2]].le(-threshold).all(axis=1)
                oe = part[effects[2:]].ge(threshold).all(axis=1)
                if threshold==0:
                    amo = part[effects[:2]].lt(0).all(axis=1);oe = part[effects[2:]].gt(0).all(axis=1)
                reverse = part[effects[:2]].ge(threshold).all(axis=1) & part[effects[2:]].le(-threshold).all(axis=1)
                if threshold==0:
                    reverse = part[effects[:2]].gt(0).all(axis=1) & part[effects[2:]].lt(0).all(axis=1)
                rows.append(dict(control_coverage_gate=gate,group=group,log2_threshold=threshold,n=len(part),
                                 AMO_both_proximal=int(amo.sum()),OE_both_distal=int(oe.sum()),reciprocal_expected_direction=int((amo&oe).sum()),
                                 reciprocal_reverse_direction=int(reverse.sum())))
    save(pd.DataFrame(rows),'external_direction_counts.tsv');save(pd.DataFrame(stats),'external_continuous_summary.tsv')
    save(pd.DataFrame(attrition),'external_coverage_attrition.tsv')
    allgene = inv.merge(z[['minimum_control_coverage','finite_all_effects','complete_all_queries']].reset_index(),on='gene_id',how='left',validate='one_to_one')
    allgene['primary_external_status'] = allgene.status
    allgene.loc[allgene.status.eq('static_probe_eligible') & allgene.minimum_control_coverage.lt(.05),'primary_external_status'] = 'control_coverage_below_0.05'
    finiteflag = allgene.finite_all_effects.fillna(False).astype(bool)
    queryflag = allgene.complete_all_queries.fillna(False).astype(bool)
    allgene.loc[allgene.status.eq('static_probe_eligible') & allgene.minimum_control_coverage.ge(.05) & ~finiteflag,'primary_external_status'] = 'treatment_or_control_nonfinite_ratio'
    allgene.loc[allgene.status.eq('static_probe_eligible') & allgene.minimum_control_coverage.ge(.05) & finiteflag,'primary_external_status'] = 'primary_external_eligible'
    allgene.loc[allgene.status.eq('static_probe_eligible') & ~queryflag,'primary_external_status'] = 'query_failed_or_missing'
    save(allgene,'external_all_751_coverage.tsv.gz')
    save(allgene[allgene.foreground_002],'external_all_27_coverage.tsv')
    print(pd.DataFrame(attrition).to_string(index=False),flush=True)
    print(pd.DataFrame(rows).query('control_coverage_gate==.05').to_string(index=False),flush=True)


if __name__ == '__main__':
    {'structure':structure, 'prepare':prepare_probes, 'collect':external_collect, 'summarize':external_summary}[sys.argv[1]]()
