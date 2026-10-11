"""Bind a verified integrated PAS reference to unchanged native count axes.

Reference BED positions are retained bases, not interbase cuts. Gene assignments
and sequence QC are reused because native coordinates and their genomes are fixed.
Every external support field is recalculated; old effects are never copied.
"""
import argparse
import json
from pathlib import Path
import subprocess
import time

import numpy as np
import pandas as pd

from polyapipe_prepare import digest, save

SOURCES = {'PolyASite': 'polyasite_broad', 'PolyA_DB4.1Main': 'polyadb_main',
           'GENCODE_polyA_site': 'gencode_polya'}
BED = ['chrom', 'start0', 'end0', 'pas_id', 'score', 'strand']
VARIANTS = ['v2_primary', 'v3_ge62_comparison']


def cuts(start, end, strand):
    start, end, strand = np.asarray(start), np.asarray(end), np.asarray(strand)
    assert np.all(end-start == 1) and np.all(start >= 0)
    assert np.isin(strand, ['+', '-']).all()
    return np.where(strand == '+', end, start).astype(np.int64)


def canonical(chrom):
    c = str(chrom)
    c = c[3:] if c.startswith('chr') else c
    return 'MT' if c in ['M', 'MT'] else c


def accept_mappings(source, mapped, fai):
    """Reject multiple mappings before any contig filtering or deduplication."""
    assert source.pas_id.is_unique
    assert set(mapped.pas_id) <= set(source.pas_id)
    duplicate = set(mapped.loc[mapped.pas_id.duplicated(False), 'pas_id'])
    bad_length = set(mapped.loc[mapped.end0.sub(mapped.start0).ne(1), 'pas_id'])
    aliases = {canonical(c): c for c in fai}
    assert len(aliases) == len(fai), 'Ambiguous genome contig aliases'
    mapped = mapped.copy()
    mapped['target_chrom'] = mapped.chrom.map(lambda c: aliases.get(canonical(c)))
    bad_contig = set(mapped.loc[mapped.target_chrom.isna(), 'pas_id'])
    outside = set(mapped.loc[mapped.end0.gt(mapped.target_chrom.map(fai)) | mapped.start0.lt(0), 'pas_id'])
    accepted = mapped[~mapped.pas_id.isin(duplicate | bad_length | bad_contig | outside)].copy()
    accepted['target_strand'] = accepted.strand
    accepted['target_cut_site'] = cuts(accepted.start0, accepted.end0, accepted.strand)
    result = source.merge(accepted[['pas_id', 'target_chrom', 'target_strand', 'target_cut_site']],
                          on='pas_id', how='left', validate='one_to_one')
    result['mapping_status'] = np.select([
        result.pas_id.isin(duplicate), result.pas_id.isin(bad_length), result.pas_id.isin(bad_contig),
        result.pas_id.isin(outside), result.target_cut_site.notna()],
        ['ambiguous', 'not_one_base', 'absent_target_contig', 'outside_target_contig', 'unique'], default='unmapped')
    return result


def nearest(frame, reference):
    distances = np.full(len(frame), -1, np.int64)
    positions = np.full(len(frame), -1, np.int64)
    groups = {k: np.unique(g.target_cut_site.to_numpy(np.int64))
              for k, g in reference.groupby(['target_chrom', 'target_strand'])}
    for key, ids in frame.groupby(['chrom', 'strand']).groups.items():
        values = groups.get(key)
        if values is None or not len(values):
            continue
        query = frame.loc[ids, 'cut_site'].to_numpy(np.int64)
        right = np.minimum(np.searchsorted(values, query), len(values)-1)
        left = np.maximum(right-1, 0)
        ix = np.where(abs(query-values[left]) <= abs(query-values[right]), left, right)
        positions[ids] = values[ix]
        distances[ids] = abs(query-values[ix])
    return distances, positions


def support_fields(frame, references):
    frame = frame.reset_index(drop=True).copy()
    for source in SOURCES.values():
        distance, position = nearest(frame, references[references.source_family.eq(source)])
        frame[source+'_distance_nt'] = distance
        frame[source+'_nearest_cut'] = position
    for window in [25, 50, 100]:
        frame[f'matched_source_count_{window}nt'] = sum(
            frame[s+'_distance_nt'].between(0, window).astype(int) for s in SOURCES.values())
    return frame


def source_records(good):
    records = good.assign(source=good.supporting_sources.str.split(',')).explode('source')
    assert records.source.isin(SOURCES).all(), records.loc[~records.source.isin(SOURCES), 'source'].unique()
    records['source_family'] = records.source.map(SOURCES)
    records['primary_contig'] = records.target_chrom.map(canonical).isin([str(i) for i in range(1, 23)]+['X', 'Y', 'MT'])
    return records.drop_duplicates(['source_family', 'target_chrom', 'target_strand', 'target_cut_site'])


def build_reference(data, variant, species, out, resources, inputs):
    folder = data/species/variant
    raw = pd.read_csv(folder/'integrated_pas.bed.gz', sep='\t', names=BED)
    metadata = pd.read_csv(folder/'site_metadata.tsv.gz', sep='\t')
    assert raw.pas_id.is_unique and metadata.pas_id.is_unique
    assert set(raw.pas_id) == set(metadata.pas_id)
    assert np.array_equal(cuts(raw.start0, raw.end0, raw.strand),
                          raw.start0+raw.strand.eq('+').astype(int))
    work = out/'reference'/species
    work.mkdir(parents=True)
    source = work/'input.bed'
    raw.to_csv(source, sep='\t', index=False, header=False)
    if species == 'mouse':
        mapped, unmapped = work/'mapped.bed', work/'unmapped.bed'
        command = [resources['liftover']['path'], '-multiple', '-minMatch=1', str(source),
                   resources['chain_mm10_mm39']['path'], str(mapped), str(unmapped)]
        with (work/'liftover.log').open('w') as log:
            subprocess.run(command, stdout=log, stderr=log, check=True)
        mapped_frame = pd.read_csv(mapped, sep='\t', names=BED) if mapped.stat().st_size else pd.DataFrame(columns=BED)
    else:
        mapped_frame = raw.copy()
        command = ['assembly unchanged; actual target FASTA contig aliases verified']
    fai_path = Path(resources[species+'_genome']['path']+'.fai')
    inputs[str(fai_path)] = digest(fai_path)
    fai = pd.read_csv(fai_path, sep='\t', header=None, usecols=[0, 1])
    mapping = accept_mappings(raw, mapped_frame, dict(zip(fai[0], fai[1])))
    mapping = mapping.merge(metadata, on='pas_id', validate='one_to_one')
    mapping['species'] = species
    mapping['target_build'] = 'GRCh38_2024_A' if species == 'human' else 'GRCm39_2024_A'
    mapping['reference_variant'] = variant
    good = mapping[mapping.mapping_status.eq('unique')].copy()
    good['target_cut_site'] = good.target_cut_site.astype(np.int64)
    good['catalog_pas_id'] = ['integrated:'+species+':'+b+':'+c+':'+s+':'+str(p)+':'+variant
        for b, c, s, p in zip(good.target_build, good.target_chrom, good.target_strand, good.target_cut_site)]
    mapping = mapping.merge(good[['pas_id', 'catalog_pas_id']], on='pas_id', how='left', validate='one_to_one')
    save(mapping, work/'mapping_audit.tsv.gz')
    # If two distinct input representatives lift to the same target cut, they have
    # one target identity. Preserve both originals in mapping_audit, never double-count.
    records = source_records(good)
    save(records, work/'representative_support_sources.tsv.gz')
    summary = dict(species=species, variant=variant, input_sites=len(raw),
                   mapping_status=mapping.mapping_status.value_counts().to_dict(),
                   unique_target_sites=good.catalog_pas_id.nunique(),
                   target_collisions=len(good)-good.catalog_pas_id.nunique(), command=command)
    (work/'mapping_summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return records, summary


def carry_membership_audit(native,out):
    """A reference change inherits cell/count evidence only after exact identity checks."""
    import shutil
    previous=json.loads((native/'intake.native.json').read_text())
    current=json.loads((out/'intake.native.json').read_text())
    assert current['previous_native_intake_sha256']==digest(native/'intake.native.json')
    assert current['studies']==previous['studies']
    evidence=json.loads((native/'frozen_native_membership_audit.json').read_text())
    assert evidence['status']=='PASS' and all(evidence['checks'].values())
    for field in ['cells','metacells','libraries']:
        assert evidence[field]==sum(row[field] for row in current['studies'])
    annotation=json.loads((out/'annotation.complete.json').read_text());checked=[]
    for row in annotation['studies']:
        study=row['study_id']
        assert {'metacells.tsv.gz','metacell_by_pas_counts.npz'}<=set(row['matrix_and_metadata_reused'])
        for name,sha in row['matrix_and_metadata_reused'].items():
            a=native/'studies'/study/name;b=out/'studies'/study/name
            assert digest(a)==sha and digest(b)==sha,(study,name)
            checked.append(dict(study_id=study,file=name,sha256=sha))
    sources={}
    for name in ['frozen_native_membership_audit.json','frozen_native_membership_audit.tsv']:
        source=native/name;sha=digest(source);target=out/name
        assert not target.exists(),target
        shutil.copy2(source,target);assert digest(target)==sha;sources[str(source)]=sha
    (out/'reference_membership_audit.json').write_text(json.dumps(dict(status='PASS',
        source_audits=sources,identical_count_and_metadata_files=checked,reference_only_change=True),indent=2)+'\n')


def main(args):
    args.out.mkdir(parents=True, exist_ok=False)
    data = args.data.resolve()
    checks = json.loads((data/'output_sha256.json').read_text())
    inputs = {str(data/p): h for p, h in checks.items()}
    inputs[str(data/'manifest.json')] = digest(data/'manifest.json')
    inputs[str(data/'output_sha256.json')] = digest(data/'output_sha256.json')
    for p, h in inputs.items():
        assert digest(p) == h, f'Reference changed: {p}'
    config = json.loads(args.resource_config.read_text())
    resources = {k: config['resources'][k] for k in ['liftover', 'chain_mm10_mm39', 'human_genome', 'mouse_genome']}
    for item in resources.values():
        assert digest(item['path']) == item['sha256']
        inputs[item['path']] = item['sha256']
    inputs[str(Path(__file__).resolve())] = digest(__file__)
    inputs[str(args.native/'annotation.complete.json')] = digest(args.native/'annotation.complete.json')
    old_annotation = json.loads((args.native/'annotation.complete.json').read_text())
    frozen = json.loads((args.native/'intake.native.json').read_text())
    inputs[str(args.native/'intake.native.json')] = digest(args.native/'intake.native.json')
    original_features = {p.name: digest(p/'features.tsv.gz') for p in (args.native/'studies').iterdir() if p.is_dir()}
    assert len(original_features) == 7
    for variant in VARIANTS:
        out = args.out/variant
        out.mkdir()
        refs, summaries = [], []
        for species in ['human', 'mouse']:
            reference, summary = build_reference(data, variant, species, out, resources, inputs)
            refs.append(reference); summaries.append(summary)
        refs = pd.concat(refs, ignore_index=True)
        save(refs, out/'reference/representative_support_sources.tsv.gz')
        rows = []
        for source in sorted((args.native/'studies').iterdir()):
            if not source.is_dir():
                continue
            dest = out/'studies'/source.name
            dest.mkdir(parents=True)
            original = pd.read_csv(source/'features.tsv.gz', sep='\t', low_memory=False)
            assert digest(source/'features.tsv.gz') == original_features[source.name]
            inputs[str(source/'features.tsv.gz')] = original_features[source.name]
            sp = original.species.iloc[0]
            updated = support_fields(original, refs[refs.species.eq(sp)])
            support_columns = [s+suffix for s in SOURCES.values() for suffix in ['_distance_nt', '_nearest_cut']]+[
                f'matched_source_count_{w}nt' for w in [25, 50, 100]]
            pd.testing.assert_frame_equal(original.drop(columns=support_columns), updated.drop(columns=support_columns))
            updated['support_reference_variant'] = variant
            updated['support_distance_basis'] = 'integrated_representative_interbase'
            save(updated, dest/'features.tsv.gz')
            reused = {}
            for path in sorted(source.iterdir()):
                if not path.is_file() or path.name in ['features.tsv.gz', 'effects.log']:
                    continue
                assert path.suffix in ['.gz', '.npz', '.tsv'], path
                (dest/path.name).symlink_to(path.resolve())
                reused[path.name] = digest(path)
            old = original.matched_source_count_50nt.ge(1)
            new = updated.matched_source_count_50nt.ge(1)
            rows.append(dict(study_id=source.name, species=sp, sites=len(updated), old_supported50=int(old.sum()),
                new_supported50=int(new.sum()), gained=int((~old & new).sum()), lost=int((old & ~new).sum()),
                unchanged_primary_eligibility=True, matrix_and_metadata_reused=reused,
                features_sha256=digest(dest/'features.tsv.gz')))
            print(json.dumps(dict(phase='support_rebound', variant=variant, study=source.name,
                                 old_support=int(old.sum()), new_support=int(new.sum()))), flush=True)
        manifest = dict(frozen, reference_variant=variant, run_id='polyapipe_integrated_20261003_'+variant,
                        support_reference='reference/representative_support_sources.tsv.gz',
                        support_coordinate_system='0-based-interbase', reference_inputs=inputs.copy(),
                        reference_policy='same-strand integrated representatives at inclusive 25/50/100 nt; full catalogue',
                        previous_native_intake_sha256=inputs[str(args.native/'intake.native.json')])
        (out/'intake.native.json').write_text(json.dumps(manifest, indent=2)+'\n')
        annotation = dict(status='complete', rows=sum(r['sites'] for r in rows), studies=rows,
                          reference_summaries=summaries, structural_annotation_reused_from=str(args.native),
                          structural_annotation_receipt_sha256=inputs[str(args.native/'annotation.complete.json')],
                          external_support_recomputed=True, bam_recount=False)
        assert annotation['rows'] == old_annotation['rows']
        (out/'annotation.complete.json').write_text(json.dumps(annotation, indent=2)+'\n')
        carry_membership_audit(args.native,out)
        (out/'reference_audit.json').write_text(json.dumps(dict(inputs=inputs, summaries=summaries,
            policy=manifest['reference_policy'], source_columns='compatibility names; PolyASite version specified by reference_variant'), indent=2)+'\n')
    assert all(digest(p) == h for p, h in inputs.items()), 'Input changed during reference replacement'
    (args.out/'reference-replacement.complete.json').write_text(json.dumps(dict(status='reference_rebound_core_and_statistics_pending',
        time_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), inputs=inputs, variants=VARIANTS,
        no_BAM_recount=True, original_results_preserved=True), indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--native', type=Path, required=True)
    parser.add_argument('--resource-config', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    main(parser.parse_args())
