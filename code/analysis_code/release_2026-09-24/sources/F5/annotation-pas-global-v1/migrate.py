"""Migrate frozen Figure 5 candidates using counts, not coordinate relabeling."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
PILOT = PROJECT.parent / '2026-09-02_mouse-time-aware-integration-pilot/data'
ROUTES = [('mouse', 'CIP', 'gse103221_pluripotency_reprogramming_core'),
          ('mouse', 'WOT', 'gse115943_pluripotency_reprogramming_core'),
          ('human', 'GSE118258', 'gse118258_human_pluripotency_reprogramming'),
          ('human', 'GSE150311', 'gse150311_human_pluripotency_reprogramming')]
SOURCES = ['polyadb_main', 'polyasite_broad', 'gencode_polya']


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def mask_for(start, end, strand, cuts):
    return (end <= cuts if strand == '+' else start >= cuts-1).astype(int)


def state_stats(model):
    out = {}
    n = len(model['pas'])
    for threshold in [20, 50]:
        for state in 'SREPM':
            libs = [l for l in model['profiles'].get(state, {}).get('libraries', [])
                    if l['total'] >= threshold and l['coverage'] >= .5]
            counts = np.array([l['counts'] for l in libs], float).reshape(-1, n)
            mean = (counts / counts.sum(1)[:, None]).mean(0) if len(libs) else np.full(n, np.nan)
            out[threshold, state] = (mean, len(libs), counts.sum(0))
    return out


def main():
    # Strand boundary and many-to-one count conservation are scientific invariants.
    assert mask_for(99, 100, '+', np.array([100, 101])).tolist() == [1, 1]
    assert mask_for(99, 100, '-', np.array([100, 101])).tolist() == [1, 0]
    assert np.array_equal(np.array([3, 7, 11]) @ np.array([[1, 0], [1, 0], [0, 1]]), [10, 11])
    helper = load_module('annotation_tracks', PROJECT/'site-recurrence-v1/figure5-publication-v14r4/prepare_annotated_tracks.py')
    added_path = PROJECT/'stage-backtracking-v1/analysis/human_GSE118258_added_RE_profiles.json'
    added = json.loads(added_path.read_text())
    human_reference = pd.read_csv(ROOT/'human_reference_source_records.tsv', sep='\t')
    assert human_reference.target_build.eq('GRCh38_2024_A').all()
    human_reference = human_reference.drop_duplicates(['source_family', 'target_chrom', 'target_strand', 'target_cut_site', 'catalog_pas_id'])
    reference_keys = ['source_family', 'target_chrom', 'target_strand', 'target_cut_site']
    assert not human_reference.duplicated(reference_keys).any(), 'Ambiguous source-to-catalog membership'
    human_lookup = human_reference.set_index(reference_keys).catalog_pas_id.to_dict()
    hashes = {str(added_path): hashlib.sha256(added_path.read_bytes()).hexdigest()}
    human_path = ROOT/'human_reference_source_records.tsv'
    hashes[str(human_path)] = hashlib.sha256(human_path.read_bytes()).hexdigest()
    ledger, summary, conservation, regions, all_models = [], [], [], [], {}
    for species, route, study in ROUTES:
        source = PROJECT/f'full-cohort-v2/data/{species}/{route}/models.json'
        match = PILOT / ('pas-annotation-support-v1/studies' if species == 'mouse' else 'apa-majority-shared-v1/human_annotation') / study / 'primary_pas_reference_matches.tsv.gz'
        for p in [source, match]:
            hashes[str(p)] = hashlib.sha256(p.read_bytes()).hexdigest()
        table = pd.read_csv(match, sep='\t', low_memory=False).set_index('atlas_pas_id')
        assert table.index.is_unique
        models = json.loads(source.read_text())
        transformed = {}
        for gid, old in models.items():
            chosen = []
            for pas in old['pas']:
                row = table.loc[pas['atlas_pas_id']]
                assert row.primary_gene_id == gid and row.cut_site == pas['cut_site']
                assert row.chrom == pas['chrom'] and row.strand == pas['strand']
                eligible = []
                for src in SOURCES:
                    distance = row[src+'_distance_nt']
                    if pd.isna(distance) or not 0 <= distance <= 50:
                        continue
                    identity = (row[src+'_catalog_pas_id'] if species == 'mouse' else
                                human_lookup[src, pas['chrom'], pas['strand'], int(row[src+'_nearest_cut'])])
                    parts = identity.split(':')
                    assert len(parts) == 7 and parts[1] == species and parts[3:5] == [pas['chrom'], pas['strand']]
                    source_cut = int(row[src+'_nearest_cut'])
                    assert abs(source_cut-pas['cut_site']) == distance
                    eligible.append(dict(annotation_id=identity, annotation_cut=int(parts[5]),
                                         source=src, distance_nt=int(distance)))
                selected = helper.choose(eligible)
                assert selected is not None, ('Unmatched PAS; coverage policy required', species, route, gid, pas['atlas_pas_id'])
                pick = selected[0]
                chosen.append(pick)
                ledger.append(dict(species=species, route=route, gene_id=gid, gene=old['gene_symbol'],
                    family_id=old['family']['family_id'], atlas_pas_id=pas['atlas_pas_id'],
                    raw_cut=pas['cut_site'], annotation_id=pick['annotation_id'], annotation_cut=pick['annotation_cut'],
                    sources=';'.join(sorted({s['source'] for s in selected})),
                    nearest_source_distance=min(s['distance_nt'] for s in selected),
                    representative_distance=abs(pick['annotation_cut']-pas['cut_site'])))
            annotations = sorted({(c['annotation_id'], c['annotation_cut']) for c in chosen},
                                 key=lambda c: (c[1], c[0]), reverse=old['family']['strand'] == '-')
            assert len({a[0] for a in annotations}) == len(annotations)
            membership = np.array([[int(c['annotation_id'] == aid) for aid, cut in annotations] for c in chosen])
            assert (membership.sum(1) == 1).all()
            model = copy.deepcopy(old)
            if route == 'GSE118258':
                model['profiles'].update(copy.deepcopy(added.get(gid, {})))
            for state, profile in model['profiles'].items():
                for kind, obj in [('pooled', profile)] + [('library', l) for l in profile.get('libraries', [])]:
                    if 'counts' not in obj:
                        continue
                    raw = np.asarray(obj['counts'])
                    assert np.all(raw >= 0) and np.array_equal(raw, raw.astype(np.int64))
                    assert np.isclose(raw.sum(), obj['total'])
                    merged = raw.astype(np.int64) @ membership
                    assert merged.sum() == raw.sum()
                    obj['counts'] = merged.tolist()
                    obj['usage'] = (merged/merged.sum()).tolist() if merged.sum() else [0.]*len(merged)
                    conservation.append(dict(species=species, route=route, gene_id=gid, state=state,
                        kind=kind, library=obj.get('library', ''), before=int(raw.sum()), after=int(merged.sum())))
            model['pas'] = [dict(atlas_pas_id=aid, annotation_id=aid, cut_site=cut,
                member_cut_site_min=cut, member_cut_site_max=cut, chrom=old['family']['chrom'],
                strand=old['family']['strand']) for aid, cut in annotations]
            model['family']['pas_ids'] = [a[0] for a in annotations]
            model['native_eligible_contrasts'] = old['eligible_contrasts']
            if len(annotations) < 2:
                model['eligible_contrasts'] = []
            transformed[gid] = model
            summary.append(dict(species=species, route=route, gene_id=gid, gene=old['gene_symbol'],
                native_pas=len(chosen), annotated_pas=len(annotations), single_pas_after_merge=len(annotations)<2))
            fam = model['family']
            cuts = np.array([a[1] for a in annotations])
            boundaries = cuts if fam['strand'] == '+' else cuts-1
            edges = sorted({fam['common_start0'], fam['common_end0']} |
                           {int(v) for v in boundaries if fam['common_start0'] < v < fam['common_end0']})
            stats = state_stats(model)
            for start, end in zip(edges[:-1], edges[1:]):
                mask = mask_for(start, end, fam['strand'], cuts)
                if mask.min() == mask.max():
                    continue
                for contrast in model['eligible_contrasts']:
                    if contrast not in ['SP', 'MP']:
                        continue
                    effects = {f'effect{t}': float((stats[t, contrast[1]][0]-stats[t, contrast[0]][0]) @ mask) for t in [20, 50]}
                    observed = all(min(float(stats[50, state][2] @ mask), float(stats[50, state][2] @ (1-mask))) >= 5
                                   for state in contrast)
                    regions.append(dict(species=species, route=route, gene_id=gid, gene_symbol=model['gene_symbol'],
                        chrom=fam['chrom'], strand=fam['strand'], start0=start, end0=end,
                        region_id=f'{route}|{gid}|annotation:{start}-{end}', mask=json.dumps(mask.tolist()),
                        contrast=contrast, observed_both5=observed, **effects))
        dest = ROOT/f'full-cohort-v2/data/{species}/{route}/models.json'
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(json.dumps(transformed))
        all_models[species, route] = transformed
        print('merged', species, route, len(transformed), flush=True)
    pd.DataFrame(ledger).to_csv(ROOT/'pas_annotation_mapping.tsv.gz', sep='\t', index=False)
    pd.DataFrame(summary).to_csv(ROOT/'family_migration.tsv', sep='\t', index=False)
    pd.DataFrame(conservation).to_csv(ROOT/'count_conservation.tsv.gz', sep='\t', index=False)
    dest = ROOT/'cis-regulatory-assessment-v1/analysis/all_coordinate_regions.tsv.gz'
    dest.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(regions).to_csv(dest, sep='\t', index=False)
    for species in ['mouse', 'human']:
        source = PROJECT/f'full-cohort-v2/data/{species}/site_effects.tsv.gz'
        hashes[str(source)] = hashlib.sha256(source.read_bytes()).hexdigest()
        frame = pd.read_csv(source, sep='\t', low_memory=False)
        parts = []
        for (route, gid, contrast), group in frame.groupby(['route', 'gene_id', 'contrast'], sort=False):
            m = all_models[species, route][gid]
            if contrast not in m['eligible_contrasts']:
                continue
            cuts = np.array([p['cut_site'] for p in m['pas']])
            masks = (group.end0.to_numpy()[:, None] <= cuts if m['family']['strand'] == '+'
                     else group.start0.to_numpy()[:, None] >= cuts-1).astype(int)
            g = group.copy()
            g['mask_json'] = [json.dumps(np.column_stack([r, r]).tolist(), separators=(',', ':')) for r in masks]
            g['exact'] = True
            g['variable'] = masks.min(1) != masks.max(1)
            stats = state_stats(m)
            for t in [20, 50]:
                ref, target = stats[t, contrast[0]], stats[t, contrast[1]]
                g[f'eq{t}_lo'] = masks @ (target[0]-ref[0])
                g[f'eq{t}_hi'] = g[f'eq{t}_lo']
                g[f'eq{t}_n_reference'] = ref[1]
                g[f'eq{t}_n_target'] = target[1]
            ref, target = [m['profiles'][s] for s in contrast]
            g['pooled_lo'] = masks @ (np.array(target['usage'])-np.array(ref['usage']))
            g['pooled_hi'] = g.pooled_lo
            parts.append(g)
        out = pd.concat(parts, ignore_index=True)
        out.to_csv(ROOT/f'full-cohort-v2/data/{species}/site_effects.tsv.gz', sep='\t', index=False)
        print('recomputed sites', species, len(out), flush=True)
    report = dict(status='PASS', count_checks=len(conservation), unmatched_pas=0,
        single_pas_families=sum(s['single_pas_after_merge'] for s in summary),
        native_pas=len(ledger), annotated_pas=sum(s['annotated_pas'] for s in summary),
        representative_distance_over50=sum(r['representative_distance']>50 for r in ledger),
        input_sha256=hashes, tests=3, missing_states_imputed=False)
    (ROOT/'migration_QA.json').write_text(json.dumps(report, indent=2))
    print(json.dumps({k:v for k,v in report.items() if k != 'input_sha256'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
