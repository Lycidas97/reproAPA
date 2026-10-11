"""Prepare the shared annotation identity map for tracks and PAS quantification."""
import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent.parent
MATCH_ROOT = PROJECT.parent / '2026-09-02_mouse-time-aware-integration-pilot/data/pas-annotation-support-v1/studies'
STUDIES = {'CIP': 'gse103221_pluripotency_reprogramming_core',
           'WOT': 'gse115943_pluripotency_reprogramming_core'}
GENES = {'ENSMUSG00000003038': 'Hmgn2', 'ENSMUSG00000020091': 'Eif4ebp2'}
SOURCES = ['polyadb_main', 'polyasite_broad', 'gencode_polya']


def choose(candidates):
    # Collapse source-specific cuts by the integrated annotation identity first.
    groups = {}
    for c in candidates:
        groups.setdefault(c['annotation_id'], []).append(c)
    if not groups:
        return None
    return min(groups.values(), key=lambda g: (-len({c['source'] for c in g}),
               min(c['distance_nt'] for c in g), g[0]['annotation_cut'], g[0]['annotation_id']))


def main():
    # Different source coordinates for one catalog ID must still collapse to one mark.
    test = [dict(annotation_id='a', annotation_cut=100, source='x', distance_nt=9),
            dict(annotation_id='a', annotation_cut=100, source='y', distance_nt=8),
            dict(annotation_id='b', annotation_cut=110, source='x', distance_nt=1)]
    assert choose(test)[0]['annotation_id'] == 'a' and choose([]) is None
    mappings, candidates, inputs = [], [], {}
    sites = pd.read_csv(ROOT / 'panel_EF_sites.tsv', sep='\t')
    columns = ['atlas_pas_id', 'primary_gene_id', 'cut_site', 'chrom', 'strand'] + [
        s + suffix for s in SOURCES for suffix in ['_nearest_cut', '_distance_nt', '_catalog_pas_id']]
    for route, study in STUDIES.items():
        path = MATCH_ROOT / study / 'primary_pas_reference_matches.tsv.gz'
        inputs[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
        matched = pd.read_csv(path, sep='\t', usecols=columns).set_index('atlas_pas_id')
        assert matched.index.is_unique
        model_path = PROJECT / f'full-cohort-v2/data/mouse/{route}/models.json'
        inputs[str(model_path)] = hashlib.sha256(model_path.read_bytes()).hexdigest()
        models = json.loads(model_path.read_text())
        for gid, gene in GENES.items():
            model = models[gid]
            for pas in model['pas']:
                row = matched.loc[pas['atlas_pas_id']]
                assert row.primary_gene_id == gid and row.cut_site == pas['cut_site']
                assert row.chrom == pas['chrom'] and row.strand == pas['strand']
                base = dict(gene=gene, gene_id=gid, route=route, atlas_pas_id=pas['atlas_pas_id'],
                            chrom=pas['chrom'], strand=pas['strand'], raw_cut=pas['cut_site'])
                eligible = []
                for source in SOURCES:
                    distance = row[source + '_distance_nt']
                    if pd.isna(distance) or not 0 <= distance <= 50:
                        continue
                    identity = row[source + '_catalog_pas_id']
                    parts = identity.split(':')
                    assert len(parts) == 7 and parts[:3] == ['ext', 'mouse', 'GRCm39_2024_A']
                    assert parts[3:5] == [pas['chrom'], pas['strand']]
                    cut = int(parts[5])  # Frozen DOC-209 integrated-cluster representative.
                    c = dict(**base, source=source, annotation_id=identity, annotation_cut=cut,
                             source_cut=int(row[source + '_nearest_cut']), distance_nt=int(distance))
                    assert abs(c['source_cut'] - base['raw_cut']) == c['distance_nt']
                    eligible.append(c)
                    candidates.append(c)
                selected = choose(eligible)
                if selected is None:
                    mappings.append(dict(**base, status='unmatched', annotation_id='', annotation_cut=None,
                                         sources='', crossed_sites=''))
                    continue
                c = selected[0]
                crossings = []
                for site in sites[sites.gene.eq(gene)].itertuples():
                    assert pas['strand'] == '-'  # These two frozen loci are both minus-strand.
                    # Mirror frozen core.masks, including both cluster bounds and its -1 convention.
                    raw_bounds = [int(site.start0) >= int(pas[k])-1
                                  for k in ['member_cut_site_max', 'member_cut_site_min']]
                    ann_included = int(site.start0) >= c['annotation_cut']-1
                    if any(raw != ann_included for raw in raw_bounds):
                        crossings.append(str(site.site))
                mappings.append(dict(**base, status='matched', annotation_id=c['annotation_id'],
                    annotation_cut=c['annotation_cut'], sources=';'.join(sorted({x['source'] for x in selected})),
                    crossed_sites=';'.join(crossings)))
    mapping = pd.DataFrame(mappings)
    assert not mapping.duplicated(['route', 'atlas_pas_id']).any()
    mapping['annotation_cut'] = pd.array(mapping.annotation_cut, dtype='Int64')
    mapping.to_csv(ROOT / 'panel_EF_pas_annotation_mapping.tsv', sep='\t', index=False)
    pd.DataFrame(candidates).to_csv(ROOT / 'panel_EF_pas_annotation_candidates.tsv', sep='\t', index=False)
    tracks = []
    for keys, g in mapping[mapping.status.eq('matched')].groupby(
            ['gene', 'route', 'chrom', 'strand', 'annotation_id', 'annotation_cut'], sort=True):
        tracks.append(dict(zip(['gene', 'route', 'chrom', 'strand', 'annotation_id', 'annotation_cut'], keys),
            raw_pas_count=len(g), raw_pas_ids=';'.join(g.atlas_pas_id),
            raw_cuts=';'.join(map(str, sorted(g.raw_cut))),
            crosses_site_boundary=g.crossed_sites.ne('').any()))
    pd.DataFrame(tracks).to_csv(ROOT / 'panel_EF_annotated_pas_tracks.tsv', sep='\t', index=False)
    assert sum(t['raw_pas_count'] for t in tracks) == mapping.status.eq('matched').sum()
    counts = []
    for (gene, route), g in mapping.groupby(['gene', 'route']):
        counts.append(dict(gene=gene, route=route, raw_pas=len(g),
            annotated_marks=g.loc[g.status.eq('matched'), 'annotation_id'].nunique(),
            unmatched=int(g.status.eq('unmatched').sum()),
            boundary_crossing=int(g.crossed_sites.ne('').sum())))
    qa = dict(status='PASS', window_nt=50, matching_unit='DOC-209 catalog_pas_id',
        representative='coordinate encoded in frozen integrated catalog ID, not a detected-cut average',
        selection='most eligible source families; nearest source cut; lower representative; ID',
        display_only=False, downstream_quantification='quantify_annotated_pas.py', counts=counts,
        boundary_crossings=mapping[mapping.crossed_sites.ne('')].to_dict('records'), input_sha256=inputs)
    (ROOT / 'annotation_tracks_QA.json').write_text(json.dumps(qa, indent=2))
    print(json.dumps({k: v for k, v in qa.items() if k != 'input_sha256'}, indent=2))


if __name__ == '__main__':
    main()
