"""Full frozen cohort, contrast-specific QC, native PAS site retention."""
import argparse
import csv
import gzip
import importlib.util
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

BASE = Path('/root/nfsdata/apadetect/analysis_projects')
spec = importlib.util.spec_from_file_location('legacy', BASE/'figure5_grouped_preview_v1/scripts/compute.py')
legacy = importlib.util.module_from_spec(spec)
spec.loader.exec_module(legacy)
from structure_consensus import genomic_retention
from core import ambiguous_overlaps, change, classify

CLAIMS = BASE/'reproapa_integration_validation_20260913/execution_v1/A05/claim_registry_v1/claim_index.tsv'


def dump(obj, path):
    def clean(v):
        if isinstance(v, dict):
            return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [clean(x) for x in v]
        if isinstance(v, np.ndarray):
            return clean(v.tolist())
        if isinstance(v, np.generic):
            return clean(v.item())
        if isinstance(v, float) and not np.isfinite(v):
            return None
        return v
    path.write_text(json.dumps(clean(obj), allow_nan=False) + '\n')


def prepare(sp, route, study, claims, models, out, hashes):
    root = legacy.MOUSE/'results_metacell_pas_v7_v1' if sp == 'mouse' else legacy.HUMAN/'results_metacell_pas_v2'
    source = root/'studies'/study
    support_path = (legacy.MOUSE/'results_pas_annotation_support_v1/studies'/study/'primary_pas_reference_matches.tsv.gz'
                    if sp == 'mouse' else legacy.PROFILES/'human_annotation'/study/'primary_pas_reference_matches.tsv.gz')
    paths = [source/'features.tsv.gz', source/'metacells.tsv.gz', source/'metacell_by_pas_counts.npz',
             legacy.PROFILES/'profiles'/study/'selected_metacells.tsv.gz', support_path]
    hashes.update({str(p): legacy.sha(p) for p in paths})
    f, axis, meta = [legacy.read(paths[i]) for i in [0, 1, 3]]
    assert set(f.species) == {sp}
    wanted = set(''.join(claims.contrast.unique()))
    meta = meta[meta.route.eq(route) & meta.node.isin(wanted)].copy()
    assert meta.matrix_row.is_unique and meta.purity.ge(.8).all()
    assert np.array_equal(axis.iloc[meta.matrix_row].metacell_id, meta.metacell_id)
    support = legacy.read(support_path).set_index('atlas_pas_id').matched_source_count_50nt
    f = f[legacy.strict_pas(f) & f.primary_gene_id.isin(claims.gene_id) & f.atlas_pas_id.map(support).ge(1)].copy()
    f['gene_id'] = f.primary_gene_id
    raw = sparse.load_npz(paths[2]).tocsr()
    selected = raw[meta.matrix_row.to_numpy(int)][:, f.feature_index_1based.to_numpy(int)-1].tocsr()
    del raw
    assert np.isfinite(selected.data).all() and (selected.data >= 0).all()
    assert np.equal(selected.data, np.floor(selected.data)).all()
    f = f.reset_index(drop=True)
    ix_map = {p: i for i, p in enumerate(f.atlas_pas_id)}
    counts = {s: np.asarray(selected[meta.node.eq(s).to_numpy()].sum(axis=0)).ravel() for s in sorted(set(meta.node))}
    libs = {(s, lib): np.asarray(selected[np.flatnonzero(meta.node.eq(s) & meta.library_key.eq(lib))].sum(axis=0)).ravel()
            for s, lib in meta[['node', 'library_key']].drop_duplicates().itertuples(index=False, name=None)}
    genes, audit = {}, []
    fg = dict(tuple(f.groupby('gene_id')))
    for gid, gc in claims.groupby('gene_id', sort=True):
        gf = fg.get(gid, f.iloc[:0])
        families, inventory = legacy.resolve_families(gf.to_dict('records'), models.get(gid, []))
        active = [fam for fam in families if len(fam['pas_ids']) >= 2]
        reason = ('known_structure_review_Nopchap1' if gc.iloc[0].gene_symbol == 'Nopchap1' else
                  'no_supported_strict_pas' if gf.empty else
                  'no_terminal_model' if not models.get(gid) else
                  'no_multi_pas_family' if not active else
                  'multiple_multi_pas_families' if len(active) > 1 else '')
        profiles = {}
        if not reason:
            fam = active[0]
            pas = sorted([p for p in inventory if p['family_id'] == fam['family_id']], key=lambda p: p['member_cut_site_min'], reverse=fam['strand'] == '-')
            ix = [ix_map[p['atlas_pas_id']] for p in pas]
            for s, count in counts.items():
                total, denom = float(count[ix].sum()), float(count[gf.index].sum())
                lib_rows = []
                for (state, lib), lc in libs.items():
                    if state != s:
                        continue
                    lt, ld = float(lc[ix].sum()), float(lc[gf.index].sum())
                    lib_rows.append(dict(library=lib, total=lt, coverage=lt/ld if ld else 0.,
                                         counts=lc[ix], usage=lc[ix]/lt if lt else np.zeros(len(ix))))
                profiles[s] = dict(total=total, coverage=total/denom if denom else 0.,
                                   counts=count[ix], usage=count[ix]/total if total else np.zeros(len(ix)), libraries=lib_rows)
            genes[gid] = dict(species=sp, route=route, gene_id=gid, gene_symbol=gc.iloc[0].gene_symbol,
                              family=fam, pas=pas, profiles=profiles)
        eligible = []
        for c in gc.itertuples():
            a, b = c.contrast
            status = reason
            if not status:
                status = ('missing_state' if a not in profiles or b not in profiles else
                          'low_pooled_UMI' if min(profiles[s]['total'] for s in [a,b]) < 50 else
                          'low_family_coverage' if min(profiles[s]['coverage'] for s in [a,b]) < .5 else 'ready')
            if status == 'ready':
                eligible.append(c.contrast)
            audit.append(dict(species=sp, route=route, study=study, gene_id=gid, gene_symbol=c.gene_symbol,
                              contrast=c.contrast, status=status, n_supported_pas=len(gf),
                              n_multi_pas_families=len(active), inventory_statuses=json.dumps(dict(Counter(p['family_status'] for p in inventory))),
                              baseline_supported=c.baseline_supported))
        if gid in genes:
            genes[gid]['eligible_contrasts'] = eligible
    dest = out/route
    dest.mkdir()
    dump(genes, dest/'models.json')
    pd.DataFrame(audit).to_csv(dest/'eligibility.tsv', sep='\t', index=False)
    meta.to_csv(dest/'selected_metacells.tsv.gz', sep='\t', index=False)
    f.to_csv(dest/'selected_features.tsv.gz', sep='\t', index=False)
    sparse.save_npz(dest/'selected_counts.npz', selected)
    print(sp, route, dict(Counter(r['status'] for r in audit)), flush=True)
    return genes, audit


def main(sp, out):
    out.mkdir(parents=True, exist_ok=False)
    hashes = {str(CLAIMS): legacy.sha(CLAIMS)}
    assert hashes[str(CLAIMS)] == '449b59569f9bfc0c4b6b187893883f481c6e484bbe40ff0efc3d5d77a3442ee3'
    claims = legacy.read(CLAIMS)
    claims = claims[claims.species.eq(sp)]
    cat, catalog, bundle = legacy.load_catalog(sp, out, hashes)
    utr_path = legacy.MOUSE/'results_tandem_apa_v1/all_reference_terminal_utr_candidates.tsv.gz' if sp == 'mouse' else catalog/'all_reference_terminal_utr_candidates.tsv.gz'
    hashes[str(utr_path)] = legacy.sha(utr_path)
    u = legacy.read(utr_path)
    c = pd.DataFrame(cat.rows('SELECT u.*,r.chrom,r.strand,r.start0,r.end0 FROM utrs u JOIN regions r USING(region_id)'))
    cg = dict(tuple(c[c.gene_id.isin(claims.gene_id)].groupby('gene_id')))
    models = {gid: legacy.models_for_gene(g, cg.get(gid, c.iloc[:0])) for gid, g in u[u.gene_id.isin(claims.gene_id)].groupby('gene_id')}
    all_genes, audits = {}, []
    for route, study in legacy.STUDIES[sp]:
        genes, audit = prepare(sp, route, study, claims, models, out, hashes)
        all_genes[route] = genes
        audits.extend(audit)
    pd.DataFrame(audits).to_csv(out/'eligibility.tsv', sep='\t', index=False)
    regions = {}
    for genes in all_genes.values():
        for g in genes.values():
            if not g['eligible_contrasts']:
                continue
            for m in g['family']['models']:
                rid = m['region_id']
                if rid not in regions:
                    regions[rid] = cat.rows('SELECT * FROM regions WHERE region_id=?', (rid,))[0]
    print(sp, 'scan_regions', len(regions), flush=True)
    transferred = {}
    if sp == 'mouse':
        sys.path.insert(0, str(legacy.MOUSE/'regulatory_query_v2/scripts'))
        from query_expanded import scan_regions
        work = out/'scan_logs'
        work.mkdir()
        transferred, scan = scan_regions(cat, regions, bundle, Path(bundle['base_bundle']).parent/'transferred_scannable.meme', work)
    else:
        cat.prepare_query_regions(list(regions.values()), [])
        scan = cat.scan_runs
    dump(scan, out/'scan_manifest.json')
    ann_path = BASE/f'regulatory_resource_revision_v1/resources/{sp}_mirna_annotations.tsv'
    ann = legacy.read(ann_path).fillna('')
    hashes[str(ann_path)] = legacy.sha(ann_path)
    mirnames = {seed: sorted({name for aliases in g.mirbase_exact_aliases for name in aliases.split(';') if name}) for seed, g in ann.groupby('seed_dna')}
    direct_names = {mid: sorted({r['RBP_Name'] for r in refs}) for mid, refs in cat.motifs.items()}
    transfer_names = defaultdict(set)
    for r in bundle.get('transfer_provenance', []):
        transfer_names[r['Motif_ID']].add(r['mouse_symbol'])
    patterns = legacy.canonical_patterns(bundle['mirna_seeds'])
    matcher = re.compile('(?=('+'|'.join(sorted(patterns, key=lambda x: (-len(x), x)))+'))')
    cached = {}
    for rid, region in regions.items():
        hits = cat.rbp_sites(region)
        for h in hits:
            cat.decorate_site(h, region)
        seq = cat.sequences.fetch(rid).upper()
        for start, length, seed, kind in legacy.seed_hits(seq, patterns, matcher):
            h = dict(kind='mirna_seed', feature_id=seed, site_type=kind, start0=start, end0=start+length)
            cat.decorate_site(h, region)
            hits.append(h)
        hits.extend(transferred.get(rid, []))
        cached[rid] = hits
    kindmap = {'mirna_seed':'miRNA', 'rbp_pwm':'RBP', 'rbp_pwm_transfer':'RBP_transfer'}
    audit = Counter()
    coverage = []
    with gzip.open(out/'site_effects.tsv.gz', 'wt') as handle:
        writer = None
        for route, genes in all_genes.items():
            for gid, g in genes.items():
                if not g['eligible_contrasts']:
                    continue
                hits = {}
                factor_spans = defaultdict(set)
                for m in g['family']['models']:
                    for h in cached[m['region_id']]:
                        if h['kind'] not in kindmap:
                            audit['non_PWM_or_seed'] += 1
                            continue
                        names = (mirnames.get(h['feature_id'], []) if h['kind']=='mirna_seed' else
                                 direct_names.get(h['feature_id'], []) if h['kind']=='rbp_pwm' else sorted(transfer_names[h['feature_id']]))
                        if not names:
                            audit['unmapped_names'] += 1
                            continue
                        result = legacy.classify_site(h, g['family'], g['pas'])
                        if not result['site_status'].startswith('model_consensus'):
                            audit['outside_consensus'] += 1
                            continue
                        key = (kindmap[h['kind']], int(h['genomic_start0']), int(h['genomic_end0']))
                        factor = h['feature_id'] if h['kind']=='mirna_seed' else ';'.join(names)
                        if key not in hits:
                            mask = np.array([result['retention_by_pas'][p['atlas_pas_id']] for p in g['pas']], dtype=int)
                            hits[key] = dict(mask=mask, factors={}, types=set())
                        hits[key]['factors'][factor] = names
                        hits[key]['types'].add(h['site_type'])
                        factor_spans[(key[0], factor)].add(key[1:])
                ambiguous = {k: ambiguous_overlaps(v) for k, v in factor_spans.items()}
                nsites = Counter()
                for (catname, start, end), h in hits.items():
                    factors = {f:n for f,n in h['factors'].items() if (start,end) not in ambiguous[(catname,f)]}
                    if not factors:
                        audit['all_factor_links_overlap_ambiguous'] += 1
                        continue
                    nsites[catname] += 1
                    mask = h['mask']
                    fam = g['family']
                    base = dict(species=sp, route=route, gene_id=gid, gene_symbol=g['gene_symbol'],
                                family_id=fam['family_id'], category=catname, chrom=fam['chrom'], strand=fam['strand'],
                                start0=start, end0=end, names=';'.join(sorted({n for ns in factors.values() for n in ns})),
                                feature_groups=';'.join(sorted(factors)), site_types=';'.join(sorted(h['types'])),
                                mask_json=json.dumps(mask.tolist(), separators=(',', ':')),
                                exact=bool(np.array_equal(mask[:,0], mask[:,1])),
                                variable=bool(mask[:,0].min()!=mask[:,1].max()))
                    for contrast in g['eligible_contrasts']:
                        a, b = contrast
                        pa, pb = g['profiles'][a], g['profiles'][b]
                        lo, hi = change(pb['usage'], pa['usage'], mask)
                        row = dict(base, contrast=contrast, pooled_lo=lo, pooled_hi=hi,
                                   reference_UMI=pa['total'], target_UMI=pb['total'])
                        for cutoff in [20,50]:
                            means, ns = [], []
                            for pr in [pa,pb]:
                                use = [x['usage'] for x in pr['libraries'] if x['total']>=cutoff and x['coverage']>=.5]
                                ns.append(len(use))
                                means.append(np.mean(use,axis=0) if use else None)
                            low, high = change(means[1],means[0],mask) if all(ns) else (float('nan'),float('nan'))
                            row.update({f'eq{cutoff}_lo':low,f'eq{cutoff}_hi':high,
                                        f'eq{cutoff}_n_reference':ns[0],f'eq{cutoff}_n_target':ns[1]})
                        if writer is None:
                            writer = csv.DictWriter(handle, fieldnames=list(row), delimiter='\t')
                            writer.writeheader()
                        writer.writerow(row)
                coverage.append(dict(species=sp,route=route,gene_id=gid,gene_symbol=g['gene_symbol'],
                                     eligible_contrasts=';'.join(g['eligible_contrasts']),
                                     **{c+'_sites':nsites[c] for c in kindmap.values()}))
    pd.DataFrame(coverage).to_csv(out/'site_coverage.tsv',sep='\t',index=False)
    cat.close()
    cat.sequences.close()
    for p, digest in hashes.items():
        assert legacy.sha(Path(p)) == digest, p
    dump(dict(status='COMPUTED_DESCRIPTIVE_NOT_INFERENTIAL',species=sp,sources=hashes,
              audit=dict(audit), candidate_genes=claims.gene_id.nunique(),
              regions_scanned=len(regions), source_state_version='frozen_original_coarse_not_A07',
              output_hashes={str(p.relative_to(out)):legacy.sha(p) for p in out.rglob('*') if p.is_file()}), out/'complete.json')
    print(sp, 'COMPLETE', flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--species',choices=['mouse','human'],required=True)
    p.add_argument('--out',type=Path,required=True)
    a = p.parse_args()
    main(a.species,a.out)
