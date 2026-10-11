"""Strict driver intervals and complete predicted seeds, with native boundary cores."""
from collections import defaultdict
import hashlib
import json
import re

import numpy as np
import pandas as pd

from analyze_sites import HERE, ROOT, TABLES, read, save, INPUTS


def core_interval(strand, proximal_cuts, distal_cuts):
    if strand == '+':
        return max(proximal_cuts), min(distal_cuts)
    return max(distal_cuts), min(proximal_cuts)


def inclusion(strand, lo, hi, cuts, usage):
    keep = np.asarray(cuts) >= hi if strand == '+' else np.asarray(cuts) <= lo
    return float(np.asarray(usage)[keep].sum())


def main():
    contexts = read(HERE / 'site_event_contexts.tsv.gz')
    cases = contexts[contexts.source_population.isin(['main48', 'hN79']) & contexts.gene_threshold.eq(.02) &
                     contexts.parent_case & contexts.context_role.eq('shared_recurrent') & contexts.strict_double_driver]
    traces = read(HERE / 'site_usage_traces.tsv.gz')
    traces = traces[traces.source_population.isin(['main48', 'hN79']) & traces.gene_threshold.eq(.02)]
    catalog = read(TABLES / 'cis/annotation_catalog.tsv.gz').set_index('id')
    mapping_frame = read(TABLES / 'cis/detected_to_annotation.tsv.gz')
    mapping = mapping_frame.set_index('native_pas_id').id.to_dict()
    cuts_by_native = mapping_frame.set_index('native_pas_id').native_cut.to_dict()
    route = read(TABLES / 'endpoints/route_effects.tsv.gz')
    route = route[route.setting.eq('supported50_p08') &
                  ((route.contrast.eq('P_S') & route.design.eq('true_source')) |
                   (route.contrast.eq('P_M') & route.design.eq('same_library')))]
    naive = read(TABLES / 'naive_state/library_pair_effects.tsv.gz')
    naive = naive[naive.version.eq('v2_primary') & naive.setting.eq('supported50_p08') &
                  naive.node_a.eq('N') & naive.node_b.eq('P') & naive.library_a.str.endswith(('6', '8')) & naive.library_b.str.endswith('7')]
    utr = read(TABLES / 'tandem/annotation_terminal_utrs.tsv.gz')
    hits = {sp: read(TABLES / f'external_mirna/{sp}_physical_seed_sites.tsv.gz') for sp in ['human', 'mouse']}
    families = {}
    for sp in ['human', 'mouse']:
        path = ROOT / f'notes/2026-10-05-apadetect-analysis-intake/completion-inputs/external/{sp}_mirna_annotations.tsv'
        annotation = read(path)
        families[sp] = annotation.groupby('seed_dna').agg(
            mirbase_aliases=('mirbase_exact_aliases', lambda x: ';'.join(sorted(set(v for value in x.dropna() for v in str(value).split(';') if v)))),
            mirgenedb_aliases=('mirgenedb_exact_aliases', lambda x: ';'.join(sorted(set(v for value in x.dropna() for v in str(value).split(';') if v)))))
    interval_rows, effect_rows, seed_rows = [], [], []
    for row in cases.itertuples():
        context = dict(set_name=row.set_name, source_population=row.source_population, tolerance=row.tolerance,
                       contrast=row.contrast, human_gene_id=row.human_gene_id, mouse_gene_id=row.mouse_gene_id,
                       human_gene_symbol=row.human_gene_symbol)
        for sp in ['human', 'mouse']:
            gene = getattr(row, sp + '_gene_id')
            gain, loss = getattr(row, sp + '_gain_id'), getattr(row, sp + '_loss_id')
            driver = catalog.loc[[gain, loss]].copy()
            assert driver.chrom.nunique() == driver.strand.nunique() == 1
            strand, chrom = driver.strand.iloc[0], driver.chrom.iloc[0]
            ordered = driver.index[np.argsort(driver.cut_site.to_numpy() * (1 if strand == '+' else -1))]
            prox, distal = ordered[0], ordered[1]
            subset = traces[traces.source_population.eq(row.source_population) & traces.tolerance.eq(row.tolerance) &
                            traces.human_gene_id.eq(row.human_gene_id) & traces.contrast.eq(row.contrast) & traces.species.eq(sp)]
            native_profiles, driver_cuts = [], defaultdict(list)
            for p in subset.itertuples():
                if sp == 'human' and row.set_name == 'hN':
                    source = naive[naive.gene_id.eq(gene) & naive.library_a.eq(p.series)]
                    assert len(source) == 1
                    source = source.iloc[0]
                    target, reference = json.loads(source.usage_a), json.loads(source.usage_b)
                else:
                    contrast = 'P_M' if row.set_name == 'hN' else row.contrast
                    source = route[route.species.eq(sp) & route.study_id.eq(p.study_id) & route.gene_id.eq(gene) & route.contrast.eq(contrast)]
                    if p.study_id == 'GSE178325':
                        source = source[source.route.str.endswith(p.series.removeprefix('178325_'))]
                    if p.study_id in ['GSE103221', 'GSE115943']:
                        source = source[source.route.eq(p.series)]
                    assert len(source) == 1
                    source = source.iloc[0]
                    target, reference = json.loads(source.target_usage), json.loads(source.reference_usage)
                ids = json.loads(source.PAS_ids)
                cuts = []
                for key in ids:
                    if key in cuts_by_native:
                        cut = int(cuts_by_native[key])
                    else:
                        cut = int(re.search(r'_(\d+)_[fr]$', key).group(1))
                    cuts.append(cut)
                    if mapping.get(key) in [prox, distal]:
                        driver_cuts[mapping[key]].append(cut)
                assert abs(sum(target) - 1) < 1e-12 and abs(sum(reference) - 1) < 1e-12
                assert abs(float(source.delta_DWUI) - p.original_delta) < 1e-12
                native_profiles.append(dict(study_id=p.study_id, series=p.series, cuts=cuts, target=target, reference=reference))
            assert driver_cuts[prox] and driver_cuts[distal]
            lo, hi = core_interval(strand, driver_cuts[prox], driver_cuts[distal])
            tx = getattr(row, sp + '_common_transcripts').split(';')
            candidates = utr[utr.species.eq(sp) & utr.gene_id.eq(gene) & utr.transcript_id.isin(tx)]
            contains = (lo < hi and len(candidates) == len(tx) and
                        candidates.chrom.astype(str).str.removeprefix('chr').eq(chrom.removeprefix('chr')).all() and
                        candidates.strand.eq(strand).all() and candidates.utr_start.sub(1).le(lo).all() and candidates.utr_end.ge(hi).all())
            assert contains, (gene, lo, hi)
            group = []
            for p in native_profiles:
                a = inclusion(strand, lo, hi, p['cuts'], p['target'])
                b = inclusion(strand, lo, hi, p['cuts'], p['reference'])
                item = dict(context, species=sp, gene_id=gene, chrom=chrom, strand=strand,
                            start0=lo, end0=hi, study_id=p['study_id'], series=p['series'],
                            target_inclusion=a, reference_inclusion=b, delta_inclusion=a-b)
                effect_rows.append(item); group.append(item)
            effects = pd.DataFrame(group).groupby('study_id')[['target_inclusion', 'reference_inclusion', 'delta_inclusion']].mean().mean()
            interval_rows.append(dict(context, species=sp, gene_id=gene, chrom=chrom, strand=strand,
                                      annotation_start0=int(driver.cut_site.min()), annotation_end0=int(driver.cut_site.max()),
                                      start0=lo, end0=hi, length=hi-lo, common_transcripts=';'.join(tx),
                                      proximal_native_cuts=json.dumps(sorted(set(driver_cuts[prox]))),
                                      distal_native_cuts=json.dumps(sorted(set(driver_cuts[distal]))),
                                      target_inclusion=effects.target_inclusion, reference_inclusion=effects.reference_inclusion,
                                      delta_inclusion=effects.delta_inclusion))
            selected = hits[sp][hits[sp].gene_id.eq(gene) & hits[sp].chrom.eq(chrom) & hits[sp].strand.eq(strand) &
                                hits[sp].start0.ge(lo) & hits[sp].end0.le(hi)]
            for hit in selected.itertuples():
                alias = families[sp].loc[hit.seed].to_dict() if hit.seed in families[sp].index else dict(mirbase_aliases='', mirgenedb_aliases='')
                assert hit.end0 - hit.start0 in [7, 8]
                seed_rows.append(dict(context, species=sp, gene_id=gene, seed=hit.seed, chrom=hit.chrom, strand=hit.strand,
                                      start0=hit.start0, end0=hit.end0, core_start0=lo, core_end0=hi, **alias))
    seeds = pd.DataFrame(seed_rows)
    keys = ['source_population', 'tolerance', 'contrast', 'human_gene_id', 'seed']
    seeds['both_species_same_seed_group'] = seeds.groupby(keys).species.transform('nunique').eq(2)
    seeds['query_id'] = seeds.apply(lambda r: '|'.join(map(str, [r.species, r.gene_id, r.seed, r.chrom, r.strand, r.start0, r.end0])), axis=1)
    query = seeds[seeds.both_species_same_seed_group][['query_id', 'species', 'gene_id', 'seed', 'chrom', 'strand', 'start0', 'end0']].drop_duplicates().sort_values('query_id')
    save(pd.DataFrame(interval_rows), 'strict_case_intervals.tsv')
    save(pd.DataFrame(effect_rows), 'strict_case_interval_inclusion.tsv')
    save(seeds, 'strict_case_predicted_seeds.tsv')
    save(query, 'seed_mapping_queries.tsv')
    save(contexts[contexts.parent_case & contexts.human_gene_symbol.eq('EIF4EBP2')], 'EIF4EBP2_qualification.tsv')
    (HERE / 'interval_input_hashes.json').write_text(json.dumps(INPUTS, indent=2) + '\n')
    print('Strict case intervals', len(interval_rows), 'seed-context rows', len(seeds), 'unique small-batch queries', len(query), flush=True)


if __name__ == '__main__':
    main()
