"""Keep shared seed groups distinct from reciprocal homologous physical positions."""
import json

import pandas as pd

from analyze_sites import HERE, read, save


def same_position(hit, seed):
    return (hit['chrom'], hit['strand'], hit['start0'], hit['end0']) == (seed.chrom, seed.strand, seed.start0, seed.end0)


def main():
    seeds = read(HERE / 'strict_case_predicted_seeds.tsv')
    intervals = read(HERE / 'strict_case_intervals.tsv')
    mapped = {r['query_id']: r for r in map(json.loads, (HERE / 'seed_chain_maps.jsonl').read_text().splitlines())}
    summaries, links = [], []
    key = ['source_population', 'tolerance', 'contrast', 'human_gene_id']
    for values, group in seeds.groupby(key, sort=True):
        context = dict(zip(key, values))
        cores = intervals[(intervals[key] == pd.Series(context)).all(axis=1)].set_index('species')
        for row in group.itertuples():
            other_sp = 'mouse' if row.species == 'human' else 'human'
            status, counterpart, mapped_core = 'single_side_seed_group_prediction', '', False
            count = 0
            if row.both_species_same_seed_group:
                record = mapped[row.query_id]
                count = len(record['hits'])
                status = record['status']
                if status == 'unique':
                    hit = record['hits'][0]
                    core = cores.loc[other_sp]
                    mapped_core = hit['chrom'] == core.chrom and hit['strand'] == core.strand and hit['start0'] >= core.start0 and hit['end0'] <= core.end0
                    candidates = group[group.species.eq(other_sp) & group.seed.eq(row.seed) &
                                       group.chrom.eq(hit['chrom']) & group.strand.eq(hit['strand']) &
                                       group.start0.eq(hit['start0']) & group.end0.eq(hit['end0'])]
                    assert len(candidates) <= 1
                    status = 'unique_without_exact_opposite_physical_seed'
                    if len(candidates):
                        other = next(candidates.itertuples())
                        back = mapped[other.query_id]
                        if back['status'] == 'unique' and same_position(back['hits'][0], row):
                            status, counterpart = 'reciprocal_homologous_predicted_position', other.query_id
                            if row.species == 'human':
                                links.append(dict(context, seed=row.seed, human_query_id=row.query_id, mouse_query_id=other.query_id,
                                                  human_start0=row.start0, human_end0=row.end0, mouse_start0=other.start0, mouse_end0=other.end0,
                                                  human_aliases=row.mirbase_aliases, mouse_aliases=other.mirbase_aliases,
                                                  status=status))
                        else:
                            status = 'opposite_physical_seed_reverse_not_unique_or_not_exact'
            summaries.append(dict(context, query_id=row.query_id, species=row.species, gene_id=row.gene_id, seed=row.seed,
                                  chrom=row.chrom, strand=row.strand, start0=row.start0, end0=row.end0,
                                  same_seed_group_both_species=row.both_species_same_seed_group,
                                  mapping_candidates=count, mapped_into_opposite_core=mapped_core,
                                  status=status, counterpart_query_id=counterpart))
    output = pd.DataFrame(summaries)
    link_columns = key + ['seed', 'human_query_id', 'mouse_query_id', 'human_start0', 'human_end0', 'mouse_start0', 'mouse_end0', 'human_aliases', 'mouse_aliases', 'status']
    link = pd.DataFrame(links, columns=link_columns)
    save(output, 'seed_position_correspondence.tsv')
    save(link, 'reciprocal_predicted_seed_pairs.tsv')
    rows = []
    for values, group in intervals.groupby(key, sort=True):
        context = dict(zip(key, values))
        subset = output[(output[key] == pd.Series(context)).all(axis=1)]
        pairs = link[(link[key] == pd.Series(context)).all(axis=1)]
        rows.append(dict(context, human_gene_symbol=group.human_gene_symbol.iloc[0],
                         human_predicted_physical_sites=int(subset.species.eq('human').sum()),
                         mouse_predicted_physical_sites=int(subset.species.eq('mouse').sum()),
                         same_seed_groups_both_species=subset[subset.same_seed_group_both_species].seed.nunique(),
                         reciprocal_physical_site_pairs=len(pairs), reciprocal_seed_groups=pairs.seed.nunique(),
                         human_delta_inclusion=group[group.species.eq('human')].delta_inclusion.iloc[0],
                         mouse_delta_inclusion=group[group.species.eq('mouse')].delta_inclusion.iloc[0]))
    save(pd.DataFrame(rows), 'strict_case_seed_summary.tsv')
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
