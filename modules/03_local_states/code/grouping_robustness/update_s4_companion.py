"""Append fixed-graph partition diagnostics using the manuscript's current family labels."""
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'notes/2026-09-14-state-grouping-robustness/data'
TABLE = ROOT / 'table/S4-local-partitions.tsv'
FIELDS = ('stability_runs', 'median_ari', 'minimum_ari',
          'median_broad_family_precision', 'minimum_broad_family_precision')


def read_tsv(path):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        return reader.fieldnames, list(reader)


original_fields, partitions = read_tsv(TABLE)
_, annotations = read_tsv(ROOT / 'table/S5-local-state-annotations.tsv')
families = {row['original_state_id']: row['display_group'] for row in annotations}
assert len(partitions) == 11 and len(families) == 111

_, runs = read_tsv(DATA / 'partition_runs.tsv')
ari_by_run = {}
for row in runs:
    if row['role'] == 'sensitivity':
        key = row['unit'], row['resolution'], row['seed']
        assert key not in ari_by_run
        ari_by_run[key] = float(row['ARI'])
assert len(runs) == 110 and len(ari_by_run) == 99

counts_by_run = defaultdict(lambda: defaultdict(int))
with gzip.open(DATA / 'partition_overlaps.tsv.gz', 'rt', newline='') as stream:
    for row in csv.DictReader(stream, delimiter='\t'):
        if row['role'] != 'sensitivity':
            continue
        key = row['unit'], row['resolution'], row['seed']
        assert key in ari_by_run and row['original_state_id'] in families
        counts_by_run[key][row['alternative_cluster'], families[row['original_state_id']]] += int(row['cells'])
assert set(counts_by_run) == set(ari_by_run)

cells_by_unit = {row['analysis_unit']: int(row['n_cells']) for row in partitions}
metrics_by_unit = defaultdict(list)
for key, counts in counts_by_run.items():
    cluster_totals = Counter()
    for (cluster, _), count in counts.items():
        cluster_totals[cluster] += count
    cells = sum(cluster_totals.values())
    assert cells == cells_by_unit[key[0]]
    precision = sum(count * count / cluster_totals[cluster]
                    for (cluster, _), count in counts.items()) / cells
    metrics_by_unit[key[0]].append((ari_by_run[key], precision))
assert set(metrics_by_unit) == set(cells_by_unit)

for row in partitions:
    values = metrics_by_unit[row['analysis_unit']]
    assert len(values) == 9
    ari = [value[0] for value in values]
    precision = [value[1] for value in values]
    row.update(stability_runs=9, median_ari=round(median(ari), 6),
               minimum_ari=round(min(ari), 6),
               median_broad_family_precision=round(median(precision), 6),
               minimum_broad_family_precision=round(min(precision), 6))

fields = list(original_fields) + [field for field in FIELDS if field not in original_fields]
with TABLE.open('w', newline='') as stream:
    writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
    writer.writeheader()
    writer.writerows(partitions)
