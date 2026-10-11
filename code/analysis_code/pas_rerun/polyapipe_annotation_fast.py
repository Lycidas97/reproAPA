"""Equivalent interval lookup without per-query pandas frame construction.

Primary assignment, neighbor policy, sequence QC, and support rules are unchanged.
The original stable interval sort and candidate order are retained because ties
in primary assignment can depend on that order.
"""
from collections import defaultdict
from pathlib import Path
import subprocess
import tempfile

import numpy as np
import pandas as pd


class CandidateIndex:
    def __init__(self, annotation, intervals):
        self.ann = annotation
        original = annotation._build_interval_index(intervals)
        self.tables = {}
        for table_name, lookup in original.items():
            groups = {}
            for key, (frame, starts) in lookup.items():
                ends = frame.end.to_numpy(np.int64)
                groups[key] = (starts, ends, np.maximum.accumulate(ends), frame.to_dict('records'))
            self.tables[table_name] = groups

    @staticmethod
    def hits(item, start, end):
        starts, ends, prefix_ends, records = item
        right = int(np.searchsorted(starts, end, side='left'))
        # Every earlier interval has end <= start and cannot overlap the query.
        left = min(right, int(np.searchsorted(prefix_ends, start, side='right')))
        for row in np.flatnonzero(ends[left:right] > start) + left:
            yield records[row]

    def generate(self, pas_map, cut_site_window_bp=25):
        ann = self.ann
        rows = []
        for pas in pas_map.itertuples(index=False):
            start = max(0, int(pas.cut_site)-cut_site_window_bp)
            end = int(pas.cut_site)+cut_site_window_bp+1
            found = False
            for same, relation in [(True, 'same_strand'), (False, 'opposite_strand')]:
                for table in ann.REGION_TABLES:
                    lookup = self.tables.get(table, {})
                    if same:
                        value = lookup.get((str(pas.chrom), str(pas.strand)))
                        groups = [] if value is None else [value]
                    else:
                        groups = [value for (c, s), value in lookup.items()
                                  if c == str(pas.chrom) and s != str(pas.strand)]
                    for group in groups:
                        for hit in self.hits(group, start, end):
                            rows.append(ann._candidate_row(pas, hit, relation, start, end))
                            found = True
            for (c, s), group in self.tables.get('gene_body', {}).items():
                if c == str(pas.chrom) and s != str(pas.strand):
                    for hit in self.hits(group, start, end):
                        hit = dict(hit, candidate_region_class='opposite_strand_gene_body_pas', region_priority=12)
                        rows.append(ann._candidate_row(pas, hit, 'opposite_strand', start, end))
                        found = True
            if not found:
                rows.append(ann._unassigned_candidate(pas))
        return ann.annotate_neighbor_gene_conflicts(pd.DataFrame(rows))


class BedtoolsCandidateIndex(CandidateIndex):
    """Batch interval intersection; restore the legacy candidate ordering."""
    def __init__(self, annotation, intervals, bedtools):
        super().__init__(annotation, intervals)
        self.bedtools = str(bedtools)
        self.temporary = tempfile.TemporaryDirectory(prefix='polyapipe-bedtools-')
        self.reference_bed = Path(self.temporary.name)/'reference.bed'
        self.records = []
        with self.reference_bed.open('w') as f:
            for table in [*annotation.REGION_TABLES, 'gene_body']:
                for (chrom, strand), item in self.tables.get(table, {}).items():
                    for row in item[3]:
                        index = len(self.records)
                        self.records.append((table, row))
                        f.write(f'{chrom}\t{row["start"]}\t{row["end"]}\t{index}\t0\t{strand}\n')

    def generate(self, pas_map, cut_site_window_bp=25):
        pas_rows = list(pas_map.itertuples(index=False))
        query = ''.join(f'{p.chrom}\t{max(0,int(p.cut_site)-cut_site_window_bp)}\t'
            f'{int(p.cut_site)+cut_site_window_bp+1}\t{i}\t0\t{p.strand}\n' for i,p in enumerate(pas_rows))
        pairs = defaultdict(list)
        if self.records:
            output = subprocess.run([self.bedtools,'intersect','-a','stdin','-b',str(self.reference_bed),'-wa','-wb'],
                                    input=query,text=True,capture_output=True,check=True)
            for line in output.stdout.splitlines():
                fields = line.split('\t')
                q, r = int(fields[3]), int(fields[9])
                table, hit = self.records[r]
                same = str(pas_rows[q].strand) == str(hit['strand'])
                if table == 'gene_body' and same:
                    continue
                phase = 2 if table == 'gene_body' else (0 if same else 1)
                pairs[q].append((phase, r))
        rows = []
        for i, pas in enumerate(pas_rows):
            start=max(0,int(pas.cut_site)-cut_site_window_bp);end=int(pas.cut_site)+cut_site_window_bp+1
            for phase,r in sorted(pairs[i]):
                _,hit=self.records[r]
                if phase==2:
                    hit=dict(hit,candidate_region_class='opposite_strand_gene_body_pas',region_priority=12)
                rows.append(self.ann._candidate_row(pas,hit,'same_strand' if phase==0 else 'opposite_strand',start,end))
            if not pairs[i]:
                rows.append(self.ann._unassigned_candidate(pas))
        return self.ann.annotate_neighbor_gene_conflicts(pd.DataFrame(rows))


def choose_with_cached_context(annotation, pas_map, candidates):
    """Compute the unchanged per-transcript best-region summary once per batch."""
    valid = candidates[candidates.strand_relation.eq('same_strand') &
        candidates.candidate_gene_id.astype(str).str.len().gt(0) &
        candidates.candidate_transcript_id.astype(str).str.len().gt(0)].copy()
    valid['_priority'] = pd.to_numeric(valid.region_priority, errors='coerce').fillna(99)
    if len(valid):
        best = valid.groupby(['pas_id','candidate_gene_id','candidate_transcript_id'],sort=False)._priority.transform('min')
        valid=valid[valid._priority.eq(best)].copy()
        valid['_family']=valid.candidate_region_class.map(annotation._transcript_region_family)
    families = defaultdict(dict)
    for (pas, gene), group in valid.groupby(['pas_id','candidate_gene_id'],sort=False):
        families[pas][str(gene)] = set(x for x in group._family if x)
    conflict, context = {}, {}
    for pas, genes in families.items():
        conflict[pas] = any(len(v)>1 for v in genes.values())
        terminal={g for g,v in genes.items() if 'terminal_3end' in v}
        intron={g for g,v in genes.items() if 'intronic' in v}
        if any({'terminal_3end','intronic'} <= v for v in genes.values()):
            context[pas]='alternative_terminal_intron'
        elif any(a!=b for a in terminal for b in intron):
            context[pas]='cross_gene_terminal_intron'
        elif conflict[pas]:
            context[pas]='alternative_transcript_context'
        else:
            context[pas]='none'
    old_conflict=annotation._has_transcript_conflict
    old_context=annotation._transcript_context_class
    def cached_conflict(frame):
        return False if frame.empty else conflict.get(frame.pas_id.iloc[0],False)
    def cached_context(frame):
        return 'none' if frame.empty else context.get(frame.pas_id.iloc[0],'none')
    # Worker-local, synchronous call; no global patch survives this function.
    annotation._has_transcript_conflict=cached_conflict
    annotation._transcript_context_class=cached_context
    try:
        return annotation.choose_primary_assignments(pas_map,candidates)
    finally:
        annotation._has_transcript_conflict=old_conflict
        annotation._transcript_context_class=old_context
