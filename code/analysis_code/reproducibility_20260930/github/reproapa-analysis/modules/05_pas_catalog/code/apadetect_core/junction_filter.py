"""Single-pass splice-junction collection and candidate filtering."""

from __future__ import annotations

import bisect
from collections import Counter
from dataclasses import dataclass, field

import pandas as pd


MATCH_OPERATIONS = {0, 7, 8}
REFERENCE_OPERATIONS = {0, 2, 3, 7, 8}


def junction_observations(
    cigartuples,
    reference_start: int,
    strand: str,
    min_intron: int = 70,
    max_intron: int = 500_000,
):
    """Yield intron coordinates and adjacent anchors from one CIGAR."""
    if not cigartuples:
        return
    reference_position = int(reference_start)
    for index, (operation, length) in enumerate(cigartuples):
        if operation == 3:
            left_anchor = (
                cigartuples[index - 1][1]
                if index > 0 and cigartuples[index - 1][0] in MATCH_OPERATIONS
                else 0
            )
            right_anchor = (
                cigartuples[index + 1][1]
                if index + 1 < len(cigartuples)
                and cigartuples[index + 1][0] in MATCH_OPERATIONS
                else 0
            )
            if min_intron <= length <= max_intron:
                yield (
                    reference_position,
                    reference_position + length,
                    strand,
                    left_anchor,
                    right_anchor,
                )
        if operation in REFERENCE_OPERATIONS:
            reference_position += length


def read_junctions(
    cigartuples,
    reference_start: int,
    strand: str,
    min_anchor: int = 8,
    min_intron: int = 70,
    max_intron: int = 500_000,
):
    """Yield introns whose two anchors pass on this alignment."""
    for start, end, current_strand, left_anchor, right_anchor in junction_observations(
        cigartuples,
        reference_start,
        strand,
        min_intron,
        max_intron,
    ):
        if left_anchor >= min_anchor and right_anchor >= min_anchor:
            yield start, end, current_strand


@dataclass
class JunctionEvidence:
    """Regtools-compatible counts and maximum anchors for one chromosome."""

    records: dict = field(default_factory=dict)

    @property
    def counts(self):
        return Counter({key: values[0] for key, values in self.records.items()})

    def add_alignment(
        self,
        cigartuples,
        reference_start: int,
        strand: str,
        min_intron: int = 70,
        max_intron: int = 500_000,
    ) -> int:
        if not cigartuples:
            return 0
        records = self.records
        reference_position = int(reference_start)
        added = 0
        cigar_count = len(cigartuples)
        for index, (operation, length) in enumerate(cigartuples):
            if operation == 3 and min_intron <= length <= max_intron:
                left = (
                    cigartuples[index - 1][1]
                    if index > 0 and cigartuples[index - 1][0] in MATCH_OPERATIONS
                    else 0
                )
                right = (
                    cigartuples[index + 1][1]
                    if index + 1 < cigar_count
                    and cigartuples[index + 1][0] in MATCH_OPERATIONS
                    else 0
                )
                key = (reference_position, reference_position + length, strand)
                values = records.get(key)
                if values is None:
                    records[key] = [1, left, right]
                else:
                    values[0] += 1
                    if left > values[1]:
                        values[1] = left
                    if right > values[2]:
                        values[2] = right
                added += 1
            if operation in REFERENCE_OPERATIONS:
                reference_position += length
        return added

    def points(
        self,
        min_count: int,
        min_anchor: int = 8,
    ) -> dict[str, list[tuple[int, int]]]:
        """Return transcript-proximal junction points and read support."""
        result = {"+": [], "-": []}
        for (intron_start, intron_end, strand), values in self.records.items():
            count, max_left_anchor, max_right_anchor = values
            if (
                count < min_count
                or max_left_anchor < min_anchor
                or max_right_anchor < min_anchor
            ):
                continue
            point = intron_start - 1 if strand == "+" else intron_end
            result[strand].append((point, int(count)))
        for strand in result:
            result[strand].sort()
        return result


def junction_overlap_mask(
    candidates: pd.DataFrame,
    points: dict[str, list[tuple[int, int]]],
    left: int,
    right: int,
) -> pd.Series:
    """Match bedtools window ``-sw -sm -v`` removal decisions."""
    positions = {
        strand: [position for position, _ in values]
        for strand, values in points.items()
    }
    removed = []
    for row in candidates.itertuples():
        current = positions.get(row.strand, ())
        if row.strand == "+":
            lower, upper = int(row.start) - left, int(row.end) + right
        else:
            lower, upper = int(row.start) - right, int(row.end) + left
        index = bisect.bisect_left(current, lower)
        removed.append(index < len(current) and current[index] < upper)
    return pd.Series(removed, index=candidates.index, dtype=bool)


def junction_frame(
    chromosome: str,
    points: dict[str, list[tuple[int, int]]],
) -> pd.DataFrame:
    rows = []
    for strand in ("+", "-"):
        for position, count in points.get(strand, ()):
            rows.append(
                {
                    "chrom": chromosome,
                    "start": position,
                    "end": position + 1,
                    "name": f"JUNC{len(rows) + 1:08d}",
                    "score": count,
                    "strand": strand,
                }
            )
    return pd.DataFrame(
        rows,
        columns=["chrom", "start", "end", "name", "score", "strand"],
    )
