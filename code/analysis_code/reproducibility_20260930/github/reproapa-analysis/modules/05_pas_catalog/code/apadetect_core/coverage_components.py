"""Coverage-defined peak components and derivative-interval membership."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class CoverageComponent:
    start: int
    end: int
    strand: str
    candidate_indices: tuple[int, ...]


def positive_intervals(
    positions: np.ndarray, deltas: np.ndarray, chrom_length: int
) -> list[tuple[int, int]]:
    """Return exact half-open intervals over which strand coverage is positive."""
    intervals = []
    coverage = 0
    start = None
    for position, delta in zip(positions, deltas):
        position = int(np.clip(position, 0, chrom_length))
        previous = coverage
        coverage += int(delta)
        if previous == 0 and coverage > 0:
            start = position
        elif previous > 0 and coverage == 0 and start is not None:
            if position > start:
                intervals.append((start, position))
            start = None
    if start is not None and start < chrom_length:
        intervals.append((start, chrom_length))
    return intervals


def assign_candidates(
    candidates: pd.DataFrame,
    intervals: list[tuple[int, int]],
    strand: str,
) -> list[CoverageComponent]:
    """Assign each derivative interval to the coverage component of max overlap."""
    groups: list[list[int]] = [[] for _ in intervals]
    current = candidates[candidates["strand"].eq(strand)]
    for index, row in current.iterrows():
        left, right = int(row["start"]), int(row["end"])
        overlaps = [
            max(0, min(right, end) - max(left, start))
            for start, end in intervals
        ]
        if overlaps and max(overlaps) > 0:
            groups[int(np.argmax(overlaps))].append(int(index))
    return [
        CoverageComponent(start, end, strand, tuple(groups[index]))
        for index, (start, end) in enumerate(intervals)
    ]


def oriented_component_profile(
    coverage: np.ndarray,
    component_start: int,
    anchor_coordinate: int,
    strand: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Orient a complete positive-coverage component around an A/T anchor."""
    genomic = np.arange(component_start, component_start + len(coverage))
    if strand == "+":
        coordinates = genomic - anchor_coordinate
        values = coverage
    elif strand == "-":
        coordinates = anchor_coordinate - genomic[::-1] - 1
        values = coverage[::-1]
    else:
        raise ValueError("strand must be '+' or '-'")
    return coordinates.astype(int), np.asarray(values, dtype=float)
