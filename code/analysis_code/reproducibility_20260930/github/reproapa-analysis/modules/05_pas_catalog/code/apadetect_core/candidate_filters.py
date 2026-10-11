"""Composable candidate-filter interface for apaDetect."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import pandas as pd


@dataclass(frozen=True)
class CandidateFilter:
    """One auditable filter step; ``evaluate`` returns a removal mask."""

    name: str
    evaluate: Callable[[pd.DataFrame], pd.Series]


def apply_candidate_filters(
    candidates: pd.DataFrame,
    steps: list[CandidateFilter],
) -> tuple[pd.DataFrame, pd.Series]:
    """Apply ordered filters and retain a separate audit column for each step."""
    result = candidates.copy()
    retained = pd.Series(True, index=result.index, dtype=bool)
    seen = set()
    for step in steps:
        if not step.name or step.name in seen:
            raise ValueError(f"Candidate filter names must be unique: {step.name!r}")
        seen.add(step.name)
        removed = step.evaluate(result)
        if not isinstance(removed, pd.Series):
            removed = pd.Series(removed, index=result.index)
        if not removed.index.equals(result.index) or len(removed) != len(result):
            raise ValueError(
                f"Candidate filter {step.name!r} returned a misaligned mask"
            )
        if removed.isna().any():
            raise ValueError(
                f"Candidate filter {step.name!r} returned missing decisions"
            )
        column = f"{step.name}_filtered"
        result[column] = removed.astype(bool)
        retained &= ~result[column]
    return result, retained
