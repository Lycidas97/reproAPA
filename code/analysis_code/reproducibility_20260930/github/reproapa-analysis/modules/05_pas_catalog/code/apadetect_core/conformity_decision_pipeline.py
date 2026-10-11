#!/usr/bin/env python3
"""Fixed single-PAS localization with BIC-selected multi-PAS mixtures."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pysam

STAGE = Path(__file__).resolve().parent
STAGES = STAGE.parent
sys.path[:0] = [
    str(STAGES / "62_multi_endpoint_mixture"),
    str(STAGES / "66_hard_interval_quantification"),
]

from hard_interval_tool import load_intervals, prepare_candidates  # noqa: E402
from multi_endpoint_mixture import (  # noqa: E402
    assign_candidates_indexed,
    component_endpoints,
    fit_mixture,
    inference,
    model_fixed_offset,
    oriented_coordinate,
    predicted_genomic,
    stream_positive_intervals,
    use_mixture_mle,
)


def bic(fit: dict, observations: int, components: int) -> float:
    return -2.0 * fit["log_likelihood"] + (2 * components - 1) * np.log(
        observations
    )


def selected_coordinate(
    selected_order: int,
    candidate_index: int,
    fixed_coordinate: float,
    fitted_coordinates: dict[int, float],
) -> float:
    """Single-candidate decisions must remain at the sample-fixed coordinate."""
    if selected_order == 1:
        return fixed_coordinate
    return fitted_coordinates.get(candidate_index, fixed_coordinate)


def select_order(
    endpoints: np.ndarray,
    candidate_indices: list[int],
    candidates: pd.DataFrame,
    strand: str,
    kernel: dict,
    threshold: float,
) -> tuple[list[int], list[dict]]:
    active = list(candidate_indices)
    history = []
    while len(active) > 1:
        p0 = np.asarray(
            [
                oriented_coordinate(
                    int(candidates.loc[index, "corrected_coordinate"]), strand
                )
                for index in active
            ]
        )
        full = fit_mixture(endpoints, p0, kernel)
        full_bic = bic(full, len(endpoints), len(active))
        reduced = []
        for drop in range(len(active)):
            fit = fit_mixture(endpoints, np.delete(p0, drop), kernel)
            reduced.append((bic(fit, len(endpoints), len(active) - 1), drop))
        best_reduced_bic, drop = min(reduced)
        improvement = best_reduced_bic - full_bic
        accepted = improvement >= threshold
        history.append(
            {
                "tested_order": len(active),
                "bic_improvement": improvement,
                "threshold": threshold,
                "accepted": accepted,
                "tested_drop_candidate_index": active[drop],
            }
        )
        if accepted:
            break
        active.pop(drop)
    return active, history


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--gate", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-reads", type=int, default=20)
    parser.add_argument("--search-radius", type=int, default=60)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    kernel = json.loads(args.kernel.read_text())
    frozen = json.loads(args.frozen.read_text())
    gate = json.loads(args.gate.read_text())
    fixed_offset = model_fixed_offset(kernel, frozen)
    threshold = float(gate["bic_improvement_threshold"])

    intervals = load_intervals(args.peaks, args.chromosome)
    candidate_path = args.output_dir / "candidates.tsv"
    candidates = prepare_candidates(
        intervals, int(kernel["read_length"]), candidate_path
    )
    base = inference(
        SimpleNamespace(
            bam=str(args.bam),
            chromosome=args.chromosome,
            candidates=candidate_path,
            kernel=args.kernel,
            frozen=args.frozen,
            search_radius=args.search_radius,
            min_reads=args.min_reads,
            mixture_mle_max_median_gap=300.0,
        )
    ).set_index("candidate_index", drop=False)

    _, positive = stream_positive_intervals(str(args.bam), args.chromosome)
    components = {
        strand: assign_candidates_indexed(candidates, positive[strand], strand)
        for strand in ("+", "-")
    }
    rows = []
    history_rows = []
    with pysam.AlignmentFile(args.bam, "rb") as bam:
        for strand in ("+", "-"):
            for component_index, component in enumerate(components[strand]):
                original = [
                    index
                    for index in component.candidate_indices
                    if index in base.index
                ]
                if not original:
                    continue
                endpoints = component_endpoints(
                    bam,
                    args.chromosome,
                    component.start,
                    component.end,
                    strand,
                )
                p0_all = np.asarray(
                    [
                        oriented_coordinate(
                            int(candidates.loc[index, "corrected_coordinate"]),
                            strand,
                        )
                        for index in original
                    ]
                )
                retained = np.any(
                    (endpoints[:, None] - p0_all[None, :] >= -220)
                    & (endpoints[:, None] - p0_all[None, :] <= 140),
                    axis=1,
                )
                endpoints = endpoints[retained]
                if len(original) == 1:
                    selected = original
                    history = []
                    decision = "fixed_single"
                    final_fit = None
                elif len(endpoints) < args.min_reads * len(original):
                    selected = original
                    history = []
                    decision = "fixed_multi_low_support"
                    final_fit = None
                else:
                    selected, history = select_order(
                        endpoints,
                        original,
                        candidates,
                        strand,
                        kernel,
                        threshold,
                    )
                    decision = (
                        "mixture_selected"
                        if len(selected) > 1
                        else "single_selected_by_bic"
                    )
                    p0_selected = np.asarray(
                        [
                            oriented_coordinate(
                                int(
                                    candidates.loc[
                                        index, "corrected_coordinate"
                                    ]
                                ),
                                strand,
                            )
                            for index in selected
                        ]
                    )
                    final_fit = (
                        fit_mixture(endpoints, p0_selected, kernel)
                        if len(selected) > 1
                        and use_mixture_mle(p0_selected, 300.0)
                        else None
                    )
                for item in history:
                    history_rows.append(
                        {
                            "component_index": component_index,
                            "strand": strand,
                            **item,
                        }
                    )
                fitted_coordinates = {}
                if final_fit is not None:
                    taus = final_fit["raw_taus"] - float(
                        kernel["endpoint_location_calibration"]
                    )
                    for offset, index in enumerate(selected):
                        p0 = int(candidates.loc[index, "corrected_coordinate"])
                        fitted_coordinates[index] = predicted_genomic(
                            p0, taus[offset], strand
                        )
                for index in original:
                    p0 = int(candidates.loc[index, "corrected_coordinate"])
                    fixed = predicted_genomic(p0, fixed_offset, strand)
                    rows.append(
                        {
                            "candidate_index": index,
                            "component_index": component_index,
                            "component_candidates": len(original),
                            "selected_order": len(selected),
                            "selected": index in selected,
                            "decision": decision,
                            "component_endpoint_reads": len(endpoints),
                            "fixed_coordinate": fixed,
                            "decision_coordinate": (
                                selected_coordinate(
                                    len(selected),
                                    index,
                                    fixed,
                                    fitted_coordinates,
                                )
                                if index in selected
                                else np.nan
                            ),
                        }
                    )

    decisions = pd.DataFrame(rows).merge(
        candidates.reset_index(names="candidate_index"),
        on="candidate_index",
        how="left",
    )
    decisions.to_csv(args.output_dir / "decision_calls.tsv", sep="\t", index=False)
    pd.DataFrame(history_rows).to_csv(
        args.output_dir / "bic_history.tsv", sep="\t", index=False
    )
    summary = {
        "input_candidates": int(len(candidates)),
        "evaluated_candidates": int(len(decisions)),
        "selected_candidates": int(decisions["selected"].sum()),
        "fixed_single_candidates": int(
            decisions["decision"].eq("fixed_single").sum()
        ),
        "mixture_selected_candidates": int(
            decisions["decision"].eq("mixture_selected").sum()
        ),
        "bic_removed_candidates": int((~decisions["selected"]).sum()),
        "fixed_offset": fixed_offset,
        "bic_improvement_threshold": threshold,
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
