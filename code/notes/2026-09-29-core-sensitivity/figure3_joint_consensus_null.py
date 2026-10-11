"""Test count-rule consensus recurrence with a joint matched null.

For each study and gene, the six count-mode event labels are permuted as one
block.  This preserves within-study cross-mode dependence, each mode's event
counts, study-specific measurability, and the union of the mode-specific
matching strata.  The primary statistic counts genes recurrent in the same
direction in every count mode.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
SOURCE_SCRIPT = HERE / "figure3_recount_recurrence.py"
INPUT = HERE / "figure3-recount-recurrence/study_events.tsv.gz"
OUT = HERE / "figure3-joint-consensus-null"
N_PERM = 10_000
SEED = 20260930


def load_source():
    spec = importlib.util.spec_from_file_location("figure3_recount", SOURCE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(values: pd.Series) -> np.ndarray:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.minimum(ranked, 1)
    return out


def consensus_mask(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return consensus, all-mode recurrence, and per-mode directions.

    Matrix dimensions are gene, study, mode.  Event value 9 denotes that a
    gene is not measurable in a study.  Direction ties follow the upstream
    recurrence code and are assigned positive.
    """
    positive = (matrix == 1).sum(axis=1)
    negative = (matrix == -1).sum(axis=1)
    recurrent = np.maximum(positive, negative) >= 2
    direction_positive = positive >= negative
    all_mode_recurrent = recurrent.all(axis=1)
    same_direction = np.logical_or(
        direction_positive.all(axis=1), (~direction_positive).all(axis=1)
    )
    return all_mode_recurrent & same_direction, all_mode_recurrent, direction_positive


def build_configuration(
    frame: pd.DataFrame, modes: list[str], source, old
) -> tuple[np.ndarray, pd.DataFrame, list[str], dict[int, list[np.ndarray]], pd.DataFrame]:
    species = frame.species.iloc[0]
    threshold = float(frame.threshold.iloc[0])
    assert frame.species.eq(species).all()
    assert frame.threshold.eq(threshold).all()

    mode_parts = []
    for mode in modes:
        current = frame[frame["mode"].eq(mode)].copy()
        assert len(current), f"No rows for mode {mode}"
        stratified = source.stratify(current, old)
        mode_parts.append(
            stratified[["unit", "gene_id", "stratum"]].rename(
                columns={"stratum": f"stratum__{mode}"}
            )
        )

    joint = mode_parts[0]
    for part in mode_parts[1:]:
        joint = joint.merge(part, on=["unit", "gene_id"], validate="one_to_one")
    stratum_columns = [f"stratum__{mode}" for mode in modes]
    joint["joint_stratum"] = joint[stratum_columns].astype(str).agg("||".join, axis=1)

    event_counts = frame.groupby(["unit", "gene_id"])["mode"].nunique()
    if not event_counts.eq(len(modes)).all():
        failed = event_counts[event_counts.ne(len(modes))].head().to_dict()
        raise AssertionError(f"Incomplete count-mode event vectors: {failed}")

    symbols = (
        frame[["gene_id", "gene_symbol"]]
        .drop_duplicates()
        .sort_values("gene_id")
        .reset_index(drop=True)
    )
    if symbols.gene_id.duplicated().any():
        raise AssertionError("A gene_id maps to multiple gene symbols")
    units = sorted(frame.unit.unique())
    gene_index = {gene: index for index, gene in enumerate(symbols.gene_id)}
    unit_index = {unit: index for index, unit in enumerate(units)}
    mode_index = {mode: index for index, mode in enumerate(modes)}
    matrix = np.full((len(symbols), len(units), len(modes)), 9, dtype=np.int8)
    for row in frame[["unit", "gene_id", "mode", "event"]].itertuples(index=False):
        matrix[
            gene_index[row.gene_id], unit_index[row.unit], mode_index[row.mode]
        ] = int(row.event)

    strata_groups: dict[int, list[np.ndarray]] = {}
    diagnostics = []
    for unit, current in joint.groupby("unit"):
        groups = []
        shuffleable = 0
        for _, group in current.groupby("joint_stratum"):
            indices = np.asarray([gene_index[gene] for gene in group.gene_id])
            groups.append(indices)
            if len(indices) > 1:
                shuffleable += len(indices)
        strata_groups[unit_index[unit]] = groups
        diagnostics.append(
            {
                "species": species,
                "threshold": threshold,
                "unit": unit,
                "eligible_gene_vectors": len(current),
                "joint_strata": current.joint_stratum.nunique(),
                "singleton_strata": int(
                    (current.groupby("joint_stratum").size() == 1).sum()
                ),
                "shuffleable_gene_vectors": shuffleable,
                "shuffleable_fraction": shuffleable / len(current),
            }
        )
    return matrix, symbols, units, strata_groups, pd.DataFrame(diagnostics)


def permute_joint(
    matrix: np.ndarray,
    strata_groups: dict[int, list[np.ndarray]],
    rng: np.random.Generator,
) -> np.ndarray:
    genes, _, modes = matrix.shape
    positive = np.zeros((N_PERM, genes, modes), dtype=np.int8)
    negative = np.zeros((N_PERM, genes, modes), dtype=np.int8)
    for unit_index, groups in strata_groups.items():
        for indices in groups:
            labels = matrix[indices, unit_index, :]
            if len(indices) == 1:
                shuffled = np.repeat(labels[None, :, :], N_PERM, axis=0)
            else:
                order = np.argsort(rng.random((N_PERM, len(indices))), axis=1)
                shuffled = labels[order, :]
            positive[:, indices, :] += shuffled == 1
            negative[:, indices, :] += shuffled == -1

    recurrent = np.maximum(positive, negative) >= 2
    direction_positive = positive >= negative
    consensus = recurrent.all(axis=2) & np.logical_or(
        direction_positive.all(axis=2), (~direction_positive).all(axis=2)
    )
    return consensus.sum(axis=1)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    assert INPUT.is_file(), f"Missing study-event input: {INPUT}"
    source = load_source()
    old = source.load_old()
    events = pd.read_csv(INPUT, sep="\t")
    modes = list(source.MODES)
    views = list(source.VIEWS)
    assert set(events["mode"]) == set(modes)
    assert set(events["view"]) == set(views)

    summary_rows = []
    member_parts = []
    diagnostic_parts = []
    null_archive = {}
    for (view, species, threshold), frame in events.groupby(
        ["view", "species", "threshold"], sort=True
    ):
        matrix, symbols, units, groups, diagnostics = build_configuration(
            frame, modes, source, old
        )
        observed_mask, all_mode_mask, direction_positive = consensus_mask(matrix)
        observed = int(observed_mask.sum())
        key = f"{view}__{species}__{threshold:.2f}"
        seed = SEED + int(hashlib.sha256(key.encode()).hexdigest()[:8], 16)
        null = permute_joint(matrix, groups, np.random.default_rng(seed))
        null_archive[key] = null
        summary_rows.append(
            {
                "view": view,
                "species": species,
                "contrast": "P_S",
                "threshold": threshold,
                "modes": len(modes),
                "studies": len(units),
                "observed_consensus_same_direction": observed,
                "observed_all_mode_recurrent_any_direction": int(all_mode_mask.sum()),
                "null_mean": float(null.mean()),
                "null_median": float(np.median(null)),
                "null_p025": float(np.quantile(null, 0.025)),
                "null_p975": float(np.quantile(null, 0.975)),
                "enrichment_vs_median": (
                    float(observed / np.median(null)) if np.median(null) > 0 else np.nan
                ),
                "empirical_p": float((1 + (null >= observed).sum()) / (N_PERM + 1)),
                "seed": seed,
                "permutations": N_PERM,
            }
        )

        positive = (matrix == 1).sum(axis=1)
        negative = (matrix == -1).sum(axis=1)
        members = symbols.loc[observed_mask].copy()
        members["direction"] = np.where(
            direction_positive[observed_mask, 0], "positive", "negative"
        )
        members["view"] = view
        members["species"] = species
        members["threshold"] = threshold
        for mode_index, mode in enumerate(modes):
            members[f"{mode}_support"] = np.maximum(
                positive[observed_mask, mode_index], negative[observed_mask, mode_index]
            )
        member_parts.append(members)
        diagnostics["view"] = view
        diagnostics["modes"] = len(modes)
        diagnostic_parts.append(diagnostics)

    summary = pd.DataFrame(summary_rows)
    summary["q_within_view_four_tests"] = np.nan
    for _, indexes in summary.groupby("view").groups.items():
        summary.loc[indexes, "q_within_view_four_tests"] = bh(
            summary.loc[indexes, "empirical_p"]
        )
    summary.to_csv(OUT / "joint_consensus_summary.tsv", sep="\t", index=False)
    pd.concat(member_parts, ignore_index=True).to_csv(
        OUT / "joint_consensus_members.tsv.gz", sep="\t", index=False
    )
    diagnostics = pd.concat(diagnostic_parts, ignore_index=True)
    diagnostics.to_csv(OUT / "joint_strata_diagnostics.tsv", sep="\t", index=False)
    np.savez_compressed(OUT / "joint_consensus_null_draws.npz", **null_archive)

    upstream_audit_path = INPUT.parent / "audit.json"
    upstream_audit = (
        json.loads(upstream_audit_path.read_text()) if upstream_audit_path.is_file() else {}
    )
    audit = {
        "status": (
            "PROVISIONAL_EXCLUDED_STUDY_ROUTE"
            if upstream_audit.get("status") == "PROVISIONAL_EXCLUDED_STUDY_ROUTE"
            else "PASS"
        ),
        "primary_statistic": (
            "genes recurrent in the same direction in every count mode"
        ),
        "permutation_unit": "six-mode event-label vector within study and joint stratum",
        "preserved": [
            "within-study cross-mode event dependence",
            "mode-specific positive, negative, and subthreshold event counts",
            "study-specific gene measurability",
            "eligible-study count, expression, APA coverage, and PAS-number strata",
            "mode-specific reference-WPP bins and directional capacities",
        ],
        "views": views,
        "modes": modes,
        "thresholds": sorted(events.threshold.unique().tolist()),
        "permutations_per_test": N_PERM,
        "tests": len(summary),
        "minimum_shuffleable_fraction": float(diagnostics.shuffleable_fraction.min()),
        "upstream_status": upstream_audit.get("status", "UNKNOWN"),
        "upstream_excluded_studies_with_incomplete_routes": upstream_audit.get(
            "excluded_studies_with_incomplete_routes", []
        ),
        "input_sha256": {
            str(INPUT): digest(INPUT),
            str(SOURCE_SCRIPT): digest(SOURCE_SCRIPT),
        },
    }
    (OUT / "audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--out", type=Path, default=OUT)
    parser.add_argument("--permutations", type=int, default=N_PERM)
    args = parser.parse_args()
    INPUT = args.input
    OUT = args.out
    N_PERM = args.permutations
    main()
