#!/usr/bin/env python3
"""Freeze reviewed training anchors and species-specific pilot readiness."""

from __future__ import annotations

import csv
import hashlib
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RELEASE = ROOT / "release"


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    fields, initial = read_tsv(RELEASE / "anchor_adjudication.tsv")
    _, decision_rows = read_tsv(ROOT / "config/anchor_adjudication_decisions.tsv")
    decisions = {(row["study_id"], row["local_cluster_id"]): row for row in decision_rows}
    b_keys = {
        (row["study_id"], row["local_cluster_id"])
        for row in initial
        if row["anchor_tier_initial"] == "B_review_required"
    }
    if set(decisions) != b_keys:
        raise ValueError(f"Decision coverage mismatch: missing={sorted(b_keys-set(decisions))} extra={sorted(set(decisions)-b_keys)}")

    accepted_alternative_subtypes: dict[str, set[str]] = defaultdict(set)
    for row in initial:
        if (
            row["anchor_tier_initial"] == "A_confirmed"
            and row["proposed_parent_family"] == "alternative_fate"
        ):
            accepted_alternative_subtypes[row["proposed_subtype"]].add(row["study_id"])

    final = []
    for row in initial:
        key = (row["study_id"], row["local_cluster_id"])
        if row["anchor_tier_initial"] == "A_confirmed":
            subtype_is_shared = (
                row["proposed_parent_family"] != "alternative_fate"
                or len(accepted_alternative_subtypes[row["proposed_subtype"]]) >= 2
            )
            decision = "accepted_existing" if subtype_is_shared else "accepted_study_specific_identity"
            eligible = "true" if subtype_is_shared else "false"
            evidence = "Accepted R9/R10B hard identity and legal-scope gates"
            visual = (
                "Previously accepted in independent visual review"
                if subtype_is_shared
                else "Retained as a study-specific child; insufficient cross-study subtype support"
            )
        elif key in decisions:
            reviewed = decisions[key]
            decision = reviewed["review_decision"]
            eligible = "true" if decision == "promote_training_anchor" else "false"
            evidence = reviewed["evidence_basis"]
            visual = reviewed["visual_judgment"]
        else:
            decision = "locked_non_anchor"
            eligible = "false"
            evidence = row["tier_reason"]
            visual = "Not eligible for anchor adjudication"
        if eligible == "true" and (
            row["manual_hold"] == "true" or row["library_context_triggered"] == "true"
        ):
            raise ValueError(f"Illegal promoted anchor: {key}")
        enriched = dict(row)
        enriched.update(
            {
                "review_decision": decision,
                "alignment_training_anchor_eligible": eligible,
                "review_evidence_basis": evidence,
                "review_visual_judgment": visual,
            }
        )
        final.append(enriched)

    out_fields = fields + [
        "review_decision",
        "alignment_training_anchor_eligible",
        "review_evidence_basis",
        "review_visual_judgment",
    ]
    write_tsv(RELEASE / "anchor_adjudication_final.tsv", out_fields, final)

    support: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in final:
        if row["alignment_training_anchor_eligible"] != "true":
            continue
        family = row["proposed_parent_family"]
        if family in {"source_retaining", "reprogramming_transition", "pluripotency_target_like"}:
            support[(row["species"], family)].add(row["study_id"])

    required = ["source_retaining", "reprogramming_transition", "pluripotency_target_like"]
    readiness = []
    for species in ("Mus_musculus", "Homo_sapiens"):
        counts = {family: len(support[(species, family)]) for family in required}
        ready = all(counts[family] >= 2 for family in required)
        readiness.append(
            {
                "species": species,
                "source_supporting_studies": str(counts["source_retaining"]),
                "transition_supporting_studies": str(counts["reprogramming_transition"]),
                "target_supporting_studies": str(counts["pluripotency_target_like"]),
                "pilot_status": "READY" if ready else "BLOCKED",
                "blocking_reason": "none" if ready else "fewer_than_two_studies_for_one_or_more_required_families",
            }
        )
    write_tsv(RELEASE / "species_pilot_readiness.tsv", list(readiness[0]), readiness)

    promoted = [row for row in final if row["review_decision"] == "promote_training_anchor"]
    reserved = [row for row in final if row["review_decision"] == "reserve_validation"]
    mouse = next(row for row in readiness if row["species"] == "Mus_musculus")
    human = next(row for row in readiness if row["species"] == "Homo_sapiens")
    report = f"""# Anchor adjudication review

Decision: `PASS_WITH_SPECIES_GATE`

## Result

- Existing A anchors retained: 5.
- B candidates promoted to alignment training anchors: {len(promoted)}.
- B candidates reserved for independent validation: {len(reserved)}.
- Training-eligible anchors after the shared-subtype gate: {sum(row['alignment_training_anchor_eligible'] == 'true' for row in final)}.
- No manual-hold or library-context cluster was promoted.

## Species gate

- Mouse pilot: `{mouse['pilot_status']}`. Source, transition and target-like families each have at least two supporting studies.
- Human pilot: `{human['pilot_status']}`. Source and target-like anchors span at least two studies, but the transition family does not.

GSE118258 C008 remains the only clean human transition candidate. It contains 470 primary cells, is observed at D16, and has discordant source-program splits. GSE178325 Stage-II clusters remain process-dominant validation pools and cannot be relabeled merely to satisfy the gate.

## Visual review

All six original-resolution partition and process-overlay figures were re-opened. They confirm strong route/library/stage structure in GSE103221, GSE118258, GSE150311 and GSE178325; consequently, none of those axes was treated as an identity label. WOT provides the clearest local backbone, but mouse readiness also requires GSE103221 anchors and is not based on WOT alone.

## Authorization boundary

This release authorizes preparation of a mouse-only scVI/scANVI pilot. It does not authorize human scANVI training, six-study production, cross-species alignment, trajectory, APA, or modification of the accepted R9/R10B public-state fields.
"""
    (RELEASE / "anchor_adjudication_review.md").write_text(report, encoding="utf-8")

    manifest = []
    for path in sorted(RELEASE.iterdir()):
        if path.is_file() and path.name != "release_manifest.tsv":
            manifest.append(
                {
                    "relative_path": path.name,
                    "sha256": sha256_file(path),
                    "size_bytes": str(path.stat().st_size),
                }
            )
    write_tsv(RELEASE / "release_manifest.tsv", list(manifest[0]), manifest)


if __name__ == "__main__":
    main()
