"""Pure, shared neighboring-gene policy contract for annotation and projection."""

from __future__ import annotations

from collections.abc import Mapping
from math import isfinite
from typing import Any


TES_REGION_CLASS = "transcript_end_zone_pas"
CANONICAL_TERMINAL_REGION_CLASSES = {"annotated_3utr_pas", "terminal_exon_pas"}
CORE_TES_MAX_DISTANCE_BP = 250
TES2KB_MAX_DISTANCE_BP = 2000

CORE_TES_POLICY_CLASS = "core_tes"
QUALIFIED_TES2KB_POLICY_CLASS = "qualified_tes2kb_rescue"
INELIGIBLE_TES_POLICY_CLASS = "policy_ineligible_tes"
NOT_SELECTED_TES_POLICY_CLASS = "not_selected_tes"

SENSITIVITY_EXONIC_CONFLICT_STATE = "sensitivity_same_strand_exonic_gene_conflict"
CANONICAL_TERMINAL_OTHER_GENE_STATE = "canonical_terminal_other_gene_over_tes2kb"
CANONICAL_MULTI_GENE_REASON = "multiple_same_priority_gene_candidates"
CANONICAL_CDS_OVERLAP_REASON = "terminal_exon_cds_overlap_retained_for_sensitivity"
CORE_TES_EXONIC_AUDIT_STATE = "core_tes_same_strand_exonic_overlap_nonblocking"
CORE_TES_STATE = "core_tes_no_blocking_neighbor_conflict"
OUT_OF_RANGE_TES_STATE = "policy_ineligible_tes_out_of_range"
INCOMPLETE_TES_STATE = "policy_ineligible_tes_incomplete_provenance"
EXONIC_CONFLICT_DECISION_RULE = "same_strand_exonic_gene_conflict_to_sensitivity"
OTHER_EXON_OVER_TES2KB_STATE = "other_exon_over_tes2kb"
INTRON_OVER_TES2KB_STATE = "intron_over_tes2kb"
OTHER_EXON_REGION_CLASSES = {"coding_exon_pas", "internal_exon_pas"}
INTRON_REGION_CLASSES = {"terminal_intron_pas", "internal_intron_pas"}

NEIGHBOR_POLICY_REQUIRED_COLUMNS = (
    "strand",
    "primary_gene_id",
    "primary_candidate_id",
    "primary_transcript_id",
    "primary_transcript_strand",
    "primary_strand_relation",
    "primary_transcript_end",
    "downstream_distance_to_primary_transcript_end",
    "primary_region_class",
    "primary_assignment_status",
    "primary_assignment_reason",
    "primary_gene_biotype",
    "primary_transcript_biotype",
    "primary_model_transcript_end",
    "primary_model_3p_distance",
    "transcript_model_conflict",
    "transcript_context_class",
    "original_tes_candidate_ids",
    "original_tes_gene_ids",
    "original_tes_transcript_ids",
    "valid_same_strand_core_tes_candidate_ids",
    "valid_same_strand_core_tes_gene_ids",
    "valid_same_strand_core_tes_transcript_ids",
    "valid_same_strand_core_tes_candidate_count",
    "valid_same_strand_core_tes_gene_count",
    "qualified_same_strand_tes_candidate_ids",
    "qualified_same_strand_tes_gene_ids",
    "qualified_same_strand_tes_transcript_ids",
    "qualified_same_strand_tes_candidate_count",
    "qualified_same_strand_tes_gene_count",
    "same_strand_other_gene_intron_overlap",
    "same_strand_exonic_gene_conflict",
    "selected_tes_candidate_id",
    "selected_tes_policy_class",
    "selected_tes_policy_reason",
    "neighbor_gene_decision_rule",
    "neighbor_gene_terminal_state",
    "neighbor_gene_exclusion_reason",
    "displaced_tes_candidate_ids",
    "displaced_tes_gene_ids",
    "displaced_tes_transcript_ids",
    "structural_priority_feature_family",
    "feature_first_candidate_id",
    "feature_first_gene_id",
    "biotype_priority_candidate_id",
    "biotype_priority_gene_id",
    "biotype_assignment_conflict",
    "biotype_tie_break_status",
)


class NeighborPolicyContractError(ValueError):
    """Raised when neighboring-gene policy provenance is missing or contradictory."""


def classify_selected_tes(
    *,
    region_class: Any,
    candidate_id: Any,
    gene_id: Any,
    transcript_id: Any,
    transcript_end: Any,
    transcript_strand: Any,
    strand_relation: Any,
    pas_strand: Any,
    downstream_distance: Any,
) -> tuple[str, str]:
    """Classify a selected feature as core TES, TES2kb rescue, or ineligible."""
    if _text(region_class) != TES_REGION_CLASS:
        return NOT_SELECTED_TES_POLICY_CLASS, "primary_region_is_not_tes"
    required = {
        "candidate": candidate_id,
        "gene": gene_id,
        "transcript": transcript_id,
        "transcript_end": transcript_end,
        "transcript_strand": transcript_strand,
        "strand_relation": strand_relation,
        "pas_strand": pas_strand,
        "downstream_distance": downstream_distance,
    }
    missing = sorted(name for name, value in required.items() if not _text(value))
    if missing:
        return INELIGIBLE_TES_POLICY_CLASS, "missing_" + "_".join(missing)
    transcript_strand = _text(transcript_strand)
    pas_strand = _text(pas_strand)
    if (
        transcript_strand not in {"+", "-"}
        or pas_strand not in {"+", "-"}
        or _text(strand_relation) != "same_strand"
        or transcript_strand != pas_strand
    ):
        return INELIGIBLE_TES_POLICY_CLASS, "invalid_same_strand_tes_provenance"
    distance = _number(downstream_distance)
    if distance is None:
        return INELIGIBLE_TES_POLICY_CLASS, "invalid_downstream_distance"
    if distance < 0:
        return INELIGIBLE_TES_POLICY_CLASS, "negative_downstream_distance"
    if distance <= CORE_TES_MAX_DISTANCE_BP:
        return CORE_TES_POLICY_CLASS, "core_tes_at_or_below_250bp"
    if distance <= TES2KB_MAX_DISTANCE_BP:
        return QUALIFIED_TES2KB_POLICY_CLASS, "qualified_tes2kb_rescue"
    return INELIGIBLE_TES_POLICY_CLASS, "distance_above_2000bp"


def validate_neighbor_policy_row(row: Any) -> tuple[bool, str]:
    """Validate one finalized annotation/map row and derive main-atlas eligibility."""
    missing = [column for column in NEIGHBOR_POLICY_REQUIRED_COLUMNS if not _has(row, column)]
    if missing:
        raise NeighborPolicyContractError(
            f"neighbor provenance missing required columns: {missing}"
        )
    intron_flag = _text(
        _get(row, "same_strand_other_gene_intron_overlap")
    ).lower()
    if intron_flag not in {"true", "false"}:
        raise NeighborPolicyContractError(
            "same_strand_other_gene_intron_overlap must be true or false"
        )
    flag = _text(_get(row, "same_strand_exonic_gene_conflict")).lower()
    if flag not in {"true", "false"}:
        raise NeighborPolicyContractError(
            "same_strand_exonic_gene_conflict must be true or false"
        )
    transcript_conflict = _text(_get(row, "transcript_model_conflict")).lower()
    if transcript_conflict not in {"true", "false"}:
        raise NeighborPolicyContractError(
            "transcript_model_conflict must be true or false"
        )
    transcript_context = _text(_get(row, "transcript_context_class"))
    if not transcript_context:
        raise NeighborPolicyContractError("transcript_context_class is required")
    if (
        transcript_context == "alternative_terminal_intron"
        and transcript_conflict != "true"
    ):
        raise NeighborPolicyContractError(
            "alternative_terminal_intron requires transcript_model_conflict=true"
        )
    biotype_conflict = _text(_get(row, "biotype_assignment_conflict")).lower()
    if biotype_conflict not in {"true", "false"}:
        raise NeighborPolicyContractError(
            "biotype_assignment_conflict must be true or false"
        )
    if not _text(_get(row, "biotype_tie_break_status")):
        raise NeighborPolicyContractError("biotype_tie_break_status is required")
    terminal = _text(_get(row, "neighbor_gene_terminal_state"))
    rule = _text(_get(row, "neighbor_gene_decision_rule"))
    exclusion = _text(_get(row, "neighbor_gene_exclusion_reason"))
    region = _text(_get(row, "primary_region_class"))
    status = _text(_get(row, "primary_assignment_status"))
    assignment_reason = _text(_get(row, "primary_assignment_reason"))
    policy_class = _text(_get(row, "selected_tes_policy_class"))
    policy_reason = _text(_get(row, "selected_tes_policy_reason"))
    if not terminal or not rule or not policy_class or not policy_reason:
        raise NeighborPolicyContractError(
            "neighbor policy terminal state, rule, selected class, and reason are required"
        )

    expected_class, expected_reason = classify_selected_tes(
        region_class=region,
        candidate_id=_get(row, "primary_candidate_id"),
        gene_id=_get(row, "primary_gene_id"),
        transcript_id=_get(row, "primary_transcript_id"),
        transcript_end=_get(row, "primary_transcript_end"),
        transcript_strand=_get(row, "primary_transcript_strand"),
        strand_relation=_get(row, "primary_strand_relation"),
        pas_strand=_get(row, "strand"),
        downstream_distance=_get(row, "downstream_distance_to_primary_transcript_end"),
    )
    if policy_class != expected_class or policy_reason != expected_reason:
        raise NeighborPolicyContractError(
            "selected TES policy classification is inconsistent: "
            f"observed=({policy_class},{policy_reason}) "
            f"expected=({expected_class},{expected_reason})"
        )

    candidate_count = _nonnegative_integer(
        _get(row, "qualified_same_strand_tes_candidate_count"),
        "qualified_same_strand_tes_candidate_count",
    )
    gene_count = _nonnegative_integer(
        _get(row, "qualified_same_strand_tes_gene_count"),
        "qualified_same_strand_tes_gene_count",
    )
    qualified_candidate_ids = _values(
        _get(row, "qualified_same_strand_tes_candidate_ids")
    )
    qualified_gene_ids = _values(_get(row, "qualified_same_strand_tes_gene_ids"))
    qualified_transcript_ids = _values(
        _get(row, "qualified_same_strand_tes_transcript_ids")
    )
    core_candidate_count = _nonnegative_integer(
        _get(row, "valid_same_strand_core_tes_candidate_count"),
        "valid_same_strand_core_tes_candidate_count",
    )
    core_gene_count = _nonnegative_integer(
        _get(row, "valid_same_strand_core_tes_gene_count"),
        "valid_same_strand_core_tes_gene_count",
    )
    core_candidate_ids = _values(
        _get(row, "valid_same_strand_core_tes_candidate_ids")
    )
    core_gene_ids = _values(_get(row, "valid_same_strand_core_tes_gene_ids"))
    core_transcript_ids = _values(
        _get(row, "valid_same_strand_core_tes_transcript_ids")
    )
    if candidate_count != len(qualified_candidate_ids) or gene_count != len(qualified_gene_ids):
        raise NeighborPolicyContractError("qualified TES counts disagree with qualified IDs")
    if candidate_count and not qualified_transcript_ids:
        raise NeighborPolicyContractError("qualified TES transcript provenance is empty")
    if core_candidate_count != len(core_candidate_ids) or core_gene_count != len(
        core_gene_ids
    ):
        raise NeighborPolicyContractError("valid core TES counts disagree with core IDs")
    if core_candidate_count and not core_transcript_ids:
        raise NeighborPolicyContractError("valid core TES transcript provenance is empty")
    if policy_class == QUALIFIED_TES2KB_POLICY_CLASS:
        if intron_flag == "true":
            raise NeighborPolicyContractError(
                "qualified TES2kb rescue is blocked by same-strand other-gene intron overlap"
            )
        if core_candidate_count:
            raise NeighborPolicyContractError(
                "qualified TES2kb rescue was selected despite an available valid core TES"
            )
        if candidate_count < 1 or gene_count < 1:
            raise NeighborPolicyContractError("qualified TES2kb selection lacks qualified candidates")
        if _text(_get(row, "selected_tes_candidate_id")) != _text(
            _get(row, "primary_candidate_id")
        ):
            raise NeighborPolicyContractError("selected qualified TES candidate is inconsistent")
        if (
            _text(_get(row, "primary_candidate_id")) not in qualified_candidate_ids
            or _text(_get(row, "primary_gene_id")) not in qualified_gene_ids
            or _text(_get(row, "primary_transcript_id"))
            not in qualified_transcript_ids
        ):
            raise NeighborPolicyContractError(
                "selected qualified TES candidate/gene/transcript is absent from its tier"
            )
    elif policy_class == CORE_TES_POLICY_CLASS:
        if core_candidate_count < 1 or core_gene_count < 1:
            raise NeighborPolicyContractError("core TES selection lacks valid core candidates")
        if _text(_get(row, "selected_tes_candidate_id")) != _text(
            _get(row, "primary_candidate_id")
        ):
            raise NeighborPolicyContractError("selected core TES candidate is inconsistent")
        if (
            _text(_get(row, "primary_candidate_id")) not in core_candidate_ids
            or _text(_get(row, "primary_gene_id")) not in core_gene_ids
            or _text(_get(row, "primary_transcript_id")) not in core_transcript_ids
        ):
            raise NeighborPolicyContractError(
                "selected core TES candidate/gene/transcript is absent from its tier"
            )
    elif policy_class == INELIGIBLE_TES_POLICY_CLASS:
        if core_candidate_count:
            raise NeighborPolicyContractError(
                "policy-ineligible TES was selected despite an available valid core TES"
            )
        if candidate_count or gene_count:
            raise NeighborPolicyContractError("policy-ineligible TES must not be counted as qualified")
    canonical_attribution_selected = region in CANONICAL_TERMINAL_REGION_CLASSES or (
        rule == "region_priority_assignment"
        and policy_class == NOT_SELECTED_TES_POLICY_CLASS
        and (
            (
                region == "multi_gene_ambiguous_pas"
                and status == "ambiguous_multi_gene"
                and assignment_reason == CANONICAL_MULTI_GENE_REASON
            )
            or (
                region == "terminal_exon_cds_overlap_pas"
                and status == "assigned_sensitivity"
                and assignment_reason == CANONICAL_CDS_OVERLAP_REASON
            )
        )
    )
    if (
        core_candidate_count
        and policy_class == NOT_SELECTED_TES_POLICY_CLASS
        and not canonical_attribution_selected
    ):
        raise NeighborPolicyContractError(
            "available valid core TES was not selected after the canonical tier"
        )

    if terminal == CANONICAL_TERMINAL_OTHER_GENE_STATE:
        primary_gene = _text(_get(row, "primary_gene_id"))
        if (
            region not in CANONICAL_TERMINAL_REGION_CLASSES
            or status != "assigned_primary"
            or rule != CANONICAL_TERMINAL_OTHER_GENE_STATE
            or policy_class != NOT_SELECTED_TES_POLICY_CLASS
            or candidate_count < 1
            or not any(gene_id != primary_gene for gene_id in qualified_gene_ids)
        ):
            raise NeighborPolicyContractError(
                "canonical terminal state lacks terminal/assigned-primary/qualified-other-gene TES evidence"
            )
        return True, ""

    if terminal in {OTHER_EXON_OVER_TES2KB_STATE, INTRON_OVER_TES2KB_STATE}:
        is_exon = terminal == OTHER_EXON_OVER_TES2KB_STATE
        expected_regions = OTHER_EXON_REGION_CLASSES if is_exon else INTRON_REGION_CLASSES
        expected_rule = (
            "same_strand_other_exon_over_qualified_tes2kb"
            if is_exon
            else "same_strand_intron_over_qualified_tes2kb"
        )
        displaced_candidate_ids = _values(
            _get(row, "displaced_tes_candidate_ids")
        )
        displaced_gene_ids = _values(_get(row, "displaced_tes_gene_ids"))
        displaced_transcript_ids = _values(
            _get(row, "displaced_tes_transcript_ids")
        )
        if (
            region not in expected_regions
            or status not in {"assigned_sensitivity", "ambiguous_multi_gene"}
            or policy_class != NOT_SELECTED_TES_POLICY_CLASS
            or rule != expected_rule
            or exclusion != terminal
            or not displaced_candidate_ids
            or not displaced_gene_ids
            or not displaced_transcript_ids
            or not qualified_candidate_ids
            or not displaced_candidate_ids.issubset(qualified_candidate_ids)
            or not displaced_gene_ids.issubset(qualified_gene_ids)
            or not displaced_transcript_ids.issubset(qualified_transcript_ids)
            or _text(_get(row, "selected_tes_candidate_id"))
        ):
            raise NeighborPolicyContractError(
                f"{terminal} state contradicts structural/displaced-TES provenance"
            )
        expected_family = "other_exon" if is_exon else "intron"
        if _text(_get(row, "structural_priority_feature_family")) != expected_family:
            raise NeighborPolicyContractError(
                f"{terminal} has an inconsistent structural feature family"
            )
        return False, terminal

    if terminal == "ambiguous_multi_gene_structural":
        if (
            region not in (OTHER_EXON_REGION_CLASSES | INTRON_REGION_CLASSES)
            or status != "ambiguous_multi_gene"
            or rule != "select_structural_feature_biotype_distance_exact_tie"
            or exclusion != "exact_feature_biotype_3p_distance_tie_across_genes"
            or _text(_get(row, "biotype_tie_break_status"))
            != "exact_cross_gene_tie"
        ):
            raise NeighborPolicyContractError(
                "ambiguous structural state contradicts tie provenance"
            )
        return False, "ambiguous_multi_gene_structural"

    if terminal == SENSITIVITY_EXONIC_CONFLICT_STATE:
        distance = _number(_get(row, "downstream_distance_to_primary_transcript_end"))
        if (
            region != TES_REGION_CLASS
            or status != "assigned_sensitivity"
            or policy_class != QUALIFIED_TES2KB_POLICY_CLASS
            or distance is None
            or not (CORE_TES_MAX_DISTANCE_BP < distance <= TES2KB_MAX_DISTANCE_BP)
            or flag != "true"
            or rule != EXONIC_CONFLICT_DECISION_RULE
            or exclusion != "same_strand_exonic_gene_conflict"
        ):
            raise NeighborPolicyContractError(
                "sensitivity exonic-conflict state contradicts region/status/distance/flag/rule"
            )
        return False, "same_strand_exonic_gene_conflict"

    if terminal in {CORE_TES_STATE, CORE_TES_EXONIC_AUDIT_STATE}:
        expected_rule = (
            "core_tes_exonic_overlap_audit_only"
            if terminal == CORE_TES_EXONIC_AUDIT_STATE
            else "core_tes_assignment_preserved"
        )
        if (
            region != TES_REGION_CLASS
            or status != "assigned_primary"
            or policy_class != CORE_TES_POLICY_CLASS
            or rule != expected_rule
            or (terminal == CORE_TES_EXONIC_AUDIT_STATE and flag != "true")
            or (terminal == CORE_TES_STATE and flag != "false")
            or exclusion
        ):
            raise NeighborPolicyContractError("core TES terminal state is inconsistent")
        return True, ""

    if terminal == OUT_OF_RANGE_TES_STATE:
        distance = _number(_get(row, "downstream_distance_to_primary_transcript_end"))
        if (
            region != TES_REGION_CLASS
            or status != "assigned_sensitivity"
            or policy_class != INELIGIBLE_TES_POLICY_CLASS
            or distance is None
            or distance <= TES2KB_MAX_DISTANCE_BP
            or rule != "tes_candidate_outside_rescue_distance"
            or exclusion != "tes_candidate_outside_rescue_distance"
        ):
            raise NeighborPolicyContractError("out-of-range TES terminal state is inconsistent")
        return False, "tes_candidate_outside_rescue_distance"

    if terminal == INCOMPLETE_TES_STATE:
        if (
            region != TES_REGION_CLASS
            or status != "assigned_sensitivity"
            or policy_class != INELIGIBLE_TES_POLICY_CLASS
            or rule != "tes_candidate_incomplete_provenance"
            or exclusion != "tes_candidate_incomplete_provenance"
        ):
            raise NeighborPolicyContractError("incomplete TES terminal state is inconsistent")
        return False, "tes_candidate_incomplete_provenance"

    if terminal.startswith("tes2kb_rescue_"):
        if region != TES_REGION_CLASS or policy_class != QUALIFIED_TES2KB_POLICY_CLASS:
            raise NeighborPolicyContractError("TES2kb rescue state lacks qualified TES2kb selection")
        if flag == "true":
            raise NeighborPolicyContractError(
                "qualified TES2kb exonic conflict lacks the sensitivity terminal state"
            )
        if intron_flag == "true":
            raise NeighborPolicyContractError(
                "qualified TES2kb intron conflict lacks the intron-over-TES2kb terminal state"
            )
        return True, ""

    if terminal == "ambiguous_multi_gene_tes":
        tier_has_exact_multi_gene_tie = (
            policy_class == CORE_TES_POLICY_CLASS
            and core_candidate_count >= 2
            and core_gene_count >= 2
        ) or (
            policy_class == QUALIFIED_TES2KB_POLICY_CLASS
            and candidate_count >= 2
            and gene_count >= 2
        )
        if (
            status != "ambiguous_multi_gene_tes"
            or rule != "select_nearest_same_strand_tes_exact_multi_gene_tie"
            or exclusion != "exact_minimum_tes_distance_tied_across_genes"
            or not tier_has_exact_multi_gene_tie
        ):
            raise NeighborPolicyContractError("ambiguous multi-gene TES state is inconsistent")
        return False, "ambiguous_multi_gene_tes"

    if policy_class == INELIGIBLE_TES_POLICY_CLASS:
        raise NeighborPolicyContractError("policy-ineligible TES lacks an explicit terminal state")
    if policy_class == CORE_TES_POLICY_CLASS:
        raise NeighborPolicyContractError("core TES lacks an explicit core terminal state")
    if policy_class == QUALIFIED_TES2KB_POLICY_CLASS:
        raise NeighborPolicyContractError("qualified TES2kb lacks an explicit TES2kb terminal state")
    return True, ""


def _has(row: Any, field: str) -> bool:
    return field in row if isinstance(row, Mapping) else hasattr(row, field)


def _get(row: Any, field: str) -> Any:
    return row.get(field, "") if isinstance(row, Mapping) else getattr(row, field, "")


def _text(value: Any) -> str:
    if value is None:
        return ""
    return str(value).strip()


def _number(value: Any) -> float | None:
    try:
        result = float(_text(value))
    except (TypeError, ValueError):
        return None
    return result if isfinite(result) else None


def _nonnegative_integer(value: Any, field: str) -> int:
    number = _number(value)
    if number is None or number < 0 or not number.is_integer():
        raise NeighborPolicyContractError(f"{field} must be a non-negative integer")
    return int(number)


def _values(value: Any) -> set[str]:
    return {item for item in _text(value).split("|") if item}
