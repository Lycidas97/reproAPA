VALID_TRANSCRIPT_STRAND_MODES = {"same", "opposite"}
VALID_REGTOOLS_STRAND_MODES = {"RF", "FR", "XS"}


def validate_transcript_strand_mode(mode: str) -> str:
    """Validate how read orientation maps to transcript strand."""
    if mode not in VALID_TRANSCRIPT_STRAND_MODES:
        valid = ", ".join(sorted(VALID_TRANSCRIPT_STRAND_MODES))
        raise ValueError(f"Invalid transcript_strand_from_read '{mode}'. Expected one of: {valid}")
    return mode


def validate_regtools_strand_mode(mode: str) -> str:
    """Validate regtools junction extraction strand mode."""
    if mode not in VALID_REGTOOLS_STRAND_MODES:
        valid = ", ".join(sorted(VALID_REGTOOLS_STRAND_MODES))
        raise ValueError(f"Invalid regtools_strand_mode '{mode}'. Expected one of: {valid}")
    return mode


def read_transcript_strand(is_reverse: bool, transcript_strand_from_read: str = "same") -> str:
    """Map a BAM read orientation to transcript strand."""
    mode = validate_transcript_strand_mode(transcript_strand_from_read)
    if mode == "same":
        return "-" if is_reverse else "+"
    return "+" if is_reverse else "-"


def read_strand_labels(transcript_strand_from_read: str = "same") -> tuple[str, str]:
    """Return transcript strand labels for forward and reverse-aligned reads."""
    return (
        read_transcript_strand(False, transcript_strand_from_read),
        read_transcript_strand(True, transcript_strand_from_read),
    )
