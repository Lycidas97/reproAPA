import os
from typing import Optional

import pandas as pd
import pysam
from Bio import SeqIO


DEFAULT_MOTIFS = ("AATAAA", "ATTAAA", "TATAAA", "AGTAAA", "AAGAAA", "AAAAAG")
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1].upper()


def load_fasta_dict(fasta_path: str) -> dict[str, str]:
    with open(fasta_path, "r") as handle:
        return {record.id: str(record.seq).upper() for record in SeqIO.parse(handle, "fasta")}


def _canonical_reference_name(name: str) -> str:
    """Canonicalize only common chromosome prefixes for safe alias lookup."""
    value = str(name).strip().casefold()
    if value.startswith("chr"):
        value = value[3:]
    if value in {"m", "mt"}:
        return "mt"
    return value


def _build_reference_aliases(reference_names) -> dict[str, str]:
    """Return only aliases that resolve to exactly one reference sequence."""
    aliases: dict[str, str] = {}
    ambiguous: set[str] = set()
    for reference_name in reference_names:
        key = _canonical_reference_name(reference_name)
        if key in aliases and aliases[key] != reference_name:
            ambiguous.add(key)
        else:
            aliases[key] = reference_name
    for key in ambiguous:
        aliases.pop(key, None)
    return aliases


def _resolve_reference_name(chrom: str, references: set[str], aliases: dict[str, str]) -> Optional[str]:
    if chrom in references:
        return chrom
    return aliases.get(_canonical_reference_name(chrom))


def build_fasta_index(fasta_path: str, index_path: str) -> str:
    """Build a derived .fai without writing beside the read-only FASTA."""
    rows = []
    current_name = None
    current_length = 0
    sequence_offset = 0
    line_bases = 0
    line_width = 0

    def finish_record() -> None:
        if current_name is not None:
            rows.append((current_name, current_length, sequence_offset, line_bases, line_width))

    with open(fasta_path, "rb") as handle:
        while True:
            line = handle.readline()
            if not line:
                break
            if line.startswith(b">"):
                finish_record()
                current_name = line[1:].strip().split(maxsplit=1)[0].decode("utf-8")
                current_length = 0
                sequence_offset = handle.tell()
                line_bases = 0
                line_width = 0
                continue
            sequence_line = line.rstrip(b"\r\n")
            if not sequence_line:
                continue
            if current_name is None:
                raise ValueError(f"Invalid FASTA sequence before first header: {fasta_path}")
            if line_bases == 0:
                line_bases = len(sequence_line)
                line_width = len(line)
            current_length += len(sequence_line)
        finish_record()

    if not rows:
        raise ValueError(f"No FASTA records found: {fasta_path}")
    os.makedirs(os.path.dirname(os.path.abspath(index_path)), exist_ok=True)
    with open(index_path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")
    return index_path


def compute_sequence_qc(
    pas_df: pd.DataFrame,
    fasta_path: Optional[str] = None,
    genome: Optional[dict[str, str]] = None,
    fasta_index_path: Optional[str] = None,
    motifs: tuple[str, ...] = DEFAULT_MOTIFS,
    motif_window_upstream_bp: tuple[int, int] = (10, 50),
    motif_search_downstream_bp: int = 200,
    downstream_bp: int = 20,
    max_a_fraction: float = 0.8,
    max_a_run: int = 8,
) -> pd.DataFrame:
    """Compute motif and internal-priming QC without changing the PAS proxy coordinate.

    Motifs in the conventional 10--50 bp upstream interval are preferred. If
    none are present, a wider strand-oriented window extending downstream of the
    observed PAS proxy is searched to tolerate upstream-biased peak boundaries.
    Wider-window hits are labelled separately so callers need not treat them as
    equivalent to conventionally positioned cleavage evidence.
    """
    fasta = None
    references: set[str] = set()
    reference_aliases: dict[str, str] = {}
    if genome is None:
        if not fasta_path:
            return _not_run(pas_df, "fasta_not_provided")
        if fasta_index_path:
            if not os.path.exists(fasta_index_path):
                build_fasta_index(fasta_path, fasta_index_path)
            fasta = pysam.FastaFile(fasta_path, filepath_index=fasta_index_path)
            references = set(fasta.references)
        else:
            genome = load_fasta_dict(fasta_path)
    if fasta is None and genome is not None:
        references = set(genome)
    reference_aliases = _build_reference_aliases(references)

    rows = []
    upstream_min, upstream_max = sorted(int(value) for value in motif_window_upstream_bp)
    max_downstream = max(int(motif_search_downstream_bp), int(downstream_bp))
    try:
        for row in pas_df.itertuples(index=False):
            chrom = str(row.chrom)
            cut_site = int(row.cut_site)
            strand = str(row.strand)
            resolved_chrom = _resolve_reference_name(chrom, references, reference_aliases)
            if fasta is not None:
                if resolved_chrom is None:
                    rows.append(_row_not_run(row.pas_id, "chrom_not_found"))
                    continue
                chrom_length = int(fasta.get_reference_length(resolved_chrom))
                if strand == "+":
                    context_start = max(0, cut_site - upstream_max)
                    context_end = min(chrom_length, cut_site + max_downstream)
                else:
                    context_start = max(0, cut_site - max_downstream)
                    context_end = min(chrom_length, cut_site + upstream_max)
                seq = fasta.fetch(resolved_chrom, context_start, context_end).upper()
                local_cut_site = cut_site - context_start
            else:
                seq = genome.get(resolved_chrom) if resolved_chrom is not None else None
                if not seq:
                    rows.append(_row_not_run(row.pas_id, "chrom_not_found"))
                    continue
                local_cut_site = cut_site

            rows.append(
                _sequence_qc_row(
                    row.pas_id,
                    seq,
                    local_cut_site,
                    strand,
                    motifs,
                    upstream_min,
                    upstream_max,
                    int(motif_search_downstream_bp),
                    int(downstream_bp),
                    float(max_a_fraction),
                    int(max_a_run),
                )
            )
    finally:
        if fasta is not None:
            fasta.close()
    return pd.DataFrame(rows)


def _sequence_qc_row(
    pas_id: str,
    seq: str,
    cut_site: int,
    strand: str,
    motifs: tuple[str, ...],
    upstream_min: int,
    upstream_max: int,
    motif_search_downstream_bp: int,
    downstream_bp: int,
    max_a_fraction: float,
    max_a_run: int,
) -> dict:
        upstream = _oriented_window(
            seq,
            cut_site,
            strand,
            left=upstream_max,
            right=upstream_min,
            upstream=True,
        )
        motif_search, first_offset = _oriented_span(
            seq,
            cut_site,
            strand,
            upstream_bp=upstream_max,
            downstream_bp=motif_search_downstream_bp,
        )
        downstream = _oriented_window(
            seq,
            cut_site,
            strand,
            left=0,
            right=downstream_bp,
            upstream=False,
        )
        motif, motif_offset, motif_class = _best_motif_in_span(
            motif_search,
            motifs,
            first_offset=first_offset,
            expected_upstream_min=upstream_min,
            expected_upstream_max=upstream_max,
        )
        a_fraction = _base_fraction(downstream, "A")
        a_run = _max_run(downstream, "A")
        internal_flag = a_fraction >= max_a_fraction or a_run >= max_a_run
        return {
                "pas_id": pas_id,
                "upstream_motif_window": upstream,
                "best_pas_motif": motif,
                "best_pas_motif_distance": str(abs(motif_offset)) if motif_offset != "" else "",
                "best_pas_motif_offset_from_proxy": motif_offset,
                "best_pas_motif_search_class": motif_class,
                "canonical_motif_present": str(motif == "AATAAA").lower(),
                "variant_motif_present": str(bool(motif and motif != "AATAAA")).lower(),
                "downstream_a_fraction_10bp": _format_float(_base_fraction(downstream[:10], "A")),
                "downstream_a_fraction_20bp": _format_float(a_fraction),
                "downstream_a_run_max": a_run,
                "internal_priming_flag": str(internal_flag).lower(),
                "internal_priming_score": _format_float(max(a_fraction, a_run / max(downstream_bp, 1))),
                "sequence_qc_status": "computed",
            }


def _not_run(pas_df: pd.DataFrame, reason: str) -> pd.DataFrame:
    return pd.DataFrame([_row_not_run(pas_id, reason) for pas_id in pas_df.get("pas_id", [])])


def _row_not_run(pas_id: str, reason: str) -> dict:
    return {
        "pas_id": pas_id,
        "upstream_motif_window": "",
        "best_pas_motif": "",
        "best_pas_motif_distance": "",
        "best_pas_motif_offset_from_proxy": "",
        "best_pas_motif_search_class": "",
        "canonical_motif_present": "false",
        "variant_motif_present": "false",
        "downstream_a_fraction_10bp": "",
        "downstream_a_fraction_20bp": "",
        "downstream_a_run_max": "",
        "internal_priming_flag": "false",
        "internal_priming_score": "",
        "sequence_qc_status": reason,
    }


def _oriented_window(seq: str, cut_site: int, strand: str, left: int, right: int, upstream: bool) -> str:
    if strand == "+":
        if upstream:
            start, end = max(0, cut_site - left), max(0, cut_site - right)
        else:
            start, end = cut_site + left, cut_site + right
        return seq[start:end].upper()

    if upstream:
        start, end = cut_site + right, cut_site + left
    else:
        start, end = max(0, cut_site - right), max(0, cut_site - left)
    return reverse_complement(seq[start:end])


def _oriented_span(
    seq: str,
    cut_site: int,
    strand: str,
    upstream_bp: int,
    downstream_bp: int,
) -> tuple[str, int]:
    """Return a transcript-oriented span and its first offset from the PAS proxy."""
    upstream_bp = max(0, int(upstream_bp))
    downstream_bp = max(0, int(downstream_bp))
    if strand == "+":
        start = max(0, cut_site - upstream_bp)
        end = min(len(seq), cut_site + downstream_bp)
        return seq[start:end].upper(), start - cut_site

    start = max(0, cut_site - downstream_bp)
    end = min(len(seq), cut_site + upstream_bp)
    return reverse_complement(seq[start:end]), -(end - cut_site)


def _best_motif_in_span(
    seq: str,
    motifs: tuple[str, ...],
    first_offset: int,
    expected_upstream_min: int,
    expected_upstream_max: int,
) -> tuple[str, int | str, str]:
    candidates = []
    for motif_rank, motif in enumerate(motifs):
        search_from = 0
        while True:
            pos = seq.find(motif, search_from)
            if pos == -1:
                break
            offset = int(first_offset) + int(pos)
            expected = -expected_upstream_max <= offset <= -expected_upstream_min
            candidates.append(
                (
                    0 if expected else 1,
                    abs(offset),
                    motif_rank,
                    offset,
                    motif,
                    "expected_upstream" if expected else "expanded_proxy_window",
                )
            )
            search_from = pos + 1
    if not candidates:
        return "", "", ""
    _, _, _, offset, motif, motif_class = min(candidates)
    return motif, offset, motif_class


def _base_fraction(seq: str, base: str) -> float:
    if not seq:
        return 0.0
    return seq.upper().count(base) / len(seq)


def _max_run(seq: str, base: str) -> int:
    longest = 0
    current = 0
    for char in seq.upper():
        if char == base:
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def _format_float(value: float) -> str:
    return f"{value:.6g}"
