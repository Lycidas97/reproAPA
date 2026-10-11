#!/usr/bin/env python3
"""Localize PAS candidates, then quantify barcode/UMI molecules by hard intervals."""

from __future__ import annotations

import argparse
import bisect
import gzip
import itertools
import json
import os
import resource
import sys
import time
from array import array
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

_blas_threads = os.environ.get("APADETECT_BLAS_THREADS", "1")
os.environ["OPENBLAS_NUM_THREADS"] = _blas_threads
os.environ["OMP_NUM_THREADS"] = _blas_threads
os.environ["MKL_NUM_THREADS"] = _blas_threads

import numpy as np
import pandas as pd
import pysam
from scipy.io import mmwrite
from scipy.sparse import coo_matrix

STAGE = Path(__file__).resolve().parent
STAGES = STAGE.parent
sys.path.insert(0, str(STAGES / "62_multi_endpoint_mixture"))
from multi_endpoint_mixture import inference, model_fixed_offset  # noqa: E402

BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


class CompactMolecules:
    """Store unique molecules without allocating one Python set per matrix cell."""

    def __init__(self, peak_names):
        self.peaks = list(dict.fromkeys(peak_names))
        self.peak_ids = {name: index for index, name in enumerate(self.peaks)}
        self.barcodes = []
        self.barcode_ids = {}
        self.molecules = set()

    def add(self, peak_name: str, barcode: str, umi: str) -> None:
        barcode_id = self.barcode_ids.get(barcode)
        if barcode_id is None:
            barcode_id = len(self.barcodes)
            self.barcode_ids[barcode] = barcode_id
            self.barcodes.append(barcode)
        self.molecules.add((self.peak_ids[peak_name], barcode_id, umi))


def load_intervals(path: Path, chromosome: str | None = None) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", header=None, names=BED_COLUMNS, usecols=range(6))
    if chromosome is not None:
        requested = str(chromosome).removeprefix("chr")
        normalized = table["chrom"].astype(str).str.removeprefix("chr")
        table = table[normalized.eq(requested)]
    table = table.astype({"start": int, "end": int, "name": str, "strand": str})
    table = table.sort_values(["chrom", "strand", "start", "end", "name"]).reset_index(drop=True)
    for (_, _), group in table.groupby(["chrom", "strand"], sort=False):
        if np.any(group["start"].to_numpy()[1:] < group["end"].to_numpy()[:-1]):
            raise ValueError("Quantification intervals must not overlap within a chromosome and strand")
    return table


def prepare_candidates(intervals: pd.DataFrame, read_length: int, output: Path) -> pd.DataFrame:
    candidates = intervals.copy()
    candidates["candidate_coordinate"] = pd.to_numeric(candidates["score"]).astype(int)
    half = read_length // 2
    candidates["corrected_coordinate"] = candidates["candidate_coordinate"] + np.where(
        candidates["strand"].eq("+"), -half, half
    )
    candidates["te_id"] = candidates["name"]
    candidates.to_csv(output, sep="\t", index=False)
    return candidates


def localize(
    bam: Path,
    chromosome: str,
    intervals: pd.DataFrame,
    kernel_path: Path,
    frozen_path: Path,
    output_dir: Path,
    min_reads: int,
) -> pd.DataFrame:
    kernel = json.loads(kernel_path.read_text())
    frozen = json.loads(frozen_path.read_text())
    candidate_path = output_dir / "candidates.tsv"
    candidates = prepare_candidates(intervals, int(kernel["read_length"]), candidate_path)
    args = SimpleNamespace(
        bam=str(bam),
        chromosome=chromosome,
        candidates=candidate_path,
        kernel=kernel_path,
        frozen=frozen_path,
        search_radius=60,
        min_reads=min_reads,
        mixture_mle_max_median_gap=300.0,
    )
    fitted = inference(args)
    keep = ["candidate_index", "model_coordinate", "component_model", "mixture_refined"]
    fitted = fitted[keep] if not fitted.empty else pd.DataFrame(columns=keep)
    calls = candidates.reset_index(names="candidate_index").merge(
        fitted, on="candidate_index", how="left"
    )
    fixed_offset = model_fixed_offset(kernel, frozen)
    fallback = calls["corrected_coordinate"] + np.where(
        calls["strand"].eq("+"), fixed_offset, -fixed_offset
    )
    calls["localization_status"] = np.where(
        calls["model_coordinate"].notna(), "model", "fixed_low_support"
    )
    calls["model_coordinate"] = calls["model_coordinate"].fillna(pd.Series(fallback))
    calls.to_csv(output_dir / "pas_calls.tsv", sep="\t", index=False)
    return calls


def read_endpoint(read: pysam.AlignedSegment, strand: str) -> int:
    return int(read.reference_end - 1 if strand == "+" else read.reference_start)


def usable_read(read: pysam.AlignedSegment) -> bool:
    return not (
        read.is_unmapped
        or read.is_duplicate
        or read.is_secondary
        or read.is_supplementary
    )


def add_read(
    counts,
    read: pysam.AlignedSegment,
    peak_name: str,
    barcode_tag: str,
    umi_tag: str,
    valid_barcodes: set[str] | None = None,
    stats=None,
) -> bool:
    try:
        barcode = str(read.get_tag(barcode_tag))
        umi = str(read.get_tag(umi_tag))
    except KeyError:
        if stats is not None:
            stats["reads_missing_tags"] = stats.get("reads_missing_tags", 0) + 1
        return False
    if valid_barcodes is not None and barcode not in valid_barcodes:
        if stats is not None:
            key = "reads_outside_barcode_whitelist"
            stats[key] = stats.get(key, 0) + 1
        return False
    if isinstance(counts, CompactMolecules):
        counts.add(peak_name, barcode, umi)
    else:
        counts[(peak_name, barcode)].add(umi)
    return True


def quantify_reference(
    bam_path: Path,
    intervals: pd.DataFrame,
    barcode_tag: str,
    umi_tag: str,
    valid_barcodes: set[str] | None = None,
) -> tuple[dict, dict]:
    counts = defaultdict(set)
    stats = defaultdict(int)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for peak in intervals.itertuples(index=False):
            for read in bam.fetch(peak.chrom, peak.start, peak.end):
                stats["reads_scanned"] += 1
                if not usable_read(read):
                    continue
                strand = "-" if read.is_reverse else "+"
                if strand != peak.strand:
                    continue
                endpoint = read_endpoint(read, strand)
                if not (peak.start <= endpoint < peak.end):
                    continue
                if add_read(
                    counts,
                    read,
                    peak.name,
                    barcode_tag,
                    umi_tag,
                    valid_barcodes,
                    stats,
                ):
                    stats["reads_assigned"] += 1
    return counts, dict(stats)


def quantify_sweep(
    bam_path: Path,
    intervals: pd.DataFrame,
    barcode_tag: str,
    umi_tag: str,
    compact: bool = False,
    valid_barcodes: set[str] | None = None,
) -> tuple[dict, dict]:
    counts = (
        CompactMolecules(intervals["name"].astype(str))
        if compact
        else defaultdict(set)
    )
    stats = defaultdict(int)
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        references = set(bam.references)
        for chrom, chrom_table in intervals.groupby("chrom", sort=False):
            fetch_chrom = chrom if chrom in references else (
                chrom[3:] if chrom.startswith("chr") and chrom[3:] in references else f"chr{chrom}"
            )
            indexes = {}
            for strand, group in chrom_table.groupby("strand", sort=False):
                group = group.sort_values("start")
                indexes[strand] = (
                    group["start"].to_list(),
                    list(group.itertuples(index=False)),
                )
            fetch_start = int(chrom_table["start"].min())
            fetch_end = int(chrom_table["end"].max())
            for read in bam.fetch(fetch_chrom, fetch_start, fetch_end):
                stats["reads_scanned"] += 1
                if not usable_read(read):
                    continue
                strand = "-" if read.is_reverse else "+"
                if strand not in indexes:
                    continue
                starts, peaks = indexes[strand]
                endpoint = read_endpoint(read, strand)
                index = bisect.bisect_right(starts, endpoint) - 1
                if index < 0:
                    continue
                peak = peaks[index]
                if endpoint >= peak.end:
                    continue
                if add_read(
                    counts,
                    read,
                    peak.name,
                    barcode_tag,
                    umi_tag,
                    valid_barcodes,
                    stats,
                ):
                    stats["reads_assigned"] += 1
    return counts, dict(stats)


def write_matrix(counts: dict, output_dir: Path) -> dict:
    if isinstance(counts, CompactMolecules):
        return write_compact_matrix(counts, output_dir)
    peaks = sorted({key[0] for key in counts})
    barcodes = sorted({key[1] for key in counts})
    peak_index = {value: index for index, value in enumerate(peaks)}
    barcode_index = {value: index for index, value in enumerate(barcodes)}
    rows, columns, values = [], [], []
    with gzip.open(output_dir / "counts.tsv.gz", "wt") as handle:
        handle.write("pas_id\tbarcode\tumi_count\n")
        for (peak, barcode), umis in sorted(counts.items()):
            rows.append(peak_index[peak])
            columns.append(barcode_index[barcode])
            values.append(len(umis))
            handle.write(f"{peak}\t{barcode}\t{len(umis)}\n")
    matrix = coo_matrix(
        (np.asarray(values, dtype=np.int32), (rows, columns)),
        shape=(len(peaks), len(barcodes)),
        dtype=np.int32,
    )
    mmwrite(output_dir / "matrix.mtx", matrix)
    (output_dir / "pas.tsv").write_text("\n".join(peaks) + ("\n" if peaks else ""))
    (output_dir / "barcodes.tsv").write_text("\n".join(barcodes) + ("\n" if barcodes else ""))
    return {
        "quantified_pas": len(peaks),
        "barcodes": len(barcodes),
        "nonzero_entries": len(values),
        "unique_molecules": int(sum(values)),
    }


def write_compact_matrix(counts: CompactMolecules, output_dir: Path) -> dict:
    used_peak_ids = {peak for peak, _, _ in counts.molecules}
    used_barcode_ids = {barcode for _, barcode, _ in counts.molecules}
    peak_ids = sorted(used_peak_ids, key=lambda value: counts.peaks[value])
    barcode_ids = sorted(
        used_barcode_ids, key=lambda value: counts.barcodes[value]
    )
    peak_rank = {value: index for index, value in enumerate(peak_ids)}
    barcode_rank = {value: index for index, value in enumerate(barcode_ids)}
    ordered = sorted(
        counts.molecules,
        key=lambda value: (
            peak_rank[value[0]], barcode_rank[value[1]], value[2]
        ),
    )
    rows, columns, values = array("i"), array("i"), array("i")
    with gzip.open(output_dir / "counts.tsv.gz", "wt") as handle:
        handle.write("pas_id\tbarcode\tumi_count\n")
        for (peak_id, barcode_id), group in itertools.groupby(
            ordered, key=lambda value: (value[0], value[1])
        ):
            molecule_count = sum(1 for _ in group)
            rows.append(peak_rank[peak_id])
            columns.append(barcode_rank[barcode_id])
            values.append(molecule_count)
            handle.write(
                f"{counts.peaks[peak_id]}\t{counts.barcodes[barcode_id]}"
                f"\t{molecule_count}\n"
            )
    matrix = coo_matrix(
        (
            np.frombuffer(values, dtype=np.int32),
            (
                np.frombuffer(rows, dtype=np.int32),
                np.frombuffer(columns, dtype=np.int32),
            ),
        ),
        shape=(len(peak_ids), len(barcode_ids)),
        dtype=np.int32,
    )
    mmwrite(output_dir / "matrix.mtx", matrix)
    peak_names = [counts.peaks[value] for value in peak_ids]
    barcode_names = [counts.barcodes[value] for value in barcode_ids]
    (output_dir / "pas.tsv").write_text(
        "\n".join(peak_names) + ("\n" if peak_names else "")
    )
    (output_dir / "barcodes.tsv").write_text(
        "\n".join(barcode_names) + ("\n" if barcode_names else "")
    )
    return {
        "quantified_pas": len(peak_ids),
        "barcodes": len(barcode_ids),
        "nonzero_entries": len(values),
        "unique_molecules": len(counts.molecules),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", type=Path, required=True)
    parser.add_argument("--chromosome", required=True)
    parser.add_argument("--peaks", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--barcode-tag", required=True)
    parser.add_argument("--umi-tag", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--backend", choices=["reference", "sweep"], default="sweep")
    parser.add_argument("--min-reads", type=int, default=20)
    args = parser.parse_args()
    total_start = time.perf_counter()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    intervals = load_intervals(args.peaks, args.chromosome)
    start = time.perf_counter()
    calls = localize(
        args.bam, args.chromosome, intervals, args.kernel, args.frozen,
        args.output_dir, args.min_reads
    )
    localization_seconds = time.perf_counter() - start
    start = time.perf_counter()
    function = quantify_sweep if args.backend == "sweep" else quantify_reference
    if args.backend == "sweep":
        counts, scan_stats = function(
            args.bam, intervals, args.barcode_tag, args.umi_tag, compact=True
        )
    else:
        counts, scan_stats = function(
            args.bam, intervals, args.barcode_tag, args.umi_tag
        )
    quantification_seconds = time.perf_counter() - start
    matrix_stats = write_matrix(counts, args.output_dir)
    summary = {
        "bam": str(args.bam),
        "chromosome": args.chromosome,
        "input_pas": len(intervals),
        "localized_pas": int(calls["model_coordinate"].notna().sum()),
        "model_localized_pas": int(calls["localization_status"].eq("model").sum()),
        "fixed_low_support_pas": int(calls["localization_status"].eq("fixed_low_support").sum()),
        "backend": args.backend,
        "barcode_tag": args.barcode_tag,
        "umi_tag": args.umi_tag,
        "localization_seconds": localization_seconds,
        "quantification_seconds": quantification_seconds,
        **scan_stats,
        **matrix_stats,
        "total_seconds": time.perf_counter() - total_start,
        "peak_rss_kb": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
