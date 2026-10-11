import os
import bisect
import time
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Set, Optional, Tuple, Dict, Any

import numpy as np
import pandas as pd
import pysam
from scipy.sparse import coo_matrix, csr_matrix
from tqdm import tqdm

from ..utils.logger import get_logger
from ..utils.io import write_mtx_matrix
from ..utils.strand import read_transcript_strand, validate_transcript_strand_mode

logger = get_logger(__name__)

CHROM_SWEEP_MIN_BLOCKS = 100


def group_contiguous_peaks(peak_bed_df: pd.DataFrame) -> pd.DataFrame:
    """
    Group contiguous peaks into blocks with unique block IDs.

    Args:
        peak_bed_df: DataFrame with peak intervals.

    Returns:
        DataFrame with added block_id column.
    """
    df = peak_bed_df.sort_values(["chrom", "strand", "start"]).copy()
    is_contiguous = (df['chrom'] == df['chrom'].shift(1)) & \
                    (df['strand'] == df['strand'].shift(1)) & \
                    (df['start'] == df['end'].shift(1))
    df['block_id'] = (~is_contiguous).cumsum()
    return df


def _nested_dd_factory():
    """Factory function for multiprocessing."""
    return defaultdict(set)


def _resolve_bam_chrom(chrom: str, references: Set[str]) -> str:
    """Resolve a peak chromosome name to the BAM reference name."""
    if chrom in references:
        return chrom
    if chrom.startswith("chr"):
        unprefixed = chrom[3:]
        if unprefixed in references:
            return unprefixed
    else:
        prefixed = f"chr{chrom}"
        if prefixed in references:
            return prefixed
    return chrom


def _assign_peak_from_overlaps(overlapping_peaks: list, strand: str) -> Optional[str]:
    if not overlapping_peaks:
        return None
    if len(overlapping_peaks) == 1:
        return overlapping_peaks[0].name
    if strand == '+':
        downstream_peak = max(overlapping_peaks, key=lambda p: p.start)
    else:
        downstream_peak = min(overlapping_peaks, key=lambda p: p.start)
    return downstream_peak.name


def _process_chunk_block_fetch(
    bam_path: str,
    chunk_df: pd.DataFrame,
    umi_tag: str,
    barcode_tag: str,
    transcript_strand_from_read: str
) -> defaultdict:
    """Process a chunk of peaks (typically all peaks from one chromosome)."""
    local_umi_counts = defaultdict(_nested_dd_factory)

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        references = set(bam.references)
        for block_id, block_df in chunk_df.groupby('block_id'):
            chrom = block_df['chrom'].iloc[0]
            block_start = block_df['start'].min()
            block_end = block_df['end'].max()
            fetch_chrom = _resolve_bam_chrom(chrom, references)

            block_peaks = list(block_df.itertuples(index=False, name='Peak'))

            if len(block_peaks) == 1:
                the_only_peak = block_peaks[0]
                peak_name = the_only_peak.name
                peak_strand = the_only_peak.strand
                for read in bam.fetch(fetch_chrom, block_start, block_end):
                    if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                        continue
                    if read_transcript_strand(read.is_reverse, transcript_strand_from_read) != peak_strand:
                        continue
                    if read.reference_start < the_only_peak.end and read.reference_end > the_only_peak.start:
                        try:
                            barcode = read.get_tag(barcode_tag)
                            umi = read.get_tag(umi_tag)
                            local_umi_counts[peak_name][barcode].add(umi)
                        except KeyError:
                            continue
            else:
                block_peak_starts = [p.start for p in block_peaks]
                strand = block_peaks[0].strand
                is_positive_strand = (strand == '+')
                for read in bam.fetch(fetch_chrom, block_start, block_end):
                    if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                        continue
                    if read_transcript_strand(read.is_reverse, transcript_strand_from_read) != strand:
                        continue
                    try:
                        barcode = read.get_tag(barcode_tag)
                        umi = read.get_tag(umi_tag)
                    except KeyError:
                        continue

                    overlapping_peaks = []
                    read_start, read_end = read.reference_start, read.reference_end
                    start_idx = bisect.bisect_right(block_peak_starts, read_start)

                    for j in range(max(0, start_idx - 1), len(block_peaks)):
                        peak = block_peaks[j]
                        if peak.start >= read_end:
                            break
                        if peak.end > read_start:
                            overlapping_peaks.append(peak)

                    assigned_peak_name = None
                    if len(overlapping_peaks) == 1:
                        assigned_peak_name = overlapping_peaks[0].name
                    elif len(overlapping_peaks) > 1:
                        if is_positive_strand:
                            downstream_peak = max(overlapping_peaks, key=lambda p: p.start)
                        else:
                            downstream_peak = min(overlapping_peaks, key=lambda p: p.start)
                        assigned_peak_name = downstream_peak.name

                    if assigned_peak_name:
                        local_umi_counts[assigned_peak_name][barcode].add(umi)

    return local_umi_counts


def _process_chunk_chrom_sweep(
    bam_path: str,
    chunk_df: pd.DataFrame,
    umi_tag: str,
    barcode_tag: str,
    transcript_strand_from_read: str
) -> defaultdict:
    """Process peaks by streaming each chromosome once per chunk."""
    local_umi_counts = defaultdict(_nested_dd_factory)

    if chunk_df.empty:
        return local_umi_counts

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        references = set(bam.references)

        for chrom, chrom_df in chunk_df.groupby('chrom', sort=False):
            fetch_chrom = _resolve_bam_chrom(chrom, references)
            fetch_start = int(chrom_df['start'].min())
            fetch_end = int(chrom_df['end'].max())
            if fetch_start >= fetch_end:
                continue

            peaks_by_strand = {}
            for strand, strand_df in chrom_df.groupby('strand', sort=False):
                peaks_by_strand[strand] = list(
                    strand_df.sort_values(["start", "end", "name"]).itertuples(index=False, name='Peak')
                )

            next_peak_idx = {strand: 0 for strand in peaks_by_strand}
            active_peaks = {strand: [] for strand in peaks_by_strand}

            for read in bam.fetch(fetch_chrom, fetch_start, fetch_end):
                if read.is_unmapped or read.is_duplicate or read.is_secondary or read.is_supplementary:
                    continue
                read_start, read_end = read.reference_start, read.reference_end
                if read_start is None or read_end is None or read_start >= read_end:
                    continue

                strand = read_transcript_strand(read.is_reverse, transcript_strand_from_read)
                if strand not in peaks_by_strand:
                    continue

                strand_peaks = peaks_by_strand[strand]
                idx = next_peak_idx[strand]
                while idx < len(strand_peaks) and strand_peaks[idx].start < read_end:
                    active_peaks[strand].append(strand_peaks[idx])
                    idx += 1
                next_peak_idx[strand] = idx

                active = [peak for peak in active_peaks[strand] if peak.end > read_start]
                active_peaks[strand] = active
                if not active:
                    continue

                try:
                    barcode = read.get_tag(barcode_tag)
                    umi = read.get_tag(umi_tag)
                except KeyError:
                    continue

                overlapping_by_block = defaultdict(list)
                for peak in active:
                    if peak.start < read_end and peak.end > read_start:
                        overlapping_by_block[peak.block_id].append(peak)

                for overlapping_peaks in overlapping_by_block.values():
                    assigned_peak_name = _assign_peak_from_overlaps(overlapping_peaks, strand)
                    if assigned_peak_name:
                        local_umi_counts[assigned_peak_name][barcode].add(umi)

    return local_umi_counts


def _process_chunk(
    bam_path: str,
    chunk_df: pd.DataFrame,
    umi_tag: str,
    barcode_tag: str,
    transcript_strand_from_read: str,
    quantification_backend: str = "auto",
) -> defaultdict:
    if quantification_backend == "auto":
        quantification_backend = _choose_quantification_backend(chunk_df, quantification_backend)
    if quantification_backend == "chrom_sweep":
        return _process_chunk_chrom_sweep(
            bam_path,
            chunk_df,
            umi_tag,
            barcode_tag,
            transcript_strand_from_read,
        )
    if quantification_backend == "block_fetch":
        return _process_chunk_block_fetch(
            bam_path,
            chunk_df,
            umi_tag,
            barcode_tag,
            transcript_strand_from_read,
        )
    raise ValueError(
        "Unsupported quantification_backend "
        f"'{quantification_backend}'. Expected auto, chrom_sweep, or block_fetch."
    )


def _choose_quantification_backend(peak_bed_df: pd.DataFrame, requested_backend: str) -> str:
    if requested_backend != "auto":
        if requested_backend not in {"chrom_sweep", "block_fetch"}:
            raise ValueError(
                "Unsupported quantification_backend "
                f"'{requested_backend}'. Expected auto, chrom_sweep, or block_fetch."
            )
        return requested_backend
    if peak_bed_df.empty or "block_id" not in peak_bed_df:
        return "block_fetch"
    n_blocks = int(peak_bed_df["block_id"].nunique())
    return "chrom_sweep" if n_blocks >= CHROM_SWEEP_MIN_BLOCKS else "block_fetch"


def _empty_qc_value(value):
    return "" if value is None else value


def _write_quantification_qc(
    peak_bed_df: pd.DataFrame,
    global_umi_counts: dict,
    all_barcodes: List[str],
    filtered_cell_tags: List[str],
    matrix_peaks: List[str],
    output_dir: str,
    min_reads_per_cell: int,
    observed_barcode_count: int,
    valid_barcode_count: Optional[int],
    overlapping_valid_barcode_count: Optional[int],
) -> Dict[str, Any]:
    """Write sample-level and PAS-level quantification QC tables."""
    os.makedirs(output_dir, exist_ok=True)

    filtered_cells = set(filtered_cell_tags)
    matrix_peak_set = set(matrix_peaks)
    filtered_pas_names = set(peak_bed_df["name"].astype(str))

    status_rows = []
    total_umi_pre_cell_filter = 0
    total_umi_post_cell_filter = 0

    for row in peak_bed_df.itertuples(index=False):
        peak_name = str(row.name)
        barcode_to_umis = global_umi_counts.get(peak_name, {})
        umi_count_pre = sum(len(umis) for umis in barcode_to_umis.values())
        cell_count_pre = sum(1 for umis in barcode_to_umis.values() if umis)
        umi_count_post = sum(
            len(umis)
            for barcode, umis in barcode_to_umis.items()
            if barcode in filtered_cells
        )
        cell_count_post = sum(
            1
            for barcode, umis in barcode_to_umis.items()
            if barcode in filtered_cells and umis
        )

        total_umi_pre_cell_filter += umi_count_pre
        total_umi_post_cell_filter += umi_count_post

        in_matrix = peak_name in matrix_peak_set
        if in_matrix:
            drop_reason = "none"
        elif umi_count_pre == 0:
            drop_reason = "no_overlapping_umi"
        else:
            drop_reason = "only_filtered_cells"

        status_rows.append({
            "chrom": row.chrom,
            "start": row.start,
            "end": row.end,
            "name": peak_name,
            "score": row.score,
            "strand": row.strand,
            "has_umi": str(umi_count_pre > 0).lower(),
            "umi_count_pre_cell_filter": umi_count_pre,
            "cell_count_pre_cell_filter": cell_count_pre,
            "umi_count_post_cell_filter": umi_count_post,
            "cell_count_post_cell_filter": cell_count_post,
            "in_matrix": str(in_matrix).lower(),
            "drop_reason": drop_reason,
        })

    filtered_pas_count = len(filtered_pas_names)
    quantified_pas_count = len(matrix_peaks)
    dropped_pas_count = filtered_pas_count - quantified_pas_count

    qc_summary = {
        "filtered_pas_count": filtered_pas_count,
        "quantified_pas_count": quantified_pas_count,
        "dropped_pas_count": dropped_pas_count,
        "observed_barcode_count": observed_barcode_count,
        "valid_barcode_count": _empty_qc_value(valid_barcode_count),
        "overlapping_valid_barcode_count": _empty_qc_value(overlapping_valid_barcode_count),
        "pre_filter_cell_count": len(all_barcodes),
        "post_filter_cell_count": len(filtered_cell_tags),
        "filtered_cell_count": len(all_barcodes) - len(filtered_cell_tags),
        "total_umi_pre_cell_filter": total_umi_pre_cell_filter,
        "total_umi_post_cell_filter": total_umi_post_cell_filter,
        "min_reads_per_cell": min_reads_per_cell,
    }

    pd.DataFrame([qc_summary]).to_csv(
        os.path.join(output_dir, "qc_summary.tsv"),
        sep="\t",
        index=False,
    )
    pd.DataFrame(status_rows).to_csv(
        os.path.join(output_dir, "pas_quantification_status.tsv"),
        sep="\t",
        index=False,
    )
    return qc_summary


def quantify_pas(
    bam_path: str,
    peak_bed_df: pd.DataFrame,
    output_dir: str,
    umi_tag: str = "UB",
    barcode_tag: str = "CB",
    n_jobs: int = 4,
    show_progress: bool = False,
    valid_barcodes: Optional[Set[str]] = None,
    min_reads_per_cell: int = 50,
    transcript_strand_from_read: str = "same",
    quantification_backend: str = "auto",
) -> Tuple[csr_matrix, List[str], List[str]]:
    """
    Quantify PAS expression using UMI deduplication and sparse matrix.

    Args:
        bam_path: Path to BAM file.
        peak_bed_df: DataFrame with peak intervals.
        output_dir: Directory to save output files.
        umi_tag: Tag for UMI in BAM.
        barcode_tag: Tag for cell barcode in BAM.
        n_jobs: Number of CPU cores to use.
        show_progress: Whether to show progress bar.
        valid_barcodes: Set of valid barcodes (optional).
        min_reads_per_cell: Minimum reads per cell to be valid.
        transcript_strand_from_read: How read orientation maps to transcript strand.
        quantification_backend: auto to pick per sample size, block_fetch for the
            legacy per-block fetcher, chrom_sweep to stream each chromosome once
            per chunk.

    Returns:
        Tuple of (sparse_matrix, peaks, barcodes).
    """
    logger.info("--- Starting PAS Quantification ---")
    transcript_strand_from_read = validate_transcript_strand_mode(transcript_strand_from_read)

    logger.info("Step 1: Identifying contiguous peak blocks...")
    peaks_with_blocks = group_contiguous_peaks(peak_bed_df)
    effective_quantification_backend = _choose_quantification_backend(peaks_with_blocks, quantification_backend)

    logger.info(
        f"Step 2: Processing reads using {n_jobs} processes "
        f"with {effective_quantification_backend} backend "
        f"(requested={quantification_backend})..."
    )

    task_chunks = []
    if effective_quantification_backend == "chrom_sweep" and n_jobs > 1 and not peaks_with_blocks.empty:
        chrom_sizes = peaks_with_blocks.groupby('chrom', sort=False).size().sort_values(ascending=False)
        chunks = [[] for _ in range(n_jobs)]
        chunk_loads = [0 for _ in range(n_jobs)]
        for chrom, size in chrom_sizes.items():
            chunk_idx = min(range(n_jobs), key=lambda idx: chunk_loads[idx])
            chunks[chunk_idx].append(peaks_with_blocks[peaks_with_blocks['chrom'] == chrom])
            chunk_loads[chunk_idx] += int(size)
        task_chunks = [pd.concat(groups) for groups in chunks if groups]
    elif n_jobs > 1 and 'block_id' in peaks_with_blocks and peaks_with_blocks['block_id'].nunique() > 1:
        block_sizes = peaks_with_blocks['block_id'].value_counts(sort=False)
        sorted_block_ids = block_sizes.sort_values(ascending=False).index
        block_to_chunk_map = {block_id: i % n_jobs for i, block_id in enumerate(sorted_block_ids)}
        peaks_with_blocks['chunk_id'] = peaks_with_blocks['block_id'].map(block_to_chunk_map)
        task_chunks = [group for _, group in peaks_with_blocks.groupby('chunk_id')]
    elif not peaks_with_blocks.empty:
        task_chunks.append(peaks_with_blocks)

    logger.info(f"Step 3: Processing {peaks_with_blocks['block_id'].nunique()} blocks in {len(task_chunks)} chunks...")

    logger.info("Step 4: Building final expression matrix using scipy.sparse...")
    start_time = time.time()

    global_umi_counts = defaultdict(lambda: defaultdict(set))

    with ProcessPoolExecutor(max_workers=n_jobs) as executor:
        futures = {
            executor.submit(
                _process_chunk,
                bam_path,
                chunk,
                umi_tag,
                barcode_tag,
                transcript_strand_from_read,
                effective_quantification_backend,
            ): i
            for i, chunk in enumerate(task_chunks)
        }

        iterable = as_completed(futures)
        if show_progress:
            iterable = tqdm(iterable, total=len(task_chunks), desc="Processing chunks")

        for future in iterable:
            local_counts = future.result()
            for peak_name, barcodes_dict in local_counts.items():
                for barcode, umis in barcodes_dict.items():
                    global_umi_counts[peak_name][barcode].update(umis)

    if not global_umi_counts:
        logger.info("No expression data found.")
        _write_quantification_qc(
            peak_bed_df=peak_bed_df,
            global_umi_counts={},
            all_barcodes=[],
            filtered_cell_tags=[],
            matrix_peaks=[],
            output_dir=output_dir,
            min_reads_per_cell=min_reads_per_cell,
            observed_barcode_count=0,
            valid_barcode_count=len(valid_barcodes) if valid_barcodes else None,
            overlapping_valid_barcode_count=0 if valid_barcodes else None,
        )
        return csr_matrix((0, 0), dtype=np.int32), [], []

    logger.info(f"Processing {len(global_umi_counts)} peaks with expression data")
    all_peaks = sorted(global_umi_counts.keys())

    start_time = time.time()

    if valid_barcodes:
        logger.info("Checking overlap between provided valid barcodes and observed data...")
        logger.info("Extracting observed barcodes from global_umi_counts...")
        obs_start_time = time.time()
        observed_barcodes = set(bc for peak_bcs in global_umi_counts.values() for bc in peak_bcs)
        observed_barcode_count = len(observed_barcodes)
        logger.info(f"Extracted {len(observed_barcodes)} unique barcodes in {time.time() - obs_start_time:.2f} seconds")
        overlapping_barcodes = valid_barcodes.intersection(observed_barcodes)

        valid_barcode_count = len(valid_barcodes)
        overlapping_valid_barcode_count = len(overlapping_barcodes)
        logger.info(f"Provided {valid_barcode_count} valid barcodes, found {overlapping_valid_barcode_count} overlapping barcodes in data.")

        if not overlapping_barcodes:
            raise ValueError("The provided valid barcode list has no overlap with barcodes found in the BAM file. Cannot build matrix.")

        logger.info("Using overlapping barcodes to build matrix.")
        all_barcodes = sorted(list(overlapping_barcodes))
    else:
        valid_barcode_count = None
        overlapping_valid_barcode_count = None
        logger.info("Deriving barcodes from observed data (this can be slow)...")
        obs_start_time = time.time()
        all_barcodes = sorted(set(bc for peak_bcs in global_umi_counts.values() for bc in peak_bcs))
        observed_barcode_count = len(all_barcodes)
        logger.info(f"Derived {len(all_barcodes)} barcodes in {time.time() - obs_start_time:.2f} seconds")

    peak_map = {peak: i for i, peak in enumerate(all_peaks)}
    barcode_map = {barcode: i for i, barcode in enumerate(all_barcodes)}

    data, row, col = [], [], []
    for peak_name, barcodes_dict in global_umi_counts.items():
        peak_idx = peak_map[peak_name]
        for barcode, umis in barcodes_dict.items():
            if barcode not in barcode_map:
                continue
            row.append(peak_idx)
            col.append(barcode_map[barcode])
            data.append(len(umis))

    sparse_matrix = coo_matrix((data, (row, col)),
                               shape=(len(all_peaks), len(all_barcodes)),
                               dtype=np.int32).tocsr()

    filtered_cell_indices = np.where(sparse_matrix.sum(axis=0).A1 >= min_reads_per_cell)[0]

    if len(filtered_cell_indices) == 0:
        _write_quantification_qc(
            peak_bed_df=peak_bed_df,
            global_umi_counts=global_umi_counts,
            all_barcodes=all_barcodes,
            filtered_cell_tags=[],
            matrix_peaks=[],
            output_dir=output_dir,
            min_reads_per_cell=min_reads_per_cell,
            observed_barcode_count=observed_barcode_count,
            valid_barcode_count=valid_barcode_count,
            overlapping_valid_barcode_count=overlapping_valid_barcode_count,
        )
        raise ValueError("No cells passed the minimum reads threshold")

    filtered_matrix = sparse_matrix[:, filtered_cell_indices]
    filtered_cell_tags = [all_barcodes[idx] for idx in filtered_cell_indices]
    retained_peak_indices = np.where(filtered_matrix.sum(axis=1).A1 > 0)[0]
    matrix_peaks = [all_peaks[idx] for idx in retained_peak_indices]
    filtered_matrix = filtered_matrix[retained_peak_indices, :]

    qc_summary = _write_quantification_qc(
        peak_bed_df=peak_bed_df,
        global_umi_counts=global_umi_counts,
        all_barcodes=all_barcodes,
        filtered_cell_tags=filtered_cell_tags,
        matrix_peaks=matrix_peaks,
        output_dir=output_dir,
        min_reads_per_cell=min_reads_per_cell,
        observed_barcode_count=observed_barcode_count,
        valid_barcode_count=valid_barcode_count,
        overlapping_valid_barcode_count=overlapping_valid_barcode_count,
    )
    logger.info(
        "Quantification QC: filtered_pas=%s, quantified_pas=%s, dropped_pas=%s, cells=%s/%s",
        qc_summary["filtered_pas_count"],
        qc_summary["quantified_pas_count"],
        qc_summary["dropped_pas_count"],
        qc_summary["post_filter_cell_count"],
        qc_summary["pre_filter_cell_count"],
    )

    logger.info(f"Matrix construction finished in {time.time() - start_time:.2f} seconds.")

    write_mtx_matrix(filtered_matrix, filtered_cell_tags, matrix_peaks, output_dir)

    logger.info(f"PAS quantification complete. Output saved to: {output_dir}")
    logger.info("--- PAS Quantification Complete ---")

    return filtered_matrix, matrix_peaks, filtered_cell_tags
