import logging
import os
import shlex
import subprocess
from typing import Optional

import pysam
import pandas as pd
from Bio import SeqIO
from pybedtools import BedTool
from tqdm.auto import tqdm

from .asnr import merge_asnr_regions
from .logger import get_logger

logger = get_logger(__name__)


def get_chromosome_names_from_bam(bam_path: str):
    """Get chromosome names from BAM file header."""
    bam = pysam.AlignmentFile(bam_path, 'rb')
    chr_list = list(bam.references)
    bam.close()
    return chr_list


def get_chromosome_names_from_fasta(fasta_path: str):
    """Get chromosome names from FASTA file."""
    chrom_names = set()
    try:
        with open(fasta_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    chrom_name = line.split()[0][1:]
                    chrom_names.add(chrom_name)
        return sorted(list(chrom_names))
    except IOError as e:
        logger.error(f"Error reading FASTA file: {e}")
        return None


def get_chromosome_names_from_gtf(gtf_path: str):
    """Get unique chromosome names from GTF file."""
    command = f"cut -f1 {shlex.quote(gtf_path)} | grep -v '^#' | sort -u"
    try:
        output = subprocess.check_output(command, shell=True, universal_newlines=True)
        chromosome_list = [line for line in output.strip().split('\n') if line]
        return chromosome_list
    except subprocess.CalledProcessError as e:
        logger.error(f"Error processing GTF file: {e}")
        return None
    except FileNotFoundError:
        logger.warning("Unix tools not available. Using Python fallback.")
        chrom_names = set()
        try:
            with open(gtf_path, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if parts:
                        chrom_names.add(parts[0])
            return sorted(list(chrom_names))
        except IOError as e:
            logger.error(f"Error reading GTF file: {e}")
            return None


def create_chrom_mapping(bam_chrom_names, fasta_chrom_names, gtf_chrom_names):
    """Create a deterministic BAM/FASTA/GTF chromosome mapping.

    The unified namespace preserves an existing ``chr`` prefix and adds it when
    absent.  Construction is O(B + F + G), preserves BAM-header order, and
    rejects ambiguous inputs such as one source containing both ``1`` and
    ``chr1``.
    """
    bam_index, bam_order = _index_chromosome_names(bam_chrom_names, "bam_chrom")
    fasta_index, _ = _index_chromosome_names(fasta_chrom_names, "fasta_chrom")
    gtf_index, _ = _index_chromosome_names(gtf_chrom_names, "gtf_chrom")

    rows = [
        {
            "unified_chrom": unified,
            "bam_chrom": bam_index[unified],
            "fasta_chrom": fasta_index[unified],
            "gtf_chrom": gtf_index[unified],
        }
        for unified in bam_order
        if unified in fasta_index and unified in gtf_index
    ]
    mapping_df = pd.DataFrame(
        rows,
        columns=["unified_chrom", "bam_chrom", "fasta_chrom", "gtf_chrom"],
    )
    if mapping_df.empty:
        raise ValueError("No chromosomes are shared across BAM, FASTA, and GTF inputs.")
    return validate_chrom_mapping(mapping_df)


def unified_chrom_name(chrom_name: str) -> str:
    """Return the canonical chromosome name used by interval outputs."""
    if chrom_name is None or pd.isna(chrom_name):
        raise ValueError("Chromosome names cannot be null.")
    name = str(chrom_name).strip()
    if not name:
        raise ValueError("Chromosome names cannot be empty.")
    return name if name.startswith("chr") else f"chr{name}"


def _index_chromosome_names(chrom_names, source_column: str):
    if chrom_names is None:
        raise ValueError(f"{source_column} names are missing.")

    by_unified = {}
    order = []
    for value in chrom_names:
        unified = unified_chrom_name(value)
        raw_name = str(value).strip()
        previous = by_unified.get(unified)
        if previous is not None and previous != raw_name:
            raise ValueError(
                f"Ambiguous {source_column} names: '{previous}' and '{raw_name}' "
                f"both normalize to '{unified}'."
            )
        if previous is None:
            by_unified[unified] = raw_name
            order.append(unified)
    return by_unified, order


def validate_chrom_mapping(mapping_df: pd.DataFrame) -> pd.DataFrame:
    """Validate and return the four-column chromosome mapping contract."""
    columns = ["unified_chrom", "bam_chrom", "fasta_chrom", "gtf_chrom"]
    missing = [column for column in columns if column not in mapping_df.columns]
    if missing:
        raise ValueError(f"Chromosome mapping is missing columns: {missing}")

    mapping = mapping_df.loc[:, columns].copy()
    if mapping.empty:
        raise ValueError("Chromosome mapping is empty.")
    if mapping.isna().any().any():
        raise ValueError("Chromosome mapping contains null values.")
    mapping = mapping.astype(str)

    for column in columns:
        duplicated = mapping[column].duplicated(keep=False)
        if duplicated.any():
            values = sorted(mapping.loc[duplicated, column].unique().tolist())
            raise ValueError(f"Chromosome mapping column '{column}' is ambiguous: {values[:5]}")

    expected = mapping["unified_chrom"].map(unified_chrom_name)
    if not expected.equals(mapping["unified_chrom"]):
        raise ValueError("unified_chrom values are not canonical chromosome names.")
    return mapping


def normalize_chromosome_column(
    table: pd.DataFrame,
    mapping_df: pd.DataFrame,
    source_column: str,
    chrom_column: str = "chrom",
    drop_unmapped: bool = True,
):
    """Map one interval table into the unified namespace in a vectorized pass.

    Already-unified values are accepted, which makes the operation idempotent
    for new outputs while retaining compatibility with legacy cached BEDs.
    Returns ``(normalized_table, unmapped_names, dropped_row_count)``.
    """
    mapping = validate_chrom_mapping(mapping_df)
    if source_column not in {"bam_chrom", "fasta_chrom", "gtf_chrom"}:
        raise ValueError(f"Unsupported chromosome source column: {source_column}")
    if chrom_column not in table.columns:
        raise ValueError(f"Interval table is missing chromosome column: {chrom_column}")

    lookup = dict(zip(mapping[source_column], mapping["unified_chrom"]))
    lookup.update({name: name for name in mapping["unified_chrom"]})

    original = table[chrom_column].astype(str)
    normalized_names = original.map(lookup)
    unmapped_mask = normalized_names.isna()
    unmapped_names = tuple(sorted(original.loc[unmapped_mask].unique().tolist()))
    dropped_count = int(unmapped_mask.sum())
    if dropped_count and not drop_unmapped:
        preview = list(unmapped_names[:5])
        raise ValueError(
            f"{dropped_count} interval rows have chromosomes absent from "
            f"'{source_column}' mapping: {preview}"
        )

    keep_mask = ~unmapped_mask if drop_unmapped else pd.Series(True, index=table.index)
    normalized = table.loc[keep_mask].copy()
    normalized[chrom_column] = normalized_names.loc[keep_mask].astype(str).to_numpy()
    return normalized, unmapped_names, dropped_count


def _find_asnr(seq: str, match_length: int, max_mismatch: int, strand: str) -> list:
    """Find A-rich (or T-rich on reverse strand) regions."""
    indices = []
    i = 0
    length = len(seq)

    if strand == "+":
        while i < length:
            if seq[i] != 'A':
                i += 1
                continue
            j = i + 1
            non_A_count = 0
            while j < len(seq) and non_A_count <= max_mismatch:
                if seq[j] != 'A':
                    non_A_count += 1
                j += 1
            while non_A_count > max_mismatch:
                j -= 1
                if seq[j] != 'A':
                    non_A_count -= 1
            if j - i >= match_length:
                indices.append(i)
            i = j
    else:
        seq = seq[::-1]
        while i < length:
            if seq[i] != 'T':
                i += 1
                continue
            j = i + 1
            non_T_count = 0
            while j < len(seq) and non_T_count <= max_mismatch:
                if seq[j] != 'T':
                    non_T_count += 1
                j += 1
            while non_T_count > max_mismatch:
                j -= 1
                if seq[j] != 'T':
                    non_T_count -= 1
            if j - i >= match_length:
                indices.append(length - i - 1)
            i = j
    return indices


def generate_asnr(
    gtf_path: str,
    fasta_path: str,
    output_path: str,
    match_length: int = 10,
    max_mismatch: int = 1,
    upstream_extend: int = 500,
    downstream_extend: int = 20,
    bedtools_path: Optional[str] = None,
    chrom_map_path: Optional[str] = None,
    show_progress: bool = False,
    progress_mininterval: float = 10.0
) -> str:
    """Generate ASNR annotation BED file."""
    logger.info("Starting ASNR annotation generation.")
    logger.info(f"Parameters: match_length={match_length}, max_mismatch={max_mismatch}, upstream_extend={upstream_extend}, downstream_extend={downstream_extend}")

    for f in [gtf_path, fasta_path]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"Input file not found: {f}")

    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        logger.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)

    if bedtools_path:
        try:
            BedTool.set_bedtools_path(bedtools_path)
            logger.info(f"PyBedTools path set to: {bedtools_path}")
        except Exception:
            logger.warning(f"Could not set bedtools path to '{bedtools_path}'")

    try:
        logger.info("Loading genome FASTA file...")
        with open(fasta_path, "r") as fasta_handle:
            genome_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))

        logger.info("Loading and parsing GTF file for 'gene' features...")
        gtf_df = pd.read_csv(gtf_path, comment='#', sep='\t', header=None,
                             names=['chrom', 'source', 'feature', 'start', 'end',
                                    'score', 'strand', 'frame', 'attribute'],
                             low_memory=False)
        gtf_df['chrom'] = gtf_df['chrom'].astype(str)
        genes_gtf = gtf_df.loc[(gtf_df['feature'] == 'gene'), ['chrom', 'start', 'end', 'feature', 'score', 'strand']].copy()

        logger.info("Loading chromosome mapping...")
        chrom_mapping_df = None
        if chrom_map_path:
            if not os.path.exists(chrom_map_path):
                raise FileNotFoundError(f"Chromosome mapping not found: {chrom_map_path}")
            chrom_mapping_df = validate_chrom_mapping(
                pd.read_csv(chrom_map_path, sep="\t", dtype=str)
            )
            gtf_to_fasta_map = dict(
                zip(chrom_mapping_df["gtf_chrom"], chrom_mapping_df["fasta_chrom"])
            )
        else:
            fasta_chroms = list(genome_dict.keys())
            gtf_chroms = genes_gtf['chrom'].unique().astype(str)
            fasta_chrom_map = {c: c for c in fasta_chroms}
            gtf_to_fasta_map = {}
            for gtf_chrom in gtf_chroms:
                if gtf_chrom in fasta_chrom_map:
                    gtf_to_fasta_map[gtf_chrom] = gtf_chrom
                elif f"chr{gtf_chrom}" in fasta_chrom_map:
                    gtf_to_fasta_map[gtf_chrom] = f"chr{gtf_chrom}"

        logger.info(f"Mapped {len(gtf_to_fasta_map)} GTF chromosomes to FASTA chromosomes")

        logger.info("Extracting gene sequences...")
        genes_gtf["fasta_chrom"] = genes_gtf['chrom'].map(gtf_to_fasta_map)
        genes_gtf = genes_gtf.dropna(subset=['fasta_chrom'])
        if genes_gtf.empty:
            raise ValueError("No GTF gene chromosomes could be mapped to FASTA chromosomes for ASNR generation.")

        gene_df = genes_gtf[['chrom', 'start', 'end', 'feature', 'score', 'strand', 'fasta_chrom']].copy()
        gene_df['start0'] = gene_df['start'] - 1

        def extract_seq(row):
            return str(genome_dict[row['fasta_chrom']].seq[row['start0']:row['end']])

        gene_df["seq"] = gene_df.apply(extract_seq, axis=1)
        logger.info(f"Found {len(gene_df)} gene entries.")

        logger.info("Finding A/T-rich regions (ASNR) within each gene...")
        if show_progress:
            tqdm.pandas(desc="Finding ASNR sites", mininterval=progress_mininterval)
            gene_df["asnr_stop"] = gene_df.progress_apply(
                lambda x: _find_asnr(x["seq"], match_length, max_mismatch, x["strand"]),
                axis=1
            )
        else:
            gene_df["asnr_stop"] = gene_df.apply(
                lambda x: _find_asnr(x["seq"], match_length, max_mismatch, x["strand"]),
                axis=1
            )

        logger.info("Calculating genomic coordinates for extended ASNR regions...")
        gene_df["asnr_expand"] = gene_df.apply(
            lambda x: [[max(0, i-upstream_extend), i+downstream_extend] for i in x["asnr_stop"]]
            if x["strand"] == "+" else
            [[i-downstream_extend, min(i+upstream_extend, len(x["seq"]))] for i in x["asnr_stop"]],
            axis=1
        )

        gene_df["asnr_coord"] = gene_df.apply(
            lambda x: [[x["start0"]+start, x["start0"]+end] for start, end in x["asnr_expand"]],
            axis=1
        )

        logger.info("Formatting data into BED format...")
        asnr_df = pd.merge(
            gene_df["asnr_coord"].explode().dropna().reset_index(),
            gene_df[['chrom', 'strand']].reset_index()
        )
        asnr_df["start"] = asnr_df["asnr_coord"].apply(lambda x: x[0])
        asnr_df["end"] = asnr_df["asnr_coord"].apply(lambda x: x[1])
        asnr_df["name"] = "intron"
        asnr_df["score"] = "."

        if asnr_df.empty:
            logger.warning("No ASNR regions were found with given parameters.")
            return None

        final_bed_df = asnr_df[["chrom", "start", "end", "name", "score", "strand"]]
        if chrom_mapping_df is not None:
            final_bed_df, unmapped_names, dropped_count = normalize_chromosome_column(
                final_bed_df,
                chrom_mapping_df,
                source_column="gtf_chrom",
            )
            if dropped_count:
                logger.info(
                    f"Dropped {dropped_count} ASNR rows on {len(unmapped_names)} "
                    "chromosomes outside the BAM/FASTA/GTF shared mapping."
                )
            if final_bed_df.empty:
                raise ValueError("No ASNR regions remain on mapped chromosomes.")

        logger.info(f"Found {len(final_bed_df)} raw ASNR regions. Sorting and merging overlapping regions...")


        merged_bed_df = merge_asnr_regions(final_bed_df)

        logger.info(f"After merging, {len(merged_bed_df)} final ASNR regions remain.")

        merged_bed_df.to_csv(output_path, sep='\t', index=False, header=False)

        logger.info(f"ASNR annotation generation complete. Output saved to: {output_path}")
        return output_path

    except Exception as e:
        logger.error(f"An error occurred during ASNR generation: {e}")
        raise
