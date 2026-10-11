import os
from typing import List, Optional

import pandas as pd
from pybedtools import BedTool


def read_bed(bed_path: str, header: bool = False) -> pd.DataFrame:
    """Read BED file to DataFrame."""
    return pd.read_csv(bed_path, sep='\t', header=header)


def write_bed(bed_path: str, df: pd.DataFrame, header: bool = False):
    """Write DataFrame to BED file."""
    os.makedirs(os.path.dirname(bed_path), exist_ok=True)
    df.to_csv(bed_path, sep='\t', index=False, header=header)


def read_bam(bam_path: str):
    """Open BAM file using pysam."""
    import pysam
    return pysam.AlignmentFile(bam_path, 'rb')


def write_mtx_matrix(matrix, barcodes: List[str], features: List[str], output_dir: str):
    """Write expression matrix in MTX format."""
    from scipy.io import mmwrite
    from scipy.sparse import csr_matrix

    os.makedirs(output_dir, exist_ok=True)

    mtx_path = os.path.join(output_dir, "matrix.mtx")
    features_path = os.path.join(output_dir, "features.tsv")
    barcodes_path = os.path.join(output_dir, "barcodes.tsv")

    mmwrite(mtx_path, matrix)

    with open(features_path, 'w') as f:
        for feature in features:
            f.write(f"{feature}\t{feature}\tpAsite\n")

    with open(barcodes_path, 'w') as f:
        for barcode in barcodes:
            f.write(f"{barcode}\n")


def bedtool_from_dataframe(df: pd.DataFrame) -> BedTool:
    """Convert DataFrame to BedTool."""
    return BedTool.from_dataframe(df)


def bedtool_from_file(bed_path: str) -> BedTool:
    """Load BED file to BedTool."""
    return BedTool(bed_path)
