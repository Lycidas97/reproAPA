import pandas as pd


def merge_asnr_regions(final_bed_df: pd.DataFrame) -> pd.DataFrame:
    """Merge overlapping ASNR regions only when chromosome and strand both match."""
    if final_bed_df.empty:
        return final_bed_df.copy()

    final_bed_df = final_bed_df.sort_values(['chrom', 'strand', 'start']).reset_index(drop=True)
    prev_row = None
    merged_rows = []

    for _, row in final_bed_df.iterrows():
        row_dict = row.to_dict()
        if prev_row is None:
            prev_row = row_dict
            continue

        if (
            row['chrom'] != prev_row['chrom']
            or row['strand'] != prev_row['strand']
            or row['start'] > prev_row['end']
        ):
            merged_rows.append(prev_row)
            prev_row = row_dict
        else:
            prev_row['end'] = max(prev_row['end'], row['end'])

    if prev_row is not None:
        merged_rows.append(prev_row)

    return pd.DataFrame(merged_rows)
