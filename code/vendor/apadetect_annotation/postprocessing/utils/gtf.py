import re
from typing import Iterable, Optional

import pandas as pd


GTF_COLUMNS = [
    "chrom",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
    "attribute",
]

ATTRIBUTE_RE = re.compile(r'(\S+)\s+"([^"]*)"')


def parse_gtf_attributes(attribute: str) -> dict:
    """Parse a GTF/GFF-style attribute string into a dictionary."""
    if not isinstance(attribute, str):
        return {}

    parsed = {key: value for key, value in ATTRIBUTE_RE.findall(attribute)}
    if parsed:
        return parsed

    fallback = {}
    for item in attribute.strip().rstrip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
        elif " " in item:
            key, value = item.split(None, 1)
        else:
            continue
        fallback[key.strip()] = value.strip().strip('"')
    return fallback


def normalize_transcript_id(transcript_id: str) -> str:
    """Normalize RSEM-style transcript ids such as ENST|symbol."""
    if not isinstance(transcript_id, str):
        return ""
    transcript_id = transcript_id.strip()
    if "|" in transcript_id:
        return transcript_id.split("|", 1)[0]
    return transcript_id


def read_gtf(path: str, features: Optional[Iterable[str]] = None) -> pd.DataFrame:
    """Read a GTF/GFF-like file and normalize core annotation attributes."""
    feature_set = set(features) if features else None
    df = pd.read_csv(
        path,
        comment="#",
        sep="\t",
        header=None,
        names=GTF_COLUMNS,
        dtype={
            "chrom": str,
            "source": str,
            "feature": str,
            "score": str,
            "strand": str,
            "frame": str,
            "attribute": str,
        },
        compression="infer",
        low_memory=False,
    )
    if feature_set is not None:
        df = df[df["feature"].isin(feature_set)].copy()
    return normalize_gtf_attributes(df)


def normalize_gtf_attributes(df: pd.DataFrame) -> pd.DataFrame:
    """Add normalized gene/transcript columns to a GTF dataframe."""
    df = df.copy()
    attrs = df["attribute"].fillna("").map(parse_gtf_attributes)

    def get_attr(keys):
        return attrs.map(lambda item: next((item.get(key, "") for key in keys if item.get(key, "")), ""))

    df["gene_id"] = get_attr(["gene_id", "gene", "ID", "Parent"])
    df["gene_symbol"] = get_attr(["gene_name", "gene", "Name", "gene_symbol"])
    df["transcript_id_raw"] = get_attr(["transcript_id", "transcript", "Parent", "ID"])
    df["transcript_id"] = df["transcript_id_raw"].map(normalize_transcript_id)
    df["gene_biotype"] = get_attr(["gene_biotype", "gene_type", "biotype"])
    df["transcript_biotype"] = get_attr(["transcript_biotype", "transcript_type", "biotype"])
    df["exon_number"] = get_attr(["exon_number"])

    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    df = df.dropna(subset=["start", "end"]).copy()
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["start0"] = df["start"] - 1
    df["end0"] = df["end"]
    df["gene_symbol"] = df["gene_symbol"].where(df["gene_symbol"].astype(str).str.len() > 0, df["gene_id"])
    return df


def build_reference_intervals(
    gtf_df: pd.DataFrame,
    transcript_end_zones_bp: Optional[list[int]] = None,
    gene_downstream_bp: int = 2000,
) -> dict[str, pd.DataFrame]:
    """Derive the interval tables needed for a first-pass PAS annotation."""
    zones = transcript_end_zones_bp or [250, 1000, 2000]
    zones = sorted(int(value) for value in zones)
    intervals: dict[str, pd.DataFrame] = {}

    genes = _build_gene_table(gtf_df)
    transcripts = _build_transcript_table(gtf_df)
    exons = _feature_table(gtf_df, "exon")
    cds = _feature_table(gtf_df, "CDS")
    utr3 = _feature_table(gtf_df, "three_prime_utr")
    if utr3.empty:
        utr3 = _feature_table(gtf_df, "3UTR")

    transcript_ends = _derive_transcript_ends(transcripts)
    exons = _attach_model_transcript_ends(exons, transcript_ends)
    cds = _attach_model_transcript_ends(cds, transcript_ends)
    utr3 = _attach_model_transcript_ends(utr3, transcript_ends)

    intervals["gene_body"] = genes
    intervals["transcript_body"] = transcripts
    intervals["exon"] = exons
    intervals["cds"] = cds
    intervals["utr3"] = _with_region(utr3, "annotated_3utr_pas", 2)

    terminal_exons, internal_exons, terminal_introns, internal_introns = _derive_exon_intervals(exons)
    intervals["terminal_exon"] = _with_region(terminal_exons, "terminal_exon_pas", 3)
    intervals["internal_exon"] = _with_region(internal_exons, "internal_exon_pas", 10)
    intervals["terminal_intron"] = _with_region(terminal_introns, "terminal_intron_pas", 7)
    intervals["internal_intron"] = _with_region(internal_introns, "internal_intron_pas", 9)
    intervals["cds"] = _with_region(cds, "coding_exon_pas", 10)

    intervals["transcript_end_site"] = transcript_ends
    prev = 0
    for zone in zones:
        name = f"tez_{_format_zone(zone)}"
        intervals[name] = _with_region(
            _derive_downstream_windows(transcript_ends, prev, zone),
            "transcript_end_zone_pas",
            4 if zone <= 250 else 5 if zone <= 1000 else 6,
        )
        prev = zone

    intervals["gene_downstream_2kb"] = _with_region(
        _derive_gene_downstream(genes, gene_downstream_bp),
        "gene_downstream_rescue_pas",
        11,
    )
    return intervals


def _format_zone(zone: int) -> str:
    return f"{zone // 1000}kb" if zone >= 1000 and zone % 1000 == 0 else f"{zone}bp"


def _base_interval_columns(df: pd.DataFrame) -> pd.DataFrame:
    columns = [
        "chrom",
        "start0",
        "end0",
        "strand",
        "gene_id",
        "gene_symbol",
        "gene_biotype",
        "transcript_id",
        "transcript_id_raw",
        "transcript_biotype",
        "model_transcript_end",
        "feature",
    ]
    for col in columns:
        if col not in df.columns:
            df[col] = ""
    out = df[columns].copy()
    return out.rename(columns={"start0": "start", "end0": "end"})


def _feature_table(gtf_df: pd.DataFrame, feature: str) -> pd.DataFrame:
    rows = gtf_df[gtf_df["feature"] == feature].copy()
    if rows.empty:
        return _empty_interval()
    return _base_interval_columns(rows)


def _empty_interval() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "chrom",
            "start",
            "end",
            "strand",
            "gene_id",
            "gene_symbol",
            "gene_biotype",
            "transcript_id",
            "transcript_id_raw",
            "transcript_biotype",
            "model_transcript_end",
            "feature",
        ]
    )


def _build_gene_table(gtf_df: pd.DataFrame) -> pd.DataFrame:
    genes = gtf_df[gtf_df["feature"] == "gene"].copy()
    if genes.empty:
        source = gtf_df[gtf_df["gene_id"].astype(str).str.len() > 0].copy()
        if source.empty:
            return _empty_interval()
        grouped = (
            source.groupby(["chrom", "strand", "gene_id"], dropna=False)
            .agg(
                start=("start0", "min"),
                end=("end0", "max"),
                gene_symbol=("gene_symbol", _first_nonempty),
                gene_biotype=("gene_biotype", _first_nonempty),
            )
            .reset_index()
        )
        grouped["transcript_id"] = ""
        grouped["transcript_id_raw"] = ""
        grouped["transcript_biotype"] = ""
        grouped["model_transcript_end"] = ""
        grouped["feature"] = "gene_reconstructed"
        grouped["gene_span_source"] = "reconstructed_from_transcripts"
        return grouped[_empty_interval().columns.tolist() + ["gene_span_source"]]

    table = _base_interval_columns(genes)
    table["gene_span_source"] = "gene_feature"
    return table[list(_empty_interval().columns) + ["gene_span_source"]]


def _build_transcript_table(gtf_df: pd.DataFrame) -> pd.DataFrame:
    transcripts = gtf_df[gtf_df["feature"] == "transcript"].copy()
    if transcripts.empty:
        source = gtf_df[gtf_df["transcript_id"].astype(str).str.len() > 0].copy()
        if source.empty:
            return _empty_interval()
        grouped = (
            source.groupby(["chrom", "strand", "gene_id", "transcript_id"], dropna=False)
            .agg(
                start=("start0", "min"),
                end=("end0", "max"),
                gene_symbol=("gene_symbol", _first_nonempty),
                gene_biotype=("gene_biotype", _first_nonempty),
                transcript_id_raw=("transcript_id_raw", _first_nonempty),
                transcript_biotype=("transcript_biotype", _first_nonempty),
            )
            .reset_index()
        )
        grouped["feature"] = "transcript_reconstructed"
        grouped["model_transcript_end"] = ""
        return grouped[_empty_interval().columns.tolist()]
    return _base_interval_columns(transcripts)


def _first_nonempty(values) -> str:
    for value in values:
        if isinstance(value, str) and value:
            return value
    return ""


def _with_region(df: pd.DataFrame, region_class: str, priority: int) -> pd.DataFrame:
    df = df.copy()
    if df.empty:
        df["candidate_region_class"] = []
        df["region_priority"] = []
        return df
    df["candidate_region_class"] = region_class
    df["region_priority"] = priority
    return df


def _derive_exon_intervals(
    exons: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if exons.empty:
        return _empty_interval(), _empty_interval(), _empty_interval(), _empty_interval()

    terminal_rows = []
    internal_exon_rows = []
    terminal_introns = []
    internal_introns = []

    for _, tx in exons.sort_values(["chrom", "strand", "transcript_id", "start"]).groupby(
        ["chrom", "strand", "transcript_id"], dropna=False
    ):
        tx = tx.sort_values("start").reset_index(drop=True)
        if tx.empty:
            continue
        terminal_idx = tx["end"].idxmax() if tx.loc[0, "strand"] == "+" else tx["start"].idxmin()
        terminal_rows.append(tx.loc[terminal_idx].to_dict())
        internal_exon_rows.extend(
            row.to_dict() for idx, row in tx.iterrows() if idx != terminal_idx
        )

        if len(tx) < 2:
            continue
        ordered = tx.sort_values("start").reset_index(drop=True)
        introns = []
        for i in range(len(ordered) - 1):
            start = int(ordered.loc[i, "end"])
            end = int(ordered.loc[i + 1, "start"])
            if end <= start:
                continue
            row = ordered.loc[i].to_dict()
            row["start"] = start
            row["end"] = end
            row["feature"] = "intron"
            introns.append(row)
        if not introns:
            continue
        terminal_intron = introns[-1] if ordered.loc[0, "strand"] == "+" else introns[0]
        terminal_introns.append(terminal_intron)
        for intron in introns:
            if intron is not terminal_intron:
                internal_introns.append(intron)

    return (
        pd.DataFrame(terminal_rows, columns=exons.columns),
        pd.DataFrame(internal_exon_rows, columns=exons.columns),
        pd.DataFrame(terminal_introns, columns=exons.columns),
        pd.DataFrame(internal_introns, columns=exons.columns),
    )


def _derive_transcript_ends(transcripts: pd.DataFrame) -> pd.DataFrame:
    if transcripts.empty:
        return _empty_interval()
    rows = []
    for row in transcripts.itertuples(index=False):
        if row.strand == "+":
            site = int(row.end)
        else:
            site = int(row.start)
        rows.append(
            {
                "chrom": row.chrom,
                "start": max(0, site - 1),
                "end": site,
                "strand": row.strand,
                "gene_id": row.gene_id,
                "gene_symbol": row.gene_symbol,
                "transcript_id": row.transcript_id,
                "transcript_id_raw": row.transcript_id_raw,
                "gene_biotype": getattr(row, "gene_biotype", ""),
                "transcript_biotype": getattr(row, "transcript_biotype", ""),
                "model_transcript_end": site,
                "feature": "transcript_end_site",
                "transcript_end": site,
            }
        )
    return pd.DataFrame(rows)


def _attach_model_transcript_ends(
    features: pd.DataFrame, transcript_ends: pd.DataFrame
) -> pd.DataFrame:
    """Attach a separate 3' model distance anchor without changing legacy distances."""
    if features.empty or transcript_ends.empty:
        return features
    lookup = (
        transcript_ends[
            [
                "chrom",
                "strand",
                "gene_id",
                "transcript_id",
                "transcript_end",
            ]
        ]
        .drop_duplicates(["chrom", "strand", "gene_id", "transcript_id"])
        .rename(columns={"transcript_end": "model_transcript_end_new"})
    )
    result = features.merge(
        lookup,
        on=["chrom", "strand", "gene_id", "transcript_id"],
        how="left",
        validate="many_to_one",
    )
    result["model_transcript_end"] = result["model_transcript_end_new"].fillna("")
    return result.drop(columns=["model_transcript_end_new"])


def _derive_downstream_windows(transcript_ends: pd.DataFrame, start_offset: int, end_offset: int) -> pd.DataFrame:
    if transcript_ends.empty:
        return _empty_interval()
    rows = []
    for row in transcript_ends.itertuples(index=False):
        site = int(row.transcript_end)
        if row.strand == "+":
            start = site + start_offset
            end = site + end_offset
        else:
            start = max(0, site - end_offset)
            end = max(0, site - start_offset)
        if end <= start:
            continue
        rows.append(
            {
                **row._asdict(),
                "start": start,
                "end": end,
                "feature": "transcript_end_zone",
                "transcript_end_zone_start_offset": start_offset,
                "transcript_end_zone_end_offset": end_offset,
            }
        )
    return pd.DataFrame(rows)


def _derive_gene_downstream(genes: pd.DataFrame, downstream_bp: int) -> pd.DataFrame:
    if genes.empty:
        return _empty_interval()
    rows = []
    for row in genes.itertuples(index=False):
        if row.strand == "+":
            start = int(row.end)
            end = int(row.end) + downstream_bp
        else:
            start = max(0, int(row.start) - downstream_bp)
            end = int(row.start)
        if end <= start:
            continue
        data = row._asdict()
        data.update({"start": start, "end": end, "feature": "gene_downstream"})
        rows.append(data)
    return pd.DataFrame(rows)
