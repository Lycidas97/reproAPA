"""One fixed Rad23b interval, one study; compare labelled archived BED sets."""
from pathlib import Path
from datetime import datetime, timezone
import csv
import hashlib
import json
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
INPUTS = HERE / "inputs"
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
PREVIOUS = Path("/root/nfsdata/apadetect_experiments/reproapa_gse99922_20261002/projects/reproapa_gse99922")
FASTA = Path("/root/nfsdata/apadetect/data/genome/mouse/refdata-gex-GRCm39-2024-A/fasta/genome.fa")
SETS = [
    ("CLIPSeq_Ago2_Leung2011a_mm9", "KO_merged", "KO1+KO2 merged; replicate membership unavailable"),
    ("CLIPSeq_Ago2_Leung2011b_mm9", "WT_merged", "WT replicates merged; replicate membership unavailable"),
    ("HITSCLIP_Ago2Leung2011a_mm9", "WT1A", "WT1"),
    ("HITSCLIP_Ago2Leung2011b_mm9", "WT1B", "WT1"),
    ("HITSCLIP_Ago2Leung2011c_mm9", "WT2", "WT2"),
]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_bed(path):
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith(("#", "track", "browser")):
            continue
        fields = line.split("\t")
        assert len(fields) == 6, (path, fields)
        chrom, start, end, name, score, strand = fields
        start, end = int(start), int(end)
        assert 0 <= start < end and strand in ("+", "-")
        rows.append(dict(chrom=chrom, start0=start, end0=end,
                         name=name, score=score, strand=strand))
    return rows


def overlap_bp(query, peak):
    if query["chrom"] != peak["chrom"] or query["strand"] != peak["strand"]:
        return 0
    return max(0, min(query["end0"], peak["end0"]) - max(query["start0"], peak["start0"]))


def write_tsv(rows, path, columns):
    with path.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    protocol = json.loads((HERE / "protocol.json").read_text())
    queries = protocol["queries"]
    assert [q["end0"] - q["start0"] for q in queries] == [415, 7]
    assert overlap_bp(queries[0], queries[1]) == 7
    tools = HERE / "tools"
    tools.mkdir(exist_ok=True)
    sources = {}
    for source, target in [
        (PREVIOUS / "tools/liftOver", tools / "liftOver"),
        (PREVIOUS / "inputs/mm9ToMm10.over.chain.gz", INPUTS / "mm9ToMm10.over.chain.gz"),
        (PREVIOUS / "inputs/mm10ToMm39.over.chain.gz", INPUTS / "mm10ToMm39.over.chain.gz"),
    ]:
        if not target.exists():
            shutil.copy2(source, target)
        assert sha(source) == sha(target)
        sources[str(target.relative_to(HERE))] = {"copied_from": str(source), "sha256": sha(source)}
    (tools / "liftOver").chmod(0o755)
    initial = OUT / "queries.mm39.bed"
    initial.write_text("".join(f'{q["chrom"]}\t{q["start0"]}\t{q["end0"]}\t{q["name"]}\t0\t{q["strand"]}\n' for q in queries))
    previous = initial
    mapping_steps = []
    for tag, chain in [
        ("mm10", "mm39ToMm10.over.chain.gz"),
        ("mm9", "mm10ToMm9.over.chain.gz"),
        ("roundtrip.mm10", "mm9ToMm10.over.chain.gz"),
        ("roundtrip.mm39", "mm10ToMm39.over.chain.gz"),
    ]:
        destination = OUT / f"queries.{tag}.bed"
        unmapped = OUT / f"queries.{tag}.unmapped.bed"
        subprocess.run([str(tools / "liftOver"), "-minMatch=1", "-multiple", str(previous),
                        str(INPUTS / chain), str(destination), str(unmapped)], check=True,
                       capture_output=True, text=True)
        rows = read_bed(destination)
        assert len(rows) == len(queries) and {r["name"] for r in rows} == {q["name"] for q in queries}
        for query in queries:
            mapped = [r for r in rows if r["name"] == query["name"]]
            assert len(mapped) == 1
            mapped = mapped[0]
            assert mapped["chrom"] == query["chrom"] and mapped["strand"] == query["strand"]
            assert mapped["end0"] - mapped["start0"] == query["end0"] - query["start0"]
        assert not unmapped.read_text().strip(), "Unmapped or multiply mapped fixed query"
        mapping_steps.append({"step": tag, "chain": chain, "mapped": rows})
        previous = destination
    coordinate_fields = ("name", "chrom", "start0", "end0", "strand")
    original_coordinates = {tuple(row[key] for key in coordinate_fields) for row in read_bed(initial)}
    roundtrip_coordinates = {tuple(row[key] for key in coordinate_fields) for row in read_bed(previous)}
    assert original_coordinates == roundtrip_coordinates, "Round trip must recover exact fixed coordinates"
    mapped_queries = read_bed(OUT / "queries.mm9.bed")
    inventory, summaries, overlaps = [], [], []
    for prefix, label, unit in SETS:
        metadata = json.loads((INPUTS / (prefix + ".json")).read_text())[0]
        assert str(metadata["PMID"]) == "21258322"
        bed = read_bed(INPUTS / (prefix + ".bed"))
        assert len(bed) == int(metadata["sites"]), "Archived metadata and BED row counts differ"
        inventory.append({"dataset": prefix, "label": label, "biological_unit": unit,
                          "sites": len(bed), "methods": metadata["methods"], "description": metadata["description"]})
        for query in mapped_queries:
            hits = [(peak, overlap_bp(query, peak)) for peak in bed if overlap_bp(query, peak) > 0]
            summaries.append({"query": query["name"], "dataset": prefix, "label": label,
                              "biological_unit": unit, "archive_sites": len(bed), "overlap_clusters": len(hits),
                              "overlap_bp_sum": sum(bp for _, bp in hits)})
            for peak, bp in hits:
                overlaps.append({"query": query["name"], "dataset": prefix, "label": label,
                                 "biological_unit": unit, "peak_chrom": peak["chrom"],
                                 "peak_start0": peak["start0"], "peak_end0": peak["end0"],
                                 "peak_strand": peak["strand"], "archive_score": peak["score"], "overlap_bp": bp})
    write_tsv(inventory, OUT / "archive_inventory.tsv", list(inventory[0]))
    write_tsv(summaries, OUT / "overlap_summary.tsv", list(summaries[0]))
    write_tsv(overlaps, OUT / "overlapping_clusters.tsv", ["query", "dataset", "label", "biological_unit",
              "peak_chrom", "peak_start0", "peak_end0", "peak_strand", "archive_score", "overlap_bp"])
    peak_mapping = []
    if overlaps:
        source = OUT / "overlapping_clusters.mm9.bed"
        source.write_text("".join(f'{r["peak_chrom"]}\t{r["peak_start0"]}\t{r["peak_end0"]}\t{r["label"]}\t0\t{r["peak_strand"]}\n' for r in overlaps))
        for assembly, chain in [("mm10", "mm9ToMm10.over.chain.gz"), ("mm39", "mm10ToMm39.over.chain.gz")]:
            destination = OUT / f"overlapping_clusters.{assembly}.bed"
            unmapped = OUT / f"overlapping_clusters.{assembly}.unmapped.bed"
            subprocess.run([str(tools / "liftOver"), "-minMatch=1", "-multiple", str(source),
                            str(INPUTS / chain), str(destination), str(unmapped)], check=True, capture_output=True)
            mapped = read_bed(destination)
            assert len(mapped) == len(overlaps) and not unmapped.read_text().strip()
            for hit in overlaps:
                records = [r for r in mapped if r["name"] == hit["label"]]
                assert len(records) == 1
                row = records[0]
                assert row["chrom"] == hit["peak_chrom"] and row["strand"] == hit["peak_strand"]
                assert row["end0"] - row["start0"] == hit["peak_end0"] - hit["peak_start0"]
            source = destination
        peak_mapping = read_bed(source)
        write_tsv(peak_mapping, OUT / "overlapping_clusters.mm39.tsv", list(peak_mapping[0]))
    unit_support = []
    for query in queries:
        counts = {r["label"]: r["overlap_clusters"] for r in summaries if r["query"] == query["name"]}
        unit_support.append({"query": query["name"], "WT1_any_technical_library": counts["WT1A"] > 0 or counts["WT1B"] > 0,
                             "WT2": counts["WT2"] > 0,
                             "WT_biological_units_supported": int(counts["WT1A"] > 0 or counts["WT1B"] > 0) + int(counts["WT2"] > 0)})
    write_tsv(unit_support, OUT / "biological_unit_support.tsv", list(unit_support[0]))
    assert Path(str(FASTA) + ".fai").is_file(), "Use an existing index; do not write shared reference inputs"
    fasta_regions = [f'{q["chrom"]}:{q["start0"] + 1}-{q["end0"]}' for q in queries]
    sequence = subprocess.check_output(["samtools", "faidx", str(FASTA), *fasta_regions], text=True)
    (OUT / "queries.mm39.fa").write_text(sequence)
    fasta_sequences = ["".join(part.splitlines()[1:]).upper() for part in sequence.split(">")[1:]]
    old_sequence = json.loads((INPUTS / "mm9_query_sequence.json").read_text())
    assert old_sequence["genome"] == "mm9" and old_sequence["chrom"] == "chr4"
    assert old_sequence["start"] == mapped_queries[0]["start0"] and old_sequence["end"] == mapped_queries[0]["end0"]
    assert old_sequence["dna"].upper() == fasta_sequences[0] and len(fasta_sequences[0]) == 415
    assert fasta_sequences[1] == "AGCACTT" and fasta_sequences[0][41:48] == "AGCACTT"
    result = {"completed_utc": datetime.now(timezone.utc).isoformat(), "server": "213docker",
              "protocol_sha256": sha(HERE / "protocol.json"), "mapping": mapping_steps,
              "archive_sets": inventory, "query_results": summaries, "biological_unit_support": unit_support,
              "overlapping_clusters_mm39": peak_mapping, "mm9_mm39_query_sequence_identical": True,
              "input_sha256": {str(p.relative_to(HERE)): sha(p) for p in sorted(INPUTS.iterdir()) if p.is_file()},
              "copied_tools_and_chains": sources, "script_sha256": sha(Path(__file__)),
              "reference_fasta": str(FASTA), "reference_fai_sha256": sha(Path(str(FASTA) + ".fai")),
              "reference_query_sequence_sha256": sha(OUT / "queries.mm39.fa"),
              "limitations": ["Archived interval boundaries include author-added flanks where preserved; original read bounds cannot be recovered from merged tracks.",
                              "Merged WT/KO tracks do not retain per-biological-replicate membership or full filter-level provenance.",
                              "Per-WT-library archived metadata states KO clusters were removed; no new statistical peak calling or WT-versus-KO differential test.",
                              "Absence of archived overlap is not evidence of absent AGO2 occupancy."]}
    (OUT / "run.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"mapped_mm9": mapped_queries, "overlap_summary": summaries, "unit_support": unit_support}, indent=2))


if __name__ == "__main__":
    main()
