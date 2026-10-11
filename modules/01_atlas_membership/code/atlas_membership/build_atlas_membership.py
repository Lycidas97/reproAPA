"""Build Supplementary Table S3 from aligned libraries and frozen atlas cells."""

import argparse
from collections import Counter, defaultdict
import csv
import gzip
from pathlib import Path


COLUMNS = [
    "study", "species", "sample_id", "library_id", "geo_sample_title",
    "called_cells_before_downstream_qc", "expression_atlas_member",
    "expression_atlas_cells", "analysis_unit", "membership_note",
]
ESC_REFERENCE = "GSM3629848"
IPSC_REFERENCES = {"GSM3195770", "GSM3195771", "GSM3195772", "GSM3195773"}


def read_atlas(path):
    counts = Counter()
    units = defaultdict(set)
    species_counts = Counter()
    cell_ids = set()
    with gzip.open(path, "rt", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            cell_id = row["cell_id"]
            if cell_id in cell_ids:
                raise ValueError(f"Duplicate atlas cell ID: {cell_id}")
            cell_ids.add(cell_id)
            species = row["species_key"]
            sample = row["sample_accession"] if species == "human" else row["library_key"]
            if not sample.startswith("GSM"):
                raise ValueError(f"Invalid sample accession: {sample}")
            unit = row["study_id"].lower() if species == "human" else row["dataset"]
            counts[sample] += 1
            units[sample].add(unit)
            species_counts[species] += 1
    if len(counts) != 161 or species_counts != {"human": 215114, "mouse": 336087}:
        raise ValueError(f"Atlas scope drift: {len(counts)} libraries, {dict(species_counts)} cells")
    if any(len(group) != 1 for group in units.values()):
        raise ValueError("A library has more than one analysis unit")
    return counts, units


def membership_reason(row, included):
    if included:
        return "included in expression atlas"
    sample = row["sample_id"]
    if sample == ESC_REFERENCE:
        return "GSE103221 mature ESC reference outside atlas scope"
    if sample in IPSC_REFERENCES:
        return "WOT mature iPSC reference outside atlas scope"
    if row["study"] == "GSE115943" and "serum" in row["geo_sample_title"].lower():
        return "WOT serum branch outside common-to-2i atlas scope"
    raise ValueError(f"Unexplained library outside atlas: {sample}")


def build(aligned_path, atlas_path):
    counts, units = read_atlas(atlas_path)
    with aligned_path.open(newline="") as stream:
        aligned = list(csv.DictReader(stream, delimiter="\t"))
    keys = [row["sample_id"] for row in aligned]
    if len(aligned) != 210 or len(set(keys)) != 210 or set(counts) - set(keys):
        raise ValueError("Aligned library inventory does not cover the frozen atlas")
    if sum(int(row["called_cells_before_downstream_qc"]) for row in aligned) != 843893:
        raise ValueError("Cell Ranger called-cell total drift")

    result = []
    for row in aligned:
        sample = row["sample_id"]
        count = counts[sample]
        if count > int(row["called_cells_before_downstream_qc"]):
            raise ValueError(f"Atlas cell count exceeds called cells: {sample}")
        result.append({
            "study": row["study"],
            "species": row["species"],
            "sample_id": sample,
            "library_id": row["library_id"],
            "geo_sample_title": row["geo_sample_title"],
            "called_cells_before_downstream_qc": row["called_cells_before_downstream_qc"],
            "expression_atlas_member": sample in counts,
            "expression_atlas_cells": count,
            "analysis_unit": ";".join(sorted(units.get(sample, set()))),
            "membership_note": membership_reason(row, sample in counts),
        })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aligned-libraries", required=True, type=Path)
    parser.add_argument("--atlas-cells", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    rows = build(args.aligned_libraries, args.atlas_cells)
    with args.output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} libraries; {sum(row['expression_atlas_member'] for row in rows)} in atlas; "
          f"{sum(row['expression_atlas_cells'] for row in rows)} atlas cells")


if __name__ == "__main__":
    main()
