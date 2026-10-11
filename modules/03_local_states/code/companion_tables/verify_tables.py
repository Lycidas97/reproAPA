"""Check that the frozen local partitions and annotations share one atlas scope."""

import argparse
from collections import Counter
import csv
from pathlib import Path


def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def state_unit(row):
    state = row["original_state_id"]
    if "::studywide::" in state:
        assert state.startswith(row["study_id"] + "::"), state
        return row["study_id"].lower()
    assert "|" in state, state
    return state.split("|", 1)[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("s3", "s4", "s5", "s7"):
        parser.add_argument("--" + name, required=True, type=Path)
    args = parser.parse_args()
    s3, s4, s5, s7 = (read(getattr(args, name)) for name in ("s3", "s4", "s5", "s7"))

    atlas_cells = sum(int(row["expression_atlas_cells"]) for row in s3)
    assert atlas_cells == 551201
    assert sum(row["expression_atlas_member"] == "True" for row in s3) == 161
    partitions = {row["analysis_unit"]: row for row in s4}
    assert len(s4) == len(partitions) == 11
    assert sum(int(row["n_cells"]) for row in s4) == atlas_cells

    assert len(s5) == len({row["original_state_id"] for row in s5}) == 111
    units = Counter(state_unit(row) for row in s5)
    cells = Counter()
    for row in s5:
        cells[state_unit(row)] += int(row["n_cells"])
    assert set(units) == set(partitions)
    for unit, row in partitions.items():
        assert units[unit] == int(row["n_local_states"]), unit
        assert cells[unit] == int(row["n_cells"]), unit
    assert sum(cells.values()) == atlas_cells
    assert Counter(row["species"] for row in s5) == {"human": 62, "mouse": 49}
    assert len({row["display_group"] for row in s5}) == 14

    program_keys = [(row["program_id"], row["mouse_gene_id"]) for row in s7]
    assert len(s7) == len(set(program_keys)) == 193
    print("PASS: S3–S5 share 551201 cells; 11 partitions, 111 local states "
          "(62 human, 49 mouse), 14 display groups; S7 has 193 unique program-gene rows")


if __name__ == "__main__":
    main()
