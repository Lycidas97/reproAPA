"""Correct only the scientific families frozen before each module's inference."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

HERE = Path(__file__).resolve().parent
SOURCES = [
    ("background/results/raw-tests.tsv", "p_value"),
    ("cross-species/corrected-tests.tsv", "raw_p"),
    ("cycle/raw-tests.tsv", None),
]
PLANNED = {
    "A1-conditional-reference-redesign": 8,
    "A1-direction-composition-redesign": 4,
    "A3-redesign-cross-species": 8,
    "A2-cycle-alignment-redesign": 16,
    "A8-cycle-PM-redesign": 10,
}


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def bh(p):
    order = np.argsort(p)
    result = np.empty(len(p))
    result[order] = np.minimum(1., np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p)+1))[::-1])[::-1])
    return result


def main():
    inputs, frames = {}, []
    for relative, p_column in SOURCES:
        path = HERE / relative
        inputs[relative] = sha(path)
        f = pd.read_csv(path, sep="\t", float_precision="round_trip", dtype={"seed": "string"})
        if p_column is None:
            p_column = next(c for c in ["raw_p", "p", "p_value"] if c in f)
        f["raw_p"] = pd.to_numeric(f[p_column])
        valid = f.raw_p.notna()
        assert np.isfinite(f.loc[valid, "raw_p"]).all()
        assert f.loc[valid, "raw_p"].between(0., 1.).all()
        assert f.family.notna().all() and f.test_id.notna().all()
        f["source_table"] = relative
        frames.append(f)
    tests = pd.concat(frames, ignore_index=True)
    assert tests.test_id.is_unique
    assert tests.groupby("family").size().to_dict() == PLANNED
    tests["q_BH_family"] = np.nan
    tests["status_for_correction"] = np.where(tests.raw_p.notna(), "estimable", "unestimable_NA")
    summary = []
    for family, f in tests.groupby("family"):
        ix = f.index[f.raw_p.notna()]
        if len(ix):
            corrected = bh(tests.loc[ix, "raw_p"].to_numpy())
            independent = false_discovery_control(tests.loc[ix, "raw_p"].to_numpy())
            np.testing.assert_allclose(corrected, independent, atol=0., rtol=1e-12)
            tests.loc[ix, "q_BH_family"] = corrected
        summary.append(dict(family=family, planned_positions=len(f), estimable=len(ix),
                            unestimable=len(f)-len(ix),
                            rejected_BH_05=int(tests.loc[ix, "q_BH_family"].lt(.05).sum())))
    assert tests.loc[tests.raw_p.isna(), "q_BH_family"].isna().all()
    cross = tests[tests.family.eq("A3-redesign-cross-species")]
    np.testing.assert_allclose(cross.q_BH_family, cross.q_BH_8, atol=0., rtol=1e-12)
    a1 = pd.read_csv(HERE / "background-root-correction.tsv", sep="\t", float_precision="round_trip")
    a1_check = tests.merge(a1[["test_id", "q_BH_family"]], on="test_id", suffixes=("", "_previous"))
    np.testing.assert_allclose(a1_check.q_BH_family, a1_check.q_BH_family_previous, atol=0., rtol=1e-12)
    assert all(sha(HERE/p) == h for p, h in inputs.items())
    tests["rejected_BH_05"] = tests.q_BH_family.lt(.05)
    tests.to_csv(HERE / "corrected-tests.tsv", sep="\t", index=False, na_rep="NA", float_format="%.17g")
    pd.DataFrame(summary).to_csv(HERE / "test-families.tsv", sep="\t", index=False)
    receipt = dict(status="PASS", code_sha256=sha(Path(__file__)), input_sha256=inputs,
                   planned_positions=len(tests), estimable_p=int(tests.raw_p.notna().sum()),
                   scipy_BH_independent_verification=True,
                   A1_previous_root_and_A3_module_correction_identical=True,
                   protocol_sha256=sha(HERE / "PROTOCOL.md"),
                   cycle_addendum_sha256=sha(HERE / "CYCLE-FREEZE-ADDENDUM.md"))
    (HERE / "correction-verification.json").write_text(json.dumps(receipt, indent=2)+"\n")
    print(pd.DataFrame(summary).to_string(index=False))


if __name__ == "__main__":
    main()
