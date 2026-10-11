"""Independent report checks; no imports from analysis implementations."""
from pathlib import Path
import hashlib
import json
import math

import numpy as np
import pandas as pd


BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def bh(values):
    order = np.argsort(values, kind="stable")
    ranked = np.asarray(values)[order] * len(values) / np.arange(1, len(values) + 1)
    corrected = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    result = np.empty(len(values))
    result[order] = corrected
    return result


def main():
    protocol = json.loads((BASE / "protocol.json").read_text())
    assert sha(BASE / "PROTOCOL.zh.md") == protocol["protocol_sha256"]
    for relative, expected in protocol["inputs_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    manifest = json.loads((BASE / "function/functional_manifest.json").read_text())
    for relative, expected in manifest["inputs_sha256"].items():
        assert sha(ROOT / relative) == expected, relative
    for relative, expected in manifest["outputs_sha256"].items():
        assert sha(BASE / "function" / relative) == expected, relative

    results = pd.read_csv(BASE / "function/category_results.tsv", sep="\t")
    margins = pd.read_csv(BASE / "function/category_stratum_margins.tsv.gz", sep="\t")
    assert len(results) == 144 and results.primary.sum() == 6
    checked = []
    for row in results[results.primary].itertuples(index=False):
        pmf = np.array([1.0])
        expected = 0.0
        for layer in margins[margins.test_id == row.test_id].itertuples(index=False):
            denominator = math.comb(layer.N, layer.n)
            distribution = np.zeros(layer.n + 1)
            for k in range(max(0, layer.n - (layer.N - layer.K)), min(layer.n, layer.K) + 1):
                distribution[k] = (
                    math.comb(layer.K, k) * math.comb(layer.N - layer.K, layer.n - k)
                    / denominator
                )
            pmf = np.convolve(pmf, distribution)
            expected += layer.n * layer.K / layer.N
        tail = min(1.0, math.fsum(pmf[row.hits:]))
        assert math.isclose(pmf.sum(), 1, abs_tol=1e-12)
        assert math.isclose(tail, row.p, abs_tol=1e-11)
        assert math.isclose(expected, row.expected, abs_tol=1e-11)
        checked.append({"test_id": row.test_id, "p": tail, "expected": expected})
    for flag, column in [(results.primary, "q_primary6"), (~results.primary, "q_sensitivity138")]:
        np.testing.assert_allclose(bh(results.loc[flag, "p"]), results.loc[flag, column], atol=1e-12)
    np.testing.assert_allclose(bh(results.p), results.q_joint144, atol=1e-12)

    primary = results[results.primary].sort_values(["set_name", "category"])
    primary.to_csv(BASE / "primary_function_results.tsv", sep="\t", index=False)
    status = {
        "status": "PASS_FUNCTIONAL_SUMMARY_AND_FROZEN_INPUTS",
        "protocol_inputs_verified": len(protocol["inputs_sha256"]),
        "functional_outputs_verified": len(manifest["outputs_sha256"]),
        "independently_recomputed_primary_tails": checked,
        "BH_families_verified": [6, 138, 144],
        "scope": "Hash closure and six combinatorial PMFs; full raw-input/GO/144-test audit is in review/",
    }
    (BASE / "root_function_validation.json").write_text(json.dumps(status, indent=2) + "\n")
    print(json.dumps({"status": status["status"], "primary_tests": len(checked)}))


if __name__ == "__main__":
    main()
