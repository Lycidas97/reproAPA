"""Independent input identity, marker recovery and native PAS-axis checks."""
from pathlib import Path
import json

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE = ROOT / "notes/2026-10-10-marker-reference-sensitivity"


def records(bundle, key):
    return pd.DataFrame(json.loads(str(bundle[key])))


def near(a, b):
    assert np.allclose(a, b, atol=1e-12, rtol=0, equal_nan=True)


def main():
    bundle = np.load(SOURCE / "inputs.npz", allow_pickle=False)
    old = np.load(ROOT / "notes/2026-09-26-state-comparability-audit/expression_inputs.npz", allow_pickle=False)
    audit = records(bundle, "audit_json").set_index("study_id")
    mapping = records(bundle, "mapping_json")
    assert mapping.status.eq("unique").all()
    assert not mapping.duplicated(["study_id", "marker"]).any()
    results = []
    for sid, inventory in audit.iterrows():
        meta = records(bundle, sid + "_meta_json")
        selected = records(bundle, sid + "_selected_json")
        f = records(bundle, sid + "_features_json")
        assert meta.metacell_id.is_unique and selected.metacell_id.is_unique
        assert len(meta) == inventory.all_native_mc
        assert len(selected) == inventory.old_selected_mc
        assert meta.in_planned_route.sum() == inventory.planned_route_mc
        assert (meta.in_planned_route & ~meta.original_selected).sum() == inventory.unselected_planned_mc
        assert np.array_equal(meta.matrix_row, np.arange(len(meta)))
        assert np.array_equal(meta.metacell_id.iloc[selected.matrix_row], selected.metacell_id)
        assert meta.original_selected.eq(meta.metacell_id.isin(selected.metacell_id)).all()
        original = meta.set_index("metacell_id").loc[selected.metacell_id]
        assert original.original_node.eq(selected.node.to_numpy()).all()
        near(original.original_purity, selected.purity)
        assert meta.loc[meta.in_planned_route, "route"].notna().all()
        markers = bundle[sid + "_markers"]
        gene_ids = bundle[sid + "_marker_gene_ids"]
        counts = bundle[sid + "_marker_counts"]
        total = bundle[sid + "_total_gex"]
        assert counts.shape == (len(meta), len(markers))
        assert np.issubdtype(counts.dtype, np.integer) and (counts >= 0).all()
        assert len(set(markers)) == len(markers) == len(set(gene_ids))
        assert np.issubdtype(total.dtype, np.integer) and (total > 0).all()
        assert (counts.sum(axis=1) <= total).all()
        near(total, meta.total_gex)
        old_values = 0
        if sid + "_marker_logcp10k" in old.files:
            old_meta = records(old, sid + "_meta_json")
            old_markers = old[inventory.species + "_marker_symbols"]
            common = sorted(set(markers) & set(old_markers))
            rows = pd.Index(meta.metacell_id).get_indexer(old_meta.metacell_id)
            assert (rows >= 0).all()
            near(total[rows], old[sid + "_total_gex"])
            cols = pd.Index(markers).get_indexer(common)
            old_cols = pd.Index(old_markers).get_indexer(common)
            calculated = np.log1p(counts[rows][:, cols] / total[rows, None] * 10000)
            near(calculated, old[sid + "_marker_logcp10k"][:, old_cols])
            old_values = calculated.size
        signed = f.cut_site.to_numpy() * np.where(f.strand.eq("+"), 1, -1)
        s = pd.Series(signed, index=f.index)
        lo = s.groupby(f.primary_gene_id).transform("min").to_numpy()
        hi = s.groupby(f.primary_gene_id).transform("max").to_numpy()
        expected = np.divide(signed - lo, hi - lo, out=np.full(len(f), np.nan), where=hi > lo)
        near(f.position, expected)
        shape = f.groupby("primary_gene_id").agg(n=("cut_site", "nunique"), c=("chrom", "nunique"), s=("strand", "nunique"))
        expected_eligible = f.primary_gene_id.isin(shape.query("n>=2 and c==1 and s==1").index)
        assert f.coordinate_eligible.eq(expected_eligible).all()
        assert f.feature_index_1based.is_unique
        near(bundle[sid + "_pas_shape"], [len(meta), len(f)])
        pointer = bundle[sid + "_pas_indptr"]
        indices = bundle[sid + "_pas_indices"]
        assert len(pointer) == len(f) + 1 and pointer[0] == 0
        assert (np.diff(pointer) >= 0).all()
        assert pointer[-1] == len(indices) == inventory.primary_sparse_nonzero
        assert indices.min() >= 0 and indices.max() < len(meta)
        results.append(dict(study_id=sid, all_native_mc=len(meta), planned_mc=int(meta.in_planned_route.sum()),
                            unselected_planned_mc=int((meta.in_planned_route & ~meta.original_selected).sum()),
                            primary_sites=len(f), old_marker_values_recovered=old_values, status="PASS"))
        del pointer, indices
    output = dict(status="PASS", studies=results, all_native_mc=sum(r["all_native_mc"] for r in results),
                  planned_mc=sum(r["planned_mc"] for r in results),
                  unselected_planned_mc=sum(r["unselected_planned_mc"] for r in results),
                  primary_sites=sum(r["primary_sites"] for r in results),
                  old_marker_values_recovered=sum(r["old_marker_values_recovered"] for r in results))
    (HERE / "input_validation.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output))


if __name__ == "__main__":
    main()
