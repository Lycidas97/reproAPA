#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pyclustree import clustree

from cluster_utils import write_marker


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--partitions", required=True)
    parser.add_argument("--selected-json", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    target = Path(args.output_dir)
    target.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.partitions, sep="\t", dtype=str)
    keys = [column for column in frame if column.startswith("res_")]
    transitions = []
    for left, right in zip(keys[:-1], keys[1:]):
        counts = frame.groupby([left, right], sort=True).size().rename("n_cells").reset_index()
        counts["source_fraction"] = counts["n_cells"] / counts.groupby(left)["n_cells"].transform("sum")
        counts.insert(0, "from_resolution", left.removeprefix("res_"))
        counts.insert(1, "to_resolution", right.removeprefix("res_"))
        counts = counts.rename(columns={left: "from_cluster", right: "to_cluster"})
        transitions.append(counts)
    pd.concat(transitions, ignore_index=True).to_csv(target / "resolution_transitions.tsv.gz", sep="\t", index=False,
                                                     compression={"method": "gzip", "mtime": 0})
    obj = ad.AnnData(X=np.zeros((len(frame), 1), dtype=np.float32), obs=frame[keys].astype("category"))
    fig = clustree(
        obj,
        cluster_keys=keys,
        title="Leiden resolution transitions",
        show_fraction=False,
        show_cluster_keys=True,
        edge_weight_threshold=0.01,
        node_size_range=(20, 450),
        edge_width_range=(0.15, 3.0),
        x_spacing=0.8,
        y_spacing=6.0,
        graph_plot_kwargs={"font_size": 6, "arrowsize": 5},
    )
    fig.set_size_inches(18, 24)
    fig.savefig(target / "resolution_tree.png", dpi=220, bbox_inches="tight")
    plt.close(fig)
    selected = json.loads(Path(args.selected_json).read_text(encoding="utf-8"))["selected_resolution"]
    values = np.asarray([float(key.removeprefix("res_")) for key in keys])
    selected_index = int(np.argmin(np.abs(values - float(selected))))
    window = keys[max(0, selected_index - 2):min(len(keys), selected_index + 3)]
    focused = clustree(
        obj,
        cluster_keys=window,
        title=f"Selected resolution window ({float(selected):g})",
        show_fraction=True,
        show_cluster_keys=True,
        edge_weight_threshold=0.01,
        node_size_range=(60, 700),
        edge_width_range=(0.2, 4.0),
        x_spacing=1.8,
        y_spacing=4.0,
        graph_plot_kwargs={"font_size": 8, "arrowsize": 7},
    )
    focused.set_size_inches(18, 11)
    focused.savefig(target / "resolution_tree_selected_window.png", dpi=220, bbox_inches="tight")
    plt.close(focused)
    write_marker(target / ".ok", status="PASS", n_cells=len(frame), n_resolutions=len(keys))


if __name__ == "__main__":
    main()
