#!/usr/bin/env python3
"""Render readable state and direct-OT views of the raw-cell coembedding.

This creates two complementary figures from frozen coordinates:
1. a state-flow-group-labelled shared map (all raw cells); and
2. a route atlas of small multiples, so every direct raw-cell OT trajectory is
   visible without stacking all 61 arrows on one panel.
"""

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import numpy as np
import pandas as pd
from PIL import Image


def configure_publication_font():
    """Require Nimbus Sans and retain editable embedded text in vector exports."""
    font_path = Path(font_manager.findfont("Nimbus Sans", fallback_to_default=False))
    if not font_path.is_file():
        raise RuntimeError("Nimbus Sans is required for this publication rendering")
    matplotlib.rcParams.update({
        "font.family": "Nimbus Sans",
        "font.sans-serif": ["Nimbus Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })
    return font_path


def validate_code_counts(observed, expected, source):
    """Require the frozen source taxonomy to equal the reviewed state contract."""
    observed = {str(code): int(count) for code, count in observed.items()}
    expected = {str(code): int(count) for code, count in expected.items()}
    if observed != expected:
        raise ValueError(f"{source} original-state counts differ from the reviewed contract: {observed}")


def species_code_to_group(groups, species):
    """Build a non-overlapping original-state-code to display-group mapping."""
    key = f"{species}_codes"
    mapping = {}
    for group in groups:
        for code in group[key]:
            if code in mapping:
                raise ValueError(f"{species} original state {code} occurs in multiple display groups")
            mapping[code] = group["id"]
    return mapping


def assign_state_groups(nodes, cfg):
    """Assign exact frozen human/mouse state-flow codes to the reviewed 14-group contract."""
    groups = cfg["state_groups"]
    human_code_to_group = species_code_to_group(groups, "human")
    mouse_code_to_group = species_code_to_group(groups, "mouse")
    human = nodes.species_key.eq("human")
    mouse = nodes.species_key.eq("mouse")

    human_map = pd.read_csv(cfg["human_coarse_state_map"], sep="\t", dtype=str)
    if not human_map.state_id.is_unique or human_map.coarse_state.isna().any():
        raise ValueError("Frozen human state-flow map is not one coarse code per original state")
    validate_code_counts(human_map.coarse_state.value_counts().to_dict(),
                         cfg["expected_human_original_state_counts"], "Human")
    human_codes = nodes.loc[human, "original_state_id"].map(
        human_map.set_index("state_id").coarse_state
    )
    if human_codes.isna().any():
        raise ValueError("A human raw cell lacks a frozen state-flow mapping")
    if set(human_codes.unique()) != set(human_code_to_group):
        raise ValueError("Human raw-cell codes do not close against the reviewed display groups")
    nodes.loc[human, "stateflow_code"] = human_codes.to_numpy()
    nodes.loc[human, "state_group"] = human_codes.map(human_code_to_group).to_numpy()

    mouse_source = pd.read_csv(
        cfg["mouse_original_state_crosswalk"], sep="\t", dtype=str,
        usecols=["dataset", "local_state_id", "coarse_from_original_majority_state"],
    )
    code_column = "coarse_from_original_majority_state"
    conflicts = (mouse_source.groupby(["dataset", "local_state_id"])[code_column].nunique() > 1)
    if conflicts.any():
        raise ValueError("Mouse original-state crosswalk has conflicting codes for a dataset/local-state pair")
    mouse_map = mouse_source.drop_duplicates(["dataset", "local_state_id"])[
        ["dataset", "local_state_id", code_column]
    ].rename(columns={code_column: "stateflow_code"})
    if len(mouse_map) != 49 or mouse_map.stateflow_code.isna().any():
        raise ValueError("Mouse frozen original-state crosswalk must contain exactly 49 coded states")
    validate_code_counts(mouse_map.stateflow_code.value_counts().to_dict(),
                         cfg["expected_mouse_original_state_counts"], "Mouse")
    mouse_rows = nodes.loc[mouse, ["dataset", "local_state_id"]].copy()
    mouse_rows["_node_index"] = mouse_rows.index
    mouse_rows = mouse_rows.merge(mouse_map, on=["dataset", "local_state_id"], how="left",
                                  sort=False, validate="many_to_one")
    mouse_rows = mouse_rows.set_index("_node_index").reindex(nodes.index[mouse])
    if mouse_rows.stateflow_code.isna().any():
        raise ValueError("A mouse raw cell lacks a frozen original-state crosswalk code")
    if set(mouse_rows.stateflow_code.unique()) != set(mouse_code_to_group):
        raise ValueError("Mouse raw-cell codes do not close against the reviewed display groups")
    nodes.loc[mouse, "stateflow_code"] = mouse_rows.stateflow_code.to_numpy()
    nodes.loc[mouse, "state_group"] = mouse_rows.stateflow_code.map(mouse_code_to_group).to_numpy()

    if nodes.state_group.isna().any() or set(nodes.state_group) != {group["id"] for group in groups}:
        raise ValueError("Raw-cell display-state assignment does not close")
    return human_map, mouse_map


def label_anchor(xy, limits, bins=90):
    """Put each broad state label on its densest occupied UMAP bin."""
    hist, x_edges, y_edges = np.histogram2d(
        xy[:, 0], xy[:, 1], bins=bins, range=[[limits[0], limits[1]], [limits[2], limits[3]]]
    )
    i, j = np.unravel_index(np.argmax(hist), hist.shape)
    return np.array([(x_edges[i] + x_edges[i + 1]) / 2, (y_edges[j] + y_edges[j + 1]) / 2])


def figure_limits(nodes):
    xy = nodes[["UMAP1", "UMAP2"]].to_numpy(float)
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    pad = .055 * (hi - lo)
    return (lo[0] - pad[0], hi[0] + pad[0], lo[1] - pad[1], hi[1] + pad[1])


def setup_axis(ax, limits):
    ax.set(xlim=limits[:2], ylim=limits[2:])
    ax.set_aspect("equal")
    ax.set_axis_off()


def save_all(fig, output, stem):
    paths = []
    for suffix in ("png", "pdf", "svg"):
        path = output / f"{stem}.{suffix}"
        fig.savefig(path, dpi=300, facecolor="white")
        paths.append(path)
    image = np.asarray(Image.open(output / f"{stem}.png").convert("RGB"))
    y, x = np.where(image.min(axis=2) < 248)
    if not len(x) or x.min() == 0 or y.min() == 0 or x.max() == image.shape[1] - 1 or y.max() == image.shape[0] - 1:
        raise ValueError(f"{stem} has clipped or empty rendered content")
    plt.close(fig)
    return {
        "files": [path.name for path in paths],
        "png_pixels": [int(image.shape[1]), int(image.shape[0])],
        "margins_px": {"left": int(x.min()), "right": int(image.shape[1] - 1 - x.max()),
                       "top": int(y.min()), "bottom": int(image.shape[0] - 1 - y.max())},
    }


def reconcile_trajectory_node_coordinates(arrows):
    """Use one display coordinate per OT node while retaining the frozen edge table.

    Consecutive frozen display edges can carry slightly different endpoint medians
    for the same source/target node.  This reconciles only the plotted endpoint
    coordinate, using the per-route median of that node's incident endpoints;
    it never changes or adds OT edges.
    """
    source_nodes = arrows[["route", "source", "source_x", "source_y"]].rename(
        columns={"source": "node", "source_x": "x", "source_y": "y"}
    )
    target_nodes = arrows[["route", "target", "target_x", "target_y"]].rename(
        columns={"target": "node", "target_x": "x", "target_y": "y"}
    )
    node_coordinates = (pd.concat([source_nodes, target_nodes], ignore_index=True)
                        .groupby(["route", "node"], as_index=False)[["x", "y"]].median())
    display = arrows.merge(
        node_coordinates.rename(columns={"node": "source", "x": "display_source_x", "y": "display_source_y"}),
        on=["route", "source"], how="left", validate="many_to_one",
    ).merge(
        node_coordinates.rename(columns={"node": "target", "x": "display_target_x", "y": "display_target_y"}),
        on=["route", "target"], how="left", validate="many_to_one",
    )
    if len(display) != len(arrows) or display.filter(like="display_").isna().any().any():
        raise ValueError("Could not reconcile every frozen OT endpoint for display")
    source_shift = np.hypot(display.source_x - display.display_source_x,
                            display.source_y - display.display_source_y)
    target_shift = np.hypot(display.target_x - display.display_target_x,
                            display.target_y - display.display_target_y)
    shifts = np.concatenate([source_shift.to_numpy(), target_shift.to_numpy()])
    return display, {
        "coordinate_mode": "per-route median of each frozen OT node's incident endpoint coordinates",
        "route_nodes": int(len(node_coordinates)),
        "endpoints_reconciled": int((shifts > 1e-9).sum()),
        "mean_endpoint_shift": float(shifts.mean()),
        "max_endpoint_shift": float(shifts.max()),
        "edges_added_or_deleted": 0,
    }


def maximal_trajectory_chains(route_edges):
    """Return non-overlapping maximal chains for one frozen route topology."""
    children = route_edges.groupby("source", sort=False).target.apply(list).to_dict()
    parents = route_edges.groupby("target", sort=False).source.apply(list).to_dict()
    route_nodes = set(route_edges.source).union(route_edges.target)
    breakpoints = {node for node in route_nodes if len(parents.get(node, [])) != 1 or len(children.get(node, [])) != 1}
    route_starts = {node for node in route_nodes if not parents.get(node)}
    terminals = {node for node in route_nodes if not children.get(node)}
    chains, visited = [], set()
    for source in sorted(breakpoints):
        for target in sorted(children.get(source, [])):
            edge = (source, target)
            if edge in visited:
                continue
            chain = [source, target]
            visited.add(edge)
            while len(parents.get(chain[-1], [])) == 1 and len(children.get(chain[-1], [])) == 1:
                child = children[chain[-1]][0]
                edge = (chain[-1], child)
                if edge in visited:
                    raise ValueError("Frozen OT route contains a cycle")
                chain.append(child)
                visited.add(edge)
            chains.append(chain)
    expected = set(zip(route_edges.source, route_edges.target))
    if visited != expected:
        raise ValueError("Could not represent every frozen OT edge as one continuous display chain")
    return chains, route_starts, terminals


def render_state_map(nodes, groups, output, seed):
    palette = {x["id"]: x["color"] for x in groups}
    labels = {x["id"]: x["display_label"] for x in groups}
    limits = figure_limits(nodes)
    fig, ax = plt.subplots(figsize=(16, 12))
    fig.subplots_adjust(left=.035, right=.73, bottom=.075, top=.925)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(nodes))
    for state in [x["id"] for x in groups]:
        for species, marker in (("human", "o"), ("mouse", "s")):
            selected = order[(nodes.state_group.to_numpy()[order] == state) &
                             (nodes.species_key.to_numpy()[order] == species)]
            ax.scatter(nodes.UMAP1.to_numpy()[selected], nodes.UMAP2.to_numpy()[selected], s=.38,
                       c=palette[state], marker=marker, alpha=.34, linewidths=0, rasterized=True, zorder=1)
    for state in [x["id"] for x in groups]:
        subset = nodes.loc[nodes.state_group.eq(state), ["UMAP1", "UMAP2"]].to_numpy(float)
        anchor = label_anchor(subset, limits)
        ax.text(*anchor, labels[state], ha="center", va="center", fontsize=8.2, weight="bold",
                color="#18242C", zorder=5,
                bbox={"boxstyle": "round,pad=.24", "facecolor": "white", "edgecolor": palette[state],
                      "linewidth": 1.15, "alpha": .92})
    setup_axis(ax, limits)
    ax.set_title("Human–mouse raw-cell coembedding: reviewed state-flow groups", fontsize=16, pad=18)
    group_handles = [Line2D([], [], marker="o", color=palette[x["id"]], linestyle="", markersize=7,
                             label=x["legend_label"]) for x in groups]
    legend = ax.legend(handles=group_handles, title="State-flow group", frameon=False, loc="upper left",
                       bbox_to_anchor=(1.01, 1.00), fontsize=8.5, title_fontsize=9.5, labelspacing=.82)
    ax.add_artist(legend)
    ax.legend(handles=[Line2D([], [], marker="o", color="#4E5961", linestyle="", label="Human raw cell"),
                       Line2D([], [], marker="s", color="#4E5961", linestyle="", label="Mouse raw cell")],
              title="Point provenance", frameon=False, loc="lower left", bbox_to_anchor=(1.01, .00),
              fontsize=8.5, title_fontsize=9.5)
    fig.text(.38, .018,
             "All 551,201 points are raw cells. State codes follow frozen human and mouse original-state mappings; "
             "shared colour marks the reviewed display group, not cross-species identity.",
             ha="center", fontsize=8.3, color="#4D5660")
    return save_all(fig, output, "human_mouse_raw_cell_stateflow_groups")


def render_route_atlas(nodes, arrows, output, seed):
    limits = figure_limits(nodes)
    routes = sorted(arrows.route.unique())
    colors = plt.get_cmap("tab20").colors
    route_colors = {route: colors[i % len(colors)] for i, route in enumerate(routes)}
    fig, axes = plt.subplots(3, 4, figsize=(17, 13.6))
    fig.subplots_adjust(left=.025, right=.985, bottom=.065, top=.925, wspace=.055, hspace=.14)
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(nodes))
    for ax, route in zip(axes.flat, routes):
        for species, color, marker in (("human", "#BE6681", "o"), ("mouse", "#43879B", "s")):
            selected = order[nodes.species_key.to_numpy()[order] == species]
            ax.scatter(nodes.UMAP1.to_numpy()[selected], nodes.UMAP2.to_numpy()[selected], s=.18, c=color,
                       marker=marker, alpha=.075, linewidths=0, rasterized=True, zorder=1)
        selected = arrows.loc[arrows.route.eq(route)]
        for edge in selected.itertuples(index=False):
            start, stop = (edge.display_source_x, edge.display_source_y), (edge.display_target_x, edge.display_target_y)
            ax.add_patch(FancyArrowPatch(start, stop, arrowstyle="-|>", mutation_scale=11,
                                         connectionstyle="arc3,rad=0.06", linewidth=4.0, color="white",
                                         alpha=.97, shrinkA=0, shrinkB=0, zorder=3))
            ax.add_patch(FancyArrowPatch(start, stop, arrowstyle="-|>", mutation_scale=10,
                                         connectionstyle="arc3,rad=0.06", linewidth=2.0,
                                         color=route_colors[route], alpha=1., shrinkA=0, shrinkB=0, zorder=4))
        setup_axis(ax, limits)
        species = selected.species.iloc[0]
        ax.set_title(f"{route}\n{species} · {len(selected)} direct OT arrows", fontsize=8.7, pad=5)
    for ax in axes.flat[len(routes):]:
        ax.set_visible(False)
    fig.suptitle("Direct raw-cell OT trajectories, separated by study route", fontsize=16, y=.972)
    fig.text(.5, .028,
             "Every panel uses the identical human–mouse coordinate system. Separating routes prevents unrelated study paths "
             "from obscuring one another; arrows are frozen within-species raw-cell OT summaries.",
             ha="center", fontsize=8.6, color="#4D5660")
    review = save_all(fig, output, "human_mouse_raw_cell_direct_ot_route_atlas")
    return review, route_colors


def render_combined_state_trajectory_map(nodes, arrows, groups, output, seed, panel):
    """Overlay direct OT chains on raw cells and keep the complete state key outside the map."""
    palette = {x["id"]: x["color"] for x in groups}
    labels = {x["id"]: x["display_label"] for x in groups}
    trajectory_colors = {"human": "#56335E", "mouse": "#1D4E89"}
    limits = figure_limits(nodes)
    fig = plt.figure(figsize=(panel["width_mm"] / 25.4, panel["height_mm"] / 25.4))
    key_ax = fig.add_axes([.025, .055, .305, .89])
    ax = fig.add_axes([.365, .03, .62, .94])
    key_ax.set_axis_off()
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(nodes))
    for state in [x["id"] for x in groups]:
        for species, marker in (("human", "o"), ("mouse", "s")):
            selected = order[(nodes.state_group.to_numpy()[order] == state) &
                             (nodes.species_key.to_numpy()[order] == species)]
            ax.scatter(nodes.UMAP1.to_numpy()[selected], nodes.UMAP2.to_numpy()[selected], s=.48,
                       c=palette[state], marker=marker, alpha=.42, linewidths=0, rasterized=True, zorder=1)
    terminal_total, chain_total = 0, 0
    route_order = sorted(arrows.route.unique())
    for route in route_order:
        current = arrows.loc[arrows.route.eq(route)]
        chains, starts, terminals = maximal_trajectory_chains(current)
        coordinates = {}
        for edge in current.itertuples(index=False):
            coordinates[edge.source] = (edge.display_source_x, edge.display_source_y)
            coordinates[edge.target] = (edge.display_target_x, edge.display_target_y)
        for chain in chains:
            chain_total += 1
            xy = np.asarray([coordinates[node] for node in chain])
            species = current.species.iloc[0]
            ax.plot(xy[:, 0], xy[:, 1], color="white", linewidth=2.15, alpha=.96, zorder=3,
                    solid_capstyle="round", solid_joinstyle="round")
            ax.plot(xy[:, 0], xy[:, 1], color=trajectory_colors[species], linewidth=.86, alpha=.96, zorder=4,
                    solid_capstyle="round", solid_joinstyle="round")
            if chain[-1] in terminals:
                terminal_total += 1
                for color, scale in (("white", 10.0), (trajectory_colors[species], 7.2)):
                    ax.add_patch(FancyArrowPatch(tuple(xy[-2]), tuple(xy[-1]), arrowstyle="-|>", mutation_scale=scale,
                                                 linewidth=0, color=color, shrinkA=0, shrinkB=0, zorder=5))
        for node in starts:
            start = coordinates[node]
            species = current.species.iloc[0]
            ax.scatter(*start, s=9, c="white", edgecolors=trajectory_colors[species],
                       linewidths=.65, zorder=5)
    for state in [x["id"] for x in groups]:
        subset = nodes.loc[nodes.state_group.eq(state), ["UMAP1", "UMAP2"]].to_numpy(float)
        anchor = label_anchor(subset, limits)
        ax.text(*anchor, labels[state], ha="center", va="center", fontsize=6.5, weight="bold",
                color="#18242C", zorder=6,
                bbox={"boxstyle": "round,pad=.18", "facecolor": "white", "edgecolor": palette[state],
                      "linewidth": .8, "alpha": .94})
    setup_axis(ax, limits)
    key_ax.text(0, 1, "State-flow group", transform=key_ax.transAxes, va="top", fontsize=7.1, weight="bold")
    y = .943
    sections = (
        ("Shared / broad-process groups", [group for group in groups if group["mouse_codes"] and group["human_codes"]]),
        ("Mouse-only groups", [group for group in groups if group["mouse_codes"] and not group["human_codes"]]),
        ("Human-only groups", [group for group in groups if group["human_codes"] and not group["mouse_codes"]]),
    )
    for title, section_groups in sections:
        key_ax.text(0, y, title, transform=key_ax.transAxes, va="top", fontsize=6.5, weight="bold", color="#36454F")
        y -= .042
        for group in section_groups:
            key_ax.plot(.025, y + .008, marker="o", markersize=4.6, markerfacecolor=group["color"],
                        markeredgewidth=0, color=group["color"], transform=key_ax.transAxes, clip_on=False)
            key_ax.text(.085, y, group["legend_label"], transform=key_ax.transAxes, va="baseline", fontsize=6.5)
            y -= .04
        y -= .01
    key_ax.text(0, .185, "Direct raw-cell OT", transform=key_ax.transAxes, va="top", fontsize=7.1, weight="bold")
    for y, species, label in ((.133, "human", "Human direct OT"), (.085, "mouse", "Mouse direct OT")):
        key_ax.plot([.0, .065], [y, y], transform=key_ax.transAxes, color="white", linewidth=2.4,
                    solid_capstyle="round", clip_on=False)
        key_ax.plot([.0, .065], [y, y], transform=key_ax.transAxes, color=trajectory_colors[species], linewidth=1.0,
                    solid_capstyle="round", clip_on=False)
        key_ax.text(.085, y - .01, label, transform=key_ax.transAxes, va="baseline", fontsize=6.5)
    key_ax.plot(.032, .029, marker="o", markersize=4.4, markerfacecolor="white", markeredgecolor="#4E5961",
                markeredgewidth=.7, color="none", transform=key_ax.transAxes, clip_on=False)
    key_ax.text(.085, .02, "Route start; circles human, squares mouse", transform=key_ax.transAxes,
                va="baseline", fontsize=6.5)
    review = save_all(fig, output, "human_mouse_raw_cell_states_with_direct_ot_trajectories")
    review["terminal_arrowheads"] = terminal_total
    review["topological_chain_polylines"] = chain_total
    review["trajectory_style"] = "maximal route-chain polylines with terminal arrowheads"
    return review


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    cfg = json.loads(args.config.read_text())
    font_path = configure_publication_font()
    result = Path(cfg["raw_cell_result"])
    output = result / cfg["output_subdirectory"]
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Refusing to overwrite non-empty output: {output}")
    output.mkdir(parents=True)
    nodes = pd.read_csv(result / "raw_cell_nodes.tsv.gz", sep="\t", low_memory=False)
    arrows = pd.read_csv(result / "raw_cell_trajectory_arrows.tsv", sep="\t")
    if len(nodes) != 551201 or not nodes.cell_id.is_unique:
        raise ValueError("Raw-cell node registry is not closed")
    if len(arrows) != 61 or arrows.route.nunique() != 12:
        raise ValueError("Direct raw-cell trajectory table is not closed")

    human_map, mouse_map = assign_state_groups(nodes, cfg)
    groups = cfg["state_groups"]
    order = [group["id"] for group in groups]
    summary = (nodes.groupby(["state_group", "species_key"]).size().unstack(fill_value=0)
               .reindex(order, fill_value=0).reindex(columns=["human", "mouse"], fill_value=0).reset_index())
    summary["display_label"] = summary.state_group.map({group["id"]: group["display_label"] for group in groups})
    summary["legend_label"] = summary.state_group.map({group["id"]: group["legend_label"] for group in groups})
    summary["total_cells"] = summary[["human", "mouse"]].sum(axis=1)
    summary.to_csv(output / "state_group_summary.tsv", sep="\t", index=False)
    contract = pd.DataFrame([{
        "state_group": group["id"],
        "display_label": group["display_label"],
        "legend_label": group["legend_label"],
        "state_label_cn": group["state_label_cn"],
        "relationship": group["relationship"],
        "mouse_original_plot_code": ",".join(group["mouse_codes"]) or "—",
        "mouse_original_state_count": sum(cfg["expected_mouse_original_state_counts"][code] for code in group["mouse_codes"]),
        "human_original_plot_code": ",".join(group["human_codes"]) or "—",
        "human_original_state_count": sum(cfg["expected_human_original_state_counts"][code] for code in group["human_codes"]),
    } for group in groups])
    contract.to_csv(output / "stateflow_cross_species_group_contract.tsv", sep="\t", index=False)
    nodes[["cell_id", "species_key", "study_id", "library_key", "original_state_id", "local_state_id",
           "stateflow_code", "state_group", "UMAP1", "UMAP2"]].to_csv(
               output / "raw_cell_stateflow_groups.tsv.gz", sep="\t", index=False
           )
    arrows.to_csv(output / "direct_raw_cell_ot_arrows.tsv", sep="\t", index=False)
    display_arrows, coordinate_review = reconcile_trajectory_node_coordinates(arrows)
    display_arrows.to_csv(output / "direct_raw_cell_ot_arrows_display_coordinates.tsv", sep="\t", index=False)
    state_review = render_state_map(nodes, groups, output, cfg["seed"])
    atlas_review, route_colors = render_route_atlas(nodes, display_arrows, output, cfg["seed"])
    combined_review = render_combined_state_trajectory_map(
        nodes, display_arrows, groups, output, cfg["seed"], cfg["combined_panel"]
    )
    result_json = {
        "status": "PASS",
        "analysis_unit": "raw primary cell",
        "raw_cells": int(len(nodes)),
        "metacells_used": False,
        "state_groups": len(groups),
        "state_group_counts": summary.set_index("state_group")[["human", "mouse", "total_cells"]].to_dict("index"),
        "direct_raw_cell_ot_arrows": int(len(arrows)),
        "direct_raw_cell_ot_routes": int(arrows.route.nunique()),
        "trajectory_coordinate_rendering": coordinate_review,
        "trajectory_display_coordinate_table": "direct_raw_cell_ot_arrows_display_coordinates.tsv",
        "publication_font": {"family": "Nimbus Sans", "file": str(font_path)},
        "combined_panel": cfg["combined_panel"],
        "human_state_source": "frozen human integrated-state-flow mapping (62 original states)",
        "mouse_state_source": "frozen mouse original-state crosswalk (49 local states; annotation only)",
        "mouse_stateflow_crosswalk_entries": int(len(mouse_map)),
        "stateflow_group_contract": "stateflow_cross_species_group_contract.tsv",
        "cross_species_state_equivalence_created": False,
        "state_map_render": state_review,
        "route_atlas_render": atlas_review,
        "combined_state_trajectory_render": combined_review,
        "route_colors": {route: "#%02x%02x%02x" % tuple(int(v * 255) for v in color[:3])
                         for route, color in route_colors.items()},
    }
    (output / "summary.json").write_text(json.dumps(result_json, indent=2) + "\n")
    print(json.dumps(result_json, indent=2))


if __name__ == "__main__":
    main()
