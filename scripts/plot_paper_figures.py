# -*- coding: utf-8 -*-

"""Generate EpiRank paper figures from the xlsx-first Taiwan dataset."""

import argparse
import os
from pathlib import Path
import sys
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "epirank-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "epirank-cache"))

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from EpiRank import additional_analysis as aa
from EpiRank import data as epirank_data
from EpiRank import epirank
from EpiRank import plotting
from run_taiwan_dataset import resolve_exfac, run_parameter_grid


def _node_positions(graph):
    return {
        node: (attrs["posx"], attrs["posy"])
        for node, attrs in graph.nodes(data=True)
        if "posx" in attrs and "posy" in attrs
    }


def _rank_from_graph(graph, town_data, key):
    return epirank_data.rank_from_town_data(graph, town_data, key)


def _inter_flows(graph, town_data):
    nodes = list(graph.nodes())
    local = np.array([graph[n][n]["commuter_type1"] if graph.has_edge(n, n) else 0 for n in nodes], dtype=float)
    total_in = np.array([town_data[graph.nodes[n]["db_ID"]][epirank_data.KEY_IN_COMMUTER_TYPE1] for n in nodes], dtype=float)
    total_out = np.array([town_data[graph.nodes[n]["db_ID"]][epirank_data.KEY_OUT_COMMUTER_TYPE1] for n in nodes], dtype=float)
    inter_in = total_in - local
    inter_out = total_out - local
    log_ratio = np.zeros(len(nodes))
    mask = (inter_in > 0) & (inter_out > 0)
    log_ratio[mask] = np.log10(inter_in[mask] / inter_out[mask])
    return nodes, local, inter_in, inter_out, log_ratio


def _save(fig, path):
    fig.savefig(path, dpi=220)
    plt.close(fig)


def _taipei_nodes(graph):
    return [
        node for node, attrs in graph.nodes(data=True)
        if attrs.get("db_ID") in epirank_data.GTAIPEI_DB_IDS
    ]


def _stacked_hist(ax, values_by_level, bins, alpha=0.88):
    bottom = np.zeros(len(bins) - 1)
    centers = bins[:-1] + np.diff(bins) / 2
    widths = np.diff(bins) * 0.92
    stack_order = ["C-I", "C-II", "C-III", "NC"]
    for level in stack_order:
        counts, _ = np.histogram(values_by_level.get(level, []), bins=bins)
        ax.bar(
            centers,
            counts,
            bottom=bottom,
            width=widths,
            color=plotting.LEVEL_COLORS[level],
            edgecolor="white",
            linewidth=0.35,
            alpha=alpha,
            label=level,
        )
        bottom += counts


def plot_figure_2_commuter_flow(graph, town_data, output_path, base_map_path):
    """Paper Figure 2: commuter-flow data summary."""
    nodes, local, inter_in, inter_out, log_ratio = _inter_flows(graph, town_data)
    pos = _node_positions(graph)
    density = np.array([
        town_data[graph.nodes[n]["db_ID"]][epirank_data.KEY_NORMALIZED_DENSITY]
        for n in nodes
    ])

    fig = plt.figure(figsize=(14, 18))
    gs = fig.add_gridspec(4, 3, hspace=0.35, wspace=0.35)
    ax = fig.add_subplot(gs[0:3, 0:2])
    plotting.draw_taiwan_base_map(ax, base_map_path)
    for label, mask, color, size in [
        ("rural", density <= 0.01, "green", 4),
        ("regular", (density > 0.01) & (density <= 0.10), "royalblue", 10),
        ("urbanized", density > 0.10, "purple", 22),
    ]:
        selected = [n for n, keep in zip(nodes, mask) if keep and n in pos]
        ax.scatter([pos[n][0] for n in selected], [pos[n][1] for n in selected], s=size, c=color, alpha=0.7, label=label)
    ax.set_title("(a)")
    ax.legend(fontsize=8)
    plotting.apply_map_extent(ax, plotting.map_extent_from_graph(graph))
    ax.set_axis_off()
    ax.set_aspect("equal")

    axes = [fig.add_subplot(gs[i, 2]) for i in range(3)]
    indeg = np.array([graph.in_degree(n) - (1 if graph.has_edge(n, n) else 0) for n in nodes])
    outdeg = np.array([graph.out_degree(n) - (1 if graph.has_edge(n, n) else 0) for n in nodes])
    axes[0].scatter(indeg, outdeg, s=8, alpha=0.6)
    axes[0].plot([0, max(indeg.max(), outdeg.max())], [0, max(indeg.max(), outdeg.max())], color="gray")
    axes[0].set_title("(b)")
    axes[0].set_xlabel("in degree")
    axes[0].set_ylabel("out degree")
    axes[1].scatter(inter_in / 1e4, inter_out / 1e4, s=8, alpha=0.6)
    axes[1].set_title("(c)")
    axes[1].set_xlabel("in flow x 10^4")
    axes[1].set_ylabel("out flow x 10^4")
    sorted_ratio = np.array(sorted(log_ratio, reverse=True))
    visible = np.abs(sorted_ratio) > 1e-12
    axes[2].bar(np.arange(len(sorted_ratio))[visible], sorted_ratio[visible], width=1.0, color="steelblue", alpha=0.8)
    axes[2].axhline(0, color="black", linewidth=0.5)
    axes[2].set_title("(d)")
    axes[2].set_ylabel("log10(in/out)")

    edge_dist = []
    for u, v, attrs in graph.edges(data=True):
        if u == v or u not in pos or v not in pos:
            continue
        ux, uy = pos[u]
        vx, vy = pos[v]
        edge_dist.append((np.hypot(ux - vx, uy - vy) / 1000.0, attrs["commuter_type1"]))
    dist = np.array([d for d, _ in edge_dist])
    weights = np.array([w for _, w in edge_dist])
    ax_e = fig.add_subplot(gs[3, 0])
    ax_e.hist(dist, bins=np.arange(0, dist.max() + 1, 1), weights=weights / 1e4)
    ax_e.set_title("(e)")
    ax_e.set_xlabel("commuting distance (km)")
    ax_f = fig.add_subplot(gs[3, 1])
    order = np.argsort(dist)
    ax_f.plot(dist[order], np.cumsum(weights[order]) / weights.sum() * 100)
    ax_f.set_title("(f)")
    ax_f.set_xlabel("commuting distance (km)")
    ax_f.set_ylabel("cumulative %")
    ax_g = fig.add_subplot(gs[3, 2])
    ax_g.hist(local, bins=10, color="#d2691e")
    ax_g.set_title("(g)")
    ax_g.set_xlabel("local flow")
    _save(fig, output_path)


def plot_figure_3_disease_frequency(graph, town_data, output_path):
    """Paper Figure 3: disease frequency and in/out-ratio distributions."""
    nodes, _, _, _, log_ratio = _inter_flows(graph, town_data)
    disease_scores = {
        "flu": _rank_from_graph(graph, town_data, epirank_data.KEY_FLU_TOTAL_CASES),
        "ev": _rank_from_graph(graph, town_data, epirank_data.KEY_EV_AVERAGE_CASES),
        "sars": _rank_from_graph(graph, town_data, epirank_data.KEY_SARS_TOTAL_CASES),
    }
    fig, axes = plt.subplots(2, len(disease_scores), figsize=(6 * len(disease_scores), 10), constrained_layout=True)
    for col, (disease, scores) in enumerate(disease_scores.items()):
        ax = axes[0, col]
        levels, breaks = aa.classify_dict(scores)
        values = np.array(list(scores.values()), dtype=float)
        bin_count = 14 if values.max() > 10 else 12
        bins = np.linspace(max(0, values.min() - 0.5), values.max() + 0.5, bin_count)
        values_by_level = {
            level: [scores[n] for n in nodes if levels[n] == level]
            for level in plotting.LEVEL_ORDER
        }
        _stacked_hist(ax, values_by_level, bins)
        for break_value in breaks:
            ax.axvline(break_value, color="black", linestyle="--")
        ax.set_title(f"({chr(97 + col)}) frequency distribution of {disease.upper() if disease == 'sars' else disease} cases")
        ax.set_ylabel("number of townships")

        ax = axes[1, col]
        levels, _ = aa.classify_dict(scores)
        values_by_level = {
            level: [lr for n, lr in zip(nodes, log_ratio) if levels[n] == level]
            for level in plotting.LEVEL_ORDER
        }
        _stacked_hist(ax, values_by_level, np.linspace(-1, 1, 11))
        ax.axvline(0, color="black", linestyle="--")
        ax.set_title(f"({chr(97 + len(disease_scores) + col)}) in/out ratio of {disease.upper() if disease == 'sars' else disease} cases")
        ax.set_xlabel("log10(in/out)")
    handles, labels = axes[0, -1].get_legend_handles_labels()
    label_order = ["C-I", "C-II", "C-III", "NC"]
    order_map = {label: idx for idx, label in enumerate(label_order)}
    pairs = sorted(zip(handles, labels), key=lambda item: order_map.get(item[1], 99))
    axes[0, -1].legend([h for h, _ in pairs], [l for _, l in pairs], fontsize=8)
    _save(fig, output_path)


def plot_figure_6_epirank_distribution(graph, town_data, scores_by_daytime, output_path):
    """Paper Figure 6: EpiRank distributions by daytime parameter."""
    nodes, _, _, _, log_ratio = _inter_flows(graph, town_data)
    fig, axes = plt.subplots(2, len(scores_by_daytime), figsize=(14, 8), constrained_layout=True)
    for col, daytime in enumerate(sorted(scores_by_daytime)):
        scores = scores_by_daytime[daytime]
        levels, breaks = aa.classify_dict(scores)
        values = np.array(list(scores.values()), dtype=float)
        bins = np.linspace(0, values.max() * 1.04, 18)
        values_by_level = {
            level: [scores[n] for n in nodes if levels[n] == level]
            for level in plotting.LEVEL_ORDER
        }
        _stacked_hist(axes[0, col], values_by_level, bins)
        for break_value in breaks:
            axes[0, col].axvline(break_value, color="black", linestyle="--")
        axes[0, col].set_title(f"({'abc'[col]}) daytime={daytime:.1f}")
        ratio_by_level = {
            level: [lr for n, lr in zip(nodes, log_ratio) if levels[n] == level]
            for level in plotting.LEVEL_ORDER
        }
        _stacked_hist(axes[1, col], ratio_by_level, np.linspace(-1, 1, 11))
        axes[1, col].axvline(0, color="black", linestyle="--")
        axes[1, col].set_title(f"({'def'[col]}) daytime={daytime:.1f}")
    _save(fig, output_path)


def plot_figure_8_overlay(graph, disease_scores, epirank_scores, output_path, base_map_path):
    """Paper Figure 8: EpiRank core overlay on disease maps."""
    er_levels, _ = aa.classify_dict(epirank_scores)
    disease_items = list(disease_scores.items())
    fig, axes = plt.subplots(2, len(disease_items), figsize=(6 * len(disease_items), 10), constrained_layout=True)
    axes = np.asarray(axes).reshape(2, len(disease_items))
    zoom_nodes = set(_taipei_nodes(graph))
    for col, (disease, scores) in enumerate(disease_items):
        title = f"({chr(97 + col)}) {disease.upper() if disease == 'sars' else disease} cases"
        levels, _ = aa.classify_dict(scores)
        for row_idx, node_set in enumerate([None, zoom_nodes]):
            ax = axes[row_idx, col]
            map_path = plotting.sibling_township_geojson(base_map_path) if row_idx == 1 else base_map_path
            plotting.draw_taiwan_base_map(ax, map_path)
            nodes_to_plot = list(graph.nodes()) if node_set is None else [n for n in graph.nodes() if n in node_set]
            size_map = plotting.MAP_SIZE_SMALL if row_idx == 0 else plotting.MAP_SIZE_ZOOM
            for level in plotting.LEVEL_ORDER:
                rows = [graph.nodes[n] for n in nodes_to_plot if levels.get(n) == level and "posx" in graph.nodes[n]]
                if rows:
                    ax.scatter(
                        [r["posx"] for r in rows],
                        [r["posy"] for r in rows],
                        s=size_map[level],
                        c=plotting.LEVEL_COLORS[level],
                        marker="o",
                        edgecolors="gray",
                        linewidths=0.25,
                        alpha=0.9,
                    )
            for level in ["C-III", "C-II", "C-I"]:
                rows = [graph.nodes[n] for n in nodes_to_plot if er_levels.get(n) == level and "posx" in graph.nodes[n]]
                if rows:
                    ax.scatter(
                        [r["posx"] for r in rows],
                        [r["posy"] for r in rows],
                        s=size_map[level] * 1.35,
                        facecolors="none",
                        edgecolors="black",
                        linewidths=0.8,
                        zorder=10 + plotting.LEVEL_ORDER.index(level),
                    )
            ax.set_title(title if row_idx == 0 else f"{title} - Taipei Metropolitan Area")
            plotting.apply_map_extent(
                ax,
                plotting.map_extent_from_graph(graph, nodes_to_plot, pad_fraction=0.06),
            )
            ax.set_axis_off()
            ax.set_aspect("equal")
    _save(fig, output_path)


def plot_figure_9_epirank_vs_disease(graph, disease_scores, epirank_scores, output_path):
    """Paper Figure 9: EpiRank classification vs observed disease levels."""
    er_levels, _ = aa.classify_dict(epirank_scores)
    order = ["C-I", "C-II", "C-III", "NC"]
    disease_items = list(disease_scores.items())
    fig, axes = plt.subplots(1, len(disease_items), figsize=(6 * len(disease_items), 6), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for idx, (ax, (disease, scores)) in enumerate(zip(axes, disease_items)):
        observed, _ = aa.classify_dict(scores)
        bottoms = np.zeros(len(order))
        x = np.arange(len(order))
        for predicted in order:
            vals = []
            for actual in order:
                group = [n for n in graph.nodes() if observed[n] == actual]
                vals.append(100 * sum(er_levels[n] == predicted for n in group) / len(group) if group else 0)
            ax.bar(x, vals, bottom=bottoms, color=plotting.LEVEL_COLORS[predicted], label=predicted)
            bottoms += np.array(vals)
        ax.set_xticks(x, order)
        ax.set_ylim(0, 100)
        ax.set_title(f"({chr(97 + idx)}) {disease.upper() if disease == 'sars' else disease} case")
        ax.set_ylabel("percentage")
    axes[-1].legend(fontsize=8)
    _save(fig, output_path)


def plot_figure_10_index_comparison(graph, town_data, metric_scores, output_path):
    """Paper Figure 10: EpiRank, PageRank, HITS-Hub, and HITS-Authority comparison."""
    nodes, _, _, _, log_ratio = _inter_flows(graph, town_data)
    fig, axes = plt.subplots(1, 4, figsize=(16, 5), constrained_layout=True)
    for ax, (title, scores) in zip(axes, metric_scores):
        levels, _ = aa.classify_dict(scores)
        values_by_level = {
            level: [lr for n, lr in zip(nodes, log_ratio) if levels[n] == level]
            for level in plotting.LEVEL_ORDER
        }
        _stacked_hist(ax, values_by_level, np.linspace(-1, 1, 11), alpha=0.86)
        ax.axvline(0, color="black", linestyle="--")
        ax.set_title(title)
        ax.set_xlabel("log10(in/out)")
    handles, labels = axes[-1].get_legend_handles_labels()
    label_order = ["C-I", "C-II", "C-III", "NC"]
    order_map = {label: idx for idx, label in enumerate(label_order)}
    pairs = sorted(zip(handles, labels), key=lambda item: order_map.get(item[1], 99))
    axes[-1].legend([h for h, _ in pairs], [l for _, l in pairs], fontsize=8)
    _save(fig, output_path)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default=REPO_ROOT / "data", type=Path)
    parser.add_argument("--output-dir", default=REPO_ROOT / "results" / "paper_figures", type=Path)
    parser.add_argument("--base-map", default=REPO_ROOT / "data_sars" / "taiwan_county.geojson", type=Path)
    parser.add_argument("--loops", default=5000, type=int)
    parser.add_argument("--grid", action="store_true")
    parser.add_argument("--daytime-step", default=0.05, type=float)
    parser.add_argument("--damping-step", default=0.05, type=float)
    parser.add_argument("--exclude-selfloop", action="store_true")
    parser.add_argument("--exfac", choices=["uniform", "population"], default="uniform")
    return parser.parse_args()


def main():
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plotting.configure_plot_style(context="paper")

    graph, town_data = epirank_data.load_taiwan_dataset(
        args.data_dir,
        exclude_selfloop=args.exclude_selfloop,
    )
    matrices = epirank.prepare_epirank_matrices(graph)
    exfac_matrix = epirank.get_exfac(resolve_exfac(args.exfac, graph, town_data), graph)
    scores_by_daytime = {
        daytime: epirank.run_epirank_prepared(
            matrices,
            d=0.95,
            daytime=daytime,
            number_of_loops=args.loops,
            exfac_matrix=exfac_matrix,
            tol=1e-12,
            verbose=False,
        )
        for daytime in [0.0, 0.5, 1.0]
    }
    disease_scores = {
        "flu": _rank_from_graph(graph, town_data, epirank_data.KEY_FLU_TOTAL_CASES),
        "ev": _rank_from_graph(graph, town_data, epirank_data.KEY_EV_AVERAGE_CASES),
        "sars": _rank_from_graph(graph, town_data, epirank_data.KEY_SARS_TOTAL_CASES),
    }
    baselines = aa.calculate_metrices(graph, d=0.95, number_of_loops=args.loops)

    # Paper Figure 2: commuter-flow data summary.
    plot_figure_2_commuter_flow(graph, town_data, args.output_dir / "figure_2_commuter_flow.png", args.base_map)
    # Paper Figure 3: disease frequency and in/out-ratio distributions.
    plot_figure_3_disease_frequency(graph, town_data, args.output_dir / "figure_3_disease_frequency.png")
    # Paper Figure 4: spatial distributions of Flu and EV case severity.
    plotting.make_disease_map_figure(graph, disease_scores, args.output_dir / "figure_4_disease_map.png", args.base_map)
    # Paper Figure 6: EpiRank distributions by daytime parameter.
    plot_figure_6_epirank_distribution(graph, town_data, scores_by_daytime, args.output_dir / "figure_6_epirank_distribution.png")
    # Paper Figure 7: EpiRank spatial distributions by daytime parameter.
    plotting.make_epirank_map_figure(
        graph,
        scores_by_daytime,
        args.output_dir / "figure_7_epirank_map.png",
        args.base_map,
        zoom_nodes=_taipei_nodes(graph),
    )
    # Paper Figure 8: EpiRank core overlay on observed disease maps.
    plot_figure_8_overlay(graph, disease_scores, scores_by_daytime[0.5], args.output_dir / "figure_8_overlay_map.png", args.base_map)
    # Paper Figure 9: cross-tabulation of observed disease vs EpiRank class.
    plot_figure_9_epirank_vs_disease(graph, disease_scores, scores_by_daytime[0.5], args.output_dir / "figure_9_epirank_vs_disease.png")
    # Paper Figure 10: EpiRank vs PageRank/HITS index comparison.
    metric_scores = [("(a) EpiRank", scores_by_daytime[0.5])] + [(f"({chr(98 + i)}) {name}", values) for i, (name, values) in enumerate(baselines)]
    plot_figure_10_index_comparison(graph, town_data, metric_scores, args.output_dir / "figure_10_index_comparison.png")

    if args.grid:
        # Paper Figure 11: parameter-grid sensitivity analysis.
        daytime_values = np.round(np.arange(0.0, 1.0 + args.daytime_step, args.daytime_step), 10)
        damping_values = np.round(np.arange(args.damping_step, 1.0 + args.damping_step, args.damping_step), 10)
        damping_values = damping_values[damping_values <= 1.0]
        grid_df = run_parameter_grid(
            graph,
            town_data,
            daytime_values,
            damping_values,
            args.loops,
            args.exfac,
        )
        grid_df.to_csv(args.output_dir / "figure_11_parameter_grid.csv", index=False)
        plotting.make_parameter_grid_figure(
            grid_df,
            args.output_dir / "figure_11_parameter_grid.png",
            disease_prefixes=("flu", "ev", "sars"),
        )

    print("wrote paper figures to", args.output_dir)


if __name__ == "__main__":
    main()
