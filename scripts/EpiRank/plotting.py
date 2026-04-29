# -*- coding: utf-8 -*-

"""Reusable plotting helpers for EpiRank analyses."""

import json
import os
from pathlib import Path
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "epirank-matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(Path(tempfile.gettempdir()) / "epirank-cache"))

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import font_manager
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from . import additional_analysis as aa


CJK_FONT_CANDIDATES = [
    "PingFang TC",
    "PingFang SC",
    "Heiti TC",
    "STHeiti",
    "Arial Unicode MS",
    "Noto Sans CJK TC",
    "Noto Sans CJK SC",
    "Microsoft JhengHei",
]
LEVEL_ORDER = ["NC", "C-III", "C-II", "C-I"]
LEVEL_COLORS = {
    "NC": "#3a8f3e",
    "C-III": "#c8b840",
    "C-II": "#e07830",
    "C-I": "#cc2020",
}
LEVEL_LABELS = {
    "NC": "non-core",
    "C-III": "core-III",
    "C-II": "core-II",
    "C-I": "core-I",
}
MAP_SIZE_FULL = {"C-I": 38, "C-II": 24, "C-III": 14, "NC": 6}
MAP_SIZE_SMALL = {"C-I": 22, "C-II": 14, "C-III": 9, "NC": 4}
MAP_SIZE_ZOOM = {"C-I": 46, "C-II": 30, "C-III": 18, "NC": 8}


def configure_plot_style(context="talk"):
    sns.set_theme(style="whitegrid", context=context)
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font_name in CJK_FONT_CANDIDATES:
        if font_name in available_fonts:
            plt.rcParams["font.family"] = font_name
            break
    plt.rcParams["axes.unicode_minus"] = False


def make_metric_barplot(df, value_col, title, output_path):
    plot_df = df.copy().sort_values(value_col, ascending=False)
    plt.figure(figsize=(8, 4.5))
    hue = "group" if "group" in plot_df.columns else None
    ax = sns.barplot(data=plot_df, x="metric", y=value_col, hue=hue, dodge=False)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(value_col.capitalize())
    ax.tick_params(axis="x", rotation=15)
    ax.set_ylim(min(0.0, plot_df[value_col].min() - 0.05), plot_df[value_col].max() + 0.05)
    if hue is None and ax.get_legend() is not None:
        ax.get_legend().remove()
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def make_sensitivity_heatmap(sensitivity_df, value_col, title, output_path):
    heatmap_df = sensitivity_df.pivot(index="daytime", columns="damping", values=value_col)
    heatmap_df = heatmap_df.sort_index(ascending=False)
    plt.figure(figsize=(9, 6))
    ax = sns.heatmap(
        heatmap_df,
        cmap="viridis",
        annot=False,
        vmin=-1.0,
        vmax=1.0,
        cbar_kws={"label": value_col.capitalize()},
    )
    _highlight_best_heatmap_cell(ax, heatmap_df)
    ax.set_title(title)
    ax.set_xlabel("Damping")
    ax.set_ylabel("Daytime")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def _highlight_best_heatmap_cell(ax, heatmap_df, color="#111111", linewidth=2.0):
    """Draw a bracket-like rectangle around the maximum non-null heatmap cell."""
    values = heatmap_df.to_numpy(dtype=float)
    if values.size == 0 or np.all(np.isnan(values)):
        return
    row, col = np.unravel_index(np.nanargmax(values), values.shape)
    ax.add_patch(
        Rectangle(
            (col, row),
            1,
            1,
            fill=False,
            edgecolor=color,
            linewidth=linewidth,
            clip_on=False,
        )
    )


def _format_disease_label(prefix):
    return {"flu": "Flu", "ev": "EV", "sars": "SARS"}.get(prefix, prefix.upper())


def make_parameter_grid_figure(grid_df, output_path, disease_prefixes=("flu", "ev", "sars")):
    """Parameter-grid sensitivity heatmaps.

    Each panel uses its own data range while still marking the best cell.
    """
    figure_specs = []
    panel_letters = iter("abcdefghijklmnopqrstuvwxyz")
    for disease in disease_prefixes:
        for metric, metric_label in [("pearson", "Pearson's R"), ("spearman", "Spearman's Rho")]:
            value_col = f"{disease}_{metric}"
            if value_col in grid_df:
                figure_specs.append(
                    (value_col, f"({next(panel_letters)}) {_format_disease_label(disease)} ({metric_label})")
                )
    if not figure_specs:
        raise ValueError("grid_df does not include any disease correlation columns")

    ncols = 2
    nrows = int(np.ceil(len(figure_specs) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4.8 * nrows), constrained_layout=True)
    axes = np.asarray(axes).reshape(nrows, ncols)
    for ax, (value_col, title) in zip(axes.ravel(), figure_specs):
        heatmap_df = grid_df.pivot(index="daytime", columns="damping", values=value_col)
        heatmap_df = heatmap_df.sort_index(ascending=True)
        sns.heatmap(
            heatmap_df,
            ax=ax,
            cmap="GnBu",
            cbar=True,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"label": "correlation"},
        )
        _highlight_best_heatmap_cell(ax, heatmap_df)
        ax.set_title(title)
        ax.set_xlabel("damping factor")
        ax.set_ylabel("daytime parameter")
        ax.tick_params(axis="x", labelrotation=90, labelsize=6)
        ax.tick_params(axis="y", labelsize=7)
    for ax in axes.ravel()[len(figure_specs):]:
        ax.set_visible(False)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def make_metric_comparison_table(correlations_df, output_path):
    """Render disease-by-method metric correlations as a PNG table."""
    display_metrics = ["pearson", "spearman", "recall", "precision"]
    table_df = correlations_df.copy()
    table_df["disease"] = table_df["disease"].map(_format_disease_label)
    table_df["metric"] = table_df["metric"].map(
        {
            "epirank": "EpiRank",
            "population": "Population",
            "pagerank": "PageRank",
            "hub_rank": "HITS-Hub",
            "authority_rank": "HITS-Authority",
        }
    ).fillna(table_df["metric"])
    table_df = table_df[["disease", "metric"] + display_metrics]
    table_df[display_metrics] = table_df[display_metrics].astype(float).round(3)

    col_labels = ["Disease", "Method", "Pearson", "Spearman", "Recall", "Precision"]
    cell_text = table_df.to_numpy().tolist()
    row_count = len(cell_text)
    fig_height = max(4.0, 0.34 * row_count + 1.1)
    fig, ax = plt.subplots(figsize=(10.5, fig_height))
    ax.axis("off")
    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        cellLoc="center",
        colLoc="center",
        loc="center",
        colWidths=[0.12, 0.22, 0.16, 0.16, 0.16, 0.16],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.25)

    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#d7d7d7")
        cell.set_linewidth(0.7)
        if row == 0:
            cell.set_facecolor("#222222")
            cell.set_text_props(color="white", weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#f6f7f8")

    for disease in table_df["disease"].unique():
        disease_rows = table_df.index[table_df["disease"] == disease]
        for metric_col in display_metrics:
            best_idx = table_df.loc[disease_rows, metric_col].astype(float).idxmax()
            table[(best_idx + 1, col_labels.index(metric_col.capitalize()))].set_facecolor("#d9ead3")
            table[(best_idx + 1, col_labels.index(metric_col.capitalize()))].set_text_props(weight="bold")

    ax.set_title("Metric Comparison by Disease", fontweight="bold", pad=16)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def make_stability_grid_figure(stability_df, output_path, reference_daytime, reference_damping):
    """Render parameter-stability heatmaps against a reference setting."""
    metric_specs = [
        ("top10_jaccard", "Top-10 overlap"),
        ("top20_jaccard", "Top-20 overlap"),
        ("core_i_jaccard", "Core-I overlap"),
        ("rank_spearman", "Rank Spearman"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), constrained_layout=True)
    for ax, (value_col, title) in zip(axes.ravel(), metric_specs):
        if value_col not in stability_df:
            ax.set_visible(False)
            continue
        heatmap_df = stability_df.pivot(index="daytime", columns="damping", values=value_col)
        heatmap_df = heatmap_df.sort_index(ascending=True)
        sns.heatmap(
            heatmap_df,
            ax=ax,
            cmap="YlGnBu",
            vmin=0.0,
            vmax=1.0,
            xticklabels=True,
            yticklabels=True,
            cbar_kws={"label": value_col.replace("_", " ")},
        )
        ref_row = np.where(np.isclose(heatmap_df.index.to_numpy(dtype=float), float(reference_daytime)))[0]
        ref_col = np.where(np.isclose(heatmap_df.columns.to_numpy(dtype=float), float(reference_damping)))[0]
        if len(ref_row) and len(ref_col):
            ax.add_patch(
                Rectangle(
                    (int(ref_col[0]), int(ref_row[0])),
                    1,
                    1,
                    fill=False,
                    edgecolor="#8b0000",
                    linewidth=2.2,
                    clip_on=False,
                )
            )
        ax.set_title(title)
        ax.set_xlabel("damping factor")
        ax.set_ylabel("daytime parameter")
        ax.tick_params(axis="x", labelrotation=90, labelsize=7)
        ax.tick_params(axis="y", labelsize=8)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def make_top_scores_plot(top_scores_df, output_path):
    plot_df = top_scores_df.copy()
    plot_df["label"] = plot_df["rank"].astype(str) + ". " + plot_df["node"].astype(str)
    plt.figure(figsize=(11, 6))
    ax = sns.barplot(
        data=plot_df,
        y="label",
        x="epirank",
        hue="setting",
        orient="h",
    )
    ax.set_title("Top EpiRank Locations by Movement Setting")
    ax.set_xlabel("EpiRank")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def make_level_count_plot(levels, title, output_path, order=None):
    order = ["NC", "C-III", "C-II", "C-I"] if order is None else order
    level_df = pd.DataFrame({"level": list(levels)})
    plt.figure(figsize=(6, 4))
    ax = sns.countplot(data=level_df, x="level", order=order)
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Townships")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def make_distribution_plot(series_map, title, output_path):
    rows = []
    for name, values in series_map.items():
        rows.extend({"series": name, "value": value} for value in values)
    plot_df = pd.DataFrame(rows)
    plt.figure(figsize=(8, 4.5))
    ax = sns.histplot(data=plot_df, x="value", hue="series", element="step", stat="density")
    ax.set_title(title)
    ax.set_xlabel("Value")
    ax.set_ylabel("Density")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def _iter_geojson_rings(geometry):
    if geometry["type"] == "Polygon":
        for ring in geometry["coordinates"]:
            yield ring
    elif geometry["type"] == "MultiPolygon":
        for polygon in geometry["coordinates"]:
            for ring in polygon:
                yield ring


def draw_taiwan_base_map(ax, geojson_path, facecolor="#f5f5f2", edgecolor="#b8b8b0", linewidth=0.35):
    """Draw a Taiwan GeoJSON base layer under point maps."""
    if geojson_path is None:
        return False
    geojson_path = Path(geojson_path)
    if not geojson_path.exists():
        return False

    with geojson_path.open(encoding="utf-8") as fh:
        feature_collection = json.load(fh)

    for feature in feature_collection.get("features", []):
        geometry = feature.get("geometry")
        if not geometry:
            continue
        for ring in _iter_geojson_rings(geometry):
            coords = np.asarray(ring, dtype=float)
            if coords.ndim == 2 and coords.shape[1] >= 2:
                ax.fill(
                    coords[:, 0],
                    coords[:, 1],
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=linewidth,
                    zorder=0,
                )
    return True


def sibling_township_geojson(geojson_path):
    if geojson_path is None:
        return None
    geojson_path = Path(geojson_path)
    candidate = geojson_path.with_name("taiwan_township.geojson")
    return candidate if candidate.exists() else geojson_path


def _node_rows(graph, values):
    rows = []
    for node, attrs in graph.nodes(data=True):
        if node in values and "posx" in attrs and "posy" in attrs:
            rows.append({"node": node, "x": attrs["posx"], "y": attrs["posy"], "score": values[node]})
    return rows


def _finish_map_axis(ax):
    ax.set_aspect("equal")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off()


def map_extent_from_graph(graph, nodes=None, pad_fraction=0.04):
    nodes = list(graph.nodes()) if nodes is None else list(nodes)
    xs = [graph.nodes[node]["posx"] for node in nodes if node in graph.nodes and "posx" in graph.nodes[node]]
    ys = [graph.nodes[node]["posy"] for node in nodes if node in graph.nodes and "posy" in graph.nodes[node]]
    if not xs or not ys:
        return None
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)
    pad_x = (xmax - xmin) * pad_fraction or 1.0
    pad_y = (ymax - ymin) * pad_fraction or 1.0
    return xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y


def apply_map_extent(ax, extent):
    if extent is None:
        return
    xmin, xmax, ymin, ymax = extent
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)


def make_network_map(
    graph,
    output_path,
    title="Commuting Network",
    edge_alpha=0.08,
    base_map_path=None,
    extent=None,
):
    """Network map: commuting graph visualization used by the GUI Tab 1."""
    pos = {
        node: (attrs["posx"], attrs["posy"])
        for node, attrs in graph.nodes(data=True)
        if "posx" in attrs and "posy" in attrs
    }
    if not pos:
        raise ValueError("graph nodes must include posx and posy attributes")

    plt.figure(figsize=(7, 9))
    ax = plt.gca()
    draw_taiwan_base_map(ax, base_map_path)
    nx.draw_networkx_edges(
        graph,
        pos,
        edge_color="#6b7280",
        alpha=edge_alpha,
        arrows=False,
        width=0.3,
    )
    nx.draw_networkx_nodes(
        graph,
        pos,
        node_size=14,
        node_color="#2563eb",
        alpha=0.85,
        linewidths=0,
    )
    ax.set_title(title)
    apply_map_extent(ax, extent or map_extent_from_graph(graph))
    _finish_map_axis(ax)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def make_score_map(graph, scores, output_path, title="EpiRank Scores", base_map_path=None, extent=None):
    """Plot township scores at their workbook coordinates."""
    rows = _node_rows(graph, scores)
    if not rows:
        raise ValueError("no score rows could be matched to graph coordinates")

    plot_df = pd.DataFrame(rows)
    plt.figure(figsize=(7, 9))
    ax = plt.gca()
    draw_taiwan_base_map(ax, base_map_path)
    ax = sns.scatterplot(
        data=plot_df,
        x="x",
        y="y",
        hue="score",
        size="score",
        palette="viridis",
        sizes=(18, 120),
        linewidth=0,
        legend="brief",
        ax=ax,
    )
    ax.set_title(title)
    apply_map_extent(ax, extent or map_extent_from_graph(graph))
    _finish_map_axis(ax)
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def make_classification_map(
    graph,
    levels,
    output_path,
    title="Risk Classification",
    base_map_path=None,
    size_map=None,
    extent=None,
):
    """Plot head/tail risk classes at township coordinates.

    Dot size follows the GUI convention by default: higher severity classes
    draw larger and above lower classes.
    """
    size_map = MAP_SIZE_FULL if size_map is None else size_map
    rows = []
    for node, attrs in graph.nodes(data=True):
        if node in levels and "posx" in attrs and "posy" in attrs:
            rows.append({"x": attrs["posx"], "y": attrs["posy"], "level": levels[node]})
    if not rows:
        raise ValueError("no level rows could be matched to graph coordinates")

    plt.figure(figsize=(7, 9))
    ax = plt.gca()
    draw_taiwan_base_map(ax, base_map_path)
    for level in LEVEL_ORDER:
        level_rows = [row for row in rows if row["level"] == level]
        if not level_rows:
            continue
        ax.scatter(
            [row["x"] for row in level_rows],
            [row["y"] for row in level_rows],
            s=size_map[level],
            c=LEVEL_COLORS[level],
            marker="o",
            edgecolors="gray",
            linewidths=0.3,
            label=LEVEL_LABELS[level],
            zorder=2 + LEVEL_ORDER.index(level),
        )
    ax.set_title(title)
    apply_map_extent(ax, extent or map_extent_from_graph(graph))
    _finish_map_axis(ax)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.85, title="risk level")
    plt.tight_layout()
    plt.savefig(output_path, dpi=220)
    plt.close()


def make_disease_map_figure(graph, disease_scores, output_path, base_map_path=None):
    """Spatial distributions of disease case severity."""
    disease_items = list(disease_scores.items())
    fig, axes = plt.subplots(1, len(disease_items), figsize=(6 * len(disease_items), 6), constrained_layout=True)
    axes = np.atleast_1d(axes)
    for idx, (ax, (disease, scores)) in enumerate(zip(axes, disease_items)):
        title = f"({chr(97 + idx)}) {_format_disease_label(disease)} case distribution"
        values = list(scores.values())
        levels, _ = aa.classify_dict(scores)
        draw_taiwan_base_map(ax, base_map_path)
        for level in LEVEL_ORDER:
            rows = [
                {"x": attrs["posx"], "y": attrs["posy"]}
                for node, attrs in graph.nodes(data=True)
                if levels.get(node) == level and "posx" in attrs and "posy" in attrs
            ]
            if rows:
                ax.scatter(
                    [row["x"] for row in rows],
                    [row["y"] for row in rows],
                    s=MAP_SIZE_FULL[level],
                    c=LEVEL_COLORS[level],
                    marker="o",
                    edgecolors="gray",
                    linewidths=0.3,
                    label=LEVEL_LABELS[level],
                    zorder=2 + LEVEL_ORDER.index(level),
                )
        ax.set_title(title)
        apply_map_extent(ax, map_extent_from_graph(graph))
        _finish_map_axis(ax)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.8)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def make_epirank_map_figure(graph, scores_by_daytime, output_path, base_map_path=None, zoom_nodes=None):
    """Paper Figure 7: EpiRank maps for daytime 0.0, 0.5, and 1.0.

    Top row is all Taiwan; bottom row is an optional zoom region when
    ``zoom_nodes`` is supplied.
    """
    daytimes = [dt for dt in [0.0, 0.5, 1.0] if dt in scores_by_daytime]
    if not daytimes:
        daytimes = sorted(scores_by_daytime)
    fig, axes = plt.subplots(2, len(daytimes), figsize=(4.8 * len(daytimes), 9), constrained_layout=True)
    if len(daytimes) == 1:
        axes = np.asarray([[axes[0]], [axes[1]]])

    for col, daytime in enumerate(daytimes):
        levels, _ = aa.classify_dict(scores_by_daytime[daytime])
        for row_idx, node_filter in enumerate([None, set(zoom_nodes or [])]):
            ax = axes[row_idx, col]
            map_path = sibling_township_geojson(base_map_path) if row_idx == 1 else base_map_path
            draw_taiwan_base_map(ax, map_path)
            nodes = list(graph.nodes()) if node_filter is None else [node for node in graph.nodes() if node in node_filter]
            size_map = MAP_SIZE_SMALL if row_idx == 0 else MAP_SIZE_ZOOM
            for level in LEVEL_ORDER:
                rows = [
                    graph.nodes[node]
                    for node in nodes
                    if levels.get(node) == level and "posx" in graph.nodes[node] and "posy" in graph.nodes[node]
                ]
                if rows:
                    ax.scatter(
                        [row["posx"] for row in rows],
                        [row["posy"] for row in rows],
                        s=size_map[level],
                        c=LEVEL_COLORS[level],
                        edgecolors="gray",
                        linewidths=0.3,
                        label=LEVEL_LABELS[level],
                        zorder=2 + LEVEL_ORDER.index(level),
                    )
            ax.set_title(f"({'abc'[col] if row_idx == 0 else 'def'[col]}) daytime={daytime:.1f}")
            apply_map_extent(
                ax,
                map_extent_from_graph(graph, nodes if node_filter else None, pad_fraction=0.06),
            )
            _finish_map_axis(ax)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
