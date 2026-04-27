# -*- coding: utf-8 -*-

"""Plot SARS replication outputs produced by run_for_sars_modern.py."""

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import font_manager
import pandas as pd
import seaborn as sns


REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DEFAULT_RESULTS_DIR = REPO_ROOT / "results" / "sars"
DEFAULT_FIGURES_DIR = DEFAULT_RESULTS_DIR / "figures"
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


def load_results(results_dir):
    results_dir = Path(results_dir)
    movement_df = pd.read_csv(results_dir / "movement_comparison_modern.csv")
    baseline_df = pd.read_csv(results_dir / "baseline_comparison_modern.csv")
    sensitivity_df = pd.read_csv(results_dir / "sensitivity_modern.csv")
    top_scores_df = pd.read_csv(results_dir / "top_scores_modern.csv")
    return movement_df, baseline_df, sensitivity_df, top_scores_df


def configure_plot_style():
    sns.set_theme(style="whitegrid", context="talk")

    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font_name in CJK_FONT_CANDIDATES:
        if font_name in available_fonts:
            plt.rcParams["font.family"] = font_name
            break
    plt.rcParams["axes.unicode_minus"] = False


def make_metric_barplot(df, value_col, title, output_path):
    plot_df = df.copy()
    plot_df = plot_df.sort_values(value_col, ascending=False)

    plt.figure(figsize=(8, 4.5))
    ax = sns.barplot(data=plot_df, x="metric", y=value_col, palette="deep")
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel(value_col.capitalize())
    ax.tick_params(axis="x", rotation=15)
    ax.set_ylim(min(0.0, plot_df[value_col].min() - 0.05), plot_df[value_col].max() + 0.05)
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def make_sensitivity_heatmap(sensitivity_df, value_col, title, output_path):
    heatmap_df = sensitivity_df.pivot(index="daytime", columns="damping", values=value_col)
    heatmap_df = heatmap_df.sort_index(ascending=False)

    plt.figure(figsize=(9, 6))
    ax = sns.heatmap(heatmap_df, cmap="viridis", annot=False, cbar_kws={"label": value_col.capitalize()})
    ax.set_title(title)
    ax.set_xlabel("Damping")
    ax.set_ylabel("Daytime")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def make_top_scores_plot(top_scores_df, output_path):
    plot_df = top_scores_df.copy()
    plot_df["label"] = plot_df["rank"].astype(str) + ". " + plot_df["node"]

    plt.figure(figsize=(11, 6))
    ax = sns.barplot(
        data=plot_df,
        y="label",
        x="epirank",
        hue="setting",
        orient="h",
        palette="deep",
    )
    ax.set_title("Top EpiRank Locations by Movement Setting")
    ax.set_xlabel("EpiRank")
    ax.set_ylabel("")
    plt.tight_layout()
    plt.savefig(output_path, dpi=200)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR, type=Path)
    parser.add_argument("--figures-dir", default=DEFAULT_FIGURES_DIR, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    movement_df, baseline_df, sensitivity_df, top_scores_df = load_results(args.results_dir)

    args.figures_dir.mkdir(parents=True, exist_ok=True)
    configure_plot_style()

    combined_df = pd.concat(
        [
            movement_df.assign(group="EpiRank movement setting"),
            baseline_df.assign(group="Baseline"),
        ],
        ignore_index=True,
    )

    make_metric_barplot(
        combined_df,
        "spearman",
        "Spearman Correlation Across EpiRank Settings and Baselines",
        args.figures_dir / "correlation_spearman_modern.png",
    )
    make_metric_barplot(
        combined_df,
        "pearson",
        "Pearson Correlation Across EpiRank Settings and Baselines",
        args.figures_dir / "correlation_pearson_modern.png",
    )
    make_sensitivity_heatmap(
        sensitivity_df,
        "spearman",
        "Sensitivity Heatmap (Spearman)",
        args.figures_dir / "sensitivity_spearman_modern.png",
    )
    make_sensitivity_heatmap(
        sensitivity_df,
        "pearson",
        "Sensitivity Heatmap (Pearson)",
        args.figures_dir / "sensitivity_pearson_modern.png",
    )
    make_top_scores_plot(
        top_scores_df,
        args.figures_dir / "top_scores_modern.png",
    )

    print("wrote figures to", args.figures_dir)


if __name__ == "__main__":
    main()
