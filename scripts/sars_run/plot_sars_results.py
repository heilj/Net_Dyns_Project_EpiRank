# -*- coding: utf-8 -*-

"""Plot SARS replication outputs produced by run_for_sars.py."""

import argparse
from pathlib import Path
import sys

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from EpiRank import plotting

DEFAULT_RESULTS_DIR = REPO_ROOT / "results" / "sars"
DEFAULT_FIGURES_DIR = DEFAULT_RESULTS_DIR / "figures"
configure_plot_style = plotting.configure_plot_style
make_metric_barplot = plotting.make_metric_barplot
make_sensitivity_heatmap = plotting.make_sensitivity_heatmap
make_top_scores_plot = plotting.make_top_scores_plot


def load_results(results_dir):
    results_dir = Path(results_dir)
    movement_df = pd.read_csv(results_dir / "movement_comparison.csv")
    baseline_df = pd.read_csv(results_dir / "baseline_comparison.csv")
    sensitivity_df = pd.read_csv(results_dir / "sensitivity.csv")
    top_scores_df = pd.read_csv(results_dir / "top_scores.csv")
    return movement_df, baseline_df, sensitivity_df, top_scores_df


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", default=DEFAULT_RESULTS_DIR, type=Path)
    parser.add_argument("--figures-dir", default=DEFAULT_FIGURES_DIR, type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    movement_df, baseline_df, sensitivity_df, top_scores_df = load_results(args.results_dir)

    args.figures_dir.mkdir(parents=True, exist_ok=True)
    plotting.configure_plot_style()

    combined_df = pd.concat(
        [
            movement_df.assign(group="EpiRank movement setting"),
            baseline_df.assign(group="Baseline"),
        ],
        ignore_index=True,
    )

    plotting.make_metric_barplot(
        combined_df,
        "spearman",
        "Spearman Correlation Across EpiRank Settings and Baselines",
        args.figures_dir / "correlation_spearman.png",
    )
    plotting.make_metric_barplot(
        combined_df,
        "pearson",
        "Pearson Correlation Across EpiRank Settings and Baselines",
        args.figures_dir / "correlation_pearson.png",
    )
    plotting.make_sensitivity_heatmap(
        sensitivity_df,
        "spearman",
        "Sensitivity Heatmap (Spearman)",
        args.figures_dir / "sensitivity_spearman.png",
    )
    plotting.make_sensitivity_heatmap(
        sensitivity_df,
        "pearson",
        "Sensitivity Heatmap (Pearson)",
        args.figures_dir / "sensitivity_pearson.png",
    )
    plotting.make_top_scores_plot(
        top_scores_df,
        args.figures_dir / "top_scores.png",
    )

    print("wrote figures to", args.figures_dir)


if __name__ == "__main__":
    main()
