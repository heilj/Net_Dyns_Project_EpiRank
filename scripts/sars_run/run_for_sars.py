# -*- coding: utf-8 -*-

"""Reproducibility runner for the EpiRank SARS case study.

This script focuses on the project-proposal requirements:
EpiRank movement settings, PageRank/HITS baselines, and sensitivity analysis.
It intentionally does not replace the original notebook.
"""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from EpiRank import additional_analysis as aa
from EpiRank import epirank


MOVEMENT_SETTINGS = [
    ("backward_only", 0.0),
    ("bidirectional", 0.5),
    ("forward_only", 1.0),
]


def _correlation_row(name, values, observed):
    pearson_cor, pearson_p = aa.get_pearson_cor(values, observed, simplify=False)
    spearman_cor, spearman_p = aa.get_spearman_cor(values, observed, simplify=False)
    return {
        "metric": name,
        "pearson": pearson_cor,
        "pearson_p": pearson_p,
        "spearman": spearman_cor,
        "spearman_p": spearman_p,
    }


def load_case_study(data_dir):
    data_dir = Path(data_dir)
    flows = pd.read_csv(data_dir / "od_flow.csv", index_col=0)
    sars_df = pd.read_csv(data_dir / "sars_data.csv")
    sars = dict(zip(sars_df["fullname"], sars_df["sars"]))
    return flows, sars


def run_movement_comparison(matrices, observed, damping):
    rows = []
    scores = {}
    for name, daytime in MOVEMENT_SETTINGS:
        values = epirank.run_epirank_prepared(
            matrices,
            daytime=daytime,
            d=damping,
            verbose=False,
        )
        scores[name] = values
        rows.append(_correlation_row(name, values, observed))
    return pd.DataFrame(rows), scores


def run_baseline_comparison(g, observed, damping):
    rows = []
    for name, values in aa.calculate_metrices(g, d=damping):
        rows.append(_correlation_row(name, values, observed))
    return pd.DataFrame(rows)


def run_sensitivity(matrices, observed, daytime_values, damping_values):
    rows = []
    for daytime in daytime_values:
        for damping in damping_values:
            values = epirank.run_epirank_prepared(
                matrices,
                daytime=float(daytime),
                d=float(damping),
                verbose=False,
            )
            pearson_cor, pearson_p = aa.get_pearson_cor(values, observed, simplify=False)
            spearman_cor, spearman_p = aa.get_spearman_cor(values, observed, simplify=False)
            rows.append(
                {
                    "daytime": float(daytime),
                    "damping": float(damping),
                    "pearson": pearson_cor,
                    "pearson_p": pearson_p,
                    "spearman": spearman_cor,
                    "spearman_p": spearman_p,
                }
            )
    return pd.DataFrame(rows)


def summarize_top_scores(scores, top_n):
    rows = []
    for setting, values in scores.items():
        for rank, (node, value) in enumerate(
            sorted(values.items(), key=lambda item: item[1], reverse=True)[:top_n],
            start=1,
        ):
            rows.append(
                {
                    "setting": setting,
                    "rank": rank,
                    "node": node,
                    "epirank": value,
                }
            )
    return pd.DataFrame(rows)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-dir",
        default=REPO_ROOT / "data_sars",
        type=Path,
    )
    parser.add_argument("--damping", default=0.95, type=float)
    parser.add_argument("--daytime-step", default=0.05, type=float)
    parser.add_argument("--damping-step", default=0.05, type=float)
    parser.add_argument("--top-n", default=10, type=int)
    parser.add_argument("--output-dir", default=REPO_ROOT / "results" / "sars", type=Path)
    return parser.parse_args()


def main():
    args = parse_args()
    flows, sars = load_case_study(args.data_dir)
    g = epirank.make_DiGraph(
        flows,
        origin_col="origin",
        destination_col="destination",
        flow_col="flow",
        largest_connected_component=False,
        exclude_selfloop=False,
    )
    print("start preparing matrices")
    matrices = epirank.prepare_epirank_matrices(g)

    movement_df, movement_scores = run_movement_comparison(matrices, sars, args.damping)
    baseline_df = run_baseline_comparison(g, sars, args.damping)

    daytime_values = np.round(np.arange(0.0, 1.0 + args.daytime_step, args.daytime_step), 10)
    damping_values = np.round(np.arange(args.damping_step, 1.0 + args.damping_step, args.damping_step), 10)
    damping_values = damping_values[damping_values <= 1.0]
    sensitivity_df = run_sensitivity(matrices, sars, daytime_values, damping_values)
    top_scores_df = summarize_top_scores(movement_scores, args.top_n)

    print("\nMovement setting correlations")
    print(movement_df.to_string(index=False))
    print("\nBaseline correlations")
    print(baseline_df.to_string(index=False))
    print("\nBest sensitivity rows by Spearman correlation")
    print(
        sensitivity_df.sort_values("spearman", ascending=False)
        .head(args.top_n)
        .to_string(index=False)
    )

    if args.output_dir is not None:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        movement_df.to_csv(args.output_dir / "movement_comparison.csv", index=False)
        baseline_df.to_csv(args.output_dir / "baseline_comparison.csv", index=False)
        sensitivity_df.to_csv(args.output_dir / "sensitivity.csv", index=False)
        top_scores_df.to_csv(args.output_dir / "top_scores.csv", index=False)


if __name__ == "__main__":
    main()
