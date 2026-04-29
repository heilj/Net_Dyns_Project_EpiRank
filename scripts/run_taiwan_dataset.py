# -*- coding: utf-8 -*-

"""Run EpiRank on the xlsx-first Taiwan township dataset."""

import argparse
from pathlib import Path
import sys

import pandas as pd
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPTS_DIR = REPO_ROOT / "scripts"
if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))

from EpiRank import additional_analysis as aa
from EpiRank import data as epirank_data
from EpiRank import epirank
from EpiRank import plotting


DEFAULT_STABILITY_TOP_K = (10, 20)


def resolve_exfac(exfac_mode, graph, town_data):
    if exfac_mode == "uniform":
        return None
    if exfac_mode == "population":
        return epirank_data.make_population_exfac(graph, town_data)
    raise ValueError(f"unsupported exfac mode: {exfac_mode}")


def _correlation_row(metric, values, observed):
    pearson, pearson_p = aa.get_pearson_cor(values, observed, simplify=False)
    spearman, spearman_p = aa.get_spearman_cor(values, observed, simplify=False)
    kendall, kendall_p = aa.get_kendalltau_cor(values, observed, simplify=False)
    pred_levels, _ = aa.classify_dict(values)
    observed_levels, _ = aa.classify_dict(observed)
    recall, precision = aa.core_precision_recall(pred_levels, observed_levels)
    return {
        "metric": metric,
        "pearson": pearson,
        "pearson_p": pearson_p,
        "spearman": spearman,
        "spearman_p": spearman_p,
        "kendall": kendall,
        "kendall_p": kendall_p,
        "recall": recall,
        "precision": precision,
    }


def build_result_tables(graph, town_data, damping, daytime, loops, exfac_mode):
    exfac = resolve_exfac(exfac_mode, graph, town_data)
    scores = epirank.run_epirank(
        graph,
        d=damping,
        daytime=daytime,
        number_of_loops=loops,
        exfac=exfac,
        tol=1e-12,
        verbose=False,
    )
    baselines = aa.calculate_metrices(graph, d=damping, number_of_loops=loops)

    flu = epirank_data.rank_from_town_data(
        graph, town_data, epirank_data.KEY_FLU_TOTAL_CASES
    )
    ev = epirank_data.rank_from_town_data(
        graph, town_data, epirank_data.KEY_EV_AVERAGE_CASES
    )
    sars = epirank_data.rank_from_town_data(
        graph, town_data, epirank_data.KEY_SARS_TOTAL_CASES
    )
    population = epirank_data.rank_from_town_data(
        graph, town_data, epirank_data.KEY_POPULATION
    )

    metric_values = [("epirank", scores), ("population", population)] + baselines
    correlation_rows = []
    for disease_name, observed in [("flu", flu), ("ev", ev), ("sars", sars)]:
        for metric_name, values in metric_values:
            row = _correlation_row(metric_name, values, observed)
            row["disease"] = disease_name
            correlation_rows.append(row)

    levels, breaks = aa.classify_dict(scores)
    ranking_rows = []
    for rank, (node, score) in enumerate(epirank_data.sorted_map(scores), start=1):
        attrs = graph.nodes[node]
        town = town_data[attrs[epirank_data.KEY_DB_ID]]
        ranking_rows.append(
            {
                "rank": rank,
                "seq_no": node,
                "post_code": attrs[epirank_data.KEY_POST_CODE],
                "db_ID": attrs[epirank_data.KEY_DB_ID],
                "county": town[epirank_data.KEY_COUNTY],
                "town": town[epirank_data.KEY_TOWN],
                "epirank": score,
                "level": levels[node],
                "population": town[epirank_data.KEY_POPULATION],
                "flu_cases": town.get(epirank_data.KEY_FLU_TOTAL_CASES, 0),
                "ev_cases": town.get(epirank_data.KEY_EV_AVERAGE_CASES, 0),
                "sars_cases": town.get(epirank_data.KEY_SARS_TOTAL_CASES, 0),
                "out_commuters": town[epirank_data.KEY_OUT_COMMUTER_TYPE1],
                "in_commuters": town[epirank_data.KEY_IN_COMMUTER_TYPE1],
            }
        )

    return (
        pd.DataFrame(ranking_rows),
        pd.DataFrame(correlation_rows),
        pd.DataFrame({"break": breaks}),
    )


def _top_k_nodes(values, k):
    return [node for node, _ in epirank_data.sorted_map(values)[:k]]


def _jaccard_similarity(left, right):
    left_set = set(left)
    right_set = set(right)
    union = left_set | right_set
    if not union:
        return 1.0
    return float(len(left_set & right_set) / len(union))


def _rank_positions(values):
    return {node: rank for rank, (node, _) in enumerate(epirank_data.sorted_map(values), start=1)}


def _mean_rank_shift(reference_values, candidate_values, k):
    ref_top = _top_k_nodes(reference_values, k)
    ref_ranks = _rank_positions(reference_values)
    candidate_ranks = _rank_positions(candidate_values)
    shifts = [abs(candidate_ranks[node] - ref_ranks[node]) for node in ref_top if node in candidate_ranks]
    return float(np.mean(shifts)) if shifts else np.nan


def run_parameter_grid(graph, town_data, daytime_values, damping_values, loops, exfac_mode):
    """Run the xlsx disease parameter grid used by paper Figure 11."""
    matrices = epirank.prepare_epirank_matrices(graph)
    exfac_matrix = epirank.get_exfac(resolve_exfac(exfac_mode, graph, town_data), graph)
    observed_by_disease = {
        "flu": epirank_data.rank_from_town_data(
            graph, town_data, epirank_data.KEY_FLU_TOTAL_CASES
        ),
        "ev": epirank_data.rank_from_town_data(
            graph, town_data, epirank_data.KEY_EV_AVERAGE_CASES
        ),
        "sars": epirank_data.rank_from_town_data(
            graph, town_data, epirank_data.KEY_SARS_TOTAL_CASES
        ),
    }

    rows = []
    for daytime in daytime_values:
        for damping in damping_values:
            scores = epirank.run_epirank_prepared(
                matrices,
                d=float(damping),
                daytime=float(daytime),
                number_of_loops=loops,
                exfac_matrix=exfac_matrix,
                tol=1e-12,
                verbose=False,
            )
            row = {"daytime": float(daytime), "damping": float(damping)}
            for disease_name, observed in observed_by_disease.items():
                pearson, pearson_p = aa.get_pearson_cor(scores, observed, simplify=False)
                spearman, spearman_p = aa.get_spearman_cor(scores, observed, simplify=False)
                row[f"{disease_name}_pearson"] = pearson
                row[f"{disease_name}_pearson_p"] = pearson_p
                row[f"{disease_name}_spearman"] = spearman
                row[f"{disease_name}_spearman_p"] = spearman_p
            rows.append(row)
    return pd.DataFrame(rows)


def run_parameter_stability(
    graph,
    town_data,
    daytime_values,
    damping_values,
    loops,
    exfac_mode,
    reference_daytime,
    reference_damping,
    top_ks=DEFAULT_STABILITY_TOP_K,
):
    """Compare EpiRank rankings across parameter settings against a reference."""
    matrices = epirank.prepare_epirank_matrices(graph)
    exfac_matrix = epirank.get_exfac(resolve_exfac(exfac_mode, graph, town_data), graph)

    score_cache = {}
    settings = []
    for daytime in daytime_values:
        for damping in damping_values:
            settings.append((float(daytime), float(damping)))
    reference_key = (float(reference_daytime), float(reference_damping))
    if reference_key not in settings:
        settings.append(reference_key)

    for daytime, damping in settings:
        score_cache[(daytime, damping)] = epirank.run_epirank_prepared(
            matrices,
            d=damping,
            daytime=daytime,
            number_of_loops=loops,
            exfac_matrix=exfac_matrix,
            tol=1e-12,
            verbose=False,
        )

    reference_scores = score_cache[reference_key]
    reference_levels, _ = aa.classify_dict(reference_scores)
    reference_core_i = {node for node, level in reference_levels.items() if level == "C-I"}

    rows = []
    for daytime, damping in settings:
        scores = score_cache[(daytime, damping)]
        row = {
            "daytime": daytime,
            "damping": damping,
            "is_reference": daytime == reference_key[0] and damping == reference_key[1],
        }
        spearman, _ = aa.get_spearman_cor(reference_scores, scores, simplify=False)
        kendall, _ = aa.get_kendalltau_cor(reference_scores, scores, simplify=False)
        row["rank_spearman"] = spearman
        row["rank_kendall"] = kendall

        levels, _ = aa.classify_dict(scores)
        candidate_core_i = {node for node, level in levels.items() if level == "C-I"}
        row["core_i_jaccard"] = _jaccard_similarity(reference_core_i, candidate_core_i)

        for k in top_ks:
            row[f"top{k}_jaccard"] = _jaccard_similarity(
                _top_k_nodes(reference_scores, k),
                _top_k_nodes(scores, k),
            )
            row[f"top{k}_mean_rank_shift"] = _mean_rank_shift(reference_scores, scores, k)
        rows.append(row)

    return pd.DataFrame(rows).sort_values(["daytime", "damping"]).reset_index(drop=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default=REPO_ROOT / "data", type=Path)
    parser.add_argument("--output-dir", default=REPO_ROOT / "results" / "taiwan", type=Path)
    parser.add_argument("--figures-dir", default=REPO_ROOT / "results" / "taiwan" / "figures_population", type=Path)
    parser.add_argument("--damping", default=0.95, type=float)
    parser.add_argument("--daytime", default=0.5, type=float)
    parser.add_argument("--loops", default=5000, type=int)
    parser.add_argument("--parameter-grid", action="store_true")
    parser.add_argument("--stability", action="store_true")
    parser.add_argument("--daytime-step", default=0.05, type=float)
    parser.add_argument("--damping-step", default=0.05, type=float)
    parser.add_argument("--base-map", default=REPO_ROOT / "data_sars" / "taiwan_county.geojson", type=Path)
    parser.add_argument("--exclude-selfloop", action="store_true")
    parser.add_argument("--exfac", choices=["uniform", "population"], default="uniform")
    return parser.parse_args()


def main():
    args = parse_args()
    graph, town_data = epirank_data.load_taiwan_dataset(
        args.data_dir,
        exclude_selfloop=args.exclude_selfloop,
    )
    ranking_df, correlations_df, breaks_df = build_result_tables(
        graph, town_data, args.damping, args.daytime, args.loops, args.exfac
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    ranking_df.to_csv(args.output_dir / "ranking.csv", index=False)
    correlations_df.to_csv(args.output_dir / "correlations.csv", index=False)
    breaks_df.to_csv(args.output_dir / "epirank_breaks.csv", index=False)

    scores = dict(zip(ranking_df["seq_no"], ranking_df["epirank"]))
    levels = dict(zip(ranking_df["seq_no"], ranking_df["level"]))
    args.figures_dir.mkdir(parents=True, exist_ok=True)
    plotting.configure_plot_style()
    plotting.make_score_map(graph, scores, args.figures_dir / "epirank_scores.png", base_map_path=args.base_map)
    plotting.make_classification_map(graph, levels, args.figures_dir / "epirank_levels.png", base_map_path=args.base_map)
    plotting.make_level_count_plot(ranking_df["level"], "EpiRank Level Counts", args.figures_dir / "level_counts.png")
    plotting.make_metric_comparison_table(
        correlations_df,
        args.figures_dir / "correlations_table.png",
    )

    if args.parameter_grid:
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
        grid_df.to_csv(args.output_dir / "parameter_grid.csv", index=False)
        plotting.make_parameter_grid_figure(
            grid_df,
            args.figures_dir / "parameter_grid.png",
            disease_prefixes=("flu", "ev", "sars"),
        )
    if args.stability:
        daytime_values = np.round(np.arange(0.0, 1.0 + args.daytime_step, args.daytime_step), 10)
        damping_values = np.round(np.arange(args.damping_step, 1.0 + args.damping_step, args.damping_step), 10)
        damping_values = damping_values[damping_values <= 1.0]
        stability_df = run_parameter_stability(
            graph,
            town_data,
            daytime_values,
            damping_values,
            args.loops,
            args.exfac,
            args.daytime,
            args.damping,
        )
        stability_df.to_csv(args.output_dir / "parameter_stability.csv", index=False)
        plotting.make_stability_grid_figure(
            stability_df,
            args.figures_dir / "parameter_stability.png",
            reference_daytime=args.daytime,
            reference_damping=args.damping,
        )

    print("loaded graph:", graph.number_of_nodes(), "nodes,", graph.number_of_edges(), "edges")
    print("wrote results to", args.output_dir)
    print("wrote figures to", args.figures_dir)
    print(ranking_df.head(10).to_string(index=False))


if __name__ == "__main__":
    main()
