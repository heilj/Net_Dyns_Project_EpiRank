# -*- encoding: utf-8 -*-

"""Analysis helpers for EpiRank."""

import copy

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats as st


LEVEL_LABELS = ["NC", "C-III", "C-II", "C-I"]


def htbreak(adic, g=4):
    alist = [v for k, v in adic.items()]
    temp = copy.copy(alist)
    breaks = []
    for i in range(g - 1):
        if not temp:
            break
        avg = sum(temp) / float(len(temp))
        breaks.append(avg)
        temp = [v for v in temp if v > avg]

    adic2 = {}
    for k, v in adic.items():
        lvl = None
        for i in range(len(breaks)):
            if v <= breaks[i]:
                lvl = i
                break
        if lvl is None:
            lvl = len(breaks)
        adic2[k] = lvl

    return adic2, breaks


def head_tail_breaks(values, n_breaks=3):
    """Return ascending head/tail break points for an array-like input."""
    head = np.asarray(values, dtype=float)
    breaks = []
    for _ in range(n_breaks):
        if len(head) < 2:
            break
        mean_value = float(head.mean())
        breaks.append(mean_value)
        head = head[head > mean_value]
    return sorted(breaks)


def classify_by_breaks(values, breaks, labels=None):
    """Classify values with head/tail breaks into four EpiRank levels."""
    labels = LEVEL_LABELS if labels is None else labels
    if len(labels) < 2:
        raise ValueError("at least two labels are required")
    if not breaks:
        return [labels[0]] * len(values)

    breaks = list(breaks)
    while len(breaks) < len(labels) - 1:
        breaks.append(breaks[-1])

    result = []
    for value in values:
        level = len(labels) - 1
        for idx, break_value in enumerate(breaks[: len(labels) - 1]):
            if value <= break_value:
                level = idx
                break
        result.append(labels[level])
    return result


def classify_dict(values, n_breaks=3, labels=None):
    """Classify a dictionary of scores and preserve its keys."""
    keys = list(values.keys())
    breaks = head_tail_breaks([values[k] for k in keys], n_breaks=n_breaks)
    classes = classify_by_breaks([values[k] for k in keys], breaks, labels=labels)
    return dict(zip(keys, classes)), breaks


def calculate_IOratio(g, exclude_selfloop=True):
    g2 = copy.deepcopy(g)
    if exclude_selfloop:
        g2.remove_edges_from(nx.selfloop_edges(g2))

    indeg = dict(g2.in_degree(weight="weight"))
    oudeg = dict(g2.out_degree(weight="weight"))
    io_ratio = {}
    for n in g2.nodes():
        if oudeg[n] > 0 and indeg[n] > 0:
            io_ratio[n] = np.log10(indeg[n] / oudeg[n])
        else:
            io_ratio[n] = 0
    return io_ratio


def calculate_metrices(g, return_dic=True, d=0.95, number_of_loops=1000, weighted=True):
    """Calculate baseline centrality metrics.

    The original project spells this function as ``calculate_metrices``; the
    name is kept for drop-in compatibility.
    """
    if weighted:
        g2 = g
        pagerank_weight = "weight"
    else:
        g2 = copy.deepcopy(g)
        for n1, n2 in g2.edges():
            g2[n1][n2]["weight"] = 1.0
        pagerank_weight = "weight"

    pagerank = nx.pagerank(
        g2,
        alpha=d,
        max_iter=number_of_loops,
        weight=pagerank_weight,
    )
    hub_rank, authority_rank = nx.hits(
        g2,
        max_iter=number_of_loops,
        normalized=True,
    )

    metrices = [pagerank, hub_rank, authority_rank]
    metrices_names = ["pagerank", "hub_rank", "authority_rank"]
    cal_res = list(zip(metrices_names, metrices))
    if return_dic:
        return cal_res

    cal_res2 = {a: b for a, b in cal_res}
    df_res = pd.DataFrame.from_dict(cal_res2)
    return df_res[metrices_names]


def _paired_values(dic1, dic2):
    keys = [k for k in dic1.keys() if k in dic2]
    if len(keys) < 2:
        raise ValueError("at least two common keys are required for correlation")
    numbers1 = [dic1[k] for k in keys]
    numbers2 = [dic2[k] for k in keys]
    return numbers1, numbers2


def _simplify_correlation(cor, pvalue, simplify):
    if not simplify:
        return cor, pvalue

    if pvalue <= 0.001:
        star = "***"
    elif pvalue <= 0.01:
        star = "**"
    elif pvalue <= 0.05:
        star = "*"
    else:
        star = ""
    return round(cor, 4), star


def get_pearson_cor(dic1, dic2, simplify=True):
    numbers1, numbers2 = _paired_values(dic1, dic2)
    cor, pvalue = st.pearsonr(numbers1, numbers2)
    return _simplify_correlation(cor, pvalue, simplify)


def get_spearman_cor(dic1, dic2, simplify=True):
    numbers1, numbers2 = _paired_values(dic1, dic2)
    cor, pvalue = st.spearmanr(numbers1, numbers2)
    return _simplify_correlation(cor, pvalue, simplify)


def get_kendalltau_cor(dic1, dic2, simplify=True):
    numbers1, numbers2 = _paired_values(dic1, dic2)
    cor, pvalue = st.kendalltau(numbers1, numbers2)
    return _simplify_correlation(cor, pvalue, simplify)


def core_precision_recall(predicted, observed, core_labels=None):
    """Precision/recall for classifying core townships.

    By default, every level except ``NC`` is treated as core, matching the
    GUI table logic.
    """
    core_labels = {"C-III", "C-II", "C-I"} if core_labels is None else set(core_labels)
    keys = [key for key in predicted if key in observed]
    tp = sum(predicted[key] in core_labels and observed[key] in core_labels for key in keys)
    fp = sum(predicted[key] in core_labels and observed[key] not in core_labels for key in keys)
    fn = sum(predicted[key] not in core_labels and observed[key] in core_labels for key in keys)
    recall = tp / (tp + fn) if tp + fn else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    return recall, precision


def validate(validate_d, epi, trad):
    validate_pearson = {}
    validate_spearman = {}
    validate_pearson["epirank"] = get_pearson_cor(epi, validate_d)
    validate_spearman["epirank"] = get_spearman_cor(epi, validate_d)
    orders = ["epirank"]
    for n, m in trad:
        validate_pearson[n] = get_pearson_cor(m, validate_d)
        validate_spearman[n] = get_spearman_cor(m, validate_d)
        orders.append(n)
    return validate_pearson, validate_spearman, orders


def prepare_all_table(listtuple_of_dicts, g=None):
    all_dic = {}
    col_order = []
    for n, d in listtuple_of_dicts:
        all_dic[n] = d
        col_order.append(n)
    df_res = pd.DataFrame.from_dict(all_dic)
    return df_res[col_order]
