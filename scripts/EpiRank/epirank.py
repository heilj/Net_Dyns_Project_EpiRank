# -*- encoding: utf-8 -*-

"""NetworkX-compatible EpiRank implementation."""

from dataclasses import dataclass

import networkx as nx
import numpy as np


@dataclass(frozen=True)
class EpiRankMatrices:
    nodes: list
    cn: np.ndarray
    cnt: np.ndarray


def make_DiGraph(
    df,
    origin_col="origin",
    destination_col="destination",
    flow_col="flow",
    largest_connected_component=True,
    exclude_selfloop=True,
):
    """Build a weighted directed graph from an origin-destination table."""
    g = nx.DiGraph()
    for o, d, f in zip(df[origin_col], df[destination_col], df[flow_col]):
        g.add_edge(o, d, weight=float(f))

    if largest_connected_component and g.number_of_nodes() > 0:
        largest_nodes = max(nx.weakly_connected_components(g), key=len)
        g = g.subgraph(largest_nodes).copy()

    if exclude_selfloop:
        g.remove_edges_from(nx.selfloop_edges(g))

    print(
        "graph construction done,"
        + " no. nodes: "
        + str(g.number_of_nodes())
        + ", no. edges: "
        + str(g.number_of_edges())
    )
    return g


def run_epirank(
    g,
    d=0.95,
    daytime=0.5,
    number_of_loops=1000,
    exfac=None,
    tol=0.0,
    weight="weight",
    verbose=True,
):
    """Calculate EpiRank values for a weighted directed graph.

    ``daytime`` controls the blend between the two movement directions, matching
    the original implementation's update equation.
    """
    if not 0.0 <= d <= 1.0:
        raise ValueError("d must be between 0 and 1")
    if not 0.0 <= daytime <= 1.0:
        raise ValueError("daytime must be between 0 and 1")

    if verbose:
        print("start preparing matrices")
    matrices = prepare_epirank_matrices(g, weight=weight)
    exfac_matrix = get_exfac(exfac, g)
    return run_epirank_prepared(
        matrices,
        d=d,
        daytime=daytime,
        number_of_loops=number_of_loops,
        exfac_matrix=exfac_matrix,
        tol=tol,
        verbose=verbose,
    )


def prepare_epirank_matrices(g, weight="weight"):
    """Prepare normalized movement matrices once for repeated EpiRank runs."""
    nodes = list(g.nodes())
    ncount = len(nodes)
    if ncount == 0:
        return EpiRankMatrices(nodes=nodes, cn=np.empty((0, 0)), cnt=np.empty((0, 0)))

    cn_t = nx.to_numpy_array(g, nodelist=nodes, weight=weight, dtype=float)

    csum = cn_t.sum(axis=0)
    cn = np.divide(
        cn_t,
        csum,
        out=np.zeros_like(cn_t, dtype=float),
        where=csum != 0,
    )

    od_t = cn_t.T
    osum = od_t.sum(axis=0)
    cnt = np.divide(
        od_t,
        osum,
        out=np.zeros_like(od_t, dtype=float),
        where=osum != 0,
    )

    return EpiRankMatrices(nodes=nodes, cn=cn, cnt=cnt)


def run_epirank_prepared(
    matrices,
    d=0.95,
    daytime=0.5,
    number_of_loops=1000,
    exfac_matrix=None,
    tol=0.0,
    verbose=True,
):
    """Calculate EpiRank values from precomputed movement matrices."""
    if not 0.0 <= d <= 1.0:
        raise ValueError("d must be between 0 and 1")
    if not 0.0 <= daytime <= 1.0:
        raise ValueError("daytime must be between 0 and 1")

    nodes = matrices.nodes
    ncount = len(nodes)
    if ncount == 0:
        return {}

    if exfac_matrix is None:
        exfac_matrix = np.ones((ncount, 1), dtype=float) / float(ncount)
    else:
        exfac_matrix = np.asarray(exfac_matrix, dtype=float)
        if exfac_matrix.shape == (ncount,):
            exfac_matrix = exfac_matrix.reshape((ncount, 1))
        if exfac_matrix.shape != (ncount, 1):
            raise ValueError("exfac_matrix must have shape (N, 1)")

    epidemic_risk = np.ones((ncount, 1), dtype=float) / float(ncount)

    if verbose:
        print("preparation done, start iterating")
    iteration = 0
    for iteration in range(number_of_loops):
        old_epidemic_risk = epidemic_risk.copy()
        epidemic_risk = (1.0 - d) * exfac_matrix + d * (
            daytime * matrices.cnt @ epidemic_risk
            + (1.0 - daytime) * matrices.cn @ epidemic_risk
        )
        if tol > 0:
            if np.allclose(epidemic_risk, old_epidemic_risk, atol=tol, rtol=0.0):
                break
        elif np.array_equal(epidemic_risk, old_epidemic_risk):
            break

    if verbose:
        print("epirank calculation done after iteration: " + str(iteration))
    vals = epidemic_risk.ravel().tolist()
    return dict(zip(nodes, vals))


def get_exfac(exfackey, g):
    """Return the external-factor vector as an ``N x 1`` NumPy array."""
    nodes = list(g.nodes())
    ncount = len(nodes)
    if ncount == 0:
        return np.empty((0, 1), dtype=float)

    if exfackey is None:
        values = np.ones(ncount, dtype=float)
    elif isinstance(exfackey, dict):
        values = []
        for n in nodes:
            if n in exfackey:
                v = exfackey[n]
            else:
                print("the dictionary do not contain info for node %s, set as 0" % (n))
                v = 0
            values.append(v)
        values = np.asarray(values, dtype=float)
    elif isinstance(exfackey, str):
        values = []
        for n, attrs in g.nodes(data=True):
            if exfackey in attrs:
                v = attrs[exfackey]
            else:
                print(
                    "the node %s do not have the attribute %s, set as 0"
                    % (n, exfackey)
                )
                v = 0
            values.append(v)
        values = np.asarray(values, dtype=float)
    elif isinstance(exfackey, list):
        values = np.asarray(exfackey[:ncount], dtype=float)
    elif isinstance(exfackey, (np.matrix, np.ndarray)):
        values = np.asarray(exfackey, dtype=float)
        if values.shape == (ncount, 1):
            values = values.ravel()
        elif values.shape == (1, ncount):
            values = values.ravel()
        elif values.shape != (ncount,):
            print(
                "the shape is neither (%s, 1) nor (1, %s), force change to default exfac"
                % (str(ncount), str(ncount))
            )
            return get_exfac(None, g)
    else:
        print("the input is not recognized, force back to default exfac")
        return get_exfac(None, g)

    if values.size != ncount:
        print("the length does not match node count, force change to default exfac")
        return get_exfac(None, g)

    total = float(values.sum())
    if total == 0.0:
        print("the external factor sums to 0, force change to default exfac")
        return get_exfac(None, g)

    return (values / total).reshape((ncount, 1))
