"""
Null model generators.

Five null models are used:

  1. Erdős–Rényi (ER)        : random graph with same N, M
  2. Configuration model     : preserves degree sequence
  3. Watts–Strogatz (WS)     : small-world reference
  4. Barabási–Albert (BA)    : scale-free / preferential attachment
  5. Geometric random graph  : nodes in space, edges by distance

Each generates an ensemble of graphs with the same number of nodes and
edges (or as close as the model allows) as a target graph.

These five span the space of common reference topologies used in network
neuroscience (Bassett & Bullmore 2017; Sporns 2011).
"""

from __future__ import annotations

from typing import List
import numpy as np
import networkx as nx


def _ensure_connected(g: nx.Graph, max_tries: int = 50) -> nx.Graph:
    """Add minimum random edges to make graph connected."""
    rng = np.random.default_rng(0)
    tries = 0
    while not nx.is_connected(g) and tries < max_tries:
        comps = list(nx.connected_components(g))
        a = rng.choice(list(comps[0]))
        b = rng.choice(list(comps[1]))
        g.add_edge(a, b)
        tries += 1
    return g


def er_ensemble(n: int, m: int, k: int = 200,
                seed: int = 0) -> List[nx.Graph]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        g = nx.gnm_random_graph(n, m, seed=int(rng.integers(1e9)))
        if nx.is_connected(g):
            out.append(g)
    return out


def config_ensemble(degree_seq: list, k: int = 200,
                    seed: int = 0) -> List[nx.Graph]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        try:
            g = nx.configuration_model(degree_seq,
                                       seed=int(rng.integers(1e9)))
            g = nx.Graph(g)             # collapse parallels
            g.remove_edges_from(nx.selfloop_edges(g))
            if nx.is_connected(g):
                out.append(g)
        except nx.NetworkXError:
            continue
    return out


def ws_ensemble(n: int, k_avg: int, p: float = 0.1, k: int = 200,
                seed: int = 0) -> List[nx.Graph]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        g = nx.watts_strogatz_graph(n, k_avg, p,
                                    seed=int(rng.integers(1e9)))
        if nx.is_connected(g):
            out.append(g)
    return out


def ba_ensemble(n: int, m_attach: int, k: int = 200,
                seed: int = 0) -> List[nx.Graph]:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(k):
        g = nx.barabasi_albert_graph(n, m_attach,
                                     seed=int(rng.integers(1e9)))
        if nx.is_connected(g):
            out.append(g)
    return out


def geo_ensemble(n: int, target_m: int, k: int = 200,
                 seed: int = 0) -> List[nx.Graph]:
    """Random geometric graphs in [0,1]^2 with radius tuned to target M."""
    rng = np.random.default_rng(seed)
    # Approximate radius via expected degree formula for 2D
    radius = float(np.sqrt(target_m / (np.pi * n * (n - 1) / 2)))
    out = []
    for _ in range(k):
        g = nx.random_geometric_graph(n, radius,
                                      seed=int(rng.integers(1e9)))
        g = _ensure_connected(g)
        out.append(g)
    return out


def all_nulls(target: nx.Graph, k_per: int = 200) -> dict:
    """Return ensembles for all five null models, matched to target."""
    n = target.number_of_nodes()
    m = target.number_of_edges()
    deg_seq = [d for _, d in target.degree()]
    if sum(deg_seq) % 2 != 0: deg_seq[0] += 1
    k_avg = max(2, int(round(2 * m / n)))
    if k_avg % 2: k_avg += 1
    m_attach = max(1, int(round(m / n)))
    return {
        "ER":  er_ensemble(n, m, k=k_per),
        "CFG": config_ensemble(deg_seq, k=k_per),
        "WS":  ws_ensemble(n, k_avg, k=k_per),
        "BA":  ba_ensemble(n, m_attach, k=k_per),
        "GEO": geo_ensemble(n, m, k=k_per),
    }
