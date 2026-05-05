"""
Graph-theoretic invariants and similarity measures.

Invariants computed:
  basic        : nodes, edges, density, mean_degree, max_degree
  paths        : diameter, characteristic path length, efficiency
  triangles    : avg_clustering, transitivity, n_triangles
  centrality   : betweenness centrality (mean, max, top-1 normalized)
  spectral     : largest 5 eigenvalues, spectral gap
  small_world  : sigma, omega (Humphries-Gurney; Telesford)
  modularity   : Louvain modularity Q, n_communities
  rich_club    : rich-club coefficient at multiple thresholds
  assortativity: degree assortativity

Similarity measures:
  d_invariant       : normalized distance over the invariants above
  d_spectral        : Frobenius norm of sorted-eigenvalue difference
  d_portrait_div    : Jensen-Shannon divergence of B-matrix portraits
                      (Bagrow & Bollt 2019)

References:
  Watts & Strogatz 1998 — small-world
  Humphries & Gurney 2008 — sigma
  Telesford et al. 2011 — omega
  Bassett & Bullmore 2017 — small-world brain networks
  van den Heuvel & Sporns 2011 — rich club
  Newman 2006 — modularity
  Blondel et al. 2008 — Louvain
  Bagrow & Bollt 2019 — network portrait divergence
"""

from __future__ import annotations

from typing import Dict, List
import warnings

import numpy as np
import networkx as nx

try:
    import community as community_louvain  # python-louvain
    HAS_LOUVAIN = True
except ImportError:
    HAS_LOUVAIN = False


# ---------------------------------------------------------------------------
# Invariants
# ---------------------------------------------------------------------------

def basic(g: nx.Graph) -> Dict[str, float]:
    deg = [d for _, d in g.degree()]
    return {
        "nodes": g.number_of_nodes(),
        "edges": g.number_of_edges(),
        "density": nx.density(g),
        "mean_degree": float(np.mean(deg)),
        "max_degree": int(max(deg)),
        "degree_var": float(np.var(deg)),
    }


def paths(g: nx.Graph) -> Dict[str, float]:
    if not nx.is_connected(g):
        return {"diameter": float("inf"), "char_path": float("inf"),
                "efficiency": 0.0}
    return {
        "diameter": nx.diameter(g),
        "char_path": nx.average_shortest_path_length(g),
        "efficiency": nx.global_efficiency(g),
    }


def triangles(g: nx.Graph) -> Dict[str, float]:
    return {
        "avg_clustering": nx.average_clustering(g),
        "transitivity": nx.transitivity(g),
        "n_triangles": sum(nx.triangles(g).values()) // 3,
    }


def centrality(g: nx.Graph) -> Dict[str, float]:
    bc = nx.betweenness_centrality(g)
    bc_vals = list(bc.values())
    return {
        "mean_betweenness": float(np.mean(bc_vals)),
        "max_betweenness": float(max(bc_vals)),
        "betweenness_centralization": (
            sum(max(bc_vals) - v for v in bc_vals) /
            ((g.number_of_nodes() - 1) * (g.number_of_nodes() - 2) / 2 + 1e-9)
        ),
    }


def spectral(g: nx.Graph) -> Dict[str, float]:
    A = nx.to_numpy_array(g)
    eigs = sorted(np.linalg.eigvals(A).real, reverse=True)
    return {
        "lambda_1": float(eigs[0]),
        "lambda_2": float(eigs[1]) if len(eigs) > 1 else 0.0,
        "spectral_gap": float(eigs[0] - eigs[1]) if len(eigs) > 1 else 0.0,
        "spectral_radius": float(max(abs(e) for e in eigs)),
    }


def assortativity(g: nx.Graph) -> Dict[str, float]:
    return {"assortativity": nx.degree_assortativity_coefficient(g)}


def _equivalent_random(g: nx.Graph, n_samples: int = 10) -> List[nx.Graph]:
    """Generate Erdős-Rényi graphs with same nodes/edges (for small-world)."""
    n = g.number_of_nodes()
    m = g.number_of_edges()
    out = []
    for s in range(n_samples):
        rg = nx.gnm_random_graph(n, m, seed=s)
        if nx.is_connected(rg):
            out.append(rg)
    return out


def _equivalent_lattice(g: nx.Graph) -> nx.Graph:
    """Build a regular ring lattice with similar mean degree (for omega)."""
    n = g.number_of_nodes()
    k = max(2, int(round(2 * g.number_of_edges() / n)))
    if k % 2: k += 1
    return nx.watts_strogatz_graph(n, k, 0.0, seed=0)


def small_world(g: nx.Graph) -> Dict[str, float]:
    """Humphries-Gurney sigma; Telesford omega.

    sigma = (C/C_rand) / (L/L_rand). >>1 indicates small-world.
    omega = L_rand/L - C/C_latt. omega in [-1, 1]; ~0 = small-world,
    -1 = lattice, +1 = random.
    """
    if not nx.is_connected(g) or g.number_of_nodes() < 4:
        return {"sigma": float("nan"), "omega": float("nan")}
    C = nx.average_clustering(g)
    L = nx.average_shortest_path_length(g)
    rands = _equivalent_random(g, n_samples=10)
    if not rands:
        return {"sigma": float("nan"), "omega": float("nan")}
    C_r = float(np.mean([nx.average_clustering(r) for r in rands]))
    L_r = float(np.mean([nx.average_shortest_path_length(r) for r in rands]))
    if C_r == 0 or L == 0:
        sigma = float("nan")
    else:
        sigma = (C / C_r) / (L / L_r)
    try:
        latt = _equivalent_lattice(g)
        C_l = nx.average_clustering(latt)
        omega = L_r / L - C / max(C_l, 1e-9)
    except Exception:
        omega = float("nan")
    return {"sigma": sigma, "omega": omega}


def modularity(g: nx.Graph) -> Dict[str, float]:
    if not HAS_LOUVAIN:
        return {"modularity_Q": float("nan"), "n_communities": -1}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        partition = community_louvain.best_partition(g, random_state=0)
    Q = community_louvain.modularity(partition, g)
    return {
        "modularity_Q": Q,
        "n_communities": len(set(partition.values())),
    }


def rich_club(g: nx.Graph) -> Dict[str, float]:
    """Rich-club coefficient at degree thresholds k = 1, 2, 5.

    phi(k) = 2 E_k / (N_k (N_k - 1))
    where E_k is the number of edges among nodes with degree > k and N_k
    is the count of those nodes (van den Heuvel & Sporns 2011).
    """
    out = {}
    try:
        rc = nx.rich_club_coefficient(g, normalized=False)
        for k in (1, 2, 5):
            out[f"rich_club_k{k}"] = rc.get(k, float("nan"))
    except Exception:
        for k in (1, 2, 5): out[f"rich_club_k{k}"] = float("nan")
    return out


def all_metrics(g: nx.Graph, label: str = "") -> Dict[str, float]:
    out = {"label": label}
    out.update(basic(g))
    out.update(paths(g))
    out.update(triangles(g))
    out.update(centrality(g))
    out.update(spectral(g))
    out.update(assortativity(g))
    out.update(small_world(g))
    out.update(modularity(g))
    out.update(rich_club(g))
    return out


# ---------------------------------------------------------------------------
# Similarity measures
# ---------------------------------------------------------------------------

INVARIANT_KEYS = (
    "density", "mean_degree", "diameter", "char_path", "efficiency",
    "avg_clustering", "transitivity",
    "mean_betweenness", "max_betweenness", "betweenness_centralization",
    "lambda_1", "spectral_gap", "spectral_radius",
    "assortativity", "sigma", "omega",
    "modularity_Q",
)


def d_invariant(a: Dict[str, float], b: Dict[str, float]) -> float:
    """Normalized euclidean over the invariant vector."""
    diffs = []
    for k in INVARIANT_KEYS:
        va, vb = a.get(k), b.get(k)
        if va is None or vb is None: continue
        if any(np.isinf([va, vb])) or any(np.isnan([va, vb])): continue
        denom = max(abs(va), abs(vb), 1e-9)
        diffs.append(((va - vb) / denom) ** 2)
    if not diffs:
        return float("nan")
    return float(np.sqrt(np.mean(diffs)))


def d_spectral(g1: nx.Graph, g2: nx.Graph) -> float:
    """Distance between sorted spectra (zero-padded to common length)."""
    s1 = sorted(np.linalg.eigvals(nx.to_numpy_array(g1)).real, reverse=True)
    s2 = sorted(np.linalg.eigvals(nx.to_numpy_array(g2)).real, reverse=True)
    n = max(len(s1), len(s2))
    s1 = np.array(s1 + [0.0] * (n - len(s1)))
    s2 = np.array(s2 + [0.0] * (n - len(s2)))
    return float(np.linalg.norm(s1 - s2) / np.sqrt(n))


def _network_portrait(g: nx.Graph) -> np.ndarray:
    """B-matrix portrait (Bagrow & Bollt 2019).

    B[l, k] = number of nodes with k neighbors at shortest-path distance l.
    """
    n = g.number_of_nodes()
    diam = nx.diameter(g) if nx.is_connected(g) else 0
    if diam == 0:
        diam = 1
    B = np.zeros((diam + 1, n + 1), dtype=int)
    for src in g.nodes():
        lengths = nx.single_source_shortest_path_length(g, src)
        # how many nodes at each distance l from src
        from collections import Counter
        dist_counts = Counter(lengths.values())
        for l in range(diam + 1):
            k = dist_counts.get(l, 0)
            B[l, k] += 1
    return B


def d_portrait(g1: nx.Graph, g2: nx.Graph) -> float:
    """Jensen-Shannon divergence between network portraits."""
    B1 = _network_portrait(g1).astype(float)
    B2 = _network_portrait(g2).astype(float)
    # Pad to common shape
    rows = max(B1.shape[0], B2.shape[0])
    cols = max(B1.shape[1], B2.shape[1])
    P1 = np.zeros((rows, cols)); P1[:B1.shape[0], :B1.shape[1]] = B1
    P2 = np.zeros((rows, cols)); P2[:B2.shape[0], :B2.shape[1]] = B2
    P1 = P1.flatten() / max(P1.sum(), 1)
    P2 = P2.flatten() / max(P2.sum(), 1)
    M = 0.5 * (P1 + P2)
    def kl(p, q):
        mask = (p > 0) & (q > 0)
        return float(np.sum(p[mask] * np.log2(p[mask] / q[mask])))
    return 0.5 * kl(P1, M) + 0.5 * kl(P2, M)


# ---------------------------------------------------------------------------
# Topological role of a single node (used for Daath-vs-CC test)
# ---------------------------------------------------------------------------

def node_role(g: nx.Graph, node) -> Dict[str, float]:
    bc = nx.betweenness_centrality(g)
    cc = nx.closeness_centrality(g)
    g_minus = g.copy()
    g_minus.remove_node(node)
    components = nx.number_connected_components(g_minus)
    largest = max(len(c) for c in nx.connected_components(g_minus))
    return {
        "node": node,
        "degree": g.degree(node),
        "degree_z": (g.degree(node) - np.mean([d for _, d in g.degree()]))
                    / (np.std([d for _, d in g.degree()]) + 1e-9),
        "betweenness": bc[node],
        "betweenness_rank": 1 + sum(1 for v in bc.values() if v > bc[node]),
        "closeness": cc[node],
        "is_articulation": node in set(nx.articulation_points(g)),
        "components_after_removal": components,
        "fraction_isolated_after": 1 - largest / max(g.number_of_nodes() - 1, 1),
    }
