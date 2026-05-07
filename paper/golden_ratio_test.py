"""
Golden ratio test: does φ ≈ 1.618 appear in graph-theoretic ratios for
the joined-trees graph and real human brain?

Tests:
  1. Spectral ratio λ_i/λ_{i+1} for the largest eigenvalues
  2. Degree ratios: max/median, max/2nd_max, etc.
  3. Modularity / clustering ratios
  4. Path-length ratios

For each, we ask: how close is the ratio to φ ≈ 1.618? Is the deviation
smaller for joined-trees and brain than for random graphs?

φ = (1 + √5) / 2 ≈ 1.6180339887
1/φ = φ - 1 ≈ 0.6180339887
φ² = φ + 1 ≈ 2.6180339887
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import joined_trees
from nulls import er_ensemble


PHI = (1 + np.sqrt(5)) / 2
print(f"φ = {PHI:.10f}")


def golden_distance(value, phi=PHI):
    """How close is `value` to a φ-related quantity (φ, 1/φ, φ², 2-φ, etc.)?"""
    candidates = [phi, 1/phi, phi**2, phi-1, 2-phi, phi/2, 2*phi]
    return min(abs(value - c) for c in candidates)


def graph_ratios(G):
    """Compute several graph ratios that could exhibit φ."""
    A = nx.to_numpy_array(G)
    eigs = sorted(np.linalg.eigvals(A).real, reverse=True)
    deg = sorted([d for _, d in G.degree()], reverse=True)
    bc = nx.betweenness_centrality(G)
    bcs = sorted(bc.values(), reverse=True)

    ratios = {}
    # Spectral
    if len(eigs) > 1 and abs(eigs[1]) > 1e-9:
        ratios["lambda_1/lambda_2"] = float(eigs[0] / eigs[1])
    if len(eigs) > 2 and abs(eigs[2]) > 1e-9:
        ratios["lambda_2/lambda_3"] = float(eigs[1] / eigs[2])
    # Degree
    if deg[0] > 0 and deg[-1] > 0:
        ratios["max_deg/min_deg"] = float(deg[0] / deg[-1])
    if len(deg) > 1 and deg[1] > 0:
        ratios["max_deg/2nd_max_deg"] = float(deg[0] / deg[1])
    if np.median(deg) > 0:
        ratios["max_deg/median_deg"] = float(deg[0] / np.median(deg))
    # Betweenness
    if len(bcs) > 1 and bcs[1] > 0:
        ratios["max_bc/2nd_max_bc"] = float(bcs[0] / bcs[1])
    # Clustering / transitivity
    C = nx.average_clustering(G)
    T = nx.transitivity(G)
    if T > 0:
        ratios["C/T"] = float(C / T)
    # Path / diameter
    if nx.is_connected(G):
        L = nx.average_shortest_path_length(G)
        D = nx.diameter(G)
        if L > 0:
            ratios["D/L"] = float(D / L)
    # Density / mean_degree (not a ratio of length, just sanity)
    return ratios


# ---------------------------------------------------------------------------
print("\n=== Joined-trees graph ===")
J = joined_trees(shared_daath=True)
ratios_J = graph_ratios(J)
for k, v in ratios_J.items():
    print(f"  {k:<22s} = {v:.4f}  golden_dist = {golden_distance(v):.4f}")

# ---------------------------------------------------------------------------
print("\n=== Loading + coarsening Budapest connectome ===")
node_info = {}
with open("data/budapest/nodes.csv", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split(",")
        if len(parts) < 4: continue
        try: idx = int(parts[0])
        except (ValueError, TypeError): continue
        node_info[idx] = {"fsname": parts[1], "hemi": parts[2]}
G_full = nx.Graph()
G_full.add_nodes_from(node_info.keys())
with open("data/budapest/edges.csv", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split(",")
        try:
            s, t = int(parts[0]), int(parts[1]); occ = int(parts[9])
        except (ValueError, TypeError, IndexError):
            continue
        if occ < 100 or s == t: continue
        G_full.add_edge(s, t)
G_full.remove_nodes_from([n for n in list(G_full.nodes()) if G_full.degree(n) == 0])

def lobe_of(name):
    s = name.lower()
    if "frontal" in s or "parsop" in s or "parstri" in s or "parsorb" in s or "rostral" in s: return "Frontal"
    if "precentral" in s or "paracentral" in s: return "Motor"
    if "postcentral" in s: return "Somatosensory"
    if "parietal" in s or "supramarg" in s or "precune" in s: return "Parietal"
    if "transversetemporal" in s: return "Auditory"
    if "temporal" in s or "fusiform" in s or "banks" in s or "entorhinal" in s or "parahippo" in s: return "Temporal"
    if "occipital" in s or "lingual" in s or "cuneus" in s or "pericalc" in s: return "Occipital"
    if "hippo" in s: return "Hippocampus"
    if "thalam" in s: return "Thalamus"
    if "cingul" in s or "isthmus" in s: return "Cingulate"
    if "insula" in s: return "Insula"
    return "Subcortical"
coarse_map = {n: f"{node_info[n]['hemi'] or 'none'}-{lobe_of(node_info[n]['fsname'])}"
              for n in G_full.nodes()}
G_brain = nx.Graph()
weights = defaultdict(int)
for u, v in G_full.edges():
    if coarse_map[u] != coarse_map[v]:
        weights[tuple(sorted([coarse_map[u], coarse_map[v]]))] += 1
for (a, b), w in weights.items():
    if w >= 5: G_brain.add_edge(a, b)

print("\n=== Real brain (Budapest coarse) ===")
ratios_B = graph_ratios(G_brain)
for k, v in ratios_B.items():
    print(f"  {k:<22s} = {v:.4f}  golden_dist = {golden_distance(v):.4f}")

# ---------------------------------------------------------------------------
print("\n=== Random graphs (ER, n=200, same N/M as joined-trees) ===")
nulls_J = er_ensemble(J.number_of_nodes(), J.number_of_edges(), k=200)
null_ratios = []
for ng in nulls_J:
    if not nx.is_connected(ng): continue
    try:
        null_ratios.append(graph_ratios(ng))
    except Exception:
        continue

# Aggregate null golden distances
print("\n--- Null distribution: golden_distance for each ratio ---")
for k in ratios_J:
    null_dists = [golden_distance(r[k]) for r in null_ratios if k in r]
    if not null_dists: continue
    null_arr = np.array(null_dists)
    j_dist = golden_distance(ratios_J[k])
    p = float((null_arr <= j_dist).mean())
    print(f"  {k:<22s}: joined φ-dist={j_dist:.4f}  "
          f"null mean={null_arr.mean():.4f}  p(null at least as φ-close)={p:.4f}")

# ---------------------------------------------------------------------------
print("\n=== Critical comparison: are joined-trees ratios closer to φ than brain's? ===")
print(f"\n{'Metric':<22s} {'joined':>10s} {'brain':>10s} "
      f"{'φ-dist(J)':>10s} {'φ-dist(B)':>10s} {'closer':>8s}")
for k in ratios_J:
    if k not in ratios_B: continue
    j = ratios_J[k]; b = ratios_B[k]
    jd = golden_distance(j); bd = golden_distance(b)
    closer = "BOTH" if jd < 0.1 and bd < 0.1 else ("J" if jd < bd else "B")
    print(f"  {k:<22s} {j:>10.4f} {b:>10.4f} {jd:>10.4f} {bd:>10.4f} {closer:>8s}")

# ---------------------------------------------------------------------------
print("\n=== Spectral signature comparison ===")
# Show top eigenvalues of each
A_J = nx.to_numpy_array(J)
A_B = nx.to_numpy_array(G_brain)
eigs_J = sorted(np.linalg.eigvals(A_J).real, reverse=True)
eigs_B = sorted(np.linalg.eigvals(A_B).real, reverse=True)
print(f"\nTop 5 eigenvalues:")
print(f"  Joined trees: {[f'{e:.3f}' for e in eigs_J[:5]]}")
print(f"  Real brain:   {[f'{e:.3f}' for e in eigs_B[:5]]}")
print(f"\nConsecutive ratios (looking for ~1.618):")
print(f"  Joined: {[f'{eigs_J[i]/eigs_J[i+1]:.3f}' for i in range(min(4, len(eigs_J)-1)) if abs(eigs_J[i+1])>1e-9]}")
print(f"  Brain:  {[f'{eigs_B[i]/eigs_B[i+1]:.3f}' for i in range(min(4, len(eigs_B)-1)) if abs(eigs_B[i+1])>1e-9]}")

with open("data/golden_ratio_results.json", "w") as f:
    json.dump({
        "phi": PHI,
        "joined_ratios": {k: ratios_J[k] for k in ratios_J},
        "brain_ratios": ratios_B,
        "joined_phi_distances": {k: golden_distance(v) for k, v in ratios_J.items()},
        "brain_phi_distances": {k: golden_distance(v) for k, v in ratios_B.items()},
        "joined_top_eigs": eigs_J[:5],
        "brain_top_eigs": eigs_B[:5],
    }, f, indent=2,
       default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/golden_ratio_results.json")
