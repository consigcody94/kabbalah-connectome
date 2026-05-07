"""
Build the Flower of Life as an actual graph and test it against the
real brain. Previously (Appendix A) we tested only its symmetry — D₆ vs
the brain's D₁ — and rejected it. But the FoL pattern produces a specific
PLANAR GRAPH (19 circles → vertex/intersection structure) that has
never been tested as a topology.

Two graph constructions:

  (a) Center-graph: each of 19 circle centers is a node. Edges connect
      centers whose circles overlap (distance ≤ 2r). Standard FoL geometry.

  (b) Tree-derived-from-FoL: Drunvalo Melchizedek's claim is that the
      Tree of Life is formed by selecting 10 specific circles from the
      Flower of Life. We construct that subgraph and test it.

Test: does either FoL graph match real brain topology better than nulls?
Does it match better than the joined-trees graph?
"""
import csv
import json
import os
import sys
from collections import defaultdict
from itertools import combinations

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import joined_trees
from metrics import all_metrics, d_invariant
from nulls import all_nulls


# ---------------------------------------------------------------------------
# Build Flower of Life as a graph
# ---------------------------------------------------------------------------
def flower_of_life_centers(r=1.0):
    """19 circle centers on hexagonal lattice for the standard Flower of Life.

    Layout:
      - 1 center at origin
      - 6 in first ring at distance √3 (hex packing) — but FoL uses r-spacing
        where adjacent circles share boundaries (distance r, overlap)
      - 12 in second ring
    """
    centers = [(0.0, 0.0)]
    # First ring: 6 hexagonal neighbors at distance r (overlapping)
    for i in range(6):
        a = i * np.pi / 3
        centers.append((r * np.cos(a), r * np.sin(a)))
    # Second ring: 12 outer
    # 6 corners at distance √3 * r (in hexagon-corner direction)
    for i in range(6):
        a = i * np.pi / 3 + np.pi / 6
        centers.append((np.sqrt(3) * r * np.cos(a), np.sqrt(3) * r * np.sin(a)))
    # 6 edge-midpoints at distance 2r
    for i in range(6):
        a = i * np.pi / 3
        centers.append((2 * r * np.cos(a), 2 * r * np.sin(a)))
    return centers


def fol_center_graph():
    """FoL center graph: nodes = circle centers, edges = pairs whose circles
    overlap (distance ≤ 2r — circles of radius r overlap when centers are
    closer than 2r)."""
    centers = flower_of_life_centers(r=1.0)
    g = nx.Graph()
    for i, c in enumerate(centers):
        g.add_node(i, pos=c)
    for i, j in combinations(range(len(centers)), 2):
        d = np.hypot(centers[i][0] - centers[j][0],
                     centers[i][1] - centers[j][1])
        if d <= 2.0 + 1e-6:    # overlap or touch
            g.add_edge(i, j)
    return g


def fol_tree_derived():
    """Drunvalo Melchizedek's Tree of Life from Flower of Life:
    10 specific circles selected from the 19. We approximate by picking
    the 10 circles whose positions match the standard Tree of Life
    arrangement (3 columns, vertical orientation).

    Layout (matches Tree of Life's three pillars):
      Center column: top, mid-up, mid, mid-down, bottom
      Right column: 2 nodes (Chokhmah, Chesed level)
      Left column: 2 nodes (Binah, Geburah level)
      Plus Netzach, Hod
    """
    centers_full = flower_of_life_centers(r=1.0)
    # Indices 0..18 in our construction:
    #   0: center
    #   1-6: first ring (i*60°: 0=E, 1=NE, 2=NW, 3=W, 4=SW, 5=SE)
    #   7-12: outer-corner ring (i*60° + 30°)
    #   13-18: outer-edge ring (i*60°, distance 2r)

    # Pick 10 nodes corresponding to the standard TOL three-pillar layout:
    # We map by position. Inspect centers_full and select by (x, y).
    selected = [
        15,   # Keter (top of outer-edge ring, north)
        7,    # Chokhmah (NE outer-corner)
        12,   # Binah (NW outer-corner)
        2,    # Chesed (NE first ring)
        3,    # Geburah (NW first ring)
        0,    # Tiferet (center)
        5,    # Netzach (SE first ring)
        4,    # Hod (SW first ring)
        9,    # Yesod (S outer-corner)
        16,   # Malkuth (bottom of outer-edge ring)
    ]
    # Build subgraph from the full FoL center graph
    g_full = fol_center_graph()
    g = g_full.subgraph(selected).copy()
    return g


# ---------------------------------------------------------------------------
# Load real brain (same as previous tests)
# ---------------------------------------------------------------------------
print("Loading + coarsening Budapest connectome...")
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
m_brain = all_metrics(G_brain, "REAL")
print(f"  Brain coarse: N={G_brain.number_of_nodes()} E={G_brain.number_of_edges()}")


# ---------------------------------------------------------------------------
# Build candidates and test
# ---------------------------------------------------------------------------
print("\nBuilding candidate structures...")
J = joined_trees(shared_daath=True)
FoL = fol_center_graph()
TOL_from_FoL = fol_tree_derived()

candidates = {
    "Joined trees (Kabbalah baseline)": J,
    "Flower of Life (center graph)":     FoL,
    "Tree of Life derived from FoL":     TOL_from_FoL,
}

for name, g in candidates.items():
    print(f"  {name:<36s} N={g.number_of_nodes():3d}  E={g.number_of_edges():3d}  "
          f"connected={nx.is_connected(g)}")

print("\n" + "=" * 80)
print("DISTANCE TO REAL BRAIN")
print("=" * 80)

results = {}
print(f"\n{'Structure':<36s} {'N':>4} {'E':>4} {'d_real':>8} {'p_ER':>8}")
print("-" * 70)
for name, g in candidates.items():
    if not nx.is_connected(g):
        print(f"  {name:<36s} disconnected, skip")
        continue
    m = all_metrics(g, name)
    d_real = d_invariant(m, m_brain)
    nulls = all_nulls(g, k_per=100)
    er_dists = [d_invariant(all_metrics(ng, "n"), m_brain)
                for ng in nulls.get("ER", [])]
    p_er = float((np.array(er_dists) <= d_real).mean()) if er_dists else float("nan")
    results[name] = {
        "N": g.number_of_nodes(), "E": g.number_of_edges(),
        "d_real": d_real, "p_ER": p_er,
    }
    print(f"  {name:<36s} {g.number_of_nodes():>4} {g.number_of_edges():>4} "
          f"{d_real:>8.4f} {p_er:>8.4f}")

# ---------------------------------------------------------------------------
# Direct comparison: Tree-derived-from-FoL vs Joined-trees
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)
J_d = results.get("Joined trees (Kabbalah baseline)", {}).get("d_real")
FoL_d = results.get("Flower of Life (center graph)", {}).get("d_real")
TFoL_d = results.get("Tree of Life derived from FoL", {}).get("d_real")

print(f"\n  Joined trees (Kabbalah):         d = {J_d:.4f}   "
      f"p_ER = {results['Joined trees (Kabbalah baseline)']['p_ER']:.4f}")
print(f"  Flower of Life (center graph):   d = {FoL_d:.4f}   "
      f"p_ER = {results['Flower of Life (center graph)']['p_ER']:.4f}")
print(f"  Tree of Life from FoL:           d = {TFoL_d:.4f}   "
      f"p_ER = {results['Tree of Life derived from FoL']['p_ER']:.4f}")

print(f"\n  Best Kabbalistic structure: " +
      ("Joined trees" if J_d <= TFoL_d else "Tree of Life derived from FoL"))
print(f"\n  Does Flower of Life graph match brain? " +
      ("YES" if results['Flower of Life (center graph)']['p_ER'] < 0.05
       else "NO — Flower of Life graph topology does NOT match brain."))

# ---------------------------------------------------------------------------
# Show degree distribution of FoL — is it consistent with hex packing?
# ---------------------------------------------------------------------------
print("\n--- FoL center graph properties ---")
degs = [d for _, d in FoL.degree()]
print(f"  Degree distribution: {sorted(degs, reverse=True)}")
print(f"  This is the hexagonal lattice graph — central node degree {FoL.degree(0)}")

with open("data/flower_of_life_results.json", "w") as f:
    json.dump(results, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/flower_of_life_results.json")
