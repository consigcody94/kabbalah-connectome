"""
Comparative mythology test: does the joined Tree of Life uniquely match
human brain topology, or do OTHER world-tree / esoteric structures match
equally well?

Tests against five mythological / esoteric structures, each constructed as
a graph from its canonical sources:

  1. Tree of Life (Kabbalistic) — joined-trees with Daath bridge (baseline)
  2. I Ching — 64 hexagrams as nodes, edges = single-line transformations
              (canonical hypercube structure: Q_6, the 6-dim hypercube)
  3. Yggdrasil — Norse 9 worlds + tree connectivity (Snorri's Prose Edda)
  4. Chakra-Nadi — Hindu yoga 7 chakras + 3 nadis (Ida, Pingala, Sushumna)
  5. Sri Yantra — Hindu/Tantric 9 interlocking triangles graph
  6. Mayan World Tree (Wacah Chan) — 4 directions × 3 levels + center

For each, compute distance to the real Budapest connectome (coarsened to
22 nodes). If joined Tree of Life uniquely beats all others → the result
is specific to Kabbalah. If multiple match equally well → the result
generalizes to a CLASS of tree-like esoteric structures.
"""
import csv
import json
import os
import sys
from collections import defaultdict
from itertools import product

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import joined_trees
from metrics import all_metrics, d_invariant
from nulls import er_ensemble


# ===========================================================================
# Build the comparative structures
# ===========================================================================

def build_iching():
    """I Ching: 64 hexagrams, edges between hexagrams differing in one line.
    This is the 6-dimensional hypercube graph Q_6 — 64 nodes, 192 edges,
    perfectly regular degree 6. The structure of the entire I Ching as a
    state-transition graph (King Wen ordering ignored)."""
    g = nx.Graph()
    nodes = list(product([0, 1], repeat=6))
    for n in nodes:
        g.add_node(n)
    for n1 in nodes:
        for i in range(6):
            n2 = list(n1); n2[i] = 1 - n2[i]
            g.add_edge(n1, tuple(n2))
    return g


def build_yggdrasil():
    """Norse Yggdrasil: 9 worlds connected through the cosmic ash tree.

    Standard topology from Snorri's Prose Edda + Eddic poems:
      - Asgard (gods) at top
      - Vanaheim, Alfheim around Asgard
      - Midgard (humans) in middle, surrounded by Jotunheim
      - Muspelheim (fire), Niflheim (ice) at bottom levels
      - Svartalfheim (dark elves) at root
      - Helheim (dead) below all

    Connections from canonical Eddic geography:
      Asgard ↔ Vanaheim ↔ Alfheim
      Asgard ↔ Midgard (Bifrost rainbow bridge)
      Midgard ↔ Jotunheim ↔ Svartalfheim ↔ Niflheim
      Niflheim ↔ Helheim ↔ Muspelheim
      Yggdrasil itself connects all 9 (additional center node)
    """
    g = nx.Graph()
    worlds = ["Asgard", "Vanaheim", "Alfheim", "Midgard", "Jotunheim",
              "Svartalfheim", "Nidavellir", "Niflheim", "Muspelheim", "Helheim"]
    g.add_nodes_from(worlds)
    edges = [
        ("Asgard", "Vanaheim"), ("Asgard", "Alfheim"),
        ("Vanaheim", "Alfheim"),
        ("Asgard", "Midgard"),         # Bifrost
        ("Vanaheim", "Midgard"),
        ("Alfheim", "Midgard"),
        ("Midgard", "Jotunheim"),
        ("Jotunheim", "Svartalfheim"),
        ("Jotunheim", "Nidavellir"),
        ("Svartalfheim", "Nidavellir"),
        ("Svartalfheim", "Niflheim"),
        ("Nidavellir", "Niflheim"),
        ("Niflheim", "Helheim"),
        ("Helheim", "Muspelheim"),
        ("Niflheim", "Muspelheim"),    # primordial cosmic gap (Ginnungagap)
        ("Asgard", "Helheim"),         # Odin's journeys
    ]
    g.add_edges_from(edges)
    return g


def build_chakra_nadi():
    """Hindu Chakra-Nadi system: 7 chakras + 3 main nadis (Ida, Pingala,
    Sushumna), modeled as a graph following Hatha Yoga Pradipika
    + Shiva Samhita.

    Nodes: 7 chakras + 3 nadi-segment terminals = 10 effective nodes.
    Sushumna runs through all chakras vertically.
    Ida (left, lunar) and Pingala (right, solar) cross at each chakra.

    Following the standard caduceus model:
    - Each chakra is a node
    - Sushumna connects all chakras in sequence (vertical axis)
    - Ida and Pingala cross at each chakra (so each chakra has L/R bridge)
    - Nadis terminate at top (Sahasrara) and bottom (Muladhara)
    """
    g = nx.Graph()
    chakras = ["Muladhara", "Svadhisthana", "Manipura", "Anahata",
               "Vishuddha", "Ajna", "Sahasrara"]
    g.add_nodes_from(chakras)
    # Sushumna: vertical chain
    for c1, c2 in zip(chakras[:-1], chakras[1:]):
        g.add_edge(c1, c2)
    # Add Ida and Pingala as separate "side" nodes per chakra
    # Following the caduceus interpretation: Ida_n and Pingala_n cross at chakra_n
    # We'll model this as: each chakra has an Ida-side and Pingala-side anchor
    # that connects up/down to the next chakra's Pingala/Ida (the crossing).
    for c in chakras:
        g.add_node(f"Ida-{c}")
        g.add_node(f"Pingala-{c}")
        g.add_edge(c, f"Ida-{c}")
        g.add_edge(c, f"Pingala-{c}")
    # Cross-connections: Ida of chakra n → Pingala of chakra n+1 (and vice versa)
    for c1, c2 in zip(chakras[:-1], chakras[1:]):
        g.add_edge(f"Ida-{c1}", f"Pingala-{c2}")
        g.add_edge(f"Pingala-{c1}", f"Ida-{c2}")
    return g


def build_sri_yantra():
    """Sri Yantra: 9 interlocking triangles (4 upward Shiva, 5 downward
    Shakti) creating 43 small triangles. The 9 base triangles share
    intersections — we model this as the intersection graph.

    Each of 9 triangles is a node; edges connect triangles that share a
    region (i.e., overlap geometrically).
    """
    g = nx.Graph()
    # 9 triangles labeled by orientation and rank (concentrically nested)
    upward = ["U1", "U2", "U3", "U4"]      # Shiva (4)
    downward = ["D1", "D2", "D3", "D4", "D5"]  # Shakti (5)
    g.add_nodes_from(upward + downward)
    # All upward and all downward triangles intersect each other (concentric)
    # Plus each upward intersects all downward (interlocking)
    for u in upward:
        for d in downward:
            g.add_edge(u, d)
    # Plus consecutive upward-upward and downward-downward share central points
    for a, b in zip(upward[:-1], upward[1:]):
        g.add_edge(a, b)
    for a, b in zip(downward[:-1], downward[1:]):
        g.add_edge(a, b)
    return g


def build_wacah_chan():
    """Mayan World Tree (Wacah Chan / Yaxche): 4 cardinal directions +
    center, on 3 levels (Upperworld / Middleworld / Underworld).

    13 nodes: 4 corners × 3 levels = 12, plus 1 axis (the world tree itself).
    Each corner connects to its level-mates and to the corresponding
    corner above/below; the axis connects all 3 center points.
    """
    g = nx.Graph()
    levels = ["Upper", "Middle", "Lower"]
    dirs = ["E", "N", "W", "S"]
    # Corner nodes
    for lev in levels:
        for d in dirs:
            g.add_node(f"{lev}-{d}")
    g.add_node("Axis")  # the World Tree itself
    # Within-level: corners connect in cardinal cycle
    for lev in levels:
        for d1, d2 in zip(dirs, dirs[1:] + dirs[:1]):
            g.add_edge(f"{lev}-{d1}", f"{lev}-{d2}")
    # Cross-level: same direction connects across levels
    for d in dirs:
        g.add_edge(f"Upper-{d}", f"Middle-{d}")
        g.add_edge(f"Middle-{d}", f"Lower-{d}")
    # Axis connects to ONE corner per level (the W direction by tradition)
    g.add_edge("Axis", "Upper-W")
    g.add_edge("Axis", "Middle-W")
    g.add_edge("Axis", "Lower-W")
    return g


# ===========================================================================
# Load real coarsened brain
# ===========================================================================
print("Loading + coarsening Budapest connectome...")
node_info = {}
with open("data/budapest/nodes.csv", encoding="utf-8") as f:
    reader = csv.DictReader(
        (l for l in f if not l.startswith("#")),
        fieldnames=["index", "dn_fsname", "dn_hemisphere",
                    "dn_name", "dn_region", "_pos"])
    for r in reader:
        node_info[int(r["index"])] = {
            "fsname": r["dn_fsname"],
            "hemi": r["dn_hemisphere"],
        }

G = nx.Graph()
G.add_nodes_from(node_info.keys())
with open("data/budapest/edges.csv", encoding="utf-8") as f:
    reader = csv.DictReader(
        (l for l in f if not l.startswith("#")),
        fieldnames=["source", "target", "fcm", "flm", "fam",
                    "ecm", "fcmd", "flmd", "famd", "occ"])
    for r in reader:
        try: occ = int(r["occ"])
        except (ValueError, TypeError): continue
        if occ < 100: continue
        s, t = int(r["source"]), int(r["target"])
        if s != t: G.add_edge(s, t)
G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])

def lobe_of(fsname):
    s = fsname.lower()
    if "frontal" in s or "parsop" in s or "parstri" in s or "parsorb" in s \
       or "rostral" in s: return "Frontal"
    if "precentral" in s or "paracentral" in s: return "Motor"
    if "postcentral" in s: return "Somatosensory"
    if "parietal" in s or "supramarg" in s or "precune" in s: return "Parietal"
    if "transversetemporal" in s: return "Auditory"
    if "temporal" in s or "fusiform" in s or "banks" in s \
       or "entorhinal" in s or "parahippo" in s: return "Temporal"
    if "occipital" in s or "lingual" in s or "cuneus" in s \
       or "pericalc" in s: return "Occipital"
    if "hippo" in s: return "Hippocampus"
    if "thalam" in s: return "Thalamus"
    if "cingul" in s or "isthmus" in s: return "Cingulate"
    if "insula" in s: return "Insula"
    return "Subcortical"

coarse_map = {n: f"{node_info[n]['hemi'] or 'none'}-{lobe_of(node_info[n]['fsname'])}"
              for n in G.nodes()}
G_coarse = nx.Graph()
weights = defaultdict(int)
for u, v in G.edges():
    if coarse_map[u] != coarse_map[v]:
        weights[tuple(sorted([coarse_map[u], coarse_map[v]]))] += 1
for (a, b), w in weights.items():
    if w >= 5: G_coarse.add_edge(a, b)

m_real = all_metrics(G_coarse, "REAL")
print(f"  Real brain coarse: N={G_coarse.number_of_nodes()} E={G_coarse.number_of_edges()}")


# ===========================================================================
# Build all candidate structures and compare
# ===========================================================================
print("\nBuilding mythological structures...")
candidates = {
    "TreeOfLife (joined+Daath)": joined_trees(shared_daath=True),
    "IChing (6-cube)":           build_iching(),
    "Yggdrasil (9 worlds)":      build_yggdrasil(),
    "Chakra-Nadi (Hindu)":       build_chakra_nadi(),
    "SriYantra (9 triangles)":   build_sri_yantra(),
    "Wacah Chan (Mayan)":        build_wacah_chan(),
}
for name, g in candidates.items():
    print(f"  {name:<32s} N={g.number_of_nodes():3d}  E={g.number_of_edges():3d}  "
          f"connected={nx.is_connected(g)}")

print("\n" + "=" * 80)
print("DISTANCE TO REAL HUMAN BRAIN (lower = better match)")
print("=" * 80)

results = {}
print(f"\n{'Structure':<32s} {'N':>4} {'E':>4} {'d_real':>8} {'p_ER':>8} {'p_GEO':>8}")
print("-" * 78)
for name, g in candidates.items():
    if not nx.is_connected(g):
        print(f"  {name:<32s} disconnected, skip")
        continue
    m = all_metrics(g, name)
    d_real = d_invariant(m, m_real)
    # Null comparison: ER + GEO at same N, M as the candidate
    er_nulls = er_ensemble(g.number_of_nodes(), g.number_of_edges(), k=100, seed=0)
    er_dists = [d_invariant(all_metrics(ng, "n"), m_real) for ng in er_nulls]
    p_er = float((np.array(er_dists) <= d_real).mean()) if er_dists else float("nan")
    from nulls import geo_ensemble
    geo_nulls = geo_ensemble(g.number_of_nodes(), g.number_of_edges(), k=100, seed=0)
    geo_dists = [d_invariant(all_metrics(ng, "n"), m_real) for ng in geo_nulls]
    p_geo = float((np.array(geo_dists) <= d_real).mean()) if geo_dists else float("nan")
    results[name] = {
        "N": g.number_of_nodes(), "E": g.number_of_edges(),
        "d_real": d_real, "p_ER": p_er, "p_GEO": p_geo,
        "metrics": {k: m.get(k) for k in m if k != "label"},
    }
    print(f"  {name:<32s} {g.number_of_nodes():>4} {g.number_of_edges():>4} "
          f"{d_real:>8.4f} {p_er:>8.4f} {p_geo:>8.4f}")

# ===========================================================================
# Verdict
# ===========================================================================
print("\n" + "=" * 80)
print("VERDICT")
print("=" * 80)
sorted_by_d = sorted(results.items(), key=lambda x: x[1]["d_real"])
print(f"\nRanked by distance to real brain (best to worst):")
for i, (name, r) in enumerate(sorted_by_d, 1):
    print(f"  {i}. {name:<32s} d={r['d_real']:.4f}  "
          f"p_ER={r['p_ER']:.4f}  p_GEO={r['p_GEO']:.4f}")

tol_d = results.get("TreeOfLife (joined+Daath)", {}).get("d_real")
print(f"\nTree of Life baseline:    d = {tol_d:.4f}")
print(f"\nNumber of mythological structures CLOSER to brain than Tree of Life: ", end="")
closer = sum(1 for n, r in results.items()
             if n != "TreeOfLife (joined+Daath)" and r["d_real"] < tol_d)
total = sum(1 for n in results if n != "TreeOfLife (joined+Daath)")
print(f"{closer}/{total}")

if closer == 0:
    print("\n=> Tree of Life UNIQUELY closest to brain among tested mythologies.")
elif closer < total / 2:
    print(f"\n=> Tree of Life among the closest, but {closer} structures match better.")
else:
    print(f"\n=> Tree of Life NOT uniquely brain-like; {closer}/{total} mythologies "
          f"match better. Result generalizes to tree-like esoteric structures.")

# Save
with open("data/comparative_mythology_results.json", "w") as f:
    json.dump(results, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/comparative_mythology_results.json")
