"""
Real connectome test using the Budapest Reference Connectome v2.0
(Szalkai et al. 2015, Neuroscience Letters; consensus from 477 HCP subjects).

Tests H1 and H2 against actual human DTI tractography data.

  H1 (aggregate structure): does the joined-trees graph match the topology
                            of the real human connectome?
  H2 (bridge node):         in the real connectome, is there a single node
                            (or small set) with the topological signature
                            we predicted for Daath / corpus callosum?

The Budapest connectome has 1015 nodes (FreeSurfer fine parcellation) with
left/right hemisphere labels. We coarsen to lobe-level (~10 regions per
hemisphere) for fair comparison with the Trees and run the same metrics.
"""

import csv
import json
import os
import sys
from collections import defaultdict, Counter

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import (
    SEPHIROT_10, QLI_NAMES, joined_trees, tree_of_life,
    brain_model_10,
)
from metrics import all_metrics, d_invariant, INVARIANT_KEYS, node_role
from nulls import all_nulls


# ---------------------------------------------------------------------------
# 1. Load the Budapest connectome
# ---------------------------------------------------------------------------
print("=" * 78)
print("REAL CONNECTOME TEST: Budapest Reference Connectome v2.0")
print("Source: Szalkai et al. 2015, NS Lett. (consensus over 477 HCP subjects)")
print("=" * 78)

NODES_CSV = "data/budapest/nodes.csv"
EDGES_CSV = "data/budapest/edges.csv"

# Load nodes with hemisphere + region labels
node_info = {}     # idx -> dict
with open(NODES_CSV, encoding="utf-8") as f:
    reader = csv.DictReader(
        (line for line in f if not line.startswith("#")),
        fieldnames=["index", "dn_fsname", "dn_hemisphere",
                    "dn_name", "dn_region", "_pos"],
    )
    for row in reader:
        idx = int(row["index"])
        node_info[idx] = {
            "fsname": row["dn_fsname"],
            "hemi": row["dn_hemisphere"],   # "left", "right", or "" for subcortical
            "name": row["dn_name"],
            "region": row["dn_region"],     # cortical / subcortical
        }

# Load edges (use occurrences > 100 as confidence threshold to thin slightly)
G_full = nx.Graph()
G_full.add_nodes_from(node_info.keys())
inter_hemi_edges = 0
intra_hemi_edges = 0
with open(EDGES_CSV, encoding="utf-8") as f:
    reader = csv.DictReader(
        (line for line in f if not line.startswith("#")),
        fieldnames=["source", "target", "fiber_count_mean",
                    "fiber_length_mean", "fractional_anisotropy_mean",
                    "electrical_connectivity_median", "fiber_count_median",
                    "fiber_length_median", "fractional_anisotropy_median",
                    "occurences"],
    )
    for row in reader:
        try:
            occ = int(row["occurences"])
        except (ValueError, TypeError):
            continue
        if occ < 100:           # require edge present in >= 100 of 477 subjects
            continue
        s, t = int(row["source"]), int(row["target"])
        if s == t:
            continue
        G_full.add_edge(s, t, weight=float(row["fiber_count_mean"]),
                        occ=occ)
        h1 = node_info[s]["hemi"]
        h2 = node_info[t]["hemi"]
        if h1 and h2:
            if h1 != h2:
                inter_hemi_edges += 1
            else:
                intra_hemi_edges += 1

# Drop isolated nodes (those with no edges meeting threshold)
G_full.remove_nodes_from([n for n in list(G_full.nodes())
                          if G_full.degree(n) == 0])

print(f"\nFull connectome (occurrence >= 100 / 477):")
print(f"  Nodes:           {G_full.number_of_nodes()}")
print(f"  Edges:           {G_full.number_of_edges()}")
print(f"  Connected:       {nx.is_connected(G_full)}")
print(f"  Intra-hemi:      {intra_hemi_edges}")
print(f"  Inter-hemi:      {inter_hemi_edges}")
print(f"  Inter/Total:     {inter_hemi_edges/(intra_hemi_edges+inter_hemi_edges):.3f}")
print(f"  Density:         {nx.density(G_full):.4f}")
print(f"  Mean degree:     {2*G_full.number_of_edges()/G_full.number_of_nodes():.2f}")

# ---------------------------------------------------------------------------
# 2. Find the topological "corpus callosum" — the strongest bridge nodes
# ---------------------------------------------------------------------------
# At fine parcellation, no single node IS the CC. Instead, the CC is
# represented as the set of inter-hemispheric edges. We look at which
# individual nodes have the highest betweenness that bridge the two
# hemispheres — the empirical analog of "Daath".

print("\n--- Bridge-node analysis (which real nodes act like Daath?) ---")
bc = nx.betweenness_centrality(G_full)
nodes_by_bc = sorted(bc.items(), key=lambda x: -x[1])
print(f"  Top 10 betweenness nodes (real connectome):")
for rank, (n, v) in enumerate(nodes_by_bc[:10], 1):
    info = node_info[n]
    deg = G_full.degree(n)
    # Count inter-hemi edges of this node
    inter = sum(1 for nb in G_full.neighbors(n)
                if node_info[nb]["hemi"] and info["hemi"]
                and node_info[nb]["hemi"] != info["hemi"])
    intra = deg - inter
    print(f"   #{rank:2d}  bc={v:.4f}  deg={deg:3d}  "
          f"inter-hemi={inter:3d}/{deg:3d}  "
          f"hemi={info['hemi']:5s}  {info['fsname']}")

# Test: does removing the top-betweenness node (or top-K nodes) significantly
# reduce inter-hemispheric connectivity? This is the empirical version of the
# articulation-point test.

def hemi_connectivity(g):
    """Number of nodes pairs (one in each hemi) connected by *some* path."""
    comps = list(nx.connected_components(g))
    out = 0
    for c in comps:
        L = sum(1 for n in c if node_info.get(n, {}).get("hemi") == "left")
        R = sum(1 for n in c if node_info.get(n, {}).get("hemi") == "right")
        out += L * R
    return out

baseline = hemi_connectivity(G_full)
print(f"\n  Connected L-R pairs in full graph: {baseline:,}")
for k in (1, 5, 10, 20, 50):
    g = G_full.copy()
    drop = [n for n, _ in nodes_by_bc[:k]]
    g.remove_nodes_from(drop)
    rem = hemi_connectivity(g)
    print(f"  Remove top {k:2d} betweenness nodes -> "
          f"L-R pairs={rem:,} ({100*rem/baseline:5.1f}% of baseline)")

# ---------------------------------------------------------------------------
# 3. Coarsen to lobe-level (~21 nodes) for direct comparison with Trees
# ---------------------------------------------------------------------------
print("\n--- Coarse-graining to 21 nodes (10 lobes/hemi + corpus callosum) ---")

def lobe_of(fsname):
    s = fsname.lower()
    if any(k in s for k in ("lateralorbito", "medialorbito", "rostral",
                            "caudalanterior", "frontalpole", "parstri",
                            "parsoper", "parsorb", "frontal")):
        return "Frontal"
    if any(k in s for k in ("precentral", "paracentral")):
        return "Motor"
    if any(k in s for k in ("postcentral",)):
        return "Somatosensory"
    if any(k in s for k in ("superiorparietal", "inferiorparietal",
                            "supramarginal", "precuneus")):
        return "Parietal"
    if any(k in s for k in ("transversetemporal",)):
        return "Auditory"
    if any(k in s for k in ("temporal", "fusiform", "banks", "entorhinal",
                            "parahippo")):
        return "Temporal"
    if any(k in s for k in ("occipital", "lingual", "cuneus", "pericalcarine")):
        return "Occipital"
    if "hippo" in s:
        return "Hippocampus"
    if "thalam" in s:
        return "Thalamus"
    if "cingul" in s or "isthmus" in s:
        return "Cingulate"
    if "insula" in s:
        return "Insula"
    return "Other"

# Build a coarsened graph: nodes = (hemisphere, lobe), or "CorpusCallosum"
# for inter-hemi. But actually the cleanest analogy is: nodes = (hemi, lobe)
# only, with inter-hemi edges represented as direct connections (no bridge
# node — because in the real brain there isn't a single bridge node either).
node_to_coarse = {}
for n, info in node_info.items():
    if n not in G_full:
        continue
    if not info["hemi"]:
        coarse = ("none", lobe_of(info["fsname"]))
    else:
        coarse = (info["hemi"], lobe_of(info["fsname"]))
    node_to_coarse[n] = coarse

G_coarse = nx.Graph()
edge_weights = defaultdict(int)
for u, v, d in G_full.edges(data=True):
    cu, cv = node_to_coarse.get(u), node_to_coarse.get(v)
    if cu is None or cv is None: continue
    if cu == cv: continue
    edge_weights[tuple(sorted([cu, cv]))] += 1
for (a, b), w in edge_weights.items():
    if w >= 5:    # require some support for coarse edge
        G_coarse.add_edge(a, b, weight=w)

# Make node names string for compatibility
G_coarse_str = nx.Graph()
for u, v, d in G_coarse.edges(data=True):
    su = f"{u[0]}-{u[1]}" if u[0] != "none" else u[1]
    sv = f"{v[0]}-{v[1]}" if v[0] != "none" else v[1]
    G_coarse_str.add_edge(su, sv, weight=d["weight"])

print(f"  Coarsened graph: N={G_coarse_str.number_of_nodes()} "
      f"E={G_coarse_str.number_of_edges()} "
      f"connected={nx.is_connected(G_coarse_str)}")

# List nodes and their inter-hemi connectivity
inter_hemi_edges_coarse = 0
for u, v in G_coarse_str.edges():
    if u.startswith("left") != v.startswith("left"):
        if u.startswith("left") or v.startswith("left"):
            if u.startswith("right") or v.startswith("right"):
                inter_hemi_edges_coarse += 1
print(f"  Inter-hemi edges (coarse): {inter_hemi_edges_coarse}")

# ---------------------------------------------------------------------------
# 4. Bridge-node test on coarse REAL connectome
# ---------------------------------------------------------------------------
print("\n--- Coarse-scale bridge-node test ---")
bc_c = nx.betweenness_centrality(G_coarse_str)
top = sorted(bc_c.items(), key=lambda x: -x[1])[:7]
print(f"  Top betweenness nodes in REAL coarsened connectome:")
for rank, (n, v) in enumerate(top, 1):
    deg = G_coarse_str.degree(n)
    print(f"   #{rank}  bc={v:.4f}  deg={deg:3d}  {n}")

is_artic = list(nx.articulation_points(G_coarse_str))
print(f"  Articulation points in real coarse connectome: {len(is_artic)}")
if is_artic:
    print("    (none expected if rich-club provides redundancy)")

# ---------------------------------------------------------------------------
# 5. Aggregate H1 test against REAL connectome
# ---------------------------------------------------------------------------
print("\n--- H1: Joined Trees vs REAL coarsened connectome ---")
m_real = all_metrics(G_coarse_str, "REAL_coarse")
m_joined = all_metrics(joined_trees(shared_daath=True), "JOINED_node")
m_synth = all_metrics(brain_model_10(), "BRAIN_10_synth")

print("\n   metric                JOINED_node  BRAIN_10_synth  REAL_coarse")
for k in INVARIANT_KEYS:
    if k in m_real and k in m_joined:
        v1, v2, v3 = m_joined.get(k), m_synth.get(k), m_real.get(k)
        if all(isinstance(v, (int, float)) and not np.isnan(v) and not np.isinf(v)
               for v in (v1, v2, v3)):
            print(f"   {k:<22s} {v1:>12.4f} {v2:>15.4f} {v3:>12.4f}")

d_joined_real = d_invariant(m_joined, m_real)
d_synth_real  = d_invariant(m_synth,  m_real)
d_joined_synth = d_invariant(m_joined, m_synth)
print(f"\n   d(JOINED_node, REAL_coarse) = {d_joined_real:.4f}")
print(f"   d(BRAIN_10_synth, REAL_coarse) = {d_synth_real:.4f}")
print(f"   d(JOINED_node, BRAIN_10_synth) = {d_joined_synth:.4f}")

if d_joined_real < d_synth_real:
    print(f"   -> Joined trees are CLOSER to real connectome than synthetic brain!")
else:
    diff = d_joined_real - d_synth_real
    print(f"   -> Synthetic brain is closer to real ({diff:+.4f} difference)")

# Null distribution: random graphs of same N, M
print("\n--- Null distribution: 200 random graphs vs REAL ---")
nulls = all_nulls(G_coarse_str, k_per=200)
null_dists = []
for null_name, ng_list in nulls.items():
    for ng in ng_list:
        nm = all_metrics(ng, "rand")
        d = d_invariant(nm, m_real)
        null_dists.append((null_name, d))
import collections
by_null = collections.defaultdict(list)
for n, d in null_dists:
    by_null[n].append(d)
for n in ("ER", "CFG", "WS", "BA", "GEO"):
    arr = np.array(by_null[n])
    if len(arr):
        p = float((arr <= d_joined_real).mean())
        print(f"   {n}: mean={arr.mean():.4f} std={arr.std():.4f}  "
              f"p(joined trees <= null)={p:.4f}")

# ---------------------------------------------------------------------------
# 6. Save results
# ---------------------------------------------------------------------------
out = {
    "data_source": "Budapest Reference Connectome v2.0 (Szalkai et al. 2015)",
    "subjects": "477 HCP subjects, consensus all_20k",
    "full_connectome": {
        "nodes": G_full.number_of_nodes(),
        "edges": G_full.number_of_edges(),
        "intra_hemi_edges": intra_hemi_edges,
        "inter_hemi_edges": inter_hemi_edges,
    },
    "coarsened_connectome": {
        "nodes": G_coarse_str.number_of_nodes(),
        "edges": G_coarse_str.number_of_edges(),
        "articulation_points": len(is_artic),
        "top_betweenness_node": top[0][0],
        "top_betweenness_value": top[0][1],
    },
    "h2_top_bc_nodes_full": [
        {"rank": i+1, "node": node_info[n]["fsname"],
         "hemi": node_info[n]["hemi"], "bc": v,
         "degree": G_full.degree(n)}
        for i, (n, v) in enumerate(nodes_by_bc[:10])
    ],
    "h2_progressive_disconnection": (lambda: {
        f"top_{k}_removed_LR_pairs_pct": (
            lambda g, drop: (g.remove_nodes_from(drop),
                             100 * hemi_connectivity(g) / baseline)[1]
        )(G_full.copy(), [n for n, _ in nodes_by_bc[:k]])
        for k in (1, 5, 10, 20, 50)
    })(),
    "h1_distances": {
        "joined_trees_to_real": d_joined_real,
        "synth_brain_to_real": d_synth_real,
        "joined_trees_to_synth_brain": d_joined_synth,
    },
    "null_p_values": {
        n: float((np.array(by_null[n]) <= d_joined_real).mean())
        for n in ("ER", "CFG", "WS", "BA", "GEO") if by_null[n]
    },
    "metrics": {
        "JOINED_node": {k: m_joined.get(k) for k in INVARIANT_KEYS},
        "BRAIN_10_synth": {k: m_synth.get(k) for k in INVARIANT_KEYS},
        "REAL_coarse": {k: m_real.get(k) for k in INVARIANT_KEYS},
    },
}
with open("data/real_connectome_results.json", "w") as f:
    json.dump(out, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else
                                 (float(o) if isinstance(o, (np.floating,)) else str(o)))
print(f"\nWrote data/real_connectome_results.json")
print("\nDONE.")
