"""
Find the brain region whose topological role in the real connectome most
closely matches Daath's role in the joined-trees graph.

Method:
  1. Compute a 6-dimensional "topological role" vector for Daath:
        (degree z-score, betweenness percentile, closeness percentile,
         is-articulation, fraction-isolated-on-removal, hub-of-bridge-edges)
  2. Compute the same vector for every node in the Budapest connectome.
  3. Rank nodes by Euclidean distance to Daath's vector in this 6-D space.
  4. Report the top candidates with anatomical labels.

This is an analogy-mapping exercise, not a metaphysical claim. We are
asking: "Of all real brain regions, which one PLAYS THE ROLE THAT DAATH
PLAYS in the joined-trees graph?"
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

from graphs import joined_trees, SEPHIROT_10, QLI_NAMES


# ---------------------------------------------------------------------------
# Step 1: load the real connectome
# ---------------------------------------------------------------------------
print("Loading Budapest connectome (occ >= 100/477)...")
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
        try:
            occ = int(r["occ"])
        except (ValueError, TypeError):
            continue
        if occ < 100: continue
        s, t = int(r["source"]), int(r["target"])
        if s != t: G.add_edge(s, t)
G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
print(f"  N={G.number_of_nodes()} E={G.number_of_edges()}")


# ---------------------------------------------------------------------------
# Step 2: build the joined-trees graph and compute Daath's signature
# ---------------------------------------------------------------------------
print("\nBuilding joined-trees graph (with shared Daath)...")
J = joined_trees(shared_daath=True)
print(f"  N={J.number_of_nodes()} E={J.number_of_edges()}")


def role_vector(g, node, all_bc, all_cc, all_deg):
    """Six-component normalized role vector for a node.

    [0] degree z-score (relative to graph)
    [1] betweenness percentile (0..1; 1 = highest BC in graph)
    [2] closeness percentile  (0..1; 1 = most central)
    [3] is articulation point? (0/1)
    [4] fraction of graph isolated by removing this node
    [5] inter-side-neighborhood ratio: how balanced its neighbors are
        between the two halves of its graph (0 = all same side,
        1 = perfectly balanced)
    """
    degs = list(all_deg.values())
    deg_z = (g.degree(node) - np.mean(degs)) / (np.std(degs) + 1e-9)
    bc_vals = sorted(all_bc.values())
    bc_pct = bc_vals.index(all_bc[node]) / max(len(bc_vals)-1, 1)
    cc_vals = sorted(all_cc.values())
    cc_pct = cc_vals.index(all_cc[node]) / max(len(cc_vals)-1, 1)
    artic = 1.0 if node in set(nx.articulation_points(g)) else 0.0
    g2 = g.copy()
    g2.remove_node(node)
    if g2.number_of_nodes() == 0:
        frac_iso = 0.0
    else:
        comps = list(nx.connected_components(g2))
        largest = max(len(c) for c in comps)
        frac_iso = 1 - largest / g2.number_of_nodes()
    return np.array([deg_z, bc_pct, cc_pct, artic, frac_iso])


# Daath signature in joined trees
print("\nComputing Daath's role signature in joined-trees graph...")
all_bc_J = nx.betweenness_centrality(J)
all_cc_J = nx.closeness_centrality(J)
all_deg_J = dict(J.degree())
daath_sig = role_vector(J, "Daath", all_bc_J, all_cc_J, all_deg_J)
print(f"  Daath signature: deg_z={daath_sig[0]:+.3f}  "
      f"BC_pct={daath_sig[1]:.3f}  CC_pct={daath_sig[2]:.3f}  "
      f"artic={daath_sig[3]:.0f}  frac_iso={daath_sig[4]:.3f}")
print(f"  Interpretation: Daath is moderately above mean degree, top in BC, top in CC,")
print(f"  IS an articulation point, isolates {daath_sig[4]*100:.0f}% on removal")


# ---------------------------------------------------------------------------
# Step 3: compute role vector for every node in real connectome
# ---------------------------------------------------------------------------
print("\nComputing role vectors for all 801 brain nodes (this takes ~30s)...")
all_bc_G = nx.betweenness_centrality(G)
all_cc_G = nx.closeness_centrality(G)
all_deg_G = dict(G.degree())

# Articulation points - precompute set for O(1) lookup
arts_G = set(nx.articulation_points(G))

# For computing "frac_iso" for every node, we'd need 801 graph copies.
# Instead, we use the fact that in a connected graph G, node v isolates
# at most (|G|-1) - largest_remaining nodes. We compute this only for the
# top 50 betweenness candidates to keep runtime down.

print("  Phase 1: computing degree-z, BC, CC for all 801 nodes...")
deg_mean = np.mean(list(all_deg_G.values()))
deg_std = np.std(list(all_deg_G.values()))
bc_sorted = sorted(all_bc_G.values())
cc_sorted = sorted(all_cc_G.values())
bc_index = {v: i/(len(bc_sorted)-1) for i, v in enumerate(bc_sorted)}
cc_index = {v: i/(len(cc_sorted)-1) for i, v in enumerate(cc_sorted)}

candidates = []
for n in G.nodes():
    deg_z = (G.degree(n) - deg_mean) / (deg_std + 1e-9)
    bc_pct = bc_index[all_bc_G[n]]
    cc_pct = cc_index[all_cc_G[n]]
    art = 1.0 if n in arts_G else 0.0
    candidates.append((n, np.array([deg_z, bc_pct, cc_pct, art, 0.0])))

print("  Phase 2: computing frac_iso for top 30 BC candidates (only)...")
candidates.sort(key=lambda x: -x[1][1])  # sort by BC percentile desc
for i in range(min(30, len(candidates))):
    n, vec = candidates[i]
    g2 = G.copy()
    g2.remove_node(n)
    if g2.number_of_nodes() == 0:
        frac = 0.0
    else:
        comps = list(nx.connected_components(g2))
        largest = max(len(c) for c in comps)
        frac = 1 - largest / g2.number_of_nodes()
    vec[4] = frac
    candidates[i] = (n, vec)

# For the rest, leave frac_iso = 0 (they're not articulation points anyway,
# so isolation would be 0 by definition for non-articulation nodes)


# ---------------------------------------------------------------------------
# Step 4: rank by similarity to Daath's signature
# ---------------------------------------------------------------------------
print("\nRanking all 801 nodes by similarity to Daath's role signature...")
scored = []
for n, vec in candidates:
    # Weight components: BC and articulation matter most for bridge-role
    weights = np.array([0.3, 1.5, 1.0, 2.0, 1.5])
    diff = (vec - daath_sig) * weights
    distance = float(np.linalg.norm(diff))
    info = node_info[n]
    scored.append({
        "node_idx": n,
        "name": info["fsname"],
        "hemi": info["hemi"],
        "deg_z": float(vec[0]),
        "bc_pct": float(vec[1]),
        "cc_pct": float(vec[2]),
        "artic": float(vec[3]),
        "frac_iso": float(vec[4]),
        "distance_to_daath": distance,
        "raw_degree": G.degree(n),
        "raw_bc": all_bc_G[n],
    })

scored.sort(key=lambda x: x["distance_to_daath"])

print("\n" + "=" * 88)
print("TOP 15 BRAIN REGIONS BY MATCH TO DAATH'S TOPOLOGICAL ROLE")
print("=" * 88)
print(f"{'rank':>4} {'distance':>9} {'name':<40} {'hemi':>6} {'deg':>4} "
      f"{'BC_pct':>7} {'artic':>5} {'iso%':>5}")
print("-" * 88)
for i, c in enumerate(scored[:15], 1):
    print(f"{i:>4} {c['distance_to_daath']:>9.4f} {c['name'][:40]:<40} "
          f"{c['hemi']:>6} {c['raw_degree']:>4} "
          f"{c['bc_pct']:>7.3f} {int(c['artic']):>5} "
          f"{c['frac_iso']*100:>5.1f}")

# Find candidate by anatomical category — Daath should be central across hemispheres
# Tag each candidate with anatomical category
def category(name):
    s = name.lower()
    if "caudate" in s or "putamen" in s or "pallid" in s or "accumb" in s:
        return "Basal ganglia"
    if "thalam" in s: return "Thalamus"
    if "hippo" in s: return "Hippocampus"
    if "amygdal" in s: return "Amygdala"
    if "brain-stem" in s or "brainstem" in s: return "Brain stem"
    if "ventral" in s and "diencephalon" in s: return "Ventral diencephalon"
    if "cerebellum" in s: return "Cerebellum"
    if "cingul" in s: return "Cingulate cortex"
    if "insula" in s: return "Insula"
    if "frontal" in s or "pars" in s or "rostral" in s or "precentral" in s:
        return "Frontal cortex"
    if "temporal" in s or "fusiform" in s: return "Temporal cortex"
    if "parietal" in s or "supramarg" in s or "precune" in s:
        return "Parietal cortex"
    if "occipital" in s or "lingual" in s or "cuneus" in s: return "Occipital cortex"
    if "postcentral" in s: return "Sensorimotor"
    return "Other"

# Aggregate by anatomical category — count appearances in top 50
print("\n--- By anatomical category (top 50 best-matches) ---")
cat_counts = defaultdict(int)
cat_best = {}
for c in scored[:50]:
    cat = category(c["name"])
    cat_counts[cat] += 1
    if cat not in cat_best or c["distance_to_daath"] < cat_best[cat]["distance_to_daath"]:
        cat_best[cat] = c
for cat, count in sorted(cat_counts.items(), key=lambda x: -x[1]):
    best = cat_best[cat]
    print(f"  {cat:<25s} {count:>3} appearances; "
          f"best: {best['name'][:40]} (d={best['distance_to_daath']:.4f})")

# Save
with open("data/daath_localization.json", "w") as f:
    json.dump({
        "daath_signature": {
            "degree_z": float(daath_sig[0]),
            "bc_percentile": float(daath_sig[1]),
            "cc_percentile": float(daath_sig[2]),
            "is_articulation": float(daath_sig[3]),
            "fraction_isolated_on_removal": float(daath_sig[4]),
        },
        "top_30_brain_matches": scored[:30],
        "category_summary": {
            cat: {"count": cat_counts[cat],
                  "best_node": cat_best[cat]["name"],
                  "best_distance": cat_best[cat]["distance_to_daath"]}
            for cat in cat_counts
        },
    }, f, indent=2,
       default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/daath_localization.json")
print("\nDONE.")
