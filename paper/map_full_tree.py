"""
Map every sephira (and qliphothic counterpart) to its best-matching
brain region in the Budapest connectome.

For each tree node:
  1. Compute topological role vector in joined-trees
  2. For every brain node, compute same vector
  3. Rank brain nodes by Euclidean distance in role-space (weighted)
  4. Report the top-K matches per sephira
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import (
    joined_trees, SEPHIROT_10, QLI_NAMES, POS_TOL,
)


# ---------------------------------------------------------------------------
# Load real connectome
# ---------------------------------------------------------------------------
print("Loading Budapest connectome...")
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
print(f"  Brain: N={G.number_of_nodes()} E={G.number_of_edges()}")


# ---------------------------------------------------------------------------
# Compute brain-node role vectors (cached from previous run idea)
# ---------------------------------------------------------------------------
print("Computing brain-node role vectors...")
all_bc_G = nx.betweenness_centrality(G)
all_cc_G = nx.closeness_centrality(G)
all_deg_G = dict(G.degree())
arts_G = set(nx.articulation_points(G))

deg_mean = np.mean(list(all_deg_G.values()))
deg_std = np.std(list(all_deg_G.values()))
bc_sorted = sorted(all_bc_G.values())
cc_sorted = sorted(all_cc_G.values())
bc_pct = {n: bc_sorted.index(all_bc_G[n])/(len(bc_sorted)-1) for n in G.nodes()}
cc_pct = {n: cc_sorted.index(all_cc_G[n])/(len(cc_sorted)-1) for n in G.nodes()}

# For each brain node, we need fraction-isolated. This is expensive (801
# graph copies), so precompute once.
print("  Computing fraction-isolated for all 801 brain nodes (slow)...")
frac_iso_G = {}
for i, n in enumerate(G.nodes()):
    if i % 100 == 0:
        print(f"    {i}/{G.number_of_nodes()}", flush=True)
    if n not in arts_G:
        frac_iso_G[n] = 0.0
        continue
    g2 = G.copy()
    g2.remove_node(n)
    if g2.number_of_nodes() == 0:
        frac_iso_G[n] = 0.0
    else:
        comps = list(nx.connected_components(g2))
        largest = max(len(c) for c in comps)
        frac_iso_G[n] = 1 - largest / g2.number_of_nodes()

brain_vectors = {}
for n in G.nodes():
    deg_z = (G.degree(n) - deg_mean) / (deg_std + 1e-9)
    art = 1.0 if n in arts_G else 0.0
    brain_vectors[n] = np.array([deg_z, bc_pct[n], cc_pct[n], art, frac_iso_G[n]])
print("  Done.")


# ---------------------------------------------------------------------------
# Compute joined-trees role vectors
# ---------------------------------------------------------------------------
print("\nComputing joined-trees role vectors...")
J = joined_trees(shared_daath=True)
all_bc_J = nx.betweenness_centrality(J)
all_cc_J = nx.closeness_centrality(J)
all_deg_J = dict(J.degree())
arts_J = set(nx.articulation_points(J))

deg_mean_J = np.mean(list(all_deg_J.values()))
deg_std_J = np.std(list(all_deg_J.values()))
bc_sorted_J = sorted(all_bc_J.values())
cc_sorted_J = sorted(all_cc_J.values())
bc_pct_J = {n: bc_sorted_J.index(all_bc_J[n])/(len(bc_sorted_J)-1) for n in J.nodes()}
cc_pct_J = {n: cc_sorted_J.index(all_cc_J[n])/(len(cc_sorted_J)-1) for n in J.nodes()}

tree_vectors = {}
for n in J.nodes():
    deg_z = (J.degree(n) - deg_mean_J) / (deg_std_J + 1e-9)
    art = 1.0 if n in arts_J else 0.0
    g2 = J.copy()
    g2.remove_node(n)
    if g2.number_of_nodes() == 0: frac_iso = 0.0
    else:
        comps = list(nx.connected_components(g2))
        frac_iso = 1 - max(len(c) for c in comps) / g2.number_of_nodes()
    tree_vectors[n] = np.array([deg_z, bc_pct_J[n], cc_pct_J[n], art, frac_iso])

print("Joined-trees node signatures:")
for n in (SEPHIROT_10 + ["Daath"] + QLI_NAMES):
    if n in tree_vectors:
        v = tree_vectors[n]
        print(f"  {n:<15s} deg_z={v[0]:+.2f}  BC={v[1]:.2f}  CC={v[2]:.2f}  "
              f"art={int(v[3])}  iso={v[4]:.2f}")


# ---------------------------------------------------------------------------
# For each sephira, find best brain matches
# ---------------------------------------------------------------------------
print("\n" + "=" * 80)
print("MAPPING EACH TREE NODE TO ITS BEST BRAIN MATCH")
print("=" * 80)

WEIGHTS = np.array([0.3, 1.5, 1.0, 2.0, 1.5])

mapping = {}
for tname, tvec in tree_vectors.items():
    scored = []
    for bname, bvec in brain_vectors.items():
        d = float(np.linalg.norm((tvec - bvec) * WEIGHTS))
        scored.append((bname, d))
    scored.sort(key=lambda x: x[1])
    top = []
    for bidx, dist in scored[:5]:
        info = node_info[bidx]
        top.append({
            "brain_idx": bidx,
            "brain_name": info["fsname"],
            "brain_hemi": info["hemi"],
            "distance": dist,
            "brain_degree": G.degree(bidx),
            "brain_bc_pct": bc_pct[bidx],
        })
    mapping[tname] = {
        "tree_signature": [float(x) for x in tvec],
        "top_matches": top,
    }

# Print canonical sephirot mapping
print(f"\n{'Sephira':<12} {'Best brain match':<35} {'hemi':>5} {'deg':>4} "
      f"{'BC%':>5} {'dist':>7}")
print("-" * 80)
order = ["Keter", "Chokhmah", "Binah", "Daath", "Chesed", "Geburah",
         "Tiferet", "Netzach", "Hod", "Yesod", "Malkuth"]
for n in order:
    best = mapping[n]["top_matches"][0]
    print(f"{n:<12} {best['brain_name'][:35]:<35} {best['brain_hemi']:>5} "
          f"{best['brain_degree']:>4} {best['brain_bc_pct']:>5.2f} "
          f"{best['distance']:>7.4f}")

print(f"\n{'Qliphah':<15} {'Best brain match':<35} {'hemi':>5} {'deg':>4} "
      f"{'BC%':>5} {'dist':>7}")
print("-" * 80)
for n in QLI_NAMES:
    best = mapping[n]["top_matches"][0]
    print(f"{n:<15} {best['brain_name'][:35]:<35} {best['brain_hemi']:>5} "
          f"{best['brain_degree']:>4} {best['brain_bc_pct']:>5.2f} "
          f"{best['distance']:>7.4f}")

# Save
with open("data/full_tree_mapping.json", "w") as f:
    json.dump(mapping, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/full_tree_mapping.json")


# ---------------------------------------------------------------------------
# Visualization: tree with each node labeled by its best brain region
# ---------------------------------------------------------------------------
print("\nGenerating fig7_full_mapping.png...")
fig, axes = plt.subplots(1, 2, figsize=(20, 10))

def short_brain_name(name):
    """Compact label"""
    return (name.replace("_", "").replace("middlefrontal", "MidFr")
                 .replace("inferiorfrontal", "InfFr")
                 .replace("superiorfrontal", "SupFr")
                 .replace("rostralmiddle", "rMid")
                 .replace("caudalmiddle", "cMid")
                 .replace("inferiorparietal", "InfPar")
                 .replace("superiorparietal", "SupPar")
                 .replace("supramarginal", "SupMrg")
                 .replace("inferiortemporal", "InfTmp")
                 .replace("superiortemporal", "SupTmp")
                 .replace("middletemporal", "MidTmp")
                 .replace("lateralorbitofrontal", "LOFC")
                 .replace("medialorbitofrontal", "MOFC")
                 .replace("parsopercularis", "pOp")
                 .replace("parstriangularis", "pTri")
                 .replace("parsorbitalis", "pOrb")
                 .replace("rostralanteriorcingulate", "rACC")
                 .replace("caudalanteriorcingulate", "cACC")
                 .replace("posteriorcingulate", "PCC")
                 .replace("isthmuscingulate", "isthCC")
                 .replace("precentral", "PreC")
                 .replace("postcentral", "PostC")
                 .replace("paracentral", "ParaC")
                 .replace("Left-", "L-")
                 .replace("Right-", "R-")
                 .replace("Hippocampus", "Hipp")
                 .replace("Amygdala", "Amyg")
                 .replace("Caudate", "Caud")
                 .replace("Putamen", "Put")
                 .replace("Thalamus-Proper", "Thal")
                 .replace("Brain-Stem", "BSt")
                 [:18])

# Left: Tree of Life with brain labels
ax = axes[0]
pos = dict(POS_TOL)
labels_tol = {}
for n in (SEPHIROT_10 + ["Daath"]):
    if n in mapping:
        bm = mapping[n]["top_matches"][0]
        labels_tol[n] = f"{n}\n=\n{short_brain_name(bm['brain_name'])}"

# Subgraph for just TOL+Daath
TOL = J.subgraph([n for n in (SEPHIROT_10 + ["Daath"]) if n in J]).copy()
nx.draw_networkx_nodes(TOL, pos, ax=ax, node_color="#f4e3a1",
                       node_size=4500, edgecolors="black", linewidths=1)
nx.draw_networkx_edges(TOL, pos, ax=ax, edge_color="#888", width=1.5)
for n, lbl in labels_tol.items():
    ax.text(pos[n][0], pos[n][1], lbl, ha="center", va="center", fontsize=7)
ax.set_title("Tree of Life → Brain mapping\n(each sephira's best-match brain region)",
             fontsize=12)
ax.set_aspect("equal")
ax.axis("off")

# Right: distance distribution showing how close each match is
ax = axes[1]
sephira_order = ["Keter", "Chokhmah", "Binah", "Daath", "Chesed", "Geburah",
                 "Tiferet", "Netzach", "Hod", "Yesod", "Malkuth"]
distances = [mapping[n]["top_matches"][0]["distance"] for n in sephira_order]
brain_names = [short_brain_name(mapping[n]["top_matches"][0]["brain_name"])
               for n in sephira_order]
y = np.arange(len(sephira_order))
colors = []
for n in sephira_order:
    bm = mapping[n]["top_matches"][0]
    h = bm["brain_hemi"]
    if h == "left": colors.append("#3366aa")
    elif h == "right": colors.append("#aa3333")
    else: colors.append("#888")
ax.barh(y, distances, color=colors, edgecolor="black", linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels([f"{s} → {b}" for s, b in zip(sephira_order, brain_names)],
                   fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("topological role distance (lower = better match)")
ax.set_title("Match quality per sephira\n(blue = left hemi, red = right hemi)")
ax.axvline(0.6, color="green", linestyle="--", lw=0.8, label="d=0.6")
ax.axvline(0.8, color="orange", linestyle="--", lw=0.8, label="d=0.8")
ax.axvline(1.0, color="red", linestyle="--", lw=0.8, label="d=1.0")
ax.legend(fontsize=8)

plt.tight_layout()
plt.savefig("figures/fig7_full_mapping.png", dpi=160, bbox_inches="tight")
print("Saved figures/fig7_full_mapping.png")
print("\nDONE.")
