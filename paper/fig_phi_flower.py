"""Figures 14 (golden ratio) and 15 (Flower of Life)."""
import json
import matplotlib.pyplot as plt
import numpy as np

# ---- Fig 14: Golden ratio ----
with open("data/golden_ratio_results.json") as f:
    g = json.load(f)

phi = g["phi"]
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: which ratios are close to phi?
ax = axes[0]
keys = list(g["joined_phi_distances"].keys())
J_dists = [g["joined_phi_distances"][k] for k in keys]
B_dists = [g["brain_phi_distances"].get(k, np.nan) for k in keys]

x = np.arange(len(keys))
w = 0.35
ax.bar(x - w/2, J_dists, w, label="Joined trees", color="#cc8833", edgecolor="black", linewidth=0.5)
ax.bar(x + w/2, B_dists, w, label="Real brain", color="#3366cc", edgecolor="black", linewidth=0.5)
ax.axhline(0.05, color="red", linestyle="--", lw=0.7, label="d = 0.05 (very close to φ)")
ax.axhline(0.1, color="orange", linestyle="--", lw=0.5, label="d = 0.1 (close)")
ax.set_xticks(x)
ax.set_xticklabels([k.replace("_", "\n") for k in keys], rotation=45, ha="right", fontsize=8)
ax.set_ylabel("Distance to nearest φ-related value (lower = more golden)")
ax.set_title(f"Golden ratio φ = {phi:.4f}: which structural ratios are close to it?")
ax.legend(fontsize=8)

# Right: spectrum bar
ax = axes[1]
J_eigs = g["joined_top_eigs"]
B_eigs = g["brain_top_eigs"]
x = np.arange(len(J_eigs))
w = 0.35
ax.bar(x - w/2, J_eigs, w, label="Joined trees", color="#cc8833",
       edgecolor="black", linewidth=0.5)
ax.bar(x + w/2, B_eigs, w, label="Real brain", color="#3366cc",
       edgecolor="black", linewidth=0.5)
ax.set_xticks(x)
ax.set_xticklabels([f"λ_{i+1}" for i in range(len(J_eigs))])
ax.set_ylabel("Eigenvalue magnitude")
ax.set_title("Top 5 eigenvalues of joined-trees vs real brain spectra")
ax.legend()
# Add ratio annotations
for i in range(len(J_eigs) - 1):
    if abs(J_eigs[i+1]) > 0.1:
        r_J = J_eigs[i] / J_eigs[i+1]
        ax.text(i + 0.05, max(J_eigs[i], J_eigs[i+1]) + 0.3,
                f"J: {r_J:.2f}", fontsize=8, color="#cc6611",
                ha="left")
    if abs(B_eigs[i+1]) > 0.1:
        r_B = B_eigs[i] / B_eigs[i+1]
        ax.text(i + 0.55, max(B_eigs[i], B_eigs[i+1]) + 1.5,
                f"B: {r_B:.2f}", fontsize=8, color="#3344aa",
                ha="left")

plt.tight_layout()
plt.savefig("figures/fig14_golden_ratio.png", dpi=160, bbox_inches="tight")
print("Saved fig14_golden_ratio.png")

# ---- Fig 15: Flower of Life ----
with open("data/flower_of_life_results.json") as f:
    fol = json.load(f)

fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: distance comparison
ax = axes[0]
labels_short = {
    "Joined trees (Kabbalah baseline)": "Joined trees\n(Kabbalah)",
    "Flower of Life (center graph)":     "Flower of Life\n(center graph)",
    "Tree of Life derived from FoL":     "Tree of Life\nderived from FoL",
}
names = list(fol.keys())
ds = [fol[n]["d_real"] for n in names]
ps = [fol[n]["p_ER"] for n in names]
colors = ["#cc3333", "#a07020", "#cc8833"]
y = np.arange(len(names))
ax.barh(y, ds, color=colors, edgecolor="black", linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels([labels_short.get(n, n) for n in names], fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("distance to real brain")
ax.axvline(0.05, color="red", linestyle="--", lw=0.5, label="p<0.05 threshold")
ax.set_title("Flower of Life vs Kabbalistic Tree of Life as graph topology")
for i, (d, p) in enumerate(zip(ds, ps)):
    ax.text(d + 0.01, i, f" d={d:.3f}, p={p:.3f}", va="center", fontsize=10)

# Right: render the FoL pattern
ax = axes[1]
import networkx as nx
sys_path = __import__('sys').path
import os
sys_path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flower_of_life_graph import flower_of_life_centers, fol_center_graph

centers = flower_of_life_centers(r=1.0)
g = fol_center_graph()
from matplotlib.patches import Circle
for cx, cy in centers:
    ax.add_patch(Circle((cx, cy), 1.0, fill=False,
                        edgecolor="#a07020", linewidth=0.8, alpha=0.5))
pos = {i: c for i, c in enumerate(centers)}
nx.draw_networkx_edges(g, pos, ax=ax, edge_color="#5555aa", width=0.4, alpha=0.6)
nx.draw_networkx_nodes(g, pos, ax=ax, node_size=80, node_color="#cc3333",
                       edgecolors="black", linewidths=0.5)
ax.set_xlim(-3.5, 3.5); ax.set_ylim(-3.5, 3.5)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title(f"Flower of Life as a graph\n"
             f"19 circles, {g.number_of_edges()} edges (overlap connections)")

plt.tight_layout()
plt.savefig("figures/fig15_flower_of_life.png", dpi=160, bbox_inches="tight")
print("Saved fig15_flower_of_life.png")
