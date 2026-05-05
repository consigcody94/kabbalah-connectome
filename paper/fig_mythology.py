"""Figure 9: comparative mythology test."""
import json
import matplotlib.pyplot as plt
import numpy as np

with open("data/comparative_mythology_results.json") as f:
    d = json.load(f)

# Order by distance
ordered = sorted(d.items(), key=lambda x: x[1]["d_real"])
names = [k for k, _ in ordered]
ds = [v["d_real"] for _, v in ordered]
p_ers = [v["p_ER"] for _, v in ordered]
p_geos = [v["p_GEO"] for _, v in ordered]

fig, axes = plt.subplots(1, 2, figsize=(18, 7))

# Left: distance bars
ax = axes[0]
colors = []
for n in names:
    if "TreeOfLife" in n: colors.append("#cc3333")
    elif "Yggdrasil" in n: colors.append("#cc8833")
    else: colors.append("#888")
y = np.arange(len(names))
ax.barh(y, ds, color=colors, edgecolor="black", linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("topological distance to real human brain (lower = better match)")
ax.axvline(0.0, color="black", lw=0.5)
ax.axvline(0.4, color="green", linestyle="--", lw=0.7, label="d = 0.4")
ax.axvline(0.6, color="orange", linestyle="--", lw=0.7, label="d = 0.6")
ax.axvline(0.8, color="red", linestyle="--", lw=0.7, label="d = 0.8")
ax.set_title("Comparative mythology test: which world structures match real brain?\n"
             "Tree of Life is uniquely closest; Yggdrasil also significantly matches",
             fontsize=11)
ax.legend(fontsize=8, loc="lower right")

# Right: significance vs nulls
ax = axes[1]
x = np.arange(len(names))
w = 0.35
ax.bar(x - w/2, p_ers, w, label="p (vs ER null)", color="#5566cc", edgecolor="black", linewidth=0.4)
ax.bar(x + w/2, p_geos, w, label="p (vs GEO null)", color="#aa66aa", edgecolor="black", linewidth=0.4)
ax.set_xticks(x)
ax.set_xticklabels([n[:18] for n in names], rotation=30, ha="right", fontsize=8)
ax.axhline(0.05, color="red", linestyle="--", lw=0.7, label="p = 0.05")
ax.set_ylabel("p (fraction of nulls at least as close to brain)")
ax.set_title("Statistical significance against random nulls")
ax.legend(fontsize=8)
ax.set_ylim(0, 1.1)

plt.tight_layout()
plt.savefig("figures/fig9_mythology_comparison.png", dpi=160, bbox_inches="tight")
print("Saved figures/fig9_mythology_comparison.png")
