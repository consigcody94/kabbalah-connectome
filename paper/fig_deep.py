"""Figures 10 (expanded mythology) and 11 (sensitivity)."""
import json
import matplotlib.pyplot as plt
import numpy as np

with open("data/deep_research_results.json", encoding="utf-8") as f:
    d = json.load(f)

# ---- Fig 10: expanded mythology -----
myth = d["expanded_mythology"]
ordered = sorted(myth.items(), key=lambda x: x[1]["d"])
names = [k for k, _ in ordered]
ds = [v["d"] for _, v in ordered]
ps = [v["p_ER"] for _, v in ordered]

fig, ax = plt.subplots(figsize=(13, 9))
y = np.arange(len(names))
colors = []
for n, p in zip(names, ps):
    if "TreeOfLife (joined+Daath)" in n: colors.append("#cc1111")
    elif "TreeOfLife" in n: colors.append("#cc6633")
    elif p < 0.05: colors.append("#33aa55")  # significant
    else: colors.append("#999")
bars = ax.barh(y, ds, color=colors, edgecolor="black", linewidth=0.5)
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=10)
ax.invert_yaxis()
ax.set_xlabel("topological distance to real human brain (lower = better)")
ax.axvline(0.4, color="green", linestyle="--", lw=0.6, label="d=0.4")
ax.axvline(0.6, color="orange", linestyle="--", lw=0.6, label="d=0.6")
for bar, p in zip(bars, ps):
    sig = "✓ p<0.05" if p < 0.05 else f"p={p:.2f}"
    ax.text(bar.get_width() + 0.005, bar.get_y() + bar.get_height()/2,
            sig, va="center", fontsize=8,
            color="green" if p < 0.05 else "gray")
ax.set_title("Expanded comparative-mythology test (12 structures)\n"
             "Red = Tree of Life family · Green = significantly matches brain · Gray = does not", fontsize=11)
ax.set_xlim(0, max(ds) * 1.15)
plt.tight_layout()
plt.savefig("figures/fig10_expanded_mythology.png", dpi=160, bbox_inches="tight")
print("Saved fig10_expanded_mythology.png")

# ---- Fig 11: sensitivity -----
sens = d["sensitivity"]
baseline = d["baseline_d"]

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for ax, kind, color, title in [
    (axes[0], "add", "#cc3333", "Add k random edges"),
    (axes[1], "remove", "#3366cc", "Remove k random edges"),
    (axes[2], "rewire", "#aa6633", "Rewire k random edges"),
]:
    pts = sens.get(kind, [])
    if not pts: continue
    ks = [p["k"] for p in pts]
    means = [p["mean"] for p in pts]
    stds = [p["std"] for p in pts]
    ax.errorbar(ks, means, yerr=stds, marker="o", color=color, lw=2, capsize=4)
    ax.axhline(baseline, color="black", linestyle="--", lw=1,
               label=f"baseline d={baseline:.3f}")
    ax.axhline(0.6, color="red", linestyle=":", lw=0.7,
               label="random-graph mean ~0.6")
    ax.set_xlabel("k (number of edges perturbed)")
    ax.set_ylabel("d to brain")
    ax.set_title(title)
    ax.set_ylim(0.1, 0.7)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

fig.suptitle("Sensitivity / perturbation analysis: how robust is the brain-match?\n"
             "Result is robust to small edge changes but degrades with rewiring",
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig("figures/fig11_sensitivity.png", dpi=160, bbox_inches="tight")
print("Saved fig11_sensitivity.png")
