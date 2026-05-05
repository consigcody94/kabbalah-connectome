"""Generate figure 6: brain regions ranked by Daath-role similarity."""
import json
import matplotlib.pyplot as plt
import numpy as np

with open("data/daath_localization.json") as f:
    d = json.load(f)

top30 = d["top_30_brain_matches"]

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# Left: bar chart of top 15 regions by distance
ax = axes[0]
names = [r["name"][:30] for r in top30[:15]]
distances = [r["distance_to_daath"] for r in top30[:15]]
colors = []
for r in top30[:15]:
    n = r["name"].lower()
    if "amygdal" in n: colors.append("#cc3333")
    elif "thalam" in n or "caudate" in n or "putamen" in n or "pallid" in n:
        colors.append("#9933aa")
    elif "hippo" in n: colors.append("#3366aa")
    elif "frontal" in n or "pars" in n or "rostral" in n or "precentral" in n:
        colors.append("#aa6633")
    elif "parietal" in n or "supramarg" in n or "precune" in n:
        colors.append("#6699aa")
    elif "temporal" in n or "fusiform" in n: colors.append("#aaaa33")
    else: colors.append("#888")
y = np.arange(len(names))
ax.barh(y, distances, color=colors, edgecolor="black", linewidth=0.4)
ax.set_yticks(y)
ax.set_yticklabels(names, fontsize=8)
ax.invert_yaxis()
ax.set_xlabel("topological distance to Daath signature (lower = better match)")
ax.set_title("Top 15 brain regions by similarity\nto Daath's role in joined-trees graph")
ax.axvline(0, color="black", lw=0.5)

# Right: category summary as donut
ax = axes[1]
cat = d["category_summary"]
labels = list(cat.keys())
sizes = [cat[k]["count"] for k in labels]
cmap = plt.cm.tab20
clist = [cmap(i / len(labels)) for i in range(len(labels))]
wedges, texts, autotexts = ax.pie(
    sizes, labels=labels, autopct="%1.0f%%",
    colors=clist, startangle=90, textprops={"fontsize": 9},
    pctdistance=0.78,
)
centre_circle = plt.Circle((0, 0), 0.45, fc="white", ec="black", lw=0.5)
ax.add_artist(centre_circle)
ax.text(0, 0, f"top 50\nDaath-like\nbrain regions",
        ha="center", va="center", fontsize=11, fontweight="bold")
ax.set_aspect("equal")
ax.set_title("Anatomical category of top-50 matches")

plt.tight_layout()
plt.savefig("figures/fig6_daath_localization.png", dpi=160, bbox_inches="tight")
print("Saved figures/fig6_daath_localization.png")
