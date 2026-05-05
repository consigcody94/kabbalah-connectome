"""
Main analysis pipeline.

Outputs all numbers, tables, and figures cited in the paper. Intended to be
run from the paper/ directory:

    python3 run_analysis.py

Will write:
    data/metrics.csv
    data/null_distances.csv
    data/node_role_comparison.csv
    figures/fig1_structures.png
    figures/fig2_metrics_radar.png
    figures/fig3_null_distribution.png
    figures/fig4_daath_vs_callosum.png
    figures/fig5_jacobs_ladder.png

The console output reproduces every number in the paper Results section.
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

THIS = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS))

from graphs import (
    SEPHIROT_10, QLI_NAMES, POS_TOL,
    tree_of_life, tree_of_death, joined_trees, jacobs_ladder,
    brain_model_10, brain_model_34, brain_ensemble,
)
from metrics import (
    all_metrics, d_invariant, d_spectral, d_portrait,
    node_role, INVARIANT_KEYS,
)
from nulls import all_nulls


DATA = THIS / "data"
FIG  = THIS / "figures"
DATA.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

K_NULLS = 200       # random graphs per null model
K_BRAIN = 100       # brain-model realizations for ensemble


def banner(s):
    print("\n" + "=" * 78)
    print(s)
    print("=" * 78)


# ---------------------------------------------------------------------------
# 1. Build all study graphs
# ---------------------------------------------------------------------------
banner("STEP 1: building study graphs")

graphs = {
    "TOL_10":      tree_of_life(include_daath=False),
    "TOL_11":      tree_of_life(include_daath=True),
    "QLI_10":      tree_of_death(include_daath=False),
    "QLI_11":      tree_of_death(include_daath=True),
    "JOINED_edge": joined_trees(shared_daath=False),   # 10 callosal edges
    "JOINED_node": joined_trees(shared_daath=True),    # Daath as bridge
    "JACOBS":      jacobs_ladder(),
    "BRAIN_10":    brain_model_10(corpus_callosum_node=True),
    "BRAIN_10nb":  brain_model_10(corpus_callosum_node=False),
    "BRAIN_34":    brain_model_34(seed=0, corpus_callosum_node=True),
}
for name, g in graphs.items():
    print(f"  {name:<14s} N={g.number_of_nodes():3d}  "
          f"M={g.number_of_edges():4d}  "
          f"connected={nx.is_connected(g)}")


# ---------------------------------------------------------------------------
# 2. Compute invariants for all graphs
# ---------------------------------------------------------------------------
banner("STEP 2: computing invariants")

rows = []
for name, g in graphs.items():
    print(f"  ... {name}")
    m = all_metrics(g, label=name)
    rows.append(m)

df = pd.DataFrame(rows).set_index("label")
df.to_csv(DATA / "metrics.csv")
print(f"\nWrote {DATA / 'metrics.csv'}")
print(df.T.round(4).to_string())


# ---------------------------------------------------------------------------
# 3. Brain ensemble (k=100 realizations of the 34-region model)
# ---------------------------------------------------------------------------
banner("STEP 3: brain ensemble (200 realizations)")

brains_34 = brain_ensemble(n=K_BRAIN, scale=34, corpus_callosum_node=True)
print(f"  built {len(brains_34)} brain-model realizations (scale=34)")

brain_metrics = [all_metrics(b, label=f"brain_{i}") for i, b in enumerate(brains_34)]
brain_df = pd.DataFrame(brain_metrics).set_index("label")
brain_df.to_csv(DATA / "brain_ensemble_metrics.csv")
print(f"  wrote {DATA / 'brain_ensemble_metrics.csv'}")
print("\nBrain ensemble (mean ± std):")
for k in INVARIANT_KEYS:
    if k in brain_df.columns:
        v = brain_df[k].dropna()
        if len(v):
            print(f"  {k:<32s} {v.mean():>8.4f} ± {v.std():.4f}")


# ---------------------------------------------------------------------------
# 4. Null model distance distributions vs brain ensemble mean
# ---------------------------------------------------------------------------
banner("STEP 4: null-model distance test")

# Reference target = ensemble-mean brain
brain_ref = {k: brain_df[k].mean() for k in brain_df.columns
             if pd.api.types.is_numeric_dtype(brain_df[k])}

candidates = ["TOL_11", "QLI_11", "JOINED_node", "JOINED_edge",
              "JACOBS", "BRAIN_10"]

# Match scale: use BRAIN_10 (21 nodes) as the comparison target for
# the small Trees, BRAIN_34 (~69 nodes) for Jacob's Ladder.
scale_match = {
    "TOL_11":      ("BRAIN_10",      brain_model_10()),
    "QLI_11":      ("BRAIN_10",      brain_model_10()),
    "JOINED_node": ("BRAIN_10",      brain_model_10()),
    "JOINED_edge": ("BRAIN_10",      brain_model_10()),
    "JACOBS":      ("BRAIN_34",      brain_model_34(seed=0)),
    "BRAIN_10":    ("BRAIN_10",      brain_model_10()),
}

null_results = []
for cand_name in candidates:
    cand = graphs[cand_name]
    target_name, target = scale_match[cand_name]
    target_metrics = all_metrics(target, label=target_name)
    cand_metrics = all_metrics(cand, label=cand_name)
    real_d = d_invariant(cand_metrics, target_metrics)
    real_dsp = d_portrait(cand, target)

    nulls = all_nulls(target, k_per=K_NULLS)
    for null_name, null_graphs in nulls.items():
        for ng in null_graphs:
            nm = all_metrics(ng, label=f"{null_name}_null")
            null_results.append({
                "candidate": cand_name,
                "target":    target_name,
                "null":      null_name,
                "d_inv":     d_invariant(nm, target_metrics),
                "d_portrait":d_portrait(ng, target),
            })
    null_results.append({
        "candidate": cand_name, "target": target_name,
        "null": "REAL", "d_inv": real_d, "d_portrait": real_dsp,
    })

null_df = pd.DataFrame(null_results)
null_df.to_csv(DATA / "null_distances.csv", index=False)
print(f"  wrote {DATA / 'null_distances.csv'} ({len(null_df)} rows)")

# Compute per-candidate p-values
banner("STEP 4b: candidate vs null-distribution p-values")
print(f"{'candidate':<14s} {'target':<12s} {'d_inv':>8s} "
      f"{'p_inv (vs ER)':>14s} {'p_inv (vs CFG)':>16s} {'p_inv (vs WS)':>14s}")
for cand_name in candidates:
    real = null_df[(null_df.candidate == cand_name) &
                   (null_df.null == "REAL")].iloc[0]
    line = [cand_name, scale_match[cand_name][0], f"{real.d_inv:.4f}"]
    for null_name in ("ER", "CFG", "WS"):
        sub = null_df[(null_df.candidate == cand_name) &
                      (null_df.null == null_name)]
        if len(sub):
            p = float((sub.d_inv <= real.d_inv).mean())
        else:
            p = float("nan")
        line.append(f"{p:>14.4f}")
    print("  " + "  ".join(line))


# ---------------------------------------------------------------------------
# 5. Node-role comparison: Daath vs Corpus Callosum
# ---------------------------------------------------------------------------
banner("STEP 5: bridge-node role comparison (Daath vs Corpus Callosum)")

daath_role = node_role(graphs["JOINED_node"], "Daath")
cc10_role  = node_role(graphs["BRAIN_10"],   "CorpusCallosum")
# Average over brain ensemble for CC role
cc_roles_ens = [node_role(b, "CorpusCallosum") for b in brains_34
                if "CorpusCallosum" in b]

role_df = pd.DataFrame([daath_role, cc10_role]).set_index("node")
role_df.to_csv(DATA / "node_role_comparison.csv")
print(role_df.T.round(4).to_string())

if cc_roles_ens:
    print("\nCorpus callosum role across brain-34 ensemble (mean ± std):")
    for k in ("degree", "betweenness", "betweenness_rank",
              "fraction_isolated_after"):
        vals = [r[k] for r in cc_roles_ens]
        print(f"  {k:<28s} {np.mean(vals):>8.4f} ± {np.std(vals):.4f}")


# ---------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------
banner("STEP 6: rendering figures")

# Figure 1 — the four structures
fig, axes = plt.subplots(1, 4, figsize=(22, 7))
for ax, (name, posfn) in zip(axes, [
    ("TOL_11", "tree"),
    ("QLI_11", "tree-mirror"),
    ("JOINED_node", "joined"),
    ("BRAIN_10", "brain"),
]):
    g = graphs[name]
    if posfn == "tree":
        pos = POS_TOL
        nc, ec = "#f4e3a1", "#b39656"
        fc = "black"
    elif posfn == "tree-mirror":
        pos = {q: (-x, y) for q, (x, y) in zip(
            QLI_NAMES + ["Daath"], list(POS_TOL.values()))}
        # Build mirrored positions specifically
        pos = {}
        for s, q in zip(SEPHIROT_10, QLI_NAMES):
            x, y = POS_TOL[s]
            pos[q] = (-x, y)
        pos["Daath"] = (0.0, 3.6)
        nc, ec, fc = "#3a3a3a", "#7a1a1a", "white"
    elif posfn == "joined":
        pos = dict(POS_TOL)
        for s, q in zip(SEPHIROT_10, QLI_NAMES):
            x, y = POS_TOL[s]
            pos[q] = (-x, y)
        pos["Daath"] = (0.0, 3.6)
        nc = ["#cc3333" if n == "Daath"
              else ("#f4e3a1" if n in SEPHIROT_10 else "#3a3a3a")
              for n in g.nodes()]
        ec = "#777"
        fc = "black"
    else:  # brain
        pos = {"CorpusCallosum": (0, 0)}
        from graphs import BRAIN10_REGIONS, _hemi
        for i, r in enumerate(BRAIN10_REGIONS):
            angle = 2 * np.pi * i / len(BRAIN10_REGIONS)
            pos[_hemi("L", r)] = (-3 + np.cos(angle), np.sin(angle))
            pos[_hemi("R", r)] = ( 3 + np.cos(angle), np.sin(angle))
        nc = ["#cc3333" if n == "CorpusCallosum" else "#ffd2b5"
              for n in g.nodes()]
        ec = "#777"
        fc = "black"

    if isinstance(nc, list):
        sizes = [1500 if (n == "Daath" or n == "CorpusCallosum") else 700
                 for n in g.nodes()]
        nx.draw_networkx_nodes(g, pos, ax=ax, node_color=nc,
                               node_size=sizes, edgecolors="black",
                               linewidths=0.5)
    else:
        nx.draw_networkx_nodes(g, pos, ax=ax, node_color=nc, node_size=900,
                               edgecolors="black", linewidths=0.5)
    nx.draw_networkx_edges(g, pos, ax=ax, edge_color=ec, width=0.9)
    nx.draw_networkx_labels(g, pos, ax=ax, font_size=5, font_color=fc)
    ax.set_title({"TOL_11": "(a) Tree of Life (Kircher, 11 sephirot)",
                  "QLI_11": "(b) Tree of Death / Qliphoth",
                  "JOINED_node": "(c) Joined trees through shared Daath",
                  "BRAIN_10": "(d) Brain (10 regions × 2 + corpus callosum)"
                  }[name], fontsize=10)
    ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout()
plt.savefig(FIG / "fig1_structures.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"  wrote {FIG / 'fig1_structures.png'}")

# Figure 2 — radar/bar of normalized invariants
fig, ax = plt.subplots(figsize=(13, 6))
keys = ["density", "mean_degree", "diameter", "char_path",
        "avg_clustering", "transitivity",
        "max_betweenness", "modularity_Q", "spectral_radius"]
labels = ["TOL_11", "JOINED_node", "BRAIN_10", "JACOBS"]
data = []
for k in keys:
    row = []
    for lbl in labels:
        v = df.loc[lbl, k]
        row.append(v)
    data.append(row)
data = np.array(data)
# Normalize each metric column-wise so bars are comparable
data_norm = data / (np.max(np.abs(data), axis=1, keepdims=True) + 1e-9)
x = np.arange(len(keys))
w = 0.2
colors = ["#d0a050", "#cc3333", "#3370b0", "#5fa05f"]
for i, lbl in enumerate(labels):
    ax.bar(x + (i - 1.5) * w, data_norm[:, i], w,
           label=lbl, color=colors[i])
ax.set_xticks(x)
ax.set_xticklabels(keys, rotation=30, ha="right")
ax.axhline(0, color="black", lw=0.5)
ax.set_ylabel("normalized value (per metric)")
ax.legend()
ax.set_title("Figure 2. Normalized invariants across structures.\n"
             "Bars within each metric are scaled to the maximum across structures.")
plt.tight_layout()
plt.savefig(FIG / "fig2_metrics_radar.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"  wrote {FIG / 'fig2_metrics_radar.png'}")

# Figure 3 — null distributions
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
plot_cands = ["TOL_11", "QLI_11", "JOINED_node",
              "JOINED_edge", "JACOBS", "BRAIN_10"]
for ax, cand_name in zip(axes.flat, plot_cands):
    sub = null_df[null_df.candidate == cand_name]
    real_d = sub[sub.null == "REAL"].d_inv.iloc[0]
    for null_name, color in [("ER", "#888"), ("CFG", "#5566cc"),
                             ("WS", "#cc6655"), ("BA", "#55aa55"),
                             ("GEO", "#aa66aa")]:
        vals = sub[sub.null == null_name].d_inv.dropna().values
        if len(vals):
            ax.hist(vals, bins=30, color=color, alpha=0.45,
                    label=f"{null_name} (n={len(vals)})")
    ax.axvline(real_d, color="red", lw=2, label=f"real d={real_d:.3f}")
    ax.set_title(cand_name)
    ax.set_xlabel("invariant distance to brain target")
    ax.set_ylabel("count")
    ax.legend(fontsize=7)
fig.suptitle("Figure 3. Distribution of invariant-distances to the brain "
             "target for each candidate structure across five null models.\n"
             "Red line = the candidate's actual distance. p-values in "
             "Table 2 of the paper.", fontsize=11)
plt.tight_layout()
plt.savefig(FIG / "fig3_null_distribution.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"  wrote {FIG / 'fig3_null_distribution.png'}")

# Figure 4 — Daath vs Corpus Callosum (recreate)
fig, axes = plt.subplots(1, 2, figsize=(18, 8))

# Joined trees with Daath highlighted
ax = axes[0]
g = graphs["JOINED_node"]
pos = dict(POS_TOL)
for s, q in zip(SEPHIROT_10, QLI_NAMES):
    x, y = POS_TOL[s]
    pos[q] = (-x, y)
pos["Daath"] = (0.0, 3.6)
nc = ["#cc3333" if n == "Daath"
      else ("#f4e3a1" if n in SEPHIROT_10 else "#3a3a3a")
      for n in g.nodes()]
sizes = [2400 if n == "Daath" else 900 for n in g.nodes()]
nx.draw_networkx_nodes(g, pos, ax=ax, node_color=nc, node_size=sizes,
                       edgecolors="black", linewidths=1.0)
nx.draw_networkx_edges(g, pos, ax=ax, edge_color="#777", width=1)
white_labels = {n: n for n in g.nodes() if n in QLI_NAMES}
black_labels = {n: n for n in g.nodes() if n not in QLI_NAMES}
nx.draw_networkx_labels(g, pos, labels=white_labels, ax=ax,
                        font_size=6, font_color="white")
nx.draw_networkx_labels(g, pos, labels=black_labels, ax=ax,
                        font_size=6, font_color="black")
ax.set_title("(a) Joined Trees through shared Daath\n"
             f"betweenness rank: {daath_role['betweenness_rank']} of "
             f"{g.number_of_nodes()}; removes into "
             f"{daath_role['components_after_removal']} components",
             fontsize=11)
ax.set_aspect("equal"); ax.axis("off")

# Brain with corpus callosum highlighted
ax = axes[1]
g = graphs["BRAIN_10"]
pos = {"CorpusCallosum": (0, 0)}
from graphs import BRAIN10_REGIONS, _hemi
for i, r in enumerate(BRAIN10_REGIONS):
    angle = 2 * np.pi * i / len(BRAIN10_REGIONS)
    pos[_hemi("L", r)] = (-3 + np.cos(angle), np.sin(angle))
    pos[_hemi("R", r)] = ( 3 + np.cos(angle), np.sin(angle))
nc = ["#cc3333" if n == "CorpusCallosum" else "#ffd2b5" for n in g.nodes()]
sizes = [2400 if n == "CorpusCallosum" else 900 for n in g.nodes()]
nx.draw_networkx_nodes(g, pos, ax=ax, node_color=nc, node_size=sizes,
                       edgecolors="black", linewidths=1.0)
nx.draw_networkx_edges(g, pos, ax=ax, edge_color="#777", width=0.8)
nx.draw_networkx_labels(g, pos, ax=ax, font_size=6)
ax.set_title("(b) Brain with corpus callosum as bridge node\n"
             f"betweenness rank: {cc10_role['betweenness_rank']} of "
             f"{g.number_of_nodes()}; removes into "
             f"{cc10_role['components_after_removal']} components",
             fontsize=11)
ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout()
plt.savefig(FIG / "fig4_daath_vs_callosum.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"  wrote {FIG / 'fig4_daath_vs_callosum.png'}")

# Figure 5 — Jacob's Ladder
fig, ax = plt.subplots(figsize=(10, 12))
g = graphs["JACOBS"]
worlds = ["Atz", "Bri", "Yet", "Ass"]
y_off = {"Atz": 18, "Bri": 12, "Yet": 6, "Ass": 0}
pos = {}
for w in worlds:
    for s, (x, y) in POS_TOL.items():
        if s == "Daath": continue
        pos[f"{w}_{s}"] = (x, y + y_off[w])
nx.draw_networkx_nodes(g, pos, ax=ax, node_color="#f4e3a1", node_size=350,
                       edgecolors="black", linewidths=0.4)
# Highlight bridge edges (Malkuth_n -> Keter_n+1)
bridge_edges = [(f"{a}_Malkuth", f"{b}_Keter")
                for a, b in zip(worlds[:-1], worlds[1:])]
nx.draw_networkx_edges(g, pos, ax=ax, edge_color="#999", width=0.7)
nx.draw_networkx_edges(g, pos, edgelist=bridge_edges, ax=ax,
                       edge_color="red", width=2.5, style="dashed")
nx.draw_networkx_labels(g, pos, ax=ax, font_size=4)
ax.set_title("Figure 5. Jacob's Ladder (Four Worlds)\n"
             "Each world is a Tree of Life; adjacent worlds share Malkuth↔Keter (red)",
             fontsize=11)
ax.set_aspect("equal"); ax.axis("off")
plt.tight_layout()
plt.savefig(FIG / "fig5_jacobs_ladder.png", dpi=160, bbox_inches="tight")
plt.close()
print(f"  wrote {FIG / 'fig5_jacobs_ladder.png'}")


# ---------------------------------------------------------------------------
# 7. Save key numbers for the paper
# ---------------------------------------------------------------------------
banner("STEP 7: writing summary JSON")
summary = {
    "candidates": candidates,
    "graphs": {name: {"N": g.number_of_nodes(), "M": g.number_of_edges()}
               for name, g in graphs.items()},
    "daath_role": daath_role,
    "corpus_callosum_role_static": cc10_role,
    "corpus_callosum_role_ensemble_mean": (
        {k: float(np.mean([r[k] for r in cc_roles_ens]))
         for k in ("degree", "betweenness", "betweenness_rank",
                   "fraction_isolated_after")}
        if cc_roles_ens else {}
    ),
    "p_values": {},
}
for cand_name in candidates:
    real = null_df[(null_df.candidate == cand_name) &
                   (null_df.null == "REAL")].iloc[0]
    summary["p_values"][cand_name] = {"d_inv_real": float(real.d_inv)}
    for null_name in ("ER", "CFG", "WS", "BA", "GEO"):
        sub = null_df[(null_df.candidate == cand_name) &
                      (null_df.null == null_name)]
        if len(sub):
            p = float((sub.d_inv <= real.d_inv).mean())
            summary["p_values"][cand_name][f"p_{null_name}"] = p
with open(DATA / "summary.json", "w") as f:
    json.dump(summary, f, indent=2, default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"  wrote {DATA / 'summary.json'}")

print("\nDONE.")
