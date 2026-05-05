"""
Critical test: is the joined-trees graph closer to real brain topology than
ARBITRARY MODULAR GRAPHS of matched size and modularity?

If yes -> the trees-are-brain-like result is non-trivial (specific to trees).
If no  -> the result reduces to "modular small-world matches modular small-
          world," which is a property of the topology class, not the trees.

Method:
  1. Generate 500 modular random graphs (stochastic block model) with the
     same N, M, and approximate modularity Q as the real coarsened brain.
  2. Compute distance-to-real-brain for each.
  3. Compare to joined-trees distance.
  4. Report p-value: fraction of modular nulls at least as close to brain
     as the joined-trees graph.

Also test: random-modular graphs at varying community counts (2, 4, 6, 8)
to see how brain-likeness scales with modularity structure.
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

from graphs import joined_trees
from metrics import all_metrics, d_invariant, INVARIANT_KEYS


# ---------------------------------------------------------------------------
# Load real coarsened connectome (same as previous tests)
# ---------------------------------------------------------------------------
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

m_real = all_metrics(G_coarse, "REAL_coarse")
N = G_coarse.number_of_nodes()
M = G_coarse.number_of_edges()
print(f"Real coarse brain: N={N} E={M} Q={m_real['modularity_Q']:.3f} "
      f"C={m_real['avg_clustering']:.3f}")

# Joined trees baseline
J = joined_trees(shared_daath=True)
m_joined = all_metrics(J, "JOINED")
d_joined_real = d_invariant(m_joined, m_real)
print(f"Joined trees baseline: d(joined, real) = {d_joined_real:.4f}")


# ---------------------------------------------------------------------------
# Generate modular random graphs via stochastic block model
# ---------------------------------------------------------------------------
def modular_random(N, M, n_blocks, p_in_to_p_out=4.0, seed=0):
    """Stochastic block model with N nodes, n_blocks communities,
    target M edges, and within-block density p_in_to_p_out × out-block density."""
    rng = np.random.default_rng(seed)
    sizes = [N // n_blocks + (1 if i < N % n_blocks else 0) for i in range(n_blocks)]

    # Solve for p_out given target edge count
    # Expected edges ≈ p_in*sum(C(s,2)) + p_out*sum(s_i*s_j for i<j)
    # With p_in = ratio * p_out:
    in_terms = sum(s*(s-1)/2 for s in sizes)
    out_terms = sum(sizes[i]*sizes[j] for i in range(n_blocks) for j in range(i+1, n_blocks))
    if out_terms == 0:
        return None
    p_out = M / (p_in_to_p_out * in_terms + out_terms)
    p_in = min(p_in_to_p_out * p_out, 1.0)
    p_out = min(p_out, 1.0)

    P = np.full((n_blocks, n_blocks), p_out)
    np.fill_diagonal(P, p_in)
    g = nx.stochastic_block_model(sizes, P, seed=int(rng.integers(1e9)))
    g = nx.Graph(g)
    g.remove_edges_from(nx.selfloop_edges(g))
    return g


print("\nGenerating modular random graphs (SBM)...")
n_per_block_count = 100
results = {}
for n_blocks in (2, 3, 4, 6, 8):
    print(f"\n  n_blocks = {n_blocks}:")
    distances = []
    for seed in range(n_per_block_count):
        g = modular_random(N, M, n_blocks, p_in_to_p_out=4.0, seed=seed)
        if g is None or not nx.is_connected(g):
            continue
        try:
            m = all_metrics(g, f"sbm_{n_blocks}_{seed}")
            d = d_invariant(m, m_real)
            distances.append(d)
        except Exception as e:
            continue
    if distances:
        arr = np.array(distances)
        p = float((arr <= d_joined_real).mean())
        results[n_blocks] = {
            "n": len(distances),
            "mean_d": float(arr.mean()),
            "std_d":  float(arr.std()),
            "min_d":  float(arr.min()),
            "max_d":  float(arr.max()),
            "p_joined_at_least_as_close": p,
        }
        print(f"    n={len(distances)} graphs  d_mean={arr.mean():.4f}  "
              f"d_std={arr.std():.4f}  d_min={arr.min():.4f}  "
              f"p={p:.4f}")
    else:
        print(f"    no connected graphs generated")

# ---------------------------------------------------------------------------
# Most-brain-like-modular: the BEST realization across all SBM trials
# ---------------------------------------------------------------------------
print("\n=== KEY COMPARISON ===")
all_dists = [d for r in results.values() for d in [r["mean_d"]] * r["n"]]
all_min_dists = [r["min_d"] for r in results.values()]
print(f"Joined-trees:           d = {d_joined_real:.4f}")
for n_blocks, r in results.items():
    print(f"SBM n_blocks={n_blocks}:        mean d = {r['mean_d']:.4f}  "
          f"min d = {r['min_d']:.4f}  p(joined<=null) = {r['p_joined_at_least_as_close']:.4f}")

# ---------------------------------------------------------------------------
# Also: tighter modular nulls — graphs with 4 blocks AND high
# within-block density (matching brain hemispheric structure)
# ---------------------------------------------------------------------------
print("\nTight modular null: 4 blocks, p_in/p_out = 6.0 (very modular)...")
distances_tight = []
for seed in range(200):
    g = modular_random(N, M, 4, p_in_to_p_out=6.0, seed=seed)
    if g is None or not nx.is_connected(g): continue
    m = all_metrics(g, f"tight_{seed}")
    distances_tight.append(d_invariant(m, m_real))
if distances_tight:
    arr = np.array(distances_tight)
    p = float((arr <= d_joined_real).mean())
    print(f"  n={len(distances_tight)} graphs  d_mean={arr.mean():.4f}  "
          f"d_min={arr.min():.4f}  p={p:.4f}")


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------
print("\nGenerating fig8_modular_null.png...")
fig, ax = plt.subplots(figsize=(10, 6))
all_data = []
labels = []
for n_blocks in sorted(results.keys()):
    distances = []
    for seed in range(n_per_block_count):
        g = modular_random(N, M, n_blocks, p_in_to_p_out=4.0, seed=seed)
        if g is None or not nx.is_connected(g): continue
        try:
            m = all_metrics(g, f"sbm_{n_blocks}_{seed}")
            distances.append(d_invariant(m, m_real))
        except Exception:
            continue
    all_data.append(distances)
    labels.append(f"SBM\n{n_blocks} blocks")
all_data.append(distances_tight)
labels.append("Tight\n4-block SBM")

bp = ax.boxplot(all_data, labels=labels, patch_artist=True,
                boxprops=dict(facecolor="#bbcce0"))
ax.axhline(d_joined_real, color="red", lw=2,
           label=f"joined-trees d = {d_joined_real:.3f}")
ax.set_ylabel("invariant distance to real brain")
ax.set_title("Critical test: do random MODULAR graphs match real brain\n"
             "as well as the joined-trees graph does?")
ax.legend()
plt.tight_layout()
plt.savefig("figures/fig8_modular_null.png", dpi=160, bbox_inches="tight")
print("Saved figures/fig8_modular_null.png")

# Save
with open("data/modular_null_results.json", "w") as f:
    json.dump({
        "joined_trees_distance": d_joined_real,
        "by_n_blocks": results,
        "tight_4block_p": float((np.array(distances_tight) <= d_joined_real).mean())
                          if distances_tight else None,
    }, f, indent=2,
       default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print("\nDONE.")
