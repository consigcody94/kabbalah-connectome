"""
Subgraph matching test: is the joined-trees graph literally embedded
as a subgraph of the real human brain?

Two questions:
  1. Is joined-trees an EXACT subgraph isomorphism of the real brain?
     (Almost certainly no for fine parcellation, but possible for coarse.)
  2. What's the MAXIMUM COMMON SUBGRAPH between joined-trees and brain?
     I.e., how much of the Tree of Life topology actually appears in brain?

For (1) we use VF2 algorithm. For (2) we use a greedy MCS approximation
since exact MCS is NP-hard.

Also: how does the count of matching subgraphs compare to random graphs
of same size?
"""
import csv
import json
import os
import sys
from collections import defaultdict
from itertools import combinations

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx
from networkx.algorithms import isomorphism

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import joined_trees, tree_of_life
from nulls import er_ensemble


# ---------------------------------------------------------------------------
# Load real coarse brain
# ---------------------------------------------------------------------------
print("Loading + coarsening Budapest connectome...")
node_info = {}
with open("data/budapest/nodes.csv", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split(",")
        if len(parts) < 4: continue
        try: idx = int(parts[0])
        except (ValueError, TypeError): continue
        node_info[idx] = {"fsname": parts[1], "hemi": parts[2]}
G_full = nx.Graph()
G_full.add_nodes_from(node_info.keys())
with open("data/budapest/edges.csv", encoding="utf-8") as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split(",")
        try:
            s, t = int(parts[0]), int(parts[1]); occ = int(parts[9])
        except (ValueError, TypeError, IndexError):
            continue
        if occ < 100 or s == t: continue
        G_full.add_edge(s, t)
G_full.remove_nodes_from([n for n in list(G_full.nodes()) if G_full.degree(n) == 0])

def lobe_of(name):
    s = name.lower()
    if "frontal" in s or "parsop" in s or "parstri" in s or "parsorb" in s or "rostral" in s: return "Frontal"
    if "precentral" in s or "paracentral" in s: return "Motor"
    if "postcentral" in s: return "Somatosensory"
    if "parietal" in s or "supramarg" in s or "precune" in s: return "Parietal"
    if "transversetemporal" in s: return "Auditory"
    if "temporal" in s or "fusiform" in s or "banks" in s or "entorhinal" in s or "parahippo" in s: return "Temporal"
    if "occipital" in s or "lingual" in s or "cuneus" in s or "pericalc" in s: return "Occipital"
    if "hippo" in s: return "Hippocampus"
    if "thalam" in s: return "Thalamus"
    if "cingul" in s or "isthmus" in s: return "Cingulate"
    if "insula" in s: return "Insula"
    return "Subcortical"

coarse_map = {n: f"{node_info[n]['hemi'] or 'none'}-{lobe_of(node_info[n]['fsname'])}"
              for n in G_full.nodes()}
G_brain = nx.Graph()
weights = defaultdict(int)
for u, v in G_full.edges():
    if coarse_map[u] != coarse_map[v]:
        weights[tuple(sorted([coarse_map[u], coarse_map[v]]))] += 1
for (a, b), w in weights.items():
    if w >= 5: G_brain.add_edge(a, b)
print(f"  Brain coarse: N={G_brain.number_of_nodes()} E={G_brain.number_of_edges()}")

# Joined-trees baseline
J = joined_trees(shared_daath=True)
print(f"  Joined trees: N={J.number_of_nodes()} E={J.number_of_edges()}")


# ---------------------------------------------------------------------------
# (1) Exact subgraph isomorphism: is J a subgraph of brain?
# ---------------------------------------------------------------------------
print("\n--- Test 1: Is joined-trees an EXACT subgraph of real brain? ---")
gm = isomorphism.GraphMatcher(G_brain, J)
exact_match = gm.subgraph_is_isomorphic()
print(f"  Subgraph isomorphism (J in Brain): {exact_match}")
if exact_match:
    print("  YES — Tree of Life literally appears as a subgraph!")
else:
    print("  NO — Tree of Life is NOT an exact subgraph at coarse scale.")
    print("  (Expected — would require N(brain) >= N(J) AND specific edge pattern.)")


# ---------------------------------------------------------------------------
# (2) Maximum Common Edge Subgraph (greedy approximation)
# ---------------------------------------------------------------------------
def mces_size(G1, G2, n_trials=200, seed=0):
    """Approximate MCES (Maximum Common Edge Subgraph) via greedy alignment.

    Tries n_trials random node mappings and reports the best alignment
    (most edges that are present in BOTH the mapped subgraph and the target).
    """
    rng = np.random.default_rng(seed)
    G1_nodes = list(G1.nodes())
    G2_nodes = list(G2.nodes())
    if len(G1_nodes) > len(G2_nodes):
        G1, G2 = G2, G1
        G1_nodes, G2_nodes = G2_nodes, G1_nodes
    best = 0
    best_mapping = None
    for trial in range(n_trials):
        sample = list(rng.choice(G2_nodes, len(G1_nodes), replace=False))
        # Map G1 -> sample
        mapping = dict(zip(G1_nodes, sample))
        # Count common edges
        match = 0
        for u, v in G1.edges():
            if G2.has_edge(mapping[u], mapping[v]):
                match += 1
        if match > best:
            best = match
            best_mapping = mapping
    return best, best_mapping


print("\n--- Test 2: Maximum common-edge subgraph (greedy, 500 trials) ---")
mces, mapping = mces_size(J, G_brain, n_trials=500, seed=0)
print(f"  Joined-trees has {J.number_of_edges()} edges total.")
print(f"  Best alignment recovers {mces} edges in real brain "
      f"({100*mces/J.number_of_edges():.1f}% of TOL edges).")

# Compare to random graphs
print("\n--- Test 3: Same MCES test for random graphs of same size as J ---")
random_mces = []
for seed in range(50):
    Rg = nx.gnm_random_graph(J.number_of_nodes(), J.number_of_edges(), seed=seed)
    if not nx.is_connected(Rg): continue
    m, _ = mces_size(Rg, G_brain, n_trials=200, seed=seed)
    random_mces.append(m)

if random_mces:
    arr = np.array(random_mces)
    p = float((arr >= mces).mean())
    print(f"  Random-graph MCES vs brain: mean={arr.mean():.1f} ± {arr.std():.1f}")
    print(f"  Joined-trees MCES = {mces}; p(random >= joined) = {p:.4f}")
    if p < 0.05:
        print(f"  -> SIGNIFICANT: joined-trees recovers more brain edges than {(1-p)*100:.0f}% of random graphs")


# ---------------------------------------------------------------------------
# (4) Show the best alignment — which sephirot map to which brain regions?
# ---------------------------------------------------------------------------
print("\n--- Test 4: Best sephirot-to-brain-region alignment from MCES ---")
print(f"  (the alignment that recovers the most TOL edges in real brain)")
if mapping:
    from graphs import SEPHIROT_10, QLI_NAMES
    # Show TOL alignments first
    print("\n  Tree of Life:")
    for s in (SEPHIROT_10 + ["Daath"]):
        if s in mapping:
            print(f"    {s:<12s} -> {mapping[s]}")
    print("\n  Tree of Death:")
    for q in QLI_NAMES:
        if q in mapping:
            print(f"    {q:<15s} -> {mapping[q]}")


# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
with open("data/subgraph_match_results.json", "w") as f:
    json.dump({
        "exact_subgraph_isomorphism": exact_match,
        "joined_trees_total_edges": J.number_of_edges(),
        "best_mces_in_brain": mces,
        "best_mces_pct": 100 * mces / J.number_of_edges(),
        "random_mces_mean": float(np.mean(random_mces)) if random_mces else None,
        "random_mces_std": float(np.std(random_mces)) if random_mces else None,
        "p_random_at_least_as_good": float((np.array(random_mces) >= mces).mean())
                                      if random_mces else None,
        "best_mapping": {k: str(v) for k, v in mapping.items()}
                        if mapping else None,
    }, f, indent=2,
       default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/subgraph_match_results.json")
