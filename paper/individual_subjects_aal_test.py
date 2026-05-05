"""
Per-subject replication test, properly coarsened using AAL atlas labels.

The previous individual_subjects_test.py compared 21-node joined-trees to
116-node full AAL connectomes — biased on size. This version uses the
canonical Tzourio-Mazoyer 2002 AAL atlas to coarsen each subject to ~21
nodes (lobe-level), making the comparison fair.

For each subject:
  1. Load full AAL connectome
  2. Coarsen to ~21 nodes via AAL_TO_LOBE map
  3. Compute distance to joined-trees graph
  4. Compute distance to ER null graphs of same size
  5. Report whether joined-trees beats random per subject

This is the proper version of the §3.10 caveat.
"""
import json
import os
import sys
import time
import urllib.request
import zipfile

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import joined_trees
from metrics import all_metrics, d_invariant
from nulls import er_ensemble
from aal_atlas import coarsen_aal_graph


def fetch_subject(name, cache_dir="data/aal_cache"):
    os.makedirs(cache_dir, exist_ok=True)
    cache = os.path.join(cache_dir, f"{name}.zip")
    if not os.path.exists(cache):
        url = f"https://networks.skewed.de/net/human_brains/files/{name}.csv.zip"
        try:
            urllib.request.urlretrieve(url, cache)
        except Exception:
            return None
    try:
        G = nx.Graph()
        with zipfile.ZipFile(cache) as z:
            with z.open("edges.csv") as f:
                text = f.read().decode("utf-8", errors="replace")
                for line in text.split("\n"):
                    if not line or line.startswith("#"): continue
                    parts = line.split(",")
                    try:
                        s, t = int(parts[0]), int(parts[1])
                    except (ValueError, TypeError, IndexError):
                        continue
                    if s != t: G.add_edge(s, t)
        if G.number_of_nodes() == 0: return None
        # Convert 0-based to 1-based to match AAL ordering
        relabel = {n: n + 1 for n in G.nodes()}
        G = nx.relabel_nodes(G, relabel)
        comps = sorted(nx.connected_components(G), key=len, reverse=True)
        return G.subgraph(comps[0]).copy()
    except Exception:
        return None


# Build joined-trees baseline
print("Joined-trees baseline...")
J = joined_trees(shared_daath=True)
m_joined = all_metrics(J, "joined")
print(f"  N={J.number_of_nodes()} E={J.number_of_edges()}")

# Get same 50 subjects as before (cached)
print("\nFetching list of subjects...")
nets_meta = json.loads(urllib.request.urlopen(
    "https://networks.skewed.de/api/net/human_brains", timeout=30
).read().decode())
all_aal = [n for n in nets_meta.get("nets", []) if n.endswith("_AAL")]
np.random.seed(42)
sample = list(np.random.choice(all_aal, size=min(50, len(all_aal)), replace=False))

results = []
sizes = []
for i, name in enumerate(sample, 1):
    G_full = fetch_subject(name)
    if G_full is None or G_full.number_of_nodes() < 50:
        continue
    Gc = coarsen_aal_graph(G_full)
    if Gc.number_of_nodes() < 10 or not nx.is_connected(Gc):
        print(f"[{i}/{len(sample)}] {name[:40]:<40} coarse N={Gc.number_of_nodes()} "
              f"connected={nx.is_connected(Gc)} -- skip")
        continue
    sizes.append((Gc.number_of_nodes(), Gc.number_of_edges()))
    try:
        m = all_metrics(Gc, name)
    except Exception as e:
        print(f"[{i}/{len(sample)}] {name[:40]:<40} metric err -- skip")
        continue
    d_j = d_invariant(m_joined, m)
    nulls = er_ensemble(Gc.number_of_nodes(), Gc.number_of_edges(), k=30, seed=i)
    if not nulls: continue
    nm_dists = [d_invariant(all_metrics(g, "n"), m) for g in nulls]
    d_n = float(np.mean(nm_dists))
    p = float((np.array(nm_dists) <= d_j).mean())
    results.append({
        "subject": name,
        "N_full": G_full.number_of_nodes(),
        "E_full": G_full.number_of_edges(),
        "N_coarse": Gc.number_of_nodes(),
        "E_coarse": Gc.number_of_edges(),
        "d_joined_to_subject": d_j,
        "d_random_to_subject_mean": d_n,
        "p_joined_at_least_as_close": p,
        "joined_beats_random": d_j < d_n,
    })
    flag = "✓" if d_j < d_n else "✗"
    print(f"[{i:2d}/{len(sample)}] {name[:40]:<40} "
          f"coarse=({Gc.number_of_nodes()},{Gc.number_of_edges()})  "
          f"d_J={d_j:.4f}  d_R={d_n:.4f}  p={p:.4f}  {flag}")

print("\n" + "=" * 80)
print("PER-SUBJECT RESULT (AAL-COARSENED, FAIR COMPARISON)")
print("=" * 80)
if results:
    djs = np.array([r["d_joined_to_subject"] for r in results])
    nms = np.array([r["d_random_to_subject_mean"] for r in results])
    ps  = np.array([r["p_joined_at_least_as_close"] for r in results])
    n_beats = sum(1 for r in results if r["joined_beats_random"])
    print(f"  Subjects analyzed:   {len(results)}")
    print(f"  d_joined  mean ± std: {djs.mean():.4f} ± {djs.std():.4f}")
    print(f"  d_random  mean ± std: {nms.mean():.4f} ± {nms.std():.4f}")
    print(f"  joined < random:     {n_beats}/{len(results)} ({100*n_beats/len(results):.1f}%)")
    print(f"  p < 0.05:            {int(sum(ps < 0.05))}/{len(results)}")
    print(f"  p = 0.000:           {int(sum(ps == 0))}/{len(results)}")
    print(f"  median p:            {float(np.median(ps)):.4f}")

with open("data/individual_subjects_aal_results.json", "w") as f:
    json.dump(results, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/individual_subjects_aal_results.json")
