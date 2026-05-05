"""
Test the joined-trees vs brain finding across INDIVIDUAL subjects.

Individual subject node CSVs lack anatomical labels (only graphml_id),
so we test the FULL graph (no coarsening) against the joined-trees
graph and random-graph nulls. This is a cleaner test than coarsening:
no choice of parcellation is involved.

For each subject:
  - load full DTI graph (~116 AAL nodes)
  - compute invariant distance to joined-trees graph (21 nodes)
  - compute invariant distance to ER null graphs of same size as subject
  - report whether joined-trees is closer than typical random graph
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
        # Restrict to largest connected component
        if G.number_of_nodes() == 0: return None
        comps = sorted(nx.connected_components(G), key=len, reverse=True)
        return G.subgraph(comps[0]).copy()
    except Exception:
        return None


# ---------------------------------------------------------------------------
print("Computing joined-trees baseline...")
J = joined_trees(shared_daath=True)
m_joined = all_metrics(J, "joined")
print(f"  N={J.number_of_nodes()} E={J.number_of_edges()}")

print("\nFetching list of individual subjects...")
nets_meta = json.loads(urllib.request.urlopen(
    "https://networks.skewed.de/api/net/human_brains", timeout=30
).read().decode())
all_aal = [n for n in nets_meta.get("nets", []) if n.endswith("_AAL")]
np.random.seed(42)
sample = list(np.random.choice(all_aal, size=min(50, len(all_aal)), replace=False))
print(f"  Sampled {len(sample)} of {len(all_aal)} AAL connectomes")

results = []
for i, name in enumerate(sample, 1):
    print(f"\n[{i}/{len(sample)}] {name}")
    G = fetch_subject(name)
    if G is None:
        print(f"  fetch failed, skip")
        continue
    if G.number_of_nodes() < 50 or G.number_of_edges() < 100:
        print(f"  too small, skip")
        continue
    print(f"  N={G.number_of_nodes()} E={G.number_of_edges()}")
    try:
        m = all_metrics(G, name)
    except Exception as e:
        print(f"  metrics error: {e}")
        continue
    d_j = d_invariant(m_joined, m)
    nulls = er_ensemble(G.number_of_nodes(), G.number_of_edges(), k=20, seed=i)
    if not nulls:
        continue
    nm_dists = [d_invariant(all_metrics(g, "n"), m) for g in nulls]
    d_n = float(np.mean(nm_dists))
    p = float((np.array(nm_dists) <= d_j).mean())
    results.append({
        "subject": name,
        "N": G.number_of_nodes(),
        "E": G.number_of_edges(),
        "d_joined_to_subject": d_j,
        "d_random_to_subject_mean": d_n,
        "p_joined_at_least_as_close": p,
        "joined_beats_random": d_j < d_n,
    })
    print(f"  d_joined={d_j:.4f}  d_random_mean={d_n:.4f}  "
          f"p={p:.4f}  joined<random: {d_j < d_n}")
    time.sleep(1)

# ---------------------------------------------------------------------------
print("\n" + "=" * 75)
print("REPLICATION ACROSS INDIVIDUAL SUBJECTS")
print("=" * 75)
if results:
    djs = np.array([r["d_joined_to_subject"] for r in results])
    nms = np.array([r["d_random_to_subject_mean"] for r in results])
    ps  = np.array([r["p_joined_at_least_as_close"] for r in results])
    n_beats = sum(1 for r in results if r["joined_beats_random"])
    print(f"  Subjects analyzed:     {len(results)}")
    print(f"  d_joined  mean ± std:  {djs.mean():.4f} ± {djs.std():.4f}")
    print(f"  d_random  mean ± std:  {nms.mean():.4f} ± {nms.std():.4f}")
    print(f"  joined < random:       {n_beats}/{len(results)} subjects "
          f"({100*n_beats/len(results):.1f}%)")
    print(f"  p < 0.05:              {int(sum(ps < 0.05))}/{len(results)} subjects")
    print(f"  p = 0.000:             {int(sum(ps == 0))}/{len(results)} subjects")

with open("data/individual_subjects_results.json", "w") as f:
    json.dump(results, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/individual_subjects_results.json")
