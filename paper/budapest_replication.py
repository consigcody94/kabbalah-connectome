"""
Replicate the joined-trees vs brain finding across all 9 Budapest
Reference Connectome variants:
  - all/female/male  ×  20k/200k/1M fiber thresholds
"""
import csv
import io
import json
import os
import sys
import time
import urllib.request
import zipfile
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import joined_trees
from metrics import all_metrics, d_invariant
from nulls import all_nulls


VARIANTS = [
    "all_20k", "all_200k", "all_1m",
    "female_20k", "female_200k", "female_1m",
    "male_20k", "male_200k", "male_1m",
]

CACHE = "data/budapest_variants"
os.makedirs(CACHE, exist_ok=True)


def fetch_variant(name):
    cache = os.path.join(CACHE, f"{name}.zip")
    if not os.path.exists(cache):
        url = f"https://networks.skewed.de/net/budapest_connectome/files/{name}.csv.zip"
        try:
            urllib.request.urlretrieve(url, cache)
            print(f"  fetched {name}: {os.path.getsize(cache)} bytes")
        except Exception as e:
            print(f"  fetch failed for {name}: {e}")
            return None, None
    # Parse
    node_info = {}
    edges = []
    try:
        with zipfile.ZipFile(cache) as z:
            with z.open("nodes.csv") as f:
                text = f.read().decode("utf-8", errors="replace")
                for line in text.split("\n"):
                    if line.startswith("#") or not line: continue
                    parts = line.split(",")
                    try:
                        idx = int(parts[0])
                    except (ValueError, TypeError):
                        continue
                    fsname = parts[1].strip().strip('"') if len(parts) > 1 else ""
                    hemi = parts[2].strip().strip('"') if len(parts) > 2 else ""
                    node_info[idx] = {"fsname": fsname, "hemi": hemi}
            with z.open("edges.csv") as f:
                text = f.read().decode("utf-8", errors="replace")
                for line in text.split("\n"):
                    if line.startswith("#") or not line: continue
                    parts = line.split(",")
                    try:
                        s, t = int(parts[0]), int(parts[1])
                        occ = int(parts[9]) if len(parts) > 9 else 1
                    except (ValueError, TypeError, IndexError):
                        continue
                    edges.append((s, t, occ))
    except Exception as e:
        print(f"  parse error for {name}: {e}")
        return None, None
    return node_info, edges


def coarsen(node_info, edges, occ_threshold=100):
    G = nx.Graph()
    G.add_nodes_from(node_info.keys())
    for s, t, occ in edges:
        if occ < occ_threshold or s == t: continue
        G.add_edge(s, t)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])

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
                  for n in G.nodes()}
    Gc = nx.Graph()
    weights = defaultdict(int)
    for u, v in G.edges():
        if coarse_map[u] != coarse_map[v]:
            weights[tuple(sorted([coarse_map[u], coarse_map[v]]))] += 1
    for (a, b), w in weights.items():
        if w >= 5: Gc.add_edge(a, b)
    return Gc


# -----------------------------------------------------------------
J = joined_trees(shared_daath=True)
m_joined = all_metrics(J, "joined")
print(f"Joined-trees: N={J.number_of_nodes()} E={J.number_of_edges()}")

print("\nFetching + testing all 9 variants...")
results = {}
for v in VARIANTS:
    print(f"\n--- {v} ---")
    info, edges_raw = fetch_variant(v)
    if info is None: continue
    print(f"  full N={len(info)} edges_raw={len(edges_raw)}")
    Gc = coarsen(info, edges_raw, occ_threshold=100)
    if Gc.number_of_nodes() < 5 or not nx.is_connected(Gc):
        print(f"  coarse N={Gc.number_of_nodes()} disconnected — skip")
        continue
    m = all_metrics(Gc, v)
    d = d_invariant(m_joined, m)
    print(f"  coarse N={Gc.number_of_nodes()} E={Gc.number_of_edges()}  "
          f"d(joined,real) = {d:.4f}")

    # Null comparison
    nulls = all_nulls(Gc, k_per=100)
    pvals = {}
    for nname, nlist in nulls.items():
        ds = []
        for ng in nlist:
            try:
                ds.append(d_invariant(all_metrics(ng, "n"), m))
            except Exception:
                continue
        if ds:
            pvals[nname] = float((np.array(ds) <= d).mean())
    print(f"  p-values: " + " ".join(f"{k}={p:.3f}" for k, p in pvals.items()))
    results[v] = {
        "coarse_N": Gc.number_of_nodes(),
        "coarse_E": Gc.number_of_edges(),
        "d_joined": d,
        "p_values": pvals,
    }
    time.sleep(0.5)

# -----------------------------------------------------------------
print("\n" + "=" * 70)
print("REPLICATION SUMMARY (9 Budapest variants)")
print("=" * 70)
print(f"\n{'Variant':<14} {'N_c':>4} {'E_c':>4} {'d':>7} {'pER':>6} "
      f"{'pCFG':>6} {'pWS':>6} {'pBA':>6} {'pGEO':>6}")
print("-" * 70)
for v in VARIANTS:
    if v not in results: continue
    r = results[v]
    p = r["p_values"]
    print(f"{v:<14} {r['coarse_N']:>4} {r['coarse_E']:>4} "
          f"{r['d_joined']:>7.4f} {p.get('ER', float('nan')):>6.3f} "
          f"{p.get('CFG', float('nan')):>6.3f} {p.get('WS', float('nan')):>6.3f} "
          f"{p.get('BA', float('nan')):>6.3f} {p.get('GEO', float('nan')):>6.3f}")

# Aggregate
ds = [r["d_joined"] for r in results.values()]
ps_er = [r["p_values"].get("ER", 1.0) for r in results.values()]
n_sig = sum(1 for p in ps_er if p < 0.05)
print(f"\n  Variants tested: {len(results)}")
print(f"  d range: [{min(ds):.4f}, {max(ds):.4f}]  median {float(np.median(ds)):.4f}")
print(f"  Variants with p_ER < 0.05: {n_sig}/{len(results)}")
print(f"  Variants with p_ER = 0.000: {sum(1 for p in ps_er if p == 0)}/{len(results)}")

with open("data/budapest_replication_results.json", "w") as f:
    json.dump(results, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/budapest_replication_results.json")
