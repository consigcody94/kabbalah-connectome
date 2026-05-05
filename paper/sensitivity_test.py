"""
Sensitivity test: does the joined-trees-closer-than-random result depend
on the connectome threshold? And: does a directed (Lightning Flash) Tree
variant match better or worse?
"""
import csv
import json
import os
import sys
from collections import defaultdict

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import SEPHIROT_10, QLI_NAMES, joined_trees, TOL_PATHS_22, DAATH_PATHS
from metrics import all_metrics, d_invariant, INVARIANT_KEYS
from nulls import all_nulls

# ---------------------------------------------------------------------------
# Reusable: load + coarsen connectome at a given threshold
# ---------------------------------------------------------------------------
def load_connectome(occ_threshold):
    node_info = {}
    with open("data/budapest/nodes.csv", encoding="utf-8") as f:
        reader = csv.DictReader(
            (line for line in f if not line.startswith("#")),
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
            (line for line in f if not line.startswith("#")),
            fieldnames=["source", "target", "fcm", "flm", "fam",
                        "ecm", "fcmd", "flmd", "famd", "occ"])
        for r in reader:
            try:
                occ = int(r["occ"])
            except (ValueError, TypeError):
                continue
            if occ < occ_threshold:
                continue
            s, t = int(r["source"]), int(r["target"])
            if s != t:
                G.add_edge(s, t)
    G.remove_nodes_from([n for n in list(G.nodes()) if G.degree(n) == 0])
    return G, node_info


def lobe_of(fsname):
    s = fsname.lower()
    if "frontal" in s or "parsop" in s or "parstri" in s or "parsorb" in s \
       or "rostral" in s or "caudal" in s and "anterior" in s:
        return "Frontal"
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
    return "Subcortical"   # caudate/putamen/pallidum/etc + brain stem


def coarsen(G, node_info):
    coarse_map = {}
    for n in G.nodes():
        info = node_info[n]
        h = info["hemi"] or "none"
        coarse_map[n] = f"{h}-{lobe_of(info['fsname'])}"
    Gc = nx.Graph()
    weights = defaultdict(int)
    for u, v in G.edges():
        cu, cv = coarse_map[u], coarse_map[v]
        if cu != cv:
            weights[tuple(sorted([cu, cv]))] += 1
    for (a, b), w in weights.items():
        if w >= 5:
            Gc.add_edge(a, b)
    return Gc


# ---------------------------------------------------------------------------
# 1. Threshold sensitivity sweep
# ---------------------------------------------------------------------------
print("=" * 78)
print("SENSITIVITY: joined-trees -> real-connectome distance vs threshold")
print("=" * 78)

joined = joined_trees(shared_daath=True)
m_joined = all_metrics(joined, "JOINED_node")

results_threshold = {}
for thresh in (50, 100, 150, 200, 250, 350):
    G, info = load_connectome(thresh)
    Gc = coarsen(G, info)
    if Gc.number_of_nodes() < 5 or not nx.is_connected(Gc):
        print(f"  thresh={thresh}: skipped (N={Gc.number_of_nodes()}, "
              f"connected={nx.is_connected(Gc) if Gc.number_of_nodes() else False})")
        continue
    m_real = all_metrics(Gc, f"REAL_t{thresh}")
    d_real = d_invariant(m_joined, m_real)
    nulls = all_nulls(Gc, k_per=100)
    null_dists = []
    for null_name, ng_list in nulls.items():
        for ng in ng_list:
            nm = all_metrics(ng, "rand")
            null_dists.append((null_name, d_invariant(nm, m_real)))
    p_per_null = {}
    for nn in ("ER", "CFG", "WS", "BA", "GEO"):
        ds = [d for n, d in null_dists if n == nn]
        if ds:
            p_per_null[nn] = float((np.array(ds) <= d_real).mean())
    results_threshold[thresh] = {
        "Nfull": G.number_of_nodes(), "Efull": G.number_of_edges(),
        "Ncoarse": Gc.number_of_nodes(), "Ecoarse": Gc.number_of_edges(),
        "d_joined_to_real": d_real,
        "p_values": p_per_null,
    }
    print(f"  thresh={thresh:3d}  full=({G.number_of_nodes()},{G.number_of_edges()})  "
          f"coarse=({Gc.number_of_nodes()},{Gc.number_of_edges()})  "
          f"d={d_real:.4f}  pER={p_per_null.get('ER','-'):.3f}  "
          f"pCFG={p_per_null.get('CFG','-'):.3f}  "
          f"pWS={p_per_null.get('WS','-'):.3f}")

# ---------------------------------------------------------------------------
# 2. Directed Lightning Flash variant
# ---------------------------------------------------------------------------
print("\n" + "=" * 78)
print("DIRECTED LIGHTNING FLASH variant (Tree as DAG, top->bottom flow)")
print("=" * 78)

# The Lightning Flash is the canonical descent order:
# Keter -> Chokhmah -> Binah -> Daath -> Chesed -> Geburah -> Tiferet
#       -> Netzach -> Hod -> Yesod -> Malkuth
LF_ORDER = ["Keter", "Chokhmah", "Binah", "Daath", "Chesed", "Geburah",
            "Tiferet", "Netzach", "Hod", "Yesod", "Malkuth"]
LF_RANK = {n: i for i, n in enumerate(LF_ORDER)}

# Build directed Tree of Life: each undirected edge oriented down-the-flash
# (lower-rank -> higher-rank). Daath included.
all_edges = list(TOL_PATHS_22) + list(DAATH_PATHS)
DTOL = nx.DiGraph()
DTOL.add_nodes_from(LF_ORDER)
for a, b in all_edges:
    if LF_RANK[a] < LF_RANK[b]:
        DTOL.add_edge(a, b)
    else:
        DTOL.add_edge(b, a)

# Same for qliphoth
DQLI = nx.DiGraph()
qmap = dict(zip(SEPHIROT_10, QLI_NAMES))
qmap["Daath"] = "Daath"   # shared
DQLI.add_nodes_from([qmap[n] for n in LF_ORDER])
for a, b in all_edges:
    if LF_RANK[a] < LF_RANK[b]:
        DQLI.add_edge(qmap[a], qmap[b])
    else:
        DQLI.add_edge(qmap[b], qmap[a])

DJOINED = nx.compose(DTOL, DQLI)

print(f"  Directed joined trees: N={DJOINED.number_of_nodes()} "
      f"E={DJOINED.number_of_edges()}")
print(f"  Is DAG: {nx.is_directed_acyclic_graph(DJOINED)}")
print(f"  In-degree distribution: {sorted([d for _, d in DJOINED.in_degree()])}")
print(f"  Out-degree distribution: {sorted([d for _, d in DJOINED.out_degree()])}")

# Daath-specific metrics in the DAG
if "Daath" in DJOINED:
    print(f"\n  Daath in directed graph:")
    print(f"    in-degree:  {DJOINED.in_degree('Daath')}")
    print(f"    out-degree: {DJOINED.out_degree('Daath')}")
    print(f"    pagerank rank: ",end="")
    pr = nx.pagerank(DJOINED)
    rank = 1 + sum(1 for v in pr.values() if v > pr["Daath"])
    print(f"{rank} of {DJOINED.number_of_nodes()}  (pr={pr['Daath']:.4f})")

# Compare to a directed analog of the brain — just the undirected brain with
# arbitrary orientation. (Brain connectomes are intrinsically directed via
# tract polarity, but bulk DTI doesn't recover this.)
# Skip for now; the directed Tree analysis is itself the new test.

# ---------------------------------------------------------------------------
# 3. Save
# ---------------------------------------------------------------------------
out = {
    "threshold_sensitivity": results_threshold,
    "directed_tree": {
        "N": DJOINED.number_of_nodes(),
        "E": DJOINED.number_of_edges(),
        "is_dag": nx.is_directed_acyclic_graph(DJOINED),
        "daath_in_degree": DJOINED.in_degree("Daath"),
        "daath_out_degree": DJOINED.out_degree("Daath"),
        "daath_pagerank_rank": rank,
    },
}
with open("data/sensitivity_results.json", "w") as f:
    json.dump(out, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/sensitivity_results.json")
