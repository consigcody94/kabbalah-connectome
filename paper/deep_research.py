"""
Deeper-research package: three tests in one script.

  A. Cross-species — joined-trees vs macaque cortex (Young 1993, 47 nodes)
  B. Expanded mythology — Tarot, Cordoveran, Lurianic, Medicine Wheel, etc.
  C. Sensitivity — perturb joined-trees with random edge add/remove, see
     how brain-match survives
"""
import csv
import json
import os
import sys
from collections import defaultdict
from itertools import product

sys.stdout.reconfigure(encoding="utf-8")

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

THIS = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, THIS)

from graphs import (
    SEPHIROT_10, QLI_NAMES, joined_trees, tree_of_life,
    TOL_PATHS_22, DAATH_PATHS,
)
from metrics import all_metrics, d_invariant, INVARIANT_KEYS
from nulls import er_ensemble, all_nulls


# ===========================================================================
# Helper: load + coarsen Budapest brain
# ===========================================================================
def load_budapest_coarse():
    node_info = {}
    with open("data/budapest/nodes.csv", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split(",")
            if len(parts) < 4: continue
            try:
                idx = int(parts[0])
            except (ValueError, TypeError):
                continue
            node_info[idx] = {"fsname": parts[1], "hemi": parts[2]}
    G = nx.Graph()
    G.add_nodes_from(node_info.keys())
    with open("data/budapest/edges.csv", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split(",")
            try:
                s, t = int(parts[0]), int(parts[1]); occ = int(parts[9])
            except (ValueError, TypeError, IndexError):
                continue
            if occ < 100 or s == t: continue
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


def load_macaque():
    """Load Young 1993 macaque cortex network."""
    G = nx.Graph()
    with open("data/macaque_neural/edges.csv", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"): continue
            parts = line.strip().split(",")
            try:
                s, t = int(parts[0]), int(parts[1])
                if s != t: G.add_edge(s, t)
            except (ValueError, TypeError, IndexError):
                continue
    return G


# ===========================================================================
# A. CROSS-SPECIES — joined-trees vs macaque cortex
# ===========================================================================
print("=" * 78)
print("A. CROSS-SPECIES: joined-trees vs Young 1993 macaque cortex")
print("=" * 78)

Mac = load_macaque()
print(f"\nMacaque cortex: N={Mac.number_of_nodes()} E={Mac.number_of_edges()} "
      f"connected={nx.is_connected(Mac)}")

J = joined_trees(shared_daath=True)
m_macaque = all_metrics(Mac, "macaque")
m_joined = all_metrics(J, "joined")

d_jm = d_invariant(m_joined, m_macaque)
print(f"\nd(joined-trees, macaque) = {d_jm:.4f}")

# Null distribution
print("\nGenerating null distribution against macaque...")
nulls = all_nulls(Mac, k_per=200)
null_dists = {n: [] for n in nulls}
for null_name, ng_list in nulls.items():
    for ng in ng_list:
        try:
            nm = all_metrics(ng, "n")
            null_dists[null_name].append(d_invariant(nm, m_macaque))
        except Exception:
            continue
for n, ds in null_dists.items():
    if not ds: continue
    arr = np.array(ds)
    p = float((arr <= d_jm).mean())
    print(f"  {n}: mean d={arr.mean():.4f}  std={arr.std():.4f}  p={p:.4f}")

cross_species = {
    "macaque_size": (Mac.number_of_nodes(), Mac.number_of_edges()),
    "d_joined_to_macaque": d_jm,
    "p_values": {n: float((np.array(ds) <= d_jm).mean()) if ds else None
                 for n, ds in null_dists.items()},
}


# ===========================================================================
# B. EXPANDED MYTHOLOGY — more structures
# ===========================================================================
print("\n" + "=" * 78)
print("B. EXPANDED MYTHOLOGY: more world structures vs human brain")
print("=" * 78)

def build_tarot_22():
    """Tarot Major Arcana (22 cards) connected per Hermetic Kabbalistic
    correspondences: each card maps to one of the 22 paths between
    sephirot in the Tree of Life. Tarot adjacency = Tree of Life path
    adjacency.

    The 22 Major Arcana: Fool, Magician, High Priestess, Empress, Emperor,
    Hierophant, Lovers, Chariot, Strength, Hermit, Wheel of Fortune,
    Justice, Hanged Man, Death, Temperance, Devil, Tower, Star, Moon,
    Sun, Judgement, World.
    """
    arcana = ["Fool", "Magician", "HighPriestess", "Empress", "Emperor",
              "Hierophant", "Lovers", "Chariot", "Strength", "Hermit",
              "WheelFortune", "Justice", "HangedMan", "Death", "Temperance",
              "Devil", "Tower", "Star", "Moon", "Sun", "Judgement", "World"]
    # Map each Tarot card to a Tree of Life path (Crowley/Golden Dawn
    # standard attribution from Crowley's 777). Here we just need that the
    # 22 cards correspond 1-to-1 with the 22 paths, so the GRAPH structure
    # of the Tarot-as-state-space is the line graph of the 22-path TOL.
    # That is: nodes = paths; edges = paths sharing a sephira endpoint.
    g = nx.Graph()
    g.add_nodes_from(arcana)
    # Pair up arcana with TOL paths in Crowley order, then connect arcana
    # whose corresponding paths share a sephira.
    for i, p1 in enumerate(TOL_PATHS_22):
        for j, p2 in enumerate(TOL_PATHS_22[i+1:], i+1):
            if set(p1) & set(p2):
                g.add_edge(arcana[i], arcana[j])
    return g

def build_cordoveran():
    """Moshe Cordovero's variant of the Tree of Life (Pardes Rimmonim, 1548).
    Cordovero retained the 10 sephirot but described 13 channels (not 22)
    based on the 13 attributes of mercy (Exodus 34:6-7). Different
    connectivity from the Kircher/Hermetic standard.

    Cordoveran channels (approximate, after Pardes Rimmonim):
    Keter→Chokhmah, Keter→Binah, Keter→Tiferet, Chokhmah→Chesed,
    Binah→Geburah, Chesed→Geburah, Chesed→Tiferet, Geburah→Tiferet,
    Tiferet→Yesod, Tiferet→Netzach, Tiferet→Hod, Yesod→Malkuth,
    Netzach→Malkuth (the 13th).
    """
    g = nx.Graph()
    g.add_nodes_from(SEPHIROT_10)
    edges = [
        ("Keter", "Chokhmah"), ("Keter", "Binah"), ("Keter", "Tiferet"),
        ("Chokhmah", "Chesed"), ("Binah", "Geburah"),
        ("Chesed", "Geburah"), ("Chesed", "Tiferet"), ("Geburah", "Tiferet"),
        ("Tiferet", "Yesod"), ("Tiferet", "Netzach"), ("Tiferet", "Hod"),
        ("Yesod", "Malkuth"), ("Netzach", "Malkuth"),
    ]
    g.add_edges_from(edges)
    return g

def build_lurianic_partzufim():
    """Lurianic Kabbalah (Isaac Luria, 16th c.): 5 Partzufim
    (Divine Personalities) replacing/parallel to the 10 sephirot.

    Five Partzufim:
      Arikh Anpin (Long Face) = Keter
      Abba (Father) = Chokhmah
      Imma (Mother) = Binah
      Zeir Anpin (Small Face) = Chesed-Geburah-Tiferet-Netzach-Hod-Yesod
      Nukva / Malkah = Malkuth

    Each partzuf is a node; connections per Lurianic doctrine:
    Arikh ↔ Abba, Arikh ↔ Imma, Abba ↔ Imma (the supernal triad)
    Abba → Zeir, Imma → Zeir (parents → child)
    Zeir ↔ Nukva (groom and bride / yichud)
    Plus inner unfoldings — each partzuf has internal sub-structure.

    Simplified: 5 nodes + the 5 partzufim each have 5 internal points
    (their own NHY of soul levels) = 30 nodes total.
    """
    g = nx.Graph()
    parts = ["Arikh", "Abba", "Imma", "Zeir", "Nukva"]
    g.add_nodes_from(parts)
    # External connections
    g.add_edges_from([
        ("Arikh", "Abba"), ("Arikh", "Imma"), ("Abba", "Imma"),
        ("Abba", "Zeir"), ("Imma", "Zeir"),
        ("Zeir", "Nukva"),
    ])
    # Each partzuf has 5 internal soul levels: Nefesh, Ruach, Neshamah,
    # Chayah, Yechidah. They're nested inside.
    levels = ["Nefesh", "Ruach", "Neshamah", "Chayah", "Yechidah"]
    for p in parts:
        for lv in levels:
            n = f"{p}-{lv}"
            g.add_node(n)
            g.add_edge(p, n)
        # Internal hierarchy
        for a, b in zip(levels[:-1], levels[1:]):
            g.add_edge(f"{p}-{a}", f"{p}-{b}")
    return g

def build_medicine_wheel():
    """Lakota / Plains Indigenous Medicine Wheel: 4 cardinal directions +
    center, on conceptual levels of Earth, Spirit, Sky.

    Lakota tradition: 4 directions (E, S, W, N) each associated with an
    element, season, and color; center represents Wakan Tanka (Great
    Spirit). Often expanded to 7-fold (4 dirs + above/below/within).
    """
    g = nx.Graph()
    dirs = ["East", "South", "West", "North"]
    levels = ["Above", "Below", "Within"]
    g.add_node("Center")
    for d in dirs: g.add_node(d)
    for l in levels: g.add_node(l)
    # Cardinal directions form a wheel
    for d1, d2 in zip(dirs, dirs[1:] + dirs[:1]):
        g.add_edge(d1, d2)
    # Center connects to all directions
    for d in dirs: g.add_edge("Center", d)
    # Above/Below/Within connect to Center (cosmic axis)
    for l in levels: g.add_edge("Center", l)
    return g

def build_buddhist_5dhyani():
    """Buddhist Five-Dhyani-Buddha mandala (Vajrayana / Tibetan).

    Center: Vairocana (white).
    Four directions: Akshobhya (E, blue), Ratnasambhava (S, yellow),
                     Amitabha (W, red), Amoghasiddhi (N, green).
    Each Buddha has a consort: Akashadhatvishvari (V), Locana (A),
    Mamaki (R), Pandara (Am), Tara (Amo).
    Plus 8 Bodhisattvas at corners and intermediate directions.

    Simplified: 5 Buddhas + 5 consorts = 10 nodes + 8 bodhisattvas = 18.
    Connections: center to 4 cardinals, each Buddha to its consort and
    flanking bodhisattvas.
    """
    g = nx.Graph()
    buddhas = ["Vairocana(C)", "Akshobhya(E)", "Ratnasambhava(S)",
               "Amitabha(W)", "Amoghasiddhi(N)"]
    g.add_nodes_from(buddhas)
    # Center → 4 directions
    for b in buddhas[1:]:
        g.add_edge(buddhas[0], b)
    # Adjacent directions connect (cardinal cycle)
    for b1, b2 in zip(buddhas[1:], buddhas[2:] + [buddhas[1]]):
        g.add_edge(b1, b2)
    # 5 consorts
    consorts = ["Akashadhatvishvari", "Locana", "Mamaki", "Pandara", "Tara"]
    for b, c in zip(buddhas, consorts):
        g.add_node(c)
        g.add_edge(b, c)
    # 8 bodhisattvas at intermediate positions
    bodhi = [f"Bodhi-{i}" for i in range(8)]
    for b in bodhi: g.add_node(b)
    # Each bodhisattva connects to two adjacent Buddhas
    for i, b in enumerate(bodhi):
        # Map bodhisattva i to two Buddhas (rough corners)
        b1 = buddhas[1 + (i // 2) % 4]
        b2 = buddhas[1 + ((i // 2) + 1) % 4]
        g.add_edge(b, b1)
        g.add_edge(b, b2)
    return g

def build_lightning_2trees_rev():
    """Two Lightning Flashes — descent (Keter→Malkuth) AND ascent
    (Malkuth→Keter), on both Tree of Life and Qliphoth.

    Adds the return path edges to the joined-trees graph. The Lightning
    Flash descent is along specific sephirot in order. The ascent traces
    the same sephirot in reverse.

    For testing purposes: same nodes as joined-trees, but with EXTRA
    direct edges Keter↔Malkuth and Thaumiel↔Lilith (the two extremes of
    each tree linked).
    """
    g = joined_trees(shared_daath=True)
    g.add_edge("Keter", "Malkuth")
    g.add_edge("Thaumiel", "Lilith")
    return g

# Earlier mythologies to include
from comparative_mythology import (build_iching, build_yggdrasil,
                                    build_chakra_nadi, build_sri_yantra,
                                    build_wacah_chan)

candidates = {
    "TreeOfLife (joined+Daath)":      joined_trees(shared_daath=True),
    "TreeOfLife + Lightning return":  build_lightning_2trees_rev(),
    "Cordoveran (13 channels)":       build_cordoveran(),
    "Lurianic Partzufim (5+inner)":   build_lurianic_partzufim(),
    "Tarot Major Arcana (22)":        build_tarot_22(),
    "Yggdrasil (9 worlds)":           build_yggdrasil(),
    "Buddhist 5-Buddha Mandala":      build_buddhist_5dhyani(),
    "Lakota Medicine Wheel":          build_medicine_wheel(),
    "Sri Yantra (9 triangles)":       build_sri_yantra(),
    "Mayan Wacah Chan":               build_wacah_chan(),
    "I Ching (6-cube)":               build_iching(),
    "Chakra-Nadi (Hindu)":            build_chakra_nadi(),
}

print("\nBuilt candidate structures:")
for name, g in candidates.items():
    print(f"  {name:<35s}  N={g.number_of_nodes():3d}  E={g.number_of_edges():4d}")

# Test all against Budapest brain
G_brain = load_budapest_coarse()
m_real = all_metrics(G_brain, "REAL")
print(f"\nReal brain: N={G_brain.number_of_nodes()} E={G_brain.number_of_edges()}")

print(f"\n{'Structure':<40s}  N    E   d_real  p_ER")
print("-" * 75)
mythology_results = {}
for name, g in candidates.items():
    if not nx.is_connected(g):
        print(f"  {name:<38s}  disconnected")
        continue
    try:
        m = all_metrics(g, name)
    except Exception as e:
        print(f"  {name:<38s}  metric err: {str(e)[:30]}")
        continue
    d = d_invariant(m, m_real)
    er = er_ensemble(g.number_of_nodes(), g.number_of_edges(), k=100, seed=0)
    er_dists = [d_invariant(all_metrics(ng, "n"), m_real) for ng in er]
    p_er = float((np.array(er_dists) <= d).mean()) if er_dists else float("nan")
    mythology_results[name] = {"N": g.number_of_nodes(), "E": g.number_of_edges(),
                                "d": d, "p_ER": p_er}
    print(f"  {name:<38s}  {g.number_of_nodes():3d}  {g.number_of_edges():4d}  "
          f"{d:.4f}  {p_er:.4f}")

# Sort
print("\nRANKED by distance to brain:")
for i, (name, r) in enumerate(sorted(mythology_results.items(),
                                     key=lambda x: x[1]["d"]), 1):
    sig = "✓" if r["p_ER"] < 0.05 else "✗"
    print(f"  {i:2d}. {name:<38s}  d={r['d']:.4f}  p_ER={r['p_ER']:.4f}  {sig}")


# ===========================================================================
# C. SENSITIVITY / PERTURBATION
# ===========================================================================
print("\n" + "=" * 78)
print("C. SENSITIVITY: perturbing joined-trees and measuring brain-match")
print("=" * 78)

J_base = joined_trees(shared_daath=True)
m_base = all_metrics(J_base, "base")
d_base = d_invariant(m_base, m_real)
print(f"\nBaseline: d(joined-trees, brain) = {d_base:.4f}")

print("\nRandom edge perturbation: add k random edges, measure new distance")
sensitivity = {"add": [], "remove": [], "rewire": []}
rng = np.random.default_rng(42)
n_trials = 30
for k in (1, 2, 3, 5, 10):
    add_dists = []
    for trial in range(n_trials):
        Jp = J_base.copy()
        # Pick k pairs of non-adjacent nodes
        nodes = list(Jp.nodes())
        added = 0
        while added < k:
            u, v = rng.choice(nodes, 2, replace=False)
            if u != v and not Jp.has_edge(u, v):
                Jp.add_edge(u, v); added += 1
        try:
            mp = all_metrics(Jp, "p")
            add_dists.append(d_invariant(mp, m_real))
        except Exception:
            continue
    arr = np.array(add_dists)
    if len(arr):
        sensitivity["add"].append({"k": k, "mean": float(arr.mean()),
                                   "std": float(arr.std()),
                                   "n": len(arr)})
        print(f"  +{k} edges:  d={arr.mean():.4f} ± {arr.std():.4f} "
              f"(baseline {d_base:.4f})")

# Edge removal
for k in (1, 2, 3, 5, 10):
    rem_dists = []
    edges = list(J_base.edges())
    for trial in range(n_trials):
        Jp = J_base.copy()
        # Remove k random edges that don't disconnect graph
        attempts = 0; removed = 0
        while removed < k and attempts < 50:
            attempts += 1
            e = edges[rng.integers(len(edges))]
            if Jp.has_edge(*e):
                Jp.remove_edge(*e)
                if nx.is_connected(Jp):
                    removed += 1
                else:
                    Jp.add_edge(*e)
        if removed != k: continue
        try:
            mp = all_metrics(Jp, "p")
            rem_dists.append(d_invariant(mp, m_real))
        except Exception:
            continue
    if rem_dists:
        arr = np.array(rem_dists)
        sensitivity["remove"].append({"k": k, "mean": float(arr.mean()),
                                      "std": float(arr.std()),
                                      "n": len(arr)})
        print(f"  -{k} edges:  d={arr.mean():.4f} ± {arr.std():.4f}")

# Random rewire
for k in (5, 10, 20):
    rew_dists = []
    for trial in range(n_trials):
        Jp = J_base.copy()
        nodes = list(Jp.nodes())
        edges = list(Jp.edges())
        rewires = 0
        while rewires < k:
            old = edges[rng.integers(len(edges))]
            if not Jp.has_edge(*old): continue
            Jp.remove_edge(*old)
            u, v = rng.choice(nodes, 2, replace=False)
            if u != v and not Jp.has_edge(u, v) and \
               nx.is_connected(Jp.copy()) if False else True:
                Jp.add_edge(u, v)
                rewires += 1
            else:
                Jp.add_edge(*old)
        if not nx.is_connected(Jp): continue
        try:
            mp = all_metrics(Jp, "p")
            rew_dists.append(d_invariant(mp, m_real))
        except Exception:
            continue
    if rew_dists:
        arr = np.array(rew_dists)
        sensitivity["rewire"].append({"k": k, "mean": float(arr.mean()),
                                      "std": float(arr.std()),
                                      "n": len(arr)})
        print(f"  rewire {k} edges:  d={arr.mean():.4f} ± {arr.std():.4f}")


# ===========================================================================
# Save everything
# ===========================================================================
output = {
    "cross_species": cross_species,
    "expanded_mythology": mythology_results,
    "sensitivity": sensitivity,
    "baseline_d": d_base,
}
with open("data/deep_research_results.json", "w") as f:
    json.dump(output, f, indent=2,
              default=lambda o: int(o) if isinstance(o, np.integer) else float(o))
print(f"\nWrote data/deep_research_results.json")

print("\nDONE.")
