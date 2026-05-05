"""
Graph constructors for the analysis.

Module exports:
  tree_of_life(include_daath=False) -> nx.Graph
  tree_of_death(include_daath=False) -> nx.Graph
  joined_trees(shared_daath=True)   -> nx.Graph     # the user's hypothesis
  jacobs_ladder()                   -> nx.Graph     # 4-Worlds variant
  brain_model(scale=20, ...)        -> nx.Graph
  brain_ensemble(n=200, ...)        -> list[nx.Graph]   # null comparison set
"""

from __future__ import annotations

import numpy as np
import networkx as nx


# ---------------------------------------------------------------------------
# Kabbalistic graphs
# ---------------------------------------------------------------------------

SEPHIROT_10 = [
    "Keter", "Chokhmah", "Binah",
    "Chesed", "Geburah", "Tiferet",
    "Netzach", "Hod", "Yesod", "Malkuth",
]

# Standard 22-path Kircher Tree of Life (the canonical Hermetic arrangement
# used by the Golden Dawn, Crowley, Regardie). Path count = 22 = number of
# Hebrew letters.
TOL_PATHS_22 = [
    ("Keter", "Chokhmah"), ("Keter", "Binah"), ("Keter", "Tiferet"),
    ("Chokhmah", "Binah"), ("Chokhmah", "Tiferet"), ("Chokhmah", "Chesed"),
    ("Binah", "Tiferet"), ("Binah", "Geburah"),
    ("Chesed", "Geburah"), ("Chesed", "Tiferet"), ("Chesed", "Netzach"),
    ("Geburah", "Tiferet"), ("Geburah", "Hod"),
    ("Tiferet", "Netzach"), ("Tiferet", "Yesod"), ("Tiferet", "Hod"),
    ("Netzach", "Hod"), ("Netzach", "Yesod"), ("Netzach", "Malkuth"),
    ("Hod", "Yesod"), ("Hod", "Malkuth"),
    ("Yesod", "Malkuth"),
]

# Daath's "hidden paths" (Crowley, 777; Regardie, Garden of Pomegranates).
# Daath does not have one of the 22 letter-paths, but is connected to the
# supernal triad and to Tiferet by veil-paths in most Hermetic accounts.
DAATH_PATHS = [
    ("Daath", "Chokhmah"),
    ("Daath", "Binah"),
    ("Daath", "Tiferet"),
]

QLI_NAMES = [
    "Thaumiel", "Ghagiel", "Sathariel",
    "Gamchicoth", "Golachab", "Thagirion",
    "Harab Serapel", "Samael", "Gamaliel", "Lilith",
]

# Canonical positions for the "three pillars" rendering
POS_TOL = {
    "Keter":    ( 0.0, 5.0),
    "Chokhmah": ( 1.0, 4.2),
    "Binah":    (-1.0, 4.2),
    "Daath":    ( 0.0, 3.6),
    "Chesed":   ( 1.0, 3.0),
    "Geburah":  (-1.0, 3.0),
    "Tiferet":  ( 0.0, 2.2),
    "Netzach":  ( 1.0, 1.2),
    "Hod":      (-1.0, 1.2),
    "Yesod":    ( 0.0, 0.4),
    "Malkuth":  ( 0.0,-0.6),
}


def tree_of_life(include_daath: bool = False) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from(SEPHIROT_10)
    g.add_edges_from(TOL_PATHS_22)
    if include_daath:
        g.add_node("Daath")
        g.add_edges_from(DAATH_PATHS)
    return g


def tree_of_death(include_daath: bool = False) -> nx.Graph:
    """Same topology as Tree of Life, with qliphothic names."""
    name_map = dict(zip(SEPHIROT_10, QLI_NAMES))
    g = nx.Graph()
    g.add_nodes_from(QLI_NAMES)
    g.add_edges_from((name_map[a], name_map[b]) for a, b in TOL_PATHS_22)
    if include_daath:
        g.add_node("Daath")
        for _, target in DAATH_PATHS:
            g.add_edge("Daath", name_map[target])
    return g


def joined_trees(shared_daath: bool = True) -> nx.Graph:
    """The user's hypothesis: Tree of Life + Tree of Death joined.

    If shared_daath: Daath is a single node belonging to both trees
    (acting as bridge / corpus-callosum analogue).
    Else: each sephira is connected to its qliphothic mirror by an edge.
    """
    tol = tree_of_life(include_daath=shared_daath)
    qli = tree_of_death(include_daath=shared_daath)
    if shared_daath:
        # Compose: Daath appears once, edges from both trees attach to it.
        g = nx.compose(tol, qli)
        return g
    else:
        g = nx.disjoint_union(tol, qli)
        # Reconstruct names since disjoint_union renumbers
        relabel = {}
        nodes = list(tol.nodes()) + list(qli.nodes())
        for i, n in enumerate(nodes):
            relabel[i] = n
        g = nx.relabel_nodes(g, relabel)
        # Add mirror edges
        for s, q in zip(SEPHIROT_10, QLI_NAMES):
            g.add_edge(s, q)
        return g


def jacobs_ladder() -> nx.Graph:
    """Four interlocking Trees of Life — the 'Four Worlds' (Atziluth, Briah,
    Yetzirah, Assiah). Adjacent worlds share Malkuth-Keter (the lower world's
    crown is the upper world's foundation). Used in Lurianic Kabbalah and
    Halevi's reconstructions.
    """
    worlds = ["Atz", "Bri", "Yet", "Ass"]
    g = nx.Graph()
    for w in worlds:
        for s in SEPHIROT_10:
            g.add_node(f"{w}_{s}")
        for a, b in TOL_PATHS_22:
            g.add_edge(f"{w}_{a}", f"{w}_{b}")
    # Bridge: Malkuth_(world n) <-> Keter_(world n+1)
    for w1, w2 in zip(worlds[:-1], worlds[1:]):
        g.add_edge(f"{w1}_Malkuth", f"{w2}_Keter")
    return g


# ---------------------------------------------------------------------------
# Brain models
# ---------------------------------------------------------------------------

# Coarse-scale parcellation (10 regions per hemisphere) — major lobes +
# subcortical landmarks. Edges follow well-documented major white-matter
# tracts (SLF, ILF, IFOF, uncinate, cingulum, arcuate, thalamic radiations).
# References: Catani & Thiebaut de Schotten (2008); Wakana et al. (2007).

BRAIN10_REGIONS = [
    "Frontal", "Prefrontal", "Motor",
    "Parietal", "Somatosensory",
    "Temporal", "Auditory",
    "Occipital",
    "Hippocampus",
    "Thalamus",
]

BRAIN10_INTRA = [
    ("Prefrontal", "Frontal"),
    ("Prefrontal", "Motor"),
    ("Frontal", "Motor"),
    ("Motor", "Somatosensory"),
    ("Somatosensory", "Parietal"),
    ("Parietal", "Occipital"),
    ("Parietal", "Temporal"),
    ("Temporal", "Auditory"),
    ("Temporal", "Occipital"),
    ("Temporal", "Hippocampus"),
    ("Frontal", "Temporal"),
    ("Thalamus", "Frontal"),
    ("Thalamus", "Parietal"),
    ("Thalamus", "Occipital"),
    ("Thalamus", "Hippocampus"),
    ("Prefrontal", "Hippocampus"),
]

# Desikan-Killiany-style ~34 regions per hemisphere, condensed for tractability
BRAIN34_REGIONS = [
    # Frontal
    "SuperiorFrontal", "RostralMiddleFrontal", "CaudalMiddleFrontal",
    "ParsOpercularis", "ParsTriangularis", "ParsOrbitalis",
    "LateralOrbitofrontal", "MedialOrbitofrontal", "Precentral", "Paracentral",
    "FrontalPole", "RostralAnteriorCingulate", "CaudalAnteriorCingulate",
    # Parietal
    "Postcentral", "SuperiorParietal", "InferiorParietal",
    "Supramarginal", "Precuneus", "PosteriorCingulate", "IsthmusCingulate",
    # Temporal
    "SuperiorTemporal", "MiddleTemporal", "InferiorTemporal",
    "BanksSTS", "Fusiform", "TransverseTemporal", "Entorhinal",
    "TemporalPole", "ParaHippocampal",
    # Occipital
    "LateralOccipital", "Lingual", "Pericalcarine", "Cuneus",
    # Insula
    "Insula",
]


def _hemi(side, name): return f"{side}-{name}"


def brain_model_10(corpus_callosum_node: bool = True) -> nx.Graph:
    """20-node brain (10 regions × 2 hemispheres), CC as bridge node."""
    g = nx.Graph()
    for side in ("L", "R"):
        for r in BRAIN10_REGIONS:
            g.add_node(_hemi(side, r))
        for a, b in BRAIN10_INTRA:
            g.add_edge(_hemi(side, a), _hemi(side, b))
    if corpus_callosum_node:
        g.add_node("CorpusCallosum")
        for side in ("L", "R"):
            for r in BRAIN10_REGIONS:
                if r in ("Hippocampus", "Thalamus"):
                    continue
                g.add_edge("CorpusCallosum", _hemi(side, r))
    else:
        for r in BRAIN10_REGIONS:
            if r != "Hippocampus":
                g.add_edge(_hemi("L", r), _hemi("R", r))
    return g


def brain_model_34(seed: int = 0,
                   corpus_callosum_node: bool = True) -> nx.Graph:
    """68-node brain (~34 regions × 2 hemispheres) with a topology built
    from anatomical priors.

    Construction:
      1. Each hemisphere has lobe communities (frontal/parietal/temporal/
         occipital/insula). Within-lobe density is high.
      2. Cross-lobe edges follow major fasciculi (SLF, ILF, IFOF, arcuate,
         uncinate, cingulum) — modelled by lobar adjacency rules.
      3. Inter-hemispheric: corpus callosum as a single bridge node attached
         to ~50% of cortical regions per side (homologous + nearby), since
         callosal projections are widespread but not exhaustive
         (Aboitiz et al. 1992).
      4. Long-range "rich club" edges added between high-degree nodes
         (van den Heuvel & Sporns 2011).
    """
    rng = np.random.default_rng(seed)

    lobe_of = {}
    for r in BRAIN34_REGIONS:
        if any(k in r for k in ("Frontal", "Cingulate", "Precentral",
                                "Paracentral", "Orbital", "Pars")):
            lobe_of[r] = "frontal"
        elif any(k in r for k in ("Parietal", "Postcentral", "Precuneus",
                                  "Supramarginal", "Isthmus")):
            lobe_of[r] = "parietal"
        elif any(k in r for k in ("Temporal", "Fusiform", "Entorhinal",
                                  "BanksSTS", "Hippo", "Insula")):
            lobe_of[r] = "temporal"
        elif any(k in r for k in ("Occipital", "Lingual", "Pericalcarine",
                                  "Cuneus")):
            lobe_of[r] = "occipital"
        else:
            lobe_of[r] = "other"

    # adjacency between lobes (anatomically-grounded)
    cross_lobe = {
        ("frontal", "parietal"):  0.45,   # SLF, central sulcus
        ("frontal", "temporal"):  0.30,   # uncinate, IFOF
        ("frontal", "occipital"): 0.05,
        ("parietal", "temporal"): 0.35,   # arcuate
        ("parietal", "occipital"):0.30,   # SLF posterior
        ("temporal", "occipital"):0.35,   # ILF
    }
    within_lobe_p = 0.55

    g = nx.Graph()
    for side in ("L", "R"):
        regions = [_hemi(side, r) for r in BRAIN34_REGIONS]
        g.add_nodes_from(regions)
        for i, r1 in enumerate(BRAIN34_REGIONS):
            for r2 in BRAIN34_REGIONS[i + 1:]:
                l1, l2 = lobe_of[r1], lobe_of[r2]
                if l1 == l2:
                    p = within_lobe_p
                else:
                    p = cross_lobe.get(tuple(sorted([l1, l2])), 0.10)
                if rng.random() < p:
                    g.add_edge(_hemi(side, r1), _hemi(side, r2))

    # Corpus callosum bridge
    if corpus_callosum_node:
        g.add_node("CorpusCallosum")
        for side in ("L", "R"):
            for r in BRAIN34_REGIONS:
                if "Hippo" in r or "Entorhinal" in r:
                    continue   # via fornix, not callosum
                if rng.random() < 0.55:
                    g.add_edge("CorpusCallosum", _hemi(side, r))
    else:
        # Direct homologous connections
        for r in BRAIN34_REGIONS:
            if "Hippo" not in r:
                g.add_edge(_hemi("L", r), _hemi("R", r))

    # Rich-club: connect top-degree nodes (van den Heuvel & Sporns 2011)
    deg = sorted(g.degree, key=lambda x: -x[1])
    hubs = [n for n, _ in deg[:int(0.10 * g.number_of_nodes())]]
    for i, h1 in enumerate(hubs):
        for h2 in hubs[i + 1:]:
            if not g.has_edge(h1, h2) and rng.random() < 0.6:
                g.add_edge(h1, h2)

    # Ensure connectedness
    while not nx.is_connected(g):
        comps = list(nx.connected_components(g))
        a = next(iter(comps[0]))
        b = next(iter(comps[1]))
        g.add_edge(a, b)

    return g


def brain_ensemble(n: int = 200,
                   scale: int = 34,
                   corpus_callosum_node: bool = True) -> list[nx.Graph]:
    """Generate an ensemble of brain-model realizations for null comparison.

    Each realization uses different random seeds for the stochastic edges
    but identical anatomical priors. This gives a distribution of brain-like
    networks against which to compare the Trees.
    """
    if scale == 10:
        # The 10-region model is deterministic; ensemble doesn't apply.
        return [brain_model_10(corpus_callosum_node)]
    return [brain_model_34(seed=s, corpus_callosum_node=corpus_callosum_node)
            for s in range(n)]
