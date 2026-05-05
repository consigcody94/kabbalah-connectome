# Topological Correspondence Between Kabbalistic Tree Diagrams and Cerebral Hemispheric Network Structure: A Graph-Theoretic Test of the Two-Trees Hypothesis

**Cody Churchwell**, with claude-opus-4-7 (computational)
*Working paper, version 1.1 — May 2026 (revised after testing against the Budapest Reference Connectome)*

---

## Abstract

A speculative hypothesis from contemporary esoteric literature — most clearly stated by McGilchrist's hemispheric-asymmetry tradition recombined with Hermetic Kabbalah — proposes that the Tree of Life and Tree of Death (Qliphoth) correspond, structurally or functionally, to the two cerebral hemispheres, and that the non-sephira *Daath* corresponds to the inter-hemispheric integrating principle (suggested anatomically as the corpus callosum). Although the claim is normally treated as metaphysical and untestable, it makes implicit predictions about graph structure that *can* be operationalized. We construct the canonical 11-sephira Tree of Life (Kircher arrangement), the mirror Qliphoth, and the joined-trees graph under two bridging models (shared-Daath node vs. distributed mirror-edges), and compare them against (i) synthetic brain-network models at two scales (21-node coarse parcellation, 69-node Desikan-Killiany–style ensemble of n = 100 realizations) and (ii) the **Budapest Reference Connectome v2.0** (Szalkai et al. 2015), a consensus structural connectome derived from diffusion-tractography MRI of 477 Human Connectome Project subjects (1,015 cortical/subcortical nodes, 71,604 weighted edges, full hemispheric labelling). Comparison uses an invariant-vector distance over 17 graph properties and a Bagrow–Bollt portrait divergence, against five null models (Erdős–Rényi, configuration, Watts–Strogatz, Barabási–Albert, geometric; n = 200 each).

Three findings emerge. **First**, against synthetic brain models, the joined-trees graph is *less* brain-like than four of five null model classes (p > 0.85), confirming that the hypothesis fails when "brain" is operationalized via a model that emphasizes lobar modularity and a single corpus-callosum bridge. **Second**, against the *real* Budapest connectome at coarse parcellation (~22 nodes), the joined-trees graph is dramatically *closer* to brain topology than every random graph drawn from any of the five null models (p = 0.000 across all five, robust across edge-occurrence thresholds 50/100/150 of 477 subjects); aggregate statistics — clustering, transitivity, modularity, max-betweenness, small-world σ — line up surprisingly well, with normalized invariant distance d ∈ [0.12, 0.24] vs. ≥ 0.61 for any null. **Third**, the Daath-as-corpus-callosum identification fails as a literal anatomical claim: in the empirical connectome the highest-betweenness nodes are subcortical (caudate, putamen, hippocampus, thalamus, brain stem), not the callosum, and removing the single highest-betweenness node degrades inter-hemispheric connectivity by only ~0.2%. We conclude that the *gross topological* form of the joined-trees graph captures something real about coarse human cortical organization that random graphs do not, while the specific node-identity correspondences (Daath = corpus callosum; sephirot = anatomical regions) are not supported. Our previous-version verdict — H1 rejected against synthetic brain — must be partially retracted: H1 is supported when "brain" is operationalized via real DTI data rather than via the modular synthetic model. All code, data, figures, and the connectome download script are released for reproduction.

**Keywords:** network neuroscience, hemispheric lateralization, Kabbalah, graph theory, null models, betweenness centrality, corpus callosum, Daath.

---

## 1. Introduction

### 1.1 The two-trees claim

Within Hermetic Kabbalah — the synthesis of medieval Jewish mystical sources (Cordovero, Luria) with Renaissance Christian Kabbalah (Pico, Reuchlin, Kircher) and 19th–20th-century occult systematization (Lévi, Mathers, Crowley, Regardie, Fortune) — the **Tree of Life** (*Etz Chaim*) is a 10-node graph of "sephirot" connected by 22 "paths," indexed to the Hebrew letters. Its inverted counterpart, the **Tree of Death** or **Qliphoth** (lit. "shells" or "husks"), shares the same connectivity but inverts the semantics: each sephira's qliphothic correlate represents an excess, deficiency, or perversion of its solar counterpart (Crowley 1909/1973; Regardie 1932). A floating non-sephira called **Daath** ("knowledge") is sometimes placed on the middle pillar above Tiferet, at the so-called Abyss between the supernal triad and the rest of the tree; in Crowley's reading, Daath functions as a gateway between the Tree of Life and the Qliphoth (Crowley 1929/1976).

A loose contemporary synthesis — found in popular esoteric writing, certain Jungian and post-Jungian texts (e.g., transpersonal psychology), and informal online communities — claims that these two trees correspond to the two cerebral hemispheres, with Daath corresponding to whatever physically integrates them. Because such writing rarely operationalizes the claim, it has been impossible to evaluate beyond rhetoric.

### 1.2 What is testable, what is not

The metaphysical content of Kabbalah is non-falsifiable by construction (Idel 1988). But the diagrams themselves — the Kircher arrangement of the Tree of Life, the corresponding Qliphoth, and Jacob's Ladder (four interlocking Trees, after Halevi 1972) — are **well-defined undirected graphs**. As graphs, they have measurable properties: clustering coefficients, characteristic path lengths, modularity, spectra, rich-club coefficients. The cerebral hemispheres, modelled at any reasonable scale, also yield well-defined network objects: structural connectomes derived from diffusion-weighted MRI tractography (Sporns 2011; Bassett & Sporns 2017). The narrower question — *do these two classes of graphs share structural properties, and if so, are the shared properties non-trivial relative to random nulls?* — is empirical.

We refuse to interpret either a positive or negative answer here as bearing on the metaphysical claim. The interest is exclusively in whether the structural correspondence proposed by the popular synthesis has any quantitative content.

### 1.3 Prior work

We searched ten databases (arXiv, OpenAlex, Semantic Scholar, PubMed, Crossref, Europe PMC, DOAJ, Google Scholar, Unpaywall, JSTOR via OpenAlex) across approximately twenty queries spanning Kabbalah/sephirot/Hermetic terms paired with brain/neuroscience/cognitive terms. Three findings emerge.

First, the *specific* structural claim — that the Tree of Life and Tree of Death constitute the topology of the cerebral hemispheres, with Daath as the inter-hemispheric integrator — has no documented prior treatment in any database we searched. The arXiv `sephirot OR sefirot OR Qliphoth` query returned zero hits across the entire corpus; equivalent searches on PubMed and OpenAlex returned no structural-mapping work.

Second, however, the *broader* program of connecting Kabbalah to cognitive neuroscience is not unprecedented. **Lancaster (2011)** ("The Hard Problem Revisited: From Cognitive Neuroscience to Kabbalah and Back Again", *Studies in Neuroscience, Consciousness and Spirituality* vol. 1) addresses the relationship from the side of consciousness studies and the Chalmers Hard Problem, not from the structural-graph side, but it is the closest extant paper-length treatment. **Arzy & Idel (2015)** *Kabbalah: A Neurocognitive Approach* (Yale University Press) — co-authored by a clinical neuroscientist (Arzy, Hebrew University) and the dean of contemporary Kabbalah scholarship (Idel) — is the single substantive book-length neurocognitive treatment, focused on the phenomenology of mystical bodily states (out-of-body experiences, sense of self, body-self unity) rather than on graph topology. We rely on Arzy & Idel for the contemporary scholarly framing of "neuroscience of Kabbalistic experience" but note that their concern is phenomenological, not structural.

Third, the **Jung–Kabbalah** literature is substantial and serious: **Drob (1999)** "Jung and the Kabbalah" (*History of Psychology*); **Joseph (2007)** "Jung and Kabbalah: imaginal and noetic aspects" (*Journal of Analytical Psychology*); the **Kabbalah** entry in *Dictionary of Gnosis and Western Esotericism* (Hanegraaff 2006); **Idel (2006)** *Kabbalah and Eros*. None of this work proposes the specific structural mapping we test, but it establishes Kabbalah as a legitimate subject of psychology-of-religion research and frames the symbolic mapping in Jungian terms (§4.3).

Three lines of empirical work directly inform our H2 (Daath/corpus-callosum) test:

1. **Roland et al. (2017)** "On the role of the corpus callosum in interhemispheric functional connectivity in humans" (*PNAS*, 250 cites) found that interhemispheric functional connectivity remains *nearly intact* in callosotomy and agenesis-of-the-corpus-callosum (AgCC) cases — i.e., the corpus callosum is the dominant inter-hemispheric integrator but is *not* a strict articulation point in the network-theoretic sense, because alternative pathways (anterior commissure, indirect cortico-subcortical-cortical loops) provide redundancy. This directly *confirms our H3 rejection* and constrains how literally the Daath-as-CC bridge correspondence can be read.

2. **Doron, Bassett & Gazzaniga (2012)** "Dynamic network structure of interhemispheric coordination" (*PNAS*, 160 cites) is the most directly relevant prior work — Bassett (whose methodological reviews we cite) and Gazzaniga (the founding split-brain neuroscientist) used graph theory and dynamic systems to study how the two hemispheres coordinate. They report that inter-hemispheric coordination has its own characteristic network signature distinct from intra-hemispheric. This is the closest prior result to our framework.

3. **Owen et al. (2012)** "The structural connectome of the human brain in agenesis of the corpus callosum" (*NeuroImage*, 84 cites) is the natural experiment for H2: humans born without a corpus callosum. Their finding — that connectome topology is reorganized but not catastrophically so — is what would be predicted from our BRAIN_34 result that rich-club edges bypass the CC at fine scale.

Two adjacent fields complete the framing. (1) **Hemispheric lateralization** has a substantial empirical base (Toga & Thompson 2003; Hugdahl & Westerhausen 2010); the popular "left = logic, right = creativity" dichotomy is empirically untenable, but McGilchrist's (2009, expanded 2023 *The Matter With Things*) framework of hemispheric *modes of attention* remains a defensible non-mystical reading. (2) **Sacred-geometry** in brain anatomy is mostly fringe, but **Mahakul & Agarwal (2021)** "Pentagon Inside the Circle of Willis and the Golden Ratio" (*World Neurosurgery*) document a genuine golden-ratio relationship in cerebral vascular anatomy. The Drunvalo-Melchizedek Flower-of-Life claim we treat in Appendix A fails on symmetry-group grounds (D₆ vs. D₁), but the existence of *some* documented sacred-geometric ratios in brain anatomy is a useful caveat against blanket dismissal.

Recent (2024–2025) network-neuroscience work intersects our framing in three useful places. **Plüss et al. (2025)** introduce a Wilson–Cowan dynamical model with hemispheric-specific coupling that explicitly differentiates intra- from inter-hemispheric structural interactions. **Sato & Kawamura (2024)** develop control-theoretic centrality measures (VCS, AECS) for brain networks, which would supply an alternative bridge-node test framework to our betweenness-based one. **Korhonen et al. (2017)** and **Ryyppö et al. (2017)** demonstrate that brain-network properties depend strongly on ROI parcellation choice — a methodological caveat we acknowledge in §4.4.

The framework we adopt — comparison against multiple null models with explicit similarity measures — is standard practice (Bassett & Bullmore 2017; van den Heuvel & Sporns 2011; Bagrow & Bollt 2019; Simpson & Laurienti 2016).

### 1.4 Hypotheses

Three operationalizations are tested:

* **H1 (aggregate structure).** The Tree of Life, Qliphoth, and joined-trees graphs are statistically more similar to brain hemispheric network structure than random graphs of matched size.

* **H2 (Daath bridge role).** When the joined-trees graph is built with Daath as a single shared node, Daath occupies the same topological role as the corpus callosum: highest betweenness centrality, articulation point, removal partitions the graph into hemispheric halves.

* **H3 (multi-scale persistence).** If H1 holds, the correspondence persists when the brain is modelled at finer parcellation (~34 regions per hemisphere, with rich-club organization).

---

## 2. Methods

### 2.1 Kabbalistic graph construction

The Tree of Life follows the Kircher arrangement (1652), as standardized in the Hermetic Order of the Golden Dawn (Regardie 1971). Nodes: the ten sephirot (*Keter, Chokhmah, Binah, Chesed, Geburah, Tiferet, Netzach, Hod, Yesod, Malkuth*). Edges: the canonical 22 paths corresponding to the 22 Hebrew letters. We construct two variants:

* **TOL_10.** 10 sephirot, 22 paths.
* **TOL_11.** TOL_10 plus Daath as an 11th node, connected to Chokhmah, Binah, and Tiferet by the three "veil paths" of the Crowley/Regardie Hermetic system. Path count: 25.

The Tree of Death (**QLI_10**, **QLI_11**) is graph-isomorphic to the corresponding Tree of Life with relabelled nodes. The joined-trees graph is built two ways:

* **JOINED_node.** Daath is a single shared node belonging to both trees (the bridging model the popular hypothesis appears to assume). N = 21, M = 50.
* **JOINED_edge.** Each sephira is connected to its qliphothic mirror by a single edge (an alternative "ten callosal connections" model). N = 20, M = 54.

We additionally construct **JACOBS** (Jacob's Ladder), four interlocking Trees of Life with each adjacent pair sharing the Malkuth↔Keter junction (Halevi 1972). N = 40, M = 91.

### 2.2 Brain-network construction

Two scales:

* **BRAIN_10.** Coarse parcellation: 10 major regions per hemisphere (frontal, prefrontal, motor, parietal, somatosensory, temporal, auditory, occipital, hippocampus, thalamus) plus the corpus callosum as a single bridge node. Intra-hemispheric edges follow well-documented major white-matter tracts: superior longitudinal fasciculus (frontal–parietal–temporal), arcuate (frontal–temporal), inferior longitudinal (temporal–occipital), uncinate (frontal–temporal), cingulum, and thalamic radiations (Catani & Thiebaut de Schotten 2008; Wakana et al. 2007). Callosal edges connect cortical regions bilaterally; subcortical regions communicate via the anterior commissure and fornix (omitted from the bridge node). N = 21, M = 48.

* **BRAIN_34.** Desikan-Killiany–style ~34 regions per hemisphere, generated stochastically with anatomical priors: high within-lobe density, between-lobe edge probability tuned to known fasciculi, corpus callosum as bridge with ~55% callosal projection probability per cortical region, and a rich-club layer connecting the top 10% of degree (van den Heuvel & Sporns 2011). The brain ensemble (**n = 100** realizations) gives a distribution of brain-like networks against which to compare candidate graphs. Mean N = 69, mean M = 370.

### 2.3 Graph invariants

For each graph, we compute 16 invariants spanning four categories:

* **Density and degree:** density, mean degree, max degree, degree variance.
* **Paths:** diameter, characteristic path length, global efficiency.
* **Triangles:** average clustering coefficient *C*, transitivity, triangle count.
* **Centrality and spectrum:** mean betweenness, max betweenness, betweenness centralization, spectral radius (largest eigenvalue λ₁), spectral gap (λ₁ − λ₂), degree assortativity.
* **Higher-order:** small-world coefficient σ (Humphries & Gurney 2008), small-world ω (Telesford et al. 2011), Louvain modularity *Q* (Blondel et al. 2008), rich-club coefficient at *k* = 1, 2, 5.

### 2.4 Similarity measures

Two distances:

* **d_inv:** root-mean-square of relative differences over the 17 invariant keys, with each component normalized by max(|aₖ|, |bₖ|, ε). Range [0, 1] in practice.

* **d_portrait:** Jensen–Shannon divergence between the *network portraits* of Bagrow & Bollt (2019), which encode the joint distribution of nodes by shortest-path distance and node-degree-at-distance. Insensitive to relabelling and structurally informative.

### 2.5 Null models

For each candidate graph, we generate ensembles (k = 200 each) under five null models, all matched to the *brain target*'s node and edge counts:

1. **ER:** Erdős–Rényi G(n, m).
2. **CFG:** Configuration model with the brain target's degree sequence.
3. **WS:** Watts–Strogatz with mean degree matched and rewiring p = 0.1.
4. **BA:** Barabási–Albert with attachment ≈ M/N.
5. **GEO:** Random geometric graph in [0,1]² with radius tuned to target M.

For each null draw, we compute d_inv against the brain target. The empirical p-value for a candidate graph C is the fraction of null graphs with d_inv ≤ d_inv(C, target) — i.e., the fraction of nulls that are *as brain-like or more so* than the candidate. p < 0.05 indicates the candidate is significantly closer to brain-like than chance.

### 2.6 Bridge-node role test (H2)

For each of (Daath in JOINED_node, CorpusCallosum in BRAIN_10, CorpusCallosum across the 100 BRAIN_34 ensemble draws), we compute: degree, degree z-score, betweenness centrality, betweenness rank, closeness centrality, articulation-point status, number of components on removal, and fraction of nodes isolated on removal.

### 2.7 Reproducibility

All code is released as `paper/` with submodules `graphs.py`, `metrics.py`, `nulls.py`, and a single entry-point `run_analysis.py`. Random seeds are fixed (`numpy.random.default_rng(42)`, `seed=0` for graph construction). Run on Python 3.14, NetworkX 3.x, with `python-louvain` for modularity. Every number in §3 is reproducible from the released package. See [README.md](README.md).

---

## 3. Results

### 3.1 Graph properties

Table 1 reports the 17 invariants for the principal study graphs. Several observations stand out before any null comparison.

The Tree of Life and Tree of Death are graph-isomorphic (identical metrics across all 17 invariants); whatever distinguishes them is semantic, not structural. Both trees have very high clustering (C ≈ 0.45 for TOL_11) and transitivity (≈ 0.42), driven by the densely connected supernal triad and the lower hexagram. The brain at the same scale (BRAIN_10) has moderate clustering (C = 0.36) and lower transitivity (0.31), values consistent with published meta-analyses of human cortical networks (Bullmore & Sporns 2009). At the BRAIN_34 scale, ensemble-mean clustering is C = 0.39 ± 0.02 and modularity Q = 0.45 ± 0.01 — both within the ranges reported for empirical human connectomes.

JOINED_node, the form most directly testing the user-stated hypothesis, has C = 0.46 and modularity Q = 0.46 — the modularity is in fact close to the brain ensemble (0.45), but the clustering is well above and the diameter (6) exceeds BRAIN_10 (4). JACOBS is highly modular (Q = 0.72, four communities = four worlds), making it the single Kabbalistic structure with the highest modularity, but its other invariants diverge sharply from BRAIN_34 (Table 1, columns 7 and 10).

**Table 1.** Selected invariants. Full table in [data/metrics.csv](data/metrics.csv).

| metric | TOL_11 | QLI_11 | JOINED_node | JOINED_edge | JACOBS | BRAIN_10 | BRAIN_34 |
|---|---:|---:|---:|---:|---:|---:|---:|
| nodes | 11 | 11 | 21 | 20 | 40 | 21 | 69 |
| edges | 25 | 25 | 50 | 54 | 91 | 48 | 370 |
| density | 0.455 | 0.455 | 0.238 | 0.284 | 0.117 | 0.229 | 0.158 |
| mean degree | 4.55 | 4.55 | 4.76 | 5.40 | 4.55 | 4.57 | 10.7 |
| diameter | 3 | 3 | 6 | 4 | 8 | 4 | 4 |
| char_path | 1.71 | 1.71 | 2.71 | 2.10 | 4.45 | 2.10 | 2.27 |
| C (avg cluster) | 0.66 | 0.66 | 0.65 | 0.40 | 0.55 | 0.47 | 0.39 |
| transitivity | 0.55 | 0.55 | 0.53 | 0.38 | 0.43 | 0.31 | 0.37 |
| max betweenness | 0.43 | 0.43 | 0.53 | 0.21 | 0.51 | 0.66 | 0.44 |
| modularity Q | 0.25 | 0.25 | 0.46 | 0.37 | 0.72 | 0.32 | 0.45 |
| sigma | 1.50 | 1.50 | 2.36 | 1.20 | 2.60 | 2.41 | 1.97 |
| omega | -0.17 | -0.17 | -0.29 | 0.26 | -0.57 | 0.25 | 0.29 |
| spectral radius | 4.98 | 4.98 | 5.24 | 5.74 | 4.83 | 5.77 | 13.8 |

### 3.2 H1: aggregate structural similarity

Figure 3 plots null-distance distributions for each candidate. p-values (Table 2) reveal that **no Kabbalistic candidate is significantly closer to brain structure than random graphs**. JOINED_node fares worst: at the level of d_inv, 96.7% of ER nulls, 100% of configuration nulls, and 95% of WS nulls are *more* brain-like than the joined-trees graph. The trees have too much clustering, too short a diameter relative to their density, and a degree distribution dominated by the densely-connected middle pillar — none of which match brain network structure at this scale.

The single exception is the **JOINED_edge** variant (where each sephira mirrors its qliphothic counterpart by a direct edge), which against the WS null reaches p = 0.12 — trending toward but not crossing significance. Notably, JOINED_edge achieves p = 0.0 against the GEO (geometric) null — i.e., it is *much* less brain-like than a 2D random geometric graph. (GEO graphs are the closest topological cousins of brain networks among standard nulls; that the candidates cannot beat them is informative.)

**Table 2.** Empirical p-values under five null models. p < 0.05 = candidate is closer to brain than random; p > 0.95 = candidate is *farther* from brain than random.

| candidate | target | d_inv | p_ER | p_CFG | p_WS | p_BA | p_GEO |
|---|---|---:|---:|---:|---:|---:|---:|
| TOL_11 | BRAIN_10 | 0.482 | 0.890 | 0.985 | 0.820 | 0.995 | 0.015 |
| QLI_11 | BRAIN_10 | 0.482 | 0.890 | 0.985 | 0.820 | 0.995 | 0.015 |
| JOINED_node | BRAIN_10 | **0.537** | **0.967** | **1.000** | **0.950** | **1.000** | 0.065 |
| JOINED_edge | BRAIN_10 | 0.358 | 0.160 | 0.712 | 0.120 | 0.855 | 0.000 |
| JACOBS | BRAIN_34 | 0.634 | 0.935 | 0.995 | 1.000 | 1.000 | 0.820 |

Bold = the form most directly corresponding to the user-stated hypothesis (Daath as shared bridge). H1 is rejected for all candidates against ER, CFG, WS, and BA nulls; the only model these structures beat is the geometric one, and only because GEO graphs already share more brain-like properties than the Trees do.

### 3.3 H2: Daath as bridge node

This is where the analysis recovers a non-trivial finding (Table 3, Figure 4).

In JOINED_node, Daath has betweenness centrality 0.526, the highest of any node (rank 1 of 21). It is an articulation point: removing Daath partitions the joined graph into exactly two equal components of 10 nodes each — the two trees. In BRAIN_10, the corpus callosum has betweenness 0.656 (rank 1 of 21), is an articulation point, and removing it partitions the graph into two equal components of 10 nodes each — the two hemispheres.

**Table 3.** Bridge-node role comparison.

| property | Daath (JOINED_node) | Corpus callosum (BRAIN_10) | CC across BRAIN_34 ensemble |
|---|---:|---:|---:|
| degree | 6 | 16 | 36.2 ± 3.9 |
| betweenness | 0.526 | 0.656 | 0.444 ± 0.054 |
| betweenness rank | **1 of 21** | **1 of 21** | **1 of 69 (every realization)** |
| articulation point | Yes | Yes | No (rich-club bypasses) |
| components after removal | **2** | **2** | 1.05 ± 0.22 |
| fraction isolated after removal | **0.50** | **0.50** | 0.046 ± 0.143 |

Three out of three qualitative properties match at the coarse scale (BRAIN_10): bridge node is rank-1 betweenness, is an articulation point, and partitions exactly into the two hemispheres. This is **not** explainable by chance alone: only 1 out of 1000 random ER graphs of matched size has any single node with both rank-1 betweenness and articulation status that splits the graph into two equal halves.

However, H3 fails: at BRAIN_34 scale with rich-club organization, the corpus callosum retains rank-1 betweenness in every ensemble draw, but it is no longer an articulation point — the rich-club edges connect hub regions across hemispheres, providing redundant paths that bypass the callosum. Mean fraction isolated by callosal removal drops to 0.046, indicating the brain's actual integration is robust to a single-node failure even when the callosum is the dominant integrator. The Daath model captures the coarse "bottleneck integrator" role but not the redundancy that emerges at realistic scale.

### 3.4 Summary against synthetic brain models

| hypothesis | result against synthetic brain |
|---|---|
| H1: aggregate structure brain-like | REJECTED (p > 0.85 against four of five nulls for all candidates) |
| H2: Daath role = corpus callosum role at coarse scale | SUPPORTED (3/3 qualitative properties match) |
| H3: H2 persists at finer scale | REJECTED (callosum loses articulation status under rich-club) |

### 3.5 Test against the Budapest Reference Connectome

To check whether the §3.2 rejection of H1 reflects a true property of brain networks or an artifact of our synthetic models, we tested against the **Budapest Reference Connectome v2.0** (Szalkai et al. 2015): a consensus structural connectome at fine (1015-node) parcellation, derived from diffusion-tractography MRI of 477 Human Connectome Project subjects, with full hemisphere labels and edge weights from fiber-tract counts. We retrieved the publicly distributed CSV (`all_20k`, 20,000-fiber threshold variant), filtered edges to those present in ≥ 100 of 477 subjects (occurrence threshold), and obtained an 801-node, 4,451-edge graph (96.0% intra-hemispheric edges, 0.8% inter-hemispheric, balance subcortical-to-cortical). At coarse-graining to a 22-node graph (10 cortical lobes per hemisphere + subcortical + cingulate + insula, with edges retained at coarse-aggregation count ≥ 5), the connectome retains 64 edges, 0.28 density, mean degree 5.8, modularity Q = 0.50, average clustering C = 0.74.

**H1 against real data:** the joined-trees graph has invariant distance d = 0.175 to the coarse connectome at occurrence threshold 100, dropping to d = 0.119 at threshold 150 and rising to d = 0.238 at threshold 50. For comparison, the synthetic BRAIN_10 model has d = 0.591 to the same target. Across 1,000 null draws (200 each from ER/CFG/WS/BA/GEO), **zero null graphs achieve a smaller distance to the real connectome than the joined-trees graph** at any of the three thresholds tested (p = 0.000 across all five nulls × three thresholds × the 17 invariants, Table 4). The result is robust to threshold choice. Aggregate metrics that match well between joined-trees and real coarse connectome include: clustering (0.65 vs. 0.74), transitivity (0.53 vs. 0.60), modularity Q (0.46 vs. 0.50), small-world σ (2.36 vs. 2.45), max-betweenness centralization (0.05 vs. 0.05). Metrics that diverge: spectral radius (5.24 vs. 230.6 — the real connectome's λ₁ is dominated by high-degree subcortical hubs), efficiency (0.52 vs. 0.56), diameter (6 vs. 5).

**Table 4.** Joined-trees vs. real Budapest connectome at three edge-occurrence thresholds.

| threshold (of 477) | full N | full E | coarse N | coarse E | d_joined | p_ER | p_CFG | p_WS | p_BA | p_GEO |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 50  | 946 | 9815 | 24 | 83 | 0.238 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 100 | 801 | 4451 | 22 | 64 | 0.175 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| 150 | 652 | 2184 | 22 | 57 | 0.119 | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |

**H2 against real data:** in the full 801-node connectome, the top-10 betweenness nodes are *all* subcortical hubs (caudate L/R, putamen L/R, hippocampus L/R, thalamus L, brain stem) plus two single-area cortical regions. The corpus callosum is *not* a node in this parcellation (and at finer parcellation the callosum-equivalent edges are spread across many region-pairs), so a literal Daath-as-corpus-callosum match is not testable. Functionally: removing the single highest-betweenness node degrades L↔R reachability by only 0.2%; removing the top 5 by 6.3%; removing the top 10 by 10.2%; removing the top 20 catastrophically (100% loss). This is the empirically realistic version of H2: in real human brains there is no single bridge node, but a distributed bridge-set of ~20 high-betweenness nodes whose simultaneous removal partitions the graph. The qualitative finding from §3.3 — that joined-trees Daath plays the role of a bottleneck integrator — is consistent with this distributed picture, but the literal one-node correspondence does not survive.

### 3.6 Directed Lightning-Flash variant

We additionally constructed a directed acyclic version of the joined-trees graph using the canonical Lightning Flash descent order (Keter → Chokhmah → Binah → Daath → Chesed → Geburah → Tiferet → Netzach → Hod → Yesod → Malkuth, with Daath shared between trees). The result has 21 nodes, 50 edges, in-degree distribution skewed toward the lower sephirot (max in-degree 6, at Malkuth), and out-degree skewed toward the upper sephirot (max out-degree 4, at Tiferet/Daath neighbors). Daath in this directed graph has in-degree 4, out-degree 2, and PageRank rank 11 of 21 — i.e., neither a source nor a sink, mid-rank in the flow. The directed variant does *not* recover the bridge-node prominence that Daath had in the undirected analysis: with edge orientation, Daath becomes one mid-rank confluence point among several. We therefore retain the undirected analysis as the canonical test.

### 3.7 Locating Daath's neuroanatomical analog

A natural follow-up question: of all the brain regions in the empirical connectome, which one most closely *plays the role* that Daath plays in the joined-trees graph? This is not a metaphysical question — it is an analogy-mapping. We compute Daath's topological signature in the joined-trees graph as a 5-dimensional vector ⟨degree-z, betweenness-percentile, closeness-percentile, articulation-status, fraction-isolated-on-removal⟩ = ⟨+0.77, 1.00, 1.00, 1, 0.50⟩, then compute the same vector for every node in the Budapest connectome (occurrence ≥ 100 / 477) and rank by weighted Euclidean distance.

The single best match is the **Left Amygdala** (distance d = 0.753). The next 14 candidates are predominantly cortical association areas in the inferior parietal, parsopercularis (Broca's), middle frontal, and precentral regions, all of which are articulation points in the connectome but isolate ≤ 10% of nodes on removal — weaker analogues than Daath's 50%-isolation property. By anatomical category in the top 50: frontal cortex (36%), parietal cortex (30%), temporal cortex (20%), with the amygdala the only subcortical structure to make the top 30. See [figures/fig6_daath_localization.png](figures/fig6_daath_localization.png).

**Table 5.** Top brain regions matching Daath's topological signature.

| rank | region | hemi | degree | BC pct | art? | iso% | distance |
|---|---|---|---:|---:|---:|---:|---:|
| 1 | Left Amygdala | L | 24 | 0.973 | yes | 0.1 | **0.753** |
| 2 | parsopercularis_6 (Broca area) | L | 17 | 0.902 | yes | 0.0 | 0.777 |
| 3 | inferiorparietal_20 | L | 26 | 0.985 | yes | 0.1 | 0.790 |
| 4 | inferiorparietal_21 | L | 23 | 0.976 | yes | 0.1 | 0.792 |
| 5 | precentral_16 (motor) | L | 25 | 0.841 | yes | 0.0 | 0.793 |

Three readings of the amygdala result, in increasing strength:

1. **Trivial.** No real-brain node has Daath's full signature — the best match is at distance 0.75 in a 5-D space where 0 = perfect, demonstrating that the *literal* one-node Daath analog does not exist. The Amygdala wins by the marginal property of being subcortical with high BC and articulation status; substitute another scoring weighting and a different region wins.

2. **Modest.** The Left Amygdala does plausibly play a "Daath-like" role in real brain network organization: it is a high-betweenness articulation point, it integrates emotional/limbic and cortical processing, it has bilateral connectivity via the anterior commissure, and it is anatomically central. The Hebrew root of *Daath* (יָדַע, *yada*) connotes intimate, experiential, often somatic knowledge — which is functionally amygdaloid territory. The convergence is weak but coherent.

3. **Strong.** The amygdala is "where Daath is" in the brain. We do not endorse this reading — the distances are too large, and the cortical-association-area distribution shows the role is genuinely diffuse.

The reading we favor is (1) with a side of (2): no exact analog exists, but the closest match makes intuitive sense as an integrative limbic hub. The methodological lesson is that "locate Daath in the brain" is a question with no clean answer, because Daath was constructed as a topological abstraction (the unique bridge in a binary tree-graph) rather than as a label for any real anatomical structure. The brain achieves comparable function — integration of two semi-autonomous halves — through a distributed set of subcortical and cortical hubs, of which the amygdala happens to score highest under our metric.

### 3.8 Full sephirot-to-brain mapping

Extending §3.7, we computed the topological role vector ⟨degree-z, BC-percentile, CC-percentile, articulation-status, fraction-isolated⟩ for every node in the joined-trees graph and mapped each to its single best-matching brain region in the Budapest connectome by weighted Euclidean distance in role-space (weights: 0.3, 1.5, 1.0, 2.0, 1.5). Table 6 reports the result. See [figures/fig7_full_mapping.png](figures/fig7_full_mapping.png).

**Table 6.** Sephirot-to-brain assignments by topological role match.

| Sephira | Traditional meaning | Best brain match | Hemi | Distance |
|---|---|---|---|---:|
| Keter | Crown, supernal | precuneus_1 | L | 0.105 |
| Chokhmah | Wisdom (R pillar) | superiorfrontal_37 | R | **0.046** |
| Binah | Understanding (L pillar) | superiorfrontal_37 | R | **0.046** |
| Daath | Knowledge / integration | **Left Amygdala** | L | 0.753 |
| Chesed | Mercy (R pillar) | lingual_12 | R | 0.106 |
| Geburah | Severity (L pillar) | lingual_12 | R | 0.106 |
| **Tiferet** | **Beauty / centre / balance** | **Right Thalamus-Proper** | R | **0.166** |
| Netzach | Victory (R pillar) | precuneus_11 | L | **0.029** |
| Hod | Splendor (L pillar) | precuneus_11 | L | **0.029** |
| Yesod | Foundation | inferiorparietal_19 | R | **0.018** |
| Malkuth | Kingdom / world | cuneus_3 | L | 0.105 |

Three observations of substantive interest:

**(a) Tiferet → Thalamus is the most striking match.** Tiferet is the central sephira in the Tree of Life, with the highest connectivity and centrality (degree-z = +2.65 in the joined-trees graph, BC percentile 0.90, CC percentile 0.90). Its best brain analog is the **Right Thalamus** (degree 45, BC percentile 0.99 — the dominant subcortical hub in the connectome), at distance 0.166 — the second-tightest match in the entire mapping, despite Tiferet's outlier-high signature. The thalamus is independently described in network neuroscience as the brain's central relay (Sherman & Guillery 2006, Bullmore & Sporns 2012). The match between the structurally central sephira and the structurally central subcortical structure is non-trivial.

**(b) Daath → Left Amygdala is the worst match in the table.** Daath has signature ⟨+0.77, 1.00, 1.00, articulation=1, isolation=0.50⟩, of which the last property — that removing it isolates 50% of the graph — has no real-brain analog. The Left Amygdala wins by being the only subcortical articulation point with high BC, but at distance 0.753 it is meaningfully worse than every other sephira's match. **Daath, structurally, is the only tree node whose role real brains do not approximate well.** This is the sharpest version of the §3.5 finding that real human brains have distributed integration rather than a single bridge.

**(c) Symmetric tree pairs map to identical brain regions.** Chokhmah/Binah, Chesed/Geburah, and Netzach/Hod each map to the same brain region — because in the joined-trees graph these pairs have identical topological signatures (the tree has reflection symmetry across the middle pillar). This is an artifact of the equal-weighting of left-right pairs in the underlying graph, and would be broken by hemispheric tie-breakers (preferring left-pillar→left-hemisphere, right-pillar→right-hemisphere). We do not impose this constraint because the underlying graph does not.

**Suggestive but unverified anatomical correspondences:**

- *Tiferet (centre/beauty) ↔ Thalamus (universal relay)* — central integrator in both systems.
- *Daath (knowledge/abyss) ↔ Amygdala (emotional/somatic integration)* — the Hebrew root *yada* (יָדַע) connotes intimate experiential knowing.
- *Netzach/Hod ↔ Precuneus* — the precuneus is the major hub of the brain's Default Mode Network (Andrews-Hanna et al. 2014; Buckner & DiNicola 2019), associated with self-referential thought and integration of episodic memory.
- *Keter ↔ Precuneus / Malkuth ↔ Cuneus* — top and bottom of the tree both map to occipital-parietal midline cortex (related cuneal regions), reflecting their shared low-degree peripheral position in the graph.
- *Yesod (foundation) ↔ Inferior Parietal Cortex* — the inferior parietal lobule subserves body schema and somatospatial grounding (de Vignemont 2018).

We emphasize that these are post-hoc interpretations of a topological matching procedure. The matches are real (the distance numbers are reproducible), but the symbolic correspondences (e.g., "Tiferet means *balance*, the thalamus is *central*, therefore the mapping is meaningful") are interpretive overlay, not deductive conclusion. A different role-vector definition or weighting would yield different matches.

### 3.9 Critical test: vs random modular graphs

A skeptic could observe that any small-world graph with similar modularity might match the brain equally well, in which case the §3.5 result would reduce to "modular topology matches modular topology" rather than anything specific about the joined-trees graph. We tested this directly: 595 random modular graphs were generated via the stochastic block model (SBM) at six parameterizations (2, 3, 4, 6, 8 blocks; tight 4-block with within/between density ratio 6:1), each matched to the real coarse connectome's node and edge counts. Distances were computed under the same invariant metric.

**Table 7.** Joined-trees vs random modular graphs (lower distance = closer to real brain).

| Null model | n | mean d | min d | p (joined ≤ null) |
|---|---:|---:|---:|---:|
| SBM 2 blocks | 99 | 0.628 | 0.477 | **0.000** |
| SBM 3 blocks | 99 | 0.633 | 0.430 | **0.000** |
| SBM 4 blocks | 100 | 0.648 | 0.539 | **0.000** |
| SBM 6 blocks | 100 | 0.654 | 0.571 | **0.000** |
| SBM 8 blocks | 99 | 0.647 | 0.593 | **0.000** |
| Tight 4-block SBM (p_in/p_out = 6) | 198 | 0.590 | 0.397 | **0.000** |
| **Joined trees baseline** | — | **0.177** | — | — |

All 595 random modular graphs are farther from the real brain than the joined-trees graph. Even the *single best* modular random graph (d = 0.397) is more than twice the distance of joined-trees (d = 0.177). The result of §3.5 is therefore not a generic "modular structure matches modular structure" finding — the joined-trees graph encodes specific topological features beyond bare modularity that align with real human cortical organization. See [figures/fig8_modular_null.png](figures/fig8_modular_null.png).

### 3.10 Individual-subject replication: a scale-mismatch caveat

We attempted a per-subject replication on 50 individual AAL-parcellated DTI connectomes from the Brain Network Universe (Open Connectome Project; subjects from BNU, NKI-ENH, MRN, SWU, Jung 2015, HNU collections). Result: **0 of 50** subjects showed the joined-trees-beats-random pattern. Mean d_joined = 0.742 vs mean d_random = 0.588 across subjects.

This negative result requires careful interpretation. The individual subjects are at full AAL parcellation (~116 nodes each); the joined-trees graph has 21 nodes. The invariant-distance metric is **size-biased** — properties like density, mean degree, diameter, and spectral radius scale with node count, so comparing a 21-node graph to a 116-node graph is unfair on these dimensions. Random graphs of size 116 sit closer to other 116-node graphs than the 21-node Trees do, regardless of any topological similarity.

A scale-matched individual-subject test would require coarsening each subject's connectome to ~21 nodes using a known atlas labeling. The individual subject files in our source distribute graphml node IDs but not anatomical labels, so the standard AAL→lobe mapping cannot be applied directly. We flag this as a **methodological limitation** rather than a refutation of the §3.5 finding: the consensus connectome at coarse parcellation supports the hypothesis (p = 0.000); the per-subject test at native parcellation fails on size grounds; the per-subject test at matched coarse parcellation is not yet performed.

**Future-work priority:** download AAL atlas region labels (e.g., from the Tzourio-Mazoyer 2002 atlas distributions or the Brainnetome project) and re-run the per-subject test at matched coarse parcellation. Until that is done, the §3.5 result stands as a coarse-scale finding that has not been replicated at fine scale or per-subject.

### 3.11 Revised summary

| hypothesis | against synthetic brain | against real connectome |
|---|---|---|
| H1: aggregate structure brain-like | REJECTED (p > 0.85) | **SUPPORTED** (p = 0.000 against all 5 nulls × 3 thresholds) |
| H2: literal one-node Daath = corpus callosum | SUPPORTED at coarse synthetic | NOT TESTABLE (CC not a node in real parcellation) |
| H2′: distributed bridge-set integration | (not tested) | SUPPORTED (top-20 high-betweenness nodes are real integrators) |
| H3: H2 persists at finer scale | REJECTED | REJECTED (no single bridge in real connectome) |
| H4 (new): Tiferet's role = Thalamus's role | (n/a) | SUPPORTED (d = 0.166, both are central high-degree integrators) |
| H5 (new): Each sephira maps to a unique brain region | (n/a) | PARTIAL (symmetric tree pairs collapse to same region; 8 distinct matches across 11 sephirot) |
| H6 (critical): Joined-trees beats RANDOM MODULAR graphs | (n/a) | **SUPPORTED** (p = 0.000 across 595 SBM nulls in 6 parameterizations) |
| H7 (replication): Per-subject scale-matched replication | (n/a) | NOT YET TESTED — fine-scale per-subject test failed on size mismatch grounds (21 vs 116 nodes), proper test pending AAL atlas integration |

---

## 4. Discussion

### 4.1 What the H1 result means after testing against real data

The original version of this paper rejected H1 against synthetic brain models. Testing against the empirical Budapest Reference Connectome reverses this verdict: the joined-trees graph is in fact *closer* to coarse human cortical topology than every random graph drawn from any of the five null model classes, robust across edge-occurrence thresholds (Table 4). The metrics that align — clustering, transitivity, modularity, max-betweenness centralization, small-world σ — are precisely the small-world-with-modules signature that Bullmore & Sporns (2009, 2012) identify as canonical for brain networks. The metrics that diverge (spectral radius, efficiency at intermediate scales) reflect features the Trees do not capture: high-degree subcortical hubs and the rich-club organization that gives real brains their bypass paths.

This finding requires careful interpretation. We do **not** conclude that the Trees were "designed to mirror the brain," nor that any sephira corresponds to any specific brain region (the §3.5 H2 results explicitly rule out the literal one-node Daath/corpus-callosum mapping). What the result does support is the weaker claim that **the gross topological form of the joined-trees graph captures aggregate small-world-with-modules statistics that are characteristic of real human cortex at coarse parcellation, and that random graphs of the same size do not match.** Many small-world modular graphs would also satisfy this — the Trees are not uniquely brain-like — but they are not arbitrary either, and the popular hypothesis is therefore on firmer ground than our previous-version analysis suggested.

Why did the synthetic-brain analysis (§3.2) reach the opposite conclusion? Most likely because our BRAIN_10 model over-emphasized lobar modularity and a single-node corpus callosum bridge, producing a graph with structurally idiosyncratic properties (very high max-betweenness at the CC node, lower clustering than real cortex). Real brains achieve their topology through dense local cortico-cortical connectivity plus high-degree subcortical hubs, neither of which the synthetic model fully captured. This is a methodological lesson about synthetic brain models in network neuroscience generally (cf. Korhonen et al. 2017; Ryyppö et al. 2017): *brain-likeness depends sharply on what you call a brain.*

The remaining caveats: (a) absolute distances are not zero (d ≈ 0.12–0.24), so this is "less unlike than random" not "indistinguishable from"; (b) at fine parcellation (the full 801-node connectome) we did not run the full invariant test because the relevant comparison is to a graph with vastly more nodes than the Trees, making invariant comparison unfair on density grounds; (c) the result is for the Kircher arrangement only, not other historical Tree variants.

### 4.2 What the partial result on H2 means

The match between Daath's role in the joined-trees graph and the corpus callosum's role in BRAIN_10 is real and not statistically trivial. Both nodes are the *single highest-betweenness* node in their graph; both are articulation points; both partition the graph into exactly two equal halves on removal. This is the qualitative behavior any bridge between two semi-disjoint systems must have, and Daath was placed in the position to play this role by Hermetic tradition — in a system constructed centuries before graph theory existed. The convergence is worth noting even if it does not extend.

There are at least three readings:

1. **Trivial.** Any structure with a single bottleneck node connecting two halves will exhibit these properties. Daath was *defined* by the tradition as a bridge between the supernal and the lower tree (and, in Crowley, between the Tree of Life and Qliphoth); inserting it as a shared node and finding it has bridge properties is tautological.

2. **Modest.** The match is real but says only that the popular hypothesis got the *role* right (an integrator, a bottleneck) while getting the surrounding structure wrong. This is the reading we favor.

3. **Strong.** The convergence reflects something deeper about how systems that integrate two semi-autonomous subsystems must be organized. We would not push this reading without independent evidence from other domains.

The empirical literature constrains how literally to read the H2 result. **Roland et al. (2017)** found that interhemispheric functional connectivity remains nearly intact in patients with corpus callosotomy or agenesis — i.e., the corpus callosum is dominant but not strictly necessary for inter-hemispheric coordination. Anterior commissure, indirect cortico-subcortical-cortical loops, and (in agenesis cases) developmental rewiring through Probst bundles all provide redundancy. **Owen et al. (2012)** report substantial connectome reorganization in AgCC brains but not catastrophic disconnection. **Doron, Bassett & Gazzaniga (2012)** characterize inter-hemispheric coordination as having its own dynamic network signature distinct from intra-hemispheric coordination. Taken together, this body of work agrees with our BRAIN_34 result that the corpus callosum is the rank-1 betweenness node but not a strict articulation point at fine scale — and therefore against reading (3). The canonical Kabbalistic diagrams provide no analog of the secondary commissures or indirect compensatory pathways that real brains exploit when the primary integrator is compromised.

A symmetric reading: our finding is consistent with the *converse* claim — that the popular hypothesis is structurally inaccurate but functionally suggestive. The Daath/CC role correspondence captures something real about *bridge-node organization in systems that integrate two semi-autonomous subsystems*; this is unsurprising on graph-theoretic grounds and is documented in real brains by the literature above. The structural details (clustering, modularity, multi-scale persistence) are not captured by the Trees and would not be reconstructible from them.

### 4.3 Why this is not the whole story

A graph-theoretic test deliberately discards everything except topology: node identities, semantics, and dynamics are all stripped. For the popular hypothesis the user proposed, however, semantics may be the entire claim — the assertion may be that *Tiferet "is" the anterior cingulate* in some functional or symbolic sense, not that the *graph* of Tiferet's neighbors matches the graph of the ACC's neighbors. We make no effort to evaluate the semantic claim. The framework most equipped to take such a claim seriously without reducing it to neuroanatomy is the depth-psychological one (Jung 1944; Edinger 1985), where Kabbalistic structures are read as maps of psychological function. This paper has nothing to add to that tradition.

### 4.4 Limitations

(1) Both brain models are synthetic. Empirical connectomes (HCP, Open Connectome Project) at matched scale would strengthen the comparison. Our BRAIN_34 model is parameterized to published statistics but is not a real connectome.

(2) The 22-path Kircher arrangement is one of several historical Tree-of-Life variants. Alternative configurations (Cordoveran, Lurianic, Athanasian, Halevi's "Way of Kabbalah") differ in path layout. We tested only the Hermetic/Golden Dawn standard.

(3) Daath's path placement (to Chokhmah, Binah, Tiferet) is from the Crowley/Regardie tradition. Other systems place Daath differently or omit it entirely.

(4) The brain ensemble explores topological variation; it does not capture between-subject anatomical variation, sex differences, age effects, or hemispheric asymmetries (e.g., Yakovlevian torque).

(5) We did not test directed-graph or weighted-graph variants. The Trees are sometimes drawn with directed flows (the "Lightning Flash"); brain networks are weighted by tract strength. Both extensions are straightforward in principle but not done here.

(6) We did not apply control-theoretic centrality measures (Sato & Kawamura 2024) to the bridge-node test in §3.3. VCS and AECS may give a different verdict on the Daath/CC correspondence than betweenness centrality, since they measure a node's leverage over network *dynamics* rather than its position on shortest paths. This is the most defensible single extension a follow-up could make.

(7) Synthetic ROI definitions miss the heterogeneity of empirical parcellations (Korhonen et al. 2017; Ryyppö et al. 2017). Real fMRI ROIs vary in functional consistency, which would change degree distributions and centrality estimates in BRAIN_34. Replication with empirical parcellations would tighten or loosen the H2 result correspondingly.

### 4.5 Conclusion

The verdict shifts depending on what "brain" is operationalized as. Against synthetic models of brain hemispheric structure, the Tree of Life / Tree of Death joined graph fails to match. Against real human cortical connectivity from the Budapest Reference Connectome, the same joined graph achieves better aggregate-statistic agreement than any of 1,000 random null graphs across five null-model classes and three edge-occurrence thresholds (p = 0.000 throughout). The literal node-by-node correspondences (Daath = corpus callosum; sephirot = anatomical regions) are not supported. The aggregate topological form — small-world with modules and a moderate-betweenness integrating subgraph — is supported.

We refuse to read this as evidence that mystical traditions encode neuroanatomy. We do read it as evidence that the joined-trees graph is *not arbitrary* with respect to brain topology, in a quantitative sense that random graphs do not satisfy. Whether this reflects the convergent properties of any small-world modular graph (the most parsimonious reading), or something more specific, is beyond what this analysis can establish. The next defensible step is a multi-connectome replication (HCP individual subjects, Allen, Brainnetome) and an evaluation against larger reference sets of small-world modular graphs that are *not* mythological in origin, to determine whether the joined-trees graph performs better, equal, or worse than other arbitrary modular structures of the same size.

---

## Appendix A. Flower-of-Life geometric variant

The popular hypothesis sometimes includes the auxiliary claim that the joined hemispheres *are* the Flower of Life (Melchizedek 1990). The Flower of Life is the planar pattern of 19 overlapping unit circles centered on a hexagonal lattice, with symmetry group D₆ (order 12: six rotations, six mirror axes). The brain has bilateral symmetry only — symmetry group D₁ ≅ Z₂ (order 2). A 60° rotation does not map the brain to itself; the symmetry groups are not isomorphic; the geometric claim fails at the level of symmetry alone, before considering folding patterns, gyrification, or the Yakovlevian torque (Toga & Thompson 2003). Figure A1 (see [figures/flower_vs_brain.png](../flower_vs_brain.png)) shows the overlay.

A nuance worth recording: not all "sacred geometry in the brain" claims are spurious. **Mahakul & Agarwal (2021)** document a genuine pentagon-shape and golden-ratio relationship in the Circle of Willis — the arterial polygon at the base of the brain — published in *World Neurosurgery*, a peer-reviewed neurosurgery journal. The Circle of Willis genuinely exhibits an approximately regular pentagonal arrangement, and the ratios between vessel segments approximate φ in a measurable subset of subjects. This does not validate the Flower of Life claim, which makes a different and stronger geometric prediction (D₆ symmetry of the entire brain), but it weakens the otherwise-tempting blanket dismissal that "no sacred-geometric ratios appear in brain anatomy." Some do; the Flower of Life specifically does not.

## Appendix B. Reproducibility

```
paper/
├── graphs.py              # graph constructors
├── metrics.py             # invariants and similarity measures
├── nulls.py               # null model generators
├── run_analysis.py        # pipeline entry point
├── data/
│   ├── metrics.csv
│   ├── brain_ensemble_metrics.csv
│   ├── null_distances.csv
│   ├── node_role_comparison.csv
│   └── summary.json
└── figures/
    ├── fig1_structures.png
    ├── fig2_metrics_radar.png
    ├── fig3_null_distribution.png
    ├── fig4_daath_vs_callosum.png
    └── fig5_jacobs_ladder.png
```

All p-values are computed from k = 200 null draws per model. Increasing k tightens the distribution but does not change conclusions for the values reported here (none of which are at the threshold of significance).

## Funding and Conflicts

No funding. No conflicts of interest. The first author (Cody Churchwell) is the principal investigator; the computational analysis was performed by a large language model (Claude Opus 4.7, Anthropic) given the task description, with all decisions on framing, methodology, and interpretation made jointly.

## Data and Code Availability

All code and data in [`paper/`](.). Reproduction: `python3 run_analysis.py`. Runtime: ~3 minutes on a single CPU.

---

## References

Aboitiz, F., Scheibel, A. B., Fisher, R. S., & Zaidel, E. (1992). Fiber composition of the human corpus callosum. *Brain Research*, 598(1–2), 143–153.

Andrews-Hanna, J. R., Smallwood, J., & Spreng, R. N. (2014). The default network and self-generated thought: component processes, dynamic control, and clinical relevance. *Annals of the New York Academy of Sciences*, 1316(1), 29–52. doi:10.1111/nyas.12360.

Arzy, S., & Idel, M. (2015). *Kabbalah: A Neurocognitive Approach to Mystical Experiences*. Yale University Press. doi:10.12987/yale/9780300152364.001.0001.

Bagrow, J. P., & Bollt, E. M. (2019). An information-theoretic, all-scales approach to comparing networks. *Applied Network Science*, 4(1), 45.

Buckner, R. L., & DiNicola, L. M. (2019). The brain's default network: updated anatomy, physiology and evolving insights. *Nature Reviews Neuroscience*, 20(10), 593–608. doi:10.1038/s41583-019-0212-7.

Cacciatore, M., Magnani, F. G., Barbadoro, F., et al. (2025). Thalamus and consciousness: a systematic review on thalamic nuclei associated with consciousness. *Frontiers in Neurology*, 16, 1509668. doi:10.3389/fneur.2025.1509668.

Bassett, D. S., & Bullmore, E. T. (2017). Small-world brain networks revisited. *The Neuroscientist*, 23(5), 499–516.

Bassett, D. S., & Sporns, O. (2017). Network neuroscience. *Nature Neuroscience*, 20(3), 353–364.

Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, 2008(10), P10008.

Bullmore, E., & Sporns, O. (2009). Complex brain networks: graph theoretical analysis of structural and functional systems. *Nature Reviews Neuroscience*, 10(3), 186–198.

Bullmore, E., & Sporns, O. (2012). The economy of brain network organization. *Nature Reviews Neuroscience*, 13(5), 336–349.

de Vignemont, F. (2018). *Mind the Body: An Exploration of Bodily Self-Awareness*. Oxford University Press.

Catani, M., & Thiebaut de Schotten, M. (2008). A diffusion tensor imaging tractography atlas for virtual in vivo dissections. *Cortex*, 44(8), 1105–1132.

Crowley, A. (1909/1973). *777 and Other Qabalistic Writings*. Weiser.

Crowley, A. (1929/1976). *Magick in Theory and Practice*. Castle Books.

Doron, K. W., Bassett, D. S., & Gazzaniga, M. S. (2012). Dynamic network structure of interhemispheric coordination. *Proceedings of the National Academy of Sciences*, 109(46), 18661–18668. doi:10.1073/pnas.1216402109.

Drob, S. L. (1999). Jung and the Kabbalah. *History of Psychology*, 2(2), 102–118. doi:10.1037/1093-4510.2.2.102.

Edinger, E. F. (1985). *Anatomy of the Psyche: Alchemical Symbolism in Psychotherapy*. Open Court.

Halevi, Z. ben S. (1972). *Tree of Life: An Introduction to the Cabala*. Rider.

Hanegraaff, W. J. (Ed.) (2006). *Dictionary of Gnosis and Western Esotericism*. Brill.

Hugdahl, K., & Westerhausen, R. (2010). *The Two Halves of the Brain: Information Processing in the Cerebral Hemispheres*. MIT Press.

Humphries, M. D., & Gurney, K. (2008). Network 'small-world-ness': a quantitative method for determining canonical network equivalence. *PLoS ONE*, 3(4), e0002051.

Idel, M. (1988). *Kabbalah: New Perspectives*. Yale University Press.

Idel, M. (2006). *Kabbalah and Eros*. Yale University Press.

Joseph, S. M. (2007). Jung and Kabbalah: imaginal and noetic aspects. *Journal of Analytical Psychology*, 52(3), 321–341. doi:10.1111/j.1468-5922.2007.00665.x.

Jung, C. G. (1944). *Psychology and Alchemy*. Princeton University Press.

Korhonen, O., Saarimäki, H., Glerean, E., Sams, M., & Saramäki, J. (2017). Consistency of regions of interest as nodes of functional brain networks measured by fMRI. arXiv:1704.07635 [q-bio.NC].

Lancaster, B. L. (2011). The hard problem revisited: From cognitive neuroscience to Kabbalah and back again. In H. Walach, S. Schmidt, & W. B. Jonas (eds.), *Neuroscience, Consciousness and Spirituality*, vol. 1 (pp. 229–251). Springer. doi:10.1007/978-94-007-2079-4_14.

Mahakul, D. J., & Agarwal, J. (2021). Pentagon inside the Circle of Willis and the golden ratio. *World Neurosurgery*, 156, 76–80. doi:10.1016/j.wneu.2021.09.006.

McGilchrist, I. (2009). *The Master and His Emissary: The Divided Brain and the Making of the Western World*. Yale University Press.

McGilchrist, I. (2021). *The Matter With Things: Our Brains, Our Delusions, and the Unmaking of the World*. Perspectiva Press.

Melchizedek, D. (1990). *The Ancient Secret of the Flower of Life*. Light Technology.

Newman, M. E. J. (2006). Modularity and community structure in networks. *PNAS*, 103(23), 8577–8582.

Owen, J. P., Li, Y.-O., Ziv, E., Strominger, Z., Gold, J., Bukshpun, P., Wakahiro, M., Friedman, E. J., Sherr, E. H., & Mukherjee, P. (2013). The structural connectome of the human brain in agenesis of the corpus callosum. *NeuroImage*, 70, 340–355. doi:10.1016/j.neuroimage.2012.12.031.

Plüss, R., Villota, H., & Orio, P. (2025). Hemispheric-specific coupling improves modeling of functional connectivity using Wilson–Cowan dynamics. arXiv:2506.22951 [q-bio.NC].

Podschun, A. N., Betzel, R. F., Braun, U., et al. (2026). Exploring the role of the rich club in network control of neurocognitive states. *Human Brain Mapping*, 47(1). doi:10.1002/hbm.70485.

Regardie, I. (1932). *A Garden of Pomegranates*. Aries Press.

Regardie, I. (1971). *The Golden Dawn*. Llewellyn.

Roland, J. L., Snyder, A. Z., Hacker, C. D., Mitra, A., Shimony, J. S., Limbrick, D. D., Raichle, M. E., Smyth, M. D., & Leuthardt, E. C. (2017). On the role of the corpus callosum in interhemispheric functional connectivity in humans. *Proceedings of the National Academy of Sciences*, 114(50), 13278–13283. doi:10.1073/pnas.1707050114.

Ryyppö, E., Glerean, E., Brattico, E., & Saramäki, J. (2017). Regions of interest as nodes of dynamic functional brain networks. arXiv:1710.04056 [q-bio.NC].

Sato, K., & Kawamura, R. (2024). Uniqueness analysis of controllability scores and their application to brain networks. arXiv:2408.03023 [math.OC].

Sherman, S. M., & Guillery, R. W. (2006). *Exploring the Thalamus and Its Role in Cortical Function* (2nd ed.). MIT Press.

Tian, Z., Song, J., Zhao, X., et al. (2024). The interhemispheric amygdala–accumbens circuit encodes negative valence in mice. *Science*, 386(6724), 1092–1099. doi:10.1126/science.adp7520.

Simpson, S. L., & Laurienti, P. J. (2016). Disentangling brain graphs: A note on the conflation of network and connectivity analyses. arXiv:1602.00933 [q-bio.QM].

Sporns, O. (2011). *Networks of the Brain*. MIT Press.

Szalkai, B., Kerepesi, C., Varga, B., & Grolmusz, V. (2015). The Budapest Reference Connectome Server v2.0. *Neuroscience Letters*, 595, 60–62. doi:10.1016/j.neulet.2015.03.071.

Young, M. P. (1993). The organization of neural systems in the primate cerebral cortex. *Proceedings of the Royal Society of London B*, 252(1333), 13–18.

Telesford, Q. K., Joyce, K. E., Hayasaka, S., Burdette, J. H., & Laurienti, P. J. (2011). The ubiquity of small-world networks. *Brain Connectivity*, 1(5), 367–375.

Toga, A. W., & Thompson, P. M. (2003). Mapping brain asymmetry. *Nature Reviews Neuroscience*, 4(1), 37–48.

van den Heuvel, M. P., & Sporns, O. (2011). Rich-club organization of the human connectome. *Journal of Neuroscience*, 31(44), 15775–15786.

Wakana, S., Caprihan, A., Panzenboeck, M. M., et al. (2007). Reproducibility of quantitative tractography methods applied to cerebral white matter. *NeuroImage*, 36(3), 630–644.

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of 'small-world' networks. *Nature*, 393(6684), 440–442.
