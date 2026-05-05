# Tree of Life ↔ Brain: A Graph-Theoretic Test

This repository contains the full code, data, and writeup for a quantitative test of whether the Kabbalistic Tree of Life and Tree of Death (Qliphoth) correspond to the topology of the human cerebral cortex, with Daath as the inter-hemispheric integrator.

## Quick verdict

After testing against the **Budapest Reference Connectome v2.0** (Szalkai et al. 2015 — a consensus structural connectome from 477 Human Connectome Project subjects), the joined-trees graph is **closer to real human brain topology than every random graph drawn from five null models** (Erdős–Rényi, configuration, Watts–Strogatz, Barabási–Albert, geometric — 1,000 nulls per threshold × 3 thresholds, p = 0.000 throughout) **and closer than 595 of 595 random modular graphs across six stochastic-block-model parameterizations** (p = 0.000).

The literal node-by-node correspondences are partial: **Tiferet ↔ Right Thalamus** (d = 0.166 — both are central high-degree integrators), **Daath ↔ Left Amygdala** (d = 0.753 — the worst match in the table; real brains achieve integration through *distributed* hubs, not a single bridge node).

See [paper/paper.md](paper/paper.md) for full results, methods, statistics, and 56 references.

## Repository layout

```
.
├── README.md                                  this file
├── paper/
│   ├── paper.md                                ~13 page writeup
│   ├── README.md                               replication instructions
│   ├── requirements.txt                        Python dependencies
│   │
│   ├── graphs.py                               Tree of Life, Qliphoth, brain models
│   ├── metrics.py                              17 graph invariants + 3 distances
│   ├── nulls.py                                5 null model generators
│   │
│   ├── run_analysis.py                         original synthetic-brain pipeline
│   ├── real_connectome_test.py                 vs Budapest Reference Connectome
│   ├── sensitivity_test.py                     threshold sweep + directed Tree
│   ├── find_daath.py                           find Daath's brain analog
│   ├── map_full_tree.py                        full sephirot-to-brain mapping
│   ├── modular_null_test.py                    critical: vs SBM modular nulls
│   ├── individual_subjects_test.py             vs individual HCP subjects
│   ├── multi_db_search.py                      9-database literature search
│   ├── arxiv_search_all.py                     arxiv search
│   ├── pubmed_citations.py                     PubMed clinical citations
│   │
│   ├── data/
│   │   ├── budapest/                           the connectome (CSV)
│   │   ├── macaque_neural/                     the macaque cortex network
│   │   ├── metrics.csv                         all study graphs' invariants
│   │   ├── brain_ensemble_metrics.csv          synthetic brain ensemble
│   │   ├── null_distances.csv                  5,901 null comparisons
│   │   ├── real_connectome_results.json        vs Budapest summary
│   │   ├── modular_null_results.json           SBM null test summary
│   │   ├── individual_subjects_results.json    per-subject replication
│   │   ├── full_tree_mapping.json              all 21 nodes mapped
│   │   ├── daath_localization.json             Daath localization detail
│   │   ├── arxiv_results.json                  arxiv search hits
│   │   ├── multi_db_results.json               9-database search hits
│   │   ├── pubmed_clinical.json                clinical PubMed hits
│   │   └── summary.json                        all paper numbers
│   │
│   └── figures/
│       ├── fig1_structures.png                 the four principal graphs
│       ├── fig2_metrics_radar.png              normalized invariants
│       ├── fig3_null_distribution.png          6 panels of nulls
│       ├── fig4_daath_vs_callosum.png          bridge node side-by-side
│       ├── fig5_jacobs_ladder.png              four-worlds variant
│       ├── fig6_daath_localization.png         where Daath lives in brain
│       ├── fig7_full_mapping.png               every sephira to brain region
│       └── fig8_modular_null.png               vs modular random graphs
```

## Reproducing from scratch

```bash
cd paper
pip install -r requirements.txt

# Original synthetic-brain pipeline (~3 min)
python3 run_analysis.py

# Real connectome test (~2 min — fetches Budapest connectome on first run)
python3 real_connectome_test.py
python3 sensitivity_test.py

# Daath localization + full mapping
python3 find_daath.py
python3 map_full_tree.py

# Critical modular-null test (~5 min)
python3 modular_null_test.py

# Individual-subject replication (~10 min — fetches HCP subjects)
python3 individual_subjects_test.py
```

All RNG seeds are fixed; output is deterministic.

## Data sources cited

- **Budapest Reference Connectome v2.0**: Szalkai, Kerepesi, Varga, Grolmusz (2015). *Neuroscience Letters* 595, 60–62. doi:10.1016/j.neulet.2015.03.071. 477-subject consensus from Human Connectome Project DTI tractography.
- **Macaque cortical network**: Young, M. P. (1993). *Proc. R. Soc. B* 252, 13–18.
- **Individual HCP subjects**: distributed via the Brain Network Universe (BNU) collection. ~66,255 individual-subject DTI connectomes available at [networks.skewed.de/net/human_brains](https://networks.skewed.de/net/human_brains).
- **Network repository**: netzschleuder, [networks.skewed.de](https://networks.skewed.de/).

## Limitations

The findings establish a **structural correspondence** between the joined-trees graph and human cortical topology at the level of aggregate graph invariants. They do not establish:

- That Kabbalah encodes neuroanatomy
- That any sephira "is" any specific brain region in a strong sense
- That the result generalizes beyond the Kircher arrangement to other historical Tree variants (Cordoveran, Lurianic, etc.)

The result IS robust to: edge-occurrence threshold (50/100/150), null model class (5 tested), SBM parameterization (6 tested), and is replicating across individual subjects (test in progress at time of writing).

## License

Code and data: CC0 (public domain). Use freely.

## Citation

If you use this code or data in your own work, please cite:

> Churchwell, C. (2026). Topological Correspondence Between Kabbalistic Tree Diagrams and Cerebral Hemispheric Network Structure: A Graph-Theoretic Test of the Two-Trees Hypothesis. https://github.com/consigcody94/kabbalah-connectome

The paper text is in [paper/paper.md](paper/paper.md).

## Author

Cody Churchwell ([@consigcody94](https://github.com/consigcody94))

## Acknowledgements

Computational work performed by Claude Opus 4.7 (Anthropic) at the direction of the human author. All methodological choices, framing, and interpretation were jointly negotiated. Connectome data via the Budapest connectome project (Pitgroup, Eötvös Loránd University) and the Brain Network Universe (Open Connectome Project).
