# Trees vs Brain — replication package

Code, data, and figures for *Topological Correspondence Between Kabbalistic Tree Diagrams and Cerebral Hemispheric Network Structure*.

## Quick start

```bash
pip install -r requirements.txt
python3 run_analysis.py
```

Runtime: about 3 minutes on a single CPU. Outputs go to `data/` (CSV, JSON) and `figures/` (PNG).

## What the analysis tests

Three claims, drawn from the popular synthesis of Hermetic Kabbalah and neuroscience:

| H | Claim | Result |
|---|---|---|
| H1 | The Tree of Life and Tree of Death together resemble the brain in structure | Rejected (p > 0.85 against four of five null models) |
| H2 | Daath, modeled as a shared bridge node, plays the role of the corpus callosum | Supported at coarse scale (3/3 qualitative properties match) |
| H3 | H2 persists at finer brain parcellation | Rejected (rich-club bypasses callosum at scale) |

Full discussion in `paper.md`.

## File layout

```
paper/
├── paper.md                     # the writeup
├── README.md                    # this file
├── requirements.txt             # Python dependencies
│
├── graphs.py                    # Tree of Life, Qliphoth, brain-model constructors
├── metrics.py                   # 17 graph invariants + 3 similarity measures
├── nulls.py                     # 5 null model generators
├── run_analysis.py              # pipeline entry point
│
├── data/
│   ├── metrics.csv              # invariants for every study graph
│   ├── brain_ensemble_metrics.csv
│   ├── null_distances.csv       # 5,901 null comparisons
│   ├── node_role_comparison.csv
│   └── summary.json             # all numbers cited in the paper
│
└── figures/
    ├── fig1_structures.png      # the four principal graphs
    ├── fig2_metrics_radar.png   # normalized invariants per structure
    ├── fig3_null_distribution.png   # 6 panels, candidate vs nulls
    ├── fig4_daath_vs_callosum.png   # bridge-node side by side
    └── fig5_jacobs_ladder.png   # the four-worlds variant
```

## Reproducing every number in the paper

`data/summary.json` is the canonical source for every numeric result cited. To
regenerate from scratch:

```bash
rm -rf data/ figures/
python3 run_analysis.py
```

All seeds are fixed; output is deterministic. The brain ensemble uses seeds 0..99 for the 34-region model; null models use seed 42 for the master RNG.

## Tuning the analysis

Two constants in `run_analysis.py` control runtime/resolution:

```python
K_NULLS = 200       # null draws per model per candidate
K_BRAIN = 100       # brain-model realizations
```

Doubling `K_NULLS` to 400 tightens the p-value distribution but does not change conclusions for any candidate (none are at the threshold of significance). Increasing `K_BRAIN` reduces the standard errors on the corpus-callosum-role estimates.

## Extending the analysis

The framework supports drop-in extensions:

* **New Kabbalistic variants:** add to `graphs.py`. Cordoveran, Lurianic, and Athanasian Tree-of-Life arrangements differ in path layout — easy to encode.
* **Empirical connectomes:** replace `brain_model_34` with a loader for HCP / Open Connectome data. Match scale and re-run.
* **Directed/weighted variants:** swap `nx.Graph` for `nx.DiGraph` and adjust the metrics module.
* **Additional null models:** add to `nulls.py`.

## Caveats

This is a structural/topological test only. The paper's negative result on H1 does not address whether Kabbalistic structures function as useful symbolic or contemplative maps of psychological function. For that question, see Jung (1944), Edinger (1985), and McGilchrist (2009).

## License

Code and data: CC0. Use freely.

## Citation

If you use this code or data, cite the paper (`paper.md`). To cite this exact replication package, include the full directory.
