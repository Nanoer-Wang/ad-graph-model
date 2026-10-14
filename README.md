# ad-graph-model

Continuous, network-resolved model of amyloid-cascade progression with an
explicit, passivatable fibril interface — code and analysis pipeline for the
manuscript

> **Fibril-interface passivation suppresses secondary nucleation and arrests
> the amyloid cascade: a continuous graph-based model of Alzheimer's disease
> progression and therapeutic design**

## Model

Each brain region (node of a weighted graph) carries six state variables:

| Symbol | Meaning |
|---|---|
| `X` | Aβ monomer burden |
| `Y` | toxic-oligomer burden |
| `Z` | fibril burden |
| `P` | passivator concentration |
| `th` | passivated-interface fraction θ ∈ [0, 1] |
| `a` | neuronal damage a ∈ [0, 1] |

Dynamics (per node `v`):

```
dX/dt = λ0(1−a)(1+βa) − cX·X − kn·X^m − ke·X·Z·(1−θ) + DX·(LX)_v
dY/dt = kn·X^m + ks·X·Z^q·(1−θ) − (cY+kf)·Y            + DY·(LY)_v
dZ/dt = kf·Y + ke·X·Z·(1−θ) − cZ·Z
dP/dt = uP − cP·P
dθ/dt = kon·P·(1−θ) − koff·θ
da/dt = ηY·Y^hY/(KY^hY+Y^hY)·(1−a) + ηN·(Wa)_v − δa·a
```

Key mechanism: **both** surface-catalysed processes (secondary nucleation and
elongation) are multiplied by the free-interface fraction (1−θ), so a
passivator occupying fibril interfaces reduces the aggregation efficiency.
Binding is Langmuir-type with KD = koff/kon. Linear stability of the
oligomer–fibril subsystem yields the analytic critical passivation fraction

θc = 1 − (cY + kf)·cZ / { Xs·[ks·kf + ke·(cY + kf)] },  Xs = λ0/cX

Three network families are supported: a 2-D cortical sheet (20×20), a
Watts–Strogatz small-world network (68 nodes), and a 68-region surrogate
structural connectome (Desikan-like labels, exponential distance rule,
hub reinforcement, spatial rewiring; density 0.215, clustering 0.54,
path length 2.13).

## Repository layout

```
src/ad_graph_model/     model, simulation, analysis and figure code
  model.py              ODE right-hand side, passivation mechanism
  graph.py              network construction (sheet / small-world / connectome)
  simulate.py           LSODA integration, intervention schedules
  analysis.py           outcome metrics
  bifurcation.py        steady states, saddle-node folds, pseudo-potential
  sensitivity.py        Latin hypercube sampling + PRCC
  fitting.py            synthetic-benchmark fitting, profile-likelihood CIs
  nondim.py             dimensionless groups and reference scales
  experiments.py        the eight experiments behind the manuscript
  figures_core.py       Figures 1–4 renderers
  figures_therapy.py    Figures 5–8 renderers
configs/                locked reference configuration
run_pipeline.py         runs all experiments (writes results_v2/data)
make_figures.py         renders all figures (reads results_v2/data)
render_equations.py     renders the manuscript equations to PNG
results_v2/             experiment data + figures (created by the pipeline)
```

## Installation

```bash
pip install -r requirements.txt
```

Python ≥ 3.10 with numpy, scipy, pandas, matplotlib. No compiled extensions.

## Reproducing the paper

```bash
# 1. run all eight experiments (~45 min on a laptop; writes results_v2/data)
python run_pipeline.py --outdir results_v2 --n-sens 400

# 2. render all figures (writes results_v2/figures, 600-dpi PNG + PDF)
python make_figures.py

# 3. render the manuscript equations (writes results_v2/equations)
python render_equations.py
```

All random procedures use fixed seeds; the pipeline is deterministic and
reproduces every figure and table from the reference configuration
(`configs/v2_baseline.json`) with the three commands above.

### Experiments

| # | Experiment | Key outputs |
|---|---|---|
| 1 | Baseline + topology robustness | `baseline_*`, `topology_*` |
| 2 | Intervention strategies | `strategy_*` |
| 3 | Passivation dose–response | `dose_response*.csv`, `dose_kinetics.csv` |
| 4 | Therapeutic window + combination | `treatment_window.csv`, `combination_grid.csv`, `clearance_sweep.csv` |
| 5 | Global sensitivity (LHS + PRCC) | `prcc_table.csv`, `lhs_samples.csv`, `lhs_ensemble.csv` |
| 6 | Parameter identifiability | `fit_result.json`, `profile_likelihood_ci.csv`, `parameter_correlation.csv` |
| 7 | Bistability and thresholds | `bifurcation_theta.csv`, `critical_curve.csv`, `potential_*.csv` |
| 8 | Spatial propagation | `spatial_region_time.csv`, `spatial_onsets.csv`, `spatial_strength_damage.csv` |

## Quick start (single simulation)

```python
import json
from ad_graph_model.simulate import run_simulation
from ad_graph_model.analysis import summary_metrics

cfg = json.load(open("configs/v2_baseline.json"))
res = run_simulation(cfg)                 # untreated baseline
print(summary_metrics(res))
```

## License

MIT License — see `LICENSE`.
