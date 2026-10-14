"""Long-horizon (t=200) rescan of dose-response, therapeutic window and
the combination grid, to separate true arrest from critical slowing down.

The t=40 horizon used in the first pass conflates durable arrest with
slow transients towards the diseased state; the full-system saddle-node
(threshold) analysis shows the true boundary is theta ~ 0.93 at the
reference clearance, so intermediate doses only *slow* the cascade.
"""
from __future__ import annotations

import copy
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from ad_graph_model.analysis import summary_metrics
from ad_graph_model.config import load_config
from ad_graph_model.experiments import theta_of_P
from ad_graph_model.simulate import run_simulation

T_END = 200.0
N_PTS = 2001


def run(cfg, interventions, t_end=T_END):
    c = copy.deepcopy(cfg)
    c["simulation"]["t_end"] = t_end
    c["simulation"]["n_points"] = N_PTS
    if interventions:
        c["interventions"] = interventions
    return run_simulation(c)


def main() -> None:
    cfg = load_config(ROOT / "configs" / "v2_baseline.json")
    p = cfg["parameters"]
    out = ROOT / "results_v2" / "data"
    t0 = time.time()

    # ---- 1. dose-response (intervention at t=8), long horizon
    rows = []
    for u in np.logspace(-3.3, -0.7, 26):
        res = run(cfg, [{"time": 8.0, "set": {"uP": float(u)}}])
        m = summary_metrics(res)
        i40 = int(np.argmin(np.abs(res.t - 40)))
        m["uP"] = float(u)
        m["P_ss"] = float(u) / p["cP"]
        m["theta_ss"] = theta_of_P(float(u) / p["cP"], p)
        m["mean_a_end40"] = float(res.a[i40].mean())
        rows.append(m)
        print(f"  dose uP={u:.5g} th={m['theta_ss']:.3f} "
              f"a40={m['mean_a_end40']:.3f} a200={m['mean_a_end']:.3f} "
              f"({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(out / "dose_response_long.csv", index=False)

    # ---- 2. therapeutic window (3 strategies x onsets), long horizon
    strategies = {
        "Passivation": {"uP": 0.012},
        "Clearance": {"cY": 0.50, "cZ": 0.20},
        "Combination": {"uP": 0.012, "cY": 0.50, "cZ": 0.20},
    }
    rows = []
    for sname, sset in strategies.items():
        for t_on in np.arange(0.0, 33.0, 2.0):
            res = run(cfg, [{"time": float(t_on), "set": dict(sset)}])
            m = summary_metrics(res)
            i40 = int(np.argmin(np.abs(res.t - 40)))
            m.update({"strategy": sname, "onset": float(t_on),
                      "mean_a_end40": float(res.a[i40].mean())})
            rows.append(m)
        print(f"  window {sname} done ({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(out / "treatment_window_long.csv", index=False)

    # ---- 3. combination grid, long horizon
    cYs = np.linspace(0.25, 0.60, 10)
    uPs = np.logspace(-3.3, -1.0, 10)
    grid40 = np.zeros((len(uPs), len(cYs)))
    grid200 = np.zeros((len(uPs), len(cYs)))
    for i, u in enumerate(uPs):
        for j, cy in enumerate(cYs):
            res = run(cfg, [{"time": 8.0,
                             "set": {"uP": float(u), "cY": float(cy)}}])
            i40 = int(np.argmin(np.abs(res.t - 40)))
            grid40[i, j] = res.a[i40].mean()
            grid200[i, j] = res.a[-1].mean()
        print(f"  grid row {i+1}/10 done ({time.time()-t0:.0f}s)", flush=True)
    pd.DataFrame(grid200, index=[f"{u:.5g}" for u in uPs],
                 columns=[f"{c:.3g}" for c in cYs]).to_csv(
                     out / "combination_grid_long.csv")
    pd.DataFrame(grid40, index=[f"{u:.5g}" for u in uPs],
                 columns=[f"{c:.3g}" for c in cYs]).to_csv(
                     out / "combination_grid_t40.csv")
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
