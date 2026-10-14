"""Rapid-infusion monotherapy scan at t = 0 onset (long horizon).

Verifies the Section 4.4 claim that monotherapy can durably arrest the
cascade only if the infusion is rapid enough for the occupancy to cross
the full-system fold (theta_SN = 0.93) before the seed region commits
to the diseased basin.  All runs: 68-region surrogate connectome,
intervention at t = 0, t_end = 200, reference tolerances.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from ad_graph_model.config import load_config
from ad_graph_model.experiments import theta_of_P
from ad_graph_model.simulate import run_simulation


def main() -> None:
    cfg = load_config(ROOT / "configs" / "v2_baseline.json")
    p = cfg["parameters"]
    out = ROOT / "results_v2" / "data"
    t0 = time.time()
    rows = []
    for u in [0.05, 0.1, 0.2, 0.3, 0.5]:
        c = json.loads(json.dumps(cfg))
        c["interventions"] = [{"time": 0.0, "set": {"uP": float(u)}}]
        c.setdefault("simulation", {})["t_end"] = 200.0
        res = run_simulation(c)
        P_ss = float(u) / p["cP"]
        rows.append({
            "uP": float(u),
            "P_ss": P_ss,
            "theta_ss": theta_of_P(P_ss, p),
            "mean_a_end40": float(res.a[int(np.argmin(np.abs(res.t - 40)))].mean()),
            "mean_a_end": float(res.a[-1].mean()),
        })
        print(f"  uP={u:.3g} theta*={rows[-1]['theta_ss']:.3f} "
              f"a40={rows[-1]['mean_a_end40']:.3f} "
              f"a200={rows[-1]['mean_a_end']:.3f} ({time.time()-t0:.0f}s)",
              flush=True)
    pd.DataFrame(rows).to_csv(out / "monotherapy_rapid_t0.csv", index=False)
    print("saved", out / "monotherapy_rapid_t0.csv")


if __name__ == "__main__":
    main()
