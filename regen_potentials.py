"""Regenerate the pseudo-potential CSVs with the segment-aware integrator.

The QSS reduction f(a) is undefined on "runaway gaps" where the fast
aggregation subsystem has no positive equilibrium; the potential is
integrated independently on each contiguous valid segment.
Run after rerun_remaining.py finishes (it used the pre-fix module).
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import numpy as np
import pandas as pd

from ad_graph_model.bifurcation import pseudo_potential, qss_damage_rate
from ad_graph_model.config import load_config
from ad_graph_model.experiments import P_of_theta
from ad_graph_model.io_utils import save_dataframe


def main() -> None:
    cfg = load_config(PROJECT_ROOT / "configs" / "v2_baseline.json")
    p = cfg["parameters"]
    data_dir = PROJECT_ROOT / "results_v2" / "data"
    a_grid = np.linspace(0.0, 1.0, 300)
    for label, th in [("none", 0.0), ("partial", 0.65), ("strong", 0.95)]:
        t0 = time.time()
        P = P_of_theta(th, p)
        f = np.array([qss_damage_rate(a, p, P) for a in a_grid])
        phi = pseudo_potential(a_grid, f)
        df = pd.DataFrame({"a": a_grid, "f": f, "phi": phi, "theta": th})
        save_dataframe(df, data_dir / f"potential_{label}.csv")
        print(f"potential_{label}: theta={th:.2f} "
              f"valid={int(np.isfinite(f).sum())}/300 ({time.time()-t0:.0f}s)",
              flush=True)


if __name__ == "__main__":
    main()
