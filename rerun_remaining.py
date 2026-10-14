"""Re-run the remaining pipeline steps (6-8) into results_v2/data.

Step 6 (identifiability) now uses the log-space fitting rewrite;
step 7 (bifurcation) includes the sustained-infusion bistability demo.
Run:  python rerun_remaining.py
"""
from __future__ import annotations

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from ad_graph_model.config import load_config
from ad_graph_model.experiments import (exp_bifurcation, exp_identifiability,
                                        exp_spatial)
from ad_graph_model.io_utils import ensure_dir


def main() -> None:
    cfg = load_config(PROJECT_ROOT / "configs" / "v2_baseline.json")
    data_dir = ensure_dir(PROJECT_ROOT / "results_v2" / "data")

    t0 = time.time()
    print("[6/8] identifiability (log-space fit + profile CI) ...", flush=True)
    exp_identifiability(cfg, data_dir)
    print(f"      identifiability done in {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    print("[7/8] bifurcation ...", flush=True)
    exp_bifurcation(cfg, data_dir)
    print(f"      bifurcation done in {time.time() - t0:.0f}s", flush=True)

    t0 = time.time()
    print("[8/8] spatial propagation ...", flush=True)
    exp_spatial(cfg, data_dir)
    print(f"      spatial done in {time.time() - t0:.0f}s", flush=True)
    print("ALL DONE", flush=True)


if __name__ == "__main__":
    main()
