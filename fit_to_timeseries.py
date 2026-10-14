from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from ad_graph_model.analysis import timeseries_dataframe
from ad_graph_model.config import load_config
from ad_graph_model.io_utils import ensure_dir, save_dataframe, save_json
from ad_graph_model.plotting import save_fit_comparison
from ad_graph_model.simulate import run_simulation


def interp_sim(sim_df: pd.DataFrame, target_time: np.ndarray, col: str) -> np.ndarray:
    return np.interp(target_time, sim_df["time"].to_numpy(), sim_df[col].to_numpy())


def main() -> None:
    parser = argparse.ArgumentParser(description="Fit selected model parameters to target timeseries CSV.")
    parser.add_argument("--config", required=True)
    parser.add_argument("--targets", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    target_df = pd.read_csv(args.targets)
    outdir = ensure_dir(args.outdir)

    fit_cfg = cfg["fit"]
    param_names: List[str] = fit_cfg["parameters"]
    lower = np.array([fit_cfg["bounds"][name][0] for name in param_names], dtype=float)
    upper = np.array([fit_cfg["bounds"][name][1] for name in param_names], dtype=float)
    x0 = np.array([cfg["parameters"][name] for name in param_names], dtype=float)

    target_cols = [c for c in target_df.columns if c != "time"]
    history = []

    def residuals(x: np.ndarray) -> np.ndarray:
        run_cfg = json.loads(json.dumps(cfg))
        for name, value in zip(param_names, x):
            run_cfg["parameters"][name] = float(value)
        res = run_simulation(run_cfg)
        sim_df = timeseries_dataframe(res)
        pieces = []
        loss = 0.0
        for col in target_cols:
            sim_col = col.replace("_target", "")
            sim_interp = interp_sim(sim_df, target_df["time"].to_numpy(), sim_col)
            target_vals = target_df[col].to_numpy(dtype=float)
            scale = max(np.max(np.abs(target_vals)), 1e-6)
            r = (sim_interp - target_vals) / scale
            loss += float(np.mean(r ** 2))
            pieces.append(r)
        history.append({"iter": len(history), "loss": loss, **{p: float(v) for p, v in zip(param_names, x)}})
        return np.concatenate(pieces)

    result = least_squares(residuals, x0=x0, bounds=(lower, upper), verbose=2, max_nfev=int(fit_cfg.get("max_nfev", 40)))

    best_cfg = json.loads(json.dumps(cfg))
    for name, value in zip(param_names, result.x):
        best_cfg["parameters"][name] = float(value)
    best_res = run_simulation(best_cfg)
    best_df = timeseries_dataframe(best_res)

    best_params = {name: float(value) for name, value in zip(param_names, result.x)}
    save_json(best_params, outdir / "best_fit_params.json")
    save_json({"success": bool(result.success), "cost": float(result.cost), "message": result.message}, outdir / "fit_report.json")
    save_dataframe(pd.DataFrame(history), outdir / "fit_loss_history.csv")
    save_dataframe(best_df, outdir / "best_fit_timeseries.csv")
    save_fit_comparison(target_df, best_df, outdir / "fit_comparison.png")

    print(f"Fit done. Outputs written to: {outdir}")
    print("Best-fit parameters:")
    for k, v in best_params.items():
        print(f"  {k}: {v:.6g}")


if __name__ == "__main__":
    main()
