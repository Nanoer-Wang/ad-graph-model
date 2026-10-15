"""Recompute profile-likelihood CIs, one parameter per process.

Run in small batches (e.g. ``--params kn ks``) because a single process
that evaluates all six parameters exhausts the memory available to the
sandbox: each profile point re-simulates the 400-node network, and the
integrator workspace accumulates faster than it is reclaimed.

Results are merged into the existing profile_likelihood_ci.csv, so
repeated invocations with different --params converge on the full table.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from ad_graph_model.config import load_config
from ad_graph_model.fitting import profile_likelihood_ci

ALL = ["kn", "ks", "ke", "cY", "etaY", "KY"]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", nargs="+", default=ALL, choices=ALL)
    ap.add_argument("--n-bisect", type=int, default=6)
    args = ap.parse_args()

    cfg = load_config(ROOT / "configs" / "v2_baseline.json")
    data_dir = ROOT / "results_v2" / "data"
    fit_json = json.load(open(data_dir / "fit_result.json"))
    target = pd.read_csv(data_dir / "synthetic_target.csv")

    p = cfg["parameters"]
    bounds = {n: (p[n] / 4.0, p[n] * 4.0) for n in ALL}
    x = np.array([fit_json["best"][n] for n in ALL])
    fit = {
        "names": ALL, "x": x,
        "sd": np.array([fit_json["sd"][n] for n in ALL]),
        "cov": np.full((6, 6), np.nan),
        "cost": fit_json["cost"],
        "columns": [c for c in target.columns if c != "time"],
        "n_obs": int(fit_json["n_obs"]),
    }

    out_csv = data_dir / "profile_likelihood_ci.csv"
    if out_csv.exists():
        table = pd.read_csv(out_csv).set_index("parameter")
    else:
        table = pd.DataFrame(
            {"best": x, "ci_lo": x, "ci_hi": x, "sd_linear": fit["sd"]},
            index=ALL)

    for name in args.params:
        i = ALL.index(name)
        sub = {k: (v if k != "x" else np.array([v[i]]))
               for k, v in fit.items()}
        sub["names"] = [name]
        # keep the per-parameter slices consistent so that the reported
        # linearised sd (and the covariance placeholder) match the single
        # free parameter of this profile
        sub["sd"] = np.array([fit["sd"][i]])
        sub["cov"] = np.full((1, 1), np.nan)
        b1 = {name: bounds[name]}
        t0 = time.time()
        # Always pass the complete fitted vector and bounds so that every
        # nuisance parameter is reoptimized. Slicing to one parameter makes
        # a conditional scan at reference nuisance values, not a profile.
        ci = profile_likelihood_ci(cfg, fit, bounds, target,
                                   n_bisect=args.n_bisect, profile_names=[name])
        row = ci.iloc[0]
        table.loc[name, ["best", "ci_lo", "ci_hi", "sd_linear"]] = [
            row["best"], row["ci_lo"], row["ci_hi"], row["sd_linear"]]
        rel = (row["ci_hi"] - row["ci_lo"]) / row["best"]
        print(f"{name}: best={row['best']:.6g} "
              f"CI=[{row['ci_lo']:.6g}, {row['ci_hi']:.6g}] "
              f"rel_width={rel:.3f} ({time.time() - t0:.0f}s)", flush=True)
        table.reindex(ALL).to_csv(out_csv, index=True)
        print("  saved (partial results preserved)", flush=True)

    print(table.reindex(ALL).to_string())


if __name__ == "__main__":
    main()
