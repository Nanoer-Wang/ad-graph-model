"""Experiment driver: computes every dataset used by the manuscript figures.

All heavy computation lives here so that figure rendering can be iterated
independently.  Results are cached as CSV/JSON under ``outdir/data``.
"""
from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd

from .analysis import (cumulative_trapezoid, front_radius_series, summary_metrics,
                       timeseries_dataframe, node_dataframe, threshold_time)
from .bifurcation import (find_steady_states, node_f4, node_jac4,
                          qss_branches, qss_damage_rate, pseudo_potential,
                          simulate_node)
from .config import load_config
from .fitting import (fit_parameters, make_synthetic_target,
                      parameter_correlation, profile_likelihood_ci)
from .graph import build_graph
from .io_utils import ensure_dir, save_dataframe, save_json
from .model import aggregation_efficiency, hill, langmuir_occupancy
from .sensitivity import default_bounds, global_sensitivity
from .simulate import SimulationResult, run_simulation

SEED = 7


def _KD(p: Dict[str, float]) -> float:
    return p["koff"] / p["kon"]


def theta_of_P(P: float, p: Dict[str, float]) -> float:
    return langmuir_occupancy(P, _KD(p))


def P_of_theta(theta: float, p: Dict[str, float]) -> float:
    return theta * _KD(p) / max(1.0 - theta, 1e-9)


def uP_of_P(P: float, p: Dict[str, float]) -> float:
    return P * p["cP"]


# ---------------------------------------------------------------------------
# 1. baseline + topology robustness
# ---------------------------------------------------------------------------
def exp_baseline(cfg: dict, outdir: Path) -> dict:
    out = {}
    res = run_simulation(cfg)
    out["sheet"] = res
    save_dataframe(timeseries_dataframe(res), outdir / "baseline_timeseries.csv")
    seed = int(cfg.get("seed_nodes", [0])[0])
    save_dataframe(node_dataframe(res, seed), outdir / "baseline_seed_node.csv")
    save_json(summary_metrics(res), outdir / "baseline_metrics.json")

    # topology robustness
    topo = {}
    for name, gcfg in {
        "sheet": {"type": "sheet", "rows": 20, "cols": 20},
        "smallworld": {"type": "smallworld", "n": 68, "k": 6, "p": 0.12, "seed": 3},
        "connectome": {"type": "connectome", "seed": 0},
    }.items():
        c = copy.deepcopy(cfg)
        c["graph"] = gcfg
        g = build_graph(gcfg)
        if name == "sheet":
            c["seed_nodes"] = [105]
        elif name == "smallworld":
            c["seed_nodes"] = [34]          # central node of the ring
        else:
            # seed the entorhinal cortex (earliest Braak stage)
            c["seed_nodes"] = [i for i, lab in enumerate(g.labels)
                               if lab.endswith("entorhinal")]
        r = run_simulation(c)
        topo[name] = {"res": r, "metrics": summary_metrics(r)}
    out["topology"] = topo
    for name, d in topo.items():
        save_dataframe(timeseries_dataframe(d["res"]),
                       outdir / f"topology_{name}_timeseries.csv")
        save_json(d["metrics"], outdir / f"topology_{name}_metrics.json")
    save_json({k: d["metrics"] for k, d in topo.items()},
              outdir / "topology_metrics.json")
    return out


# ---------------------------------------------------------------------------
# 2. intervention strategies
# ---------------------------------------------------------------------------
STRATEGIES: Dict[str, list] = {
    "Baseline": [],
    "Passivation": [{"time": 8.0, "set": {"uP": 0.012}}],
    "Clearance": [{"time": 8.0, "set": {"cY": 0.50, "cZ": 0.20}}],
    "Nucleation": [{"time": 8.0, "set": {"kn": 0.002}}],
    "Propagation": [{"time": 8.0, "set": {"DY": 0.02, "etaN": 0.05}}],
    "Combination": [{"time": 8.0, "set": {"uP": 0.012, "cY": 0.50,
                                          "cZ": 0.20}}],
}


def exp_interventions(cfg: dict, outdir: Path) -> dict:
    out = {"metrics": {}, "curves": {}}
    for name, iv in STRATEGIES.items():
        c = copy.deepcopy(cfg)
        c["interventions"] = iv
        res = run_simulation(c)
        out["metrics"][name] = summary_metrics(res)
        out["curves"][name] = timeseries_dataframe(res)
        save_dataframe(out["curves"][name],
                       outdir / f"strategy_{name}_timeseries.csv")
    df = pd.DataFrame(out["metrics"]).T
    df.index.name = "strategy"
    save_dataframe(df.reset_index(), outdir / "strategy_metrics.csv")
    save_json(out["metrics"], outdir / "strategy_metrics.json")
    return out


# ---------------------------------------------------------------------------
# 3. dose-response of passivation
# ---------------------------------------------------------------------------
def exp_dose_response(cfg: dict, outdir: Path, n: int = 26) -> dict:
    uPs = np.logspace(-3.3, -0.7, n)   # 5e-4 .. 0.2
    rows = []
    curves = []
    for u in uPs:
        c = copy.deepcopy(cfg)
        c["interventions"] = [{"time": 8.0, "set": {"uP": float(u)}}]
        res = run_simulation(c)
        m = summary_metrics(res)
        p = c["parameters"]
        m["uP"] = float(u)
        m["P_ss"] = float(u) / p["cP"]
        m["theta_ss"] = theta_of_P(float(u) / p["cP"], p)
        rows.append(m)
        curves.append(timeseries_dataframe(res).assign(uP=float(u)))
    df = pd.DataFrame(rows)
    save_dataframe(df, outdir / "dose_response.csv")
    pd.concat(curves, ignore_index=True).to_csv(
        outdir / "dose_response_curves.csv", index=False)

    # binding kinetics: theta(t) for a few doses
    kin = {}
    for u in [0.003, 0.006, 0.012, 0.025]:
        c = copy.deepcopy(cfg)
        c["interventions"] = [{"time": 8.0, "set": {"uP": float(u)}}]
        res = run_simulation(c)
        kin[f"{u}"] = timeseries_dataframe(res)
    save_json({k: 1 for k in kin}, outdir / "dose_kinetics_keys.json")
    kin_df = pd.concat([d.assign(uP=float(k)) for k, d in kin.items()],
                       ignore_index=True)
    kin_df.to_csv(outdir / "dose_kinetics.csv", index=False)

    # aggregation-efficiency landscape over the untreated baseline
    res0 = run_simulation(copy.deepcopy(cfg))
    eff = aggregation_efficiency(res0.th)
    save_dataframe(pd.DataFrame({"time": res0.t,
                                 "mean_eff": eff.mean(axis=1)}),
                   outdir / "baseline_efficiency.csv")
    return {"table": df, "curves": curves, "kinetics": kin, "baseline": res0}


# ---------------------------------------------------------------------------
# 4. therapeutic window + combination therapy
# ---------------------------------------------------------------------------
def exp_therapy(cfg: dict, outdir: Path) -> dict:
    onsets = np.arange(0.0, 33.0, 2.0)
    strategies = {
        "Passivation": {"uP": 0.012},
        "Clearance": {"cY": 0.50, "cZ": 0.20},
        "Combination": {"uP": 0.012, "cY": 0.50, "cZ": 0.20},
    }
    rows = []
    traj_store = {}
    for sname, sset in strategies.items():
        for t0 in onsets:
            c = copy.deepcopy(cfg)
            c["interventions"] = [{"time": float(t0), "set": dict(sset)}]
            res = run_simulation(c)
            m = summary_metrics(res)
            m.update({"strategy": sname, "onset": float(t0)})
            rows.append(m)
            if sname == "Passivation" and float(t0) in (0.0, 8.0, 16.0, 24.0):
                traj_store[float(t0)] = timeseries_dataframe(res)
    df = pd.DataFrame(rows)
    save_dataframe(df, outdir / "treatment_window.csv")
    for t0, d in traj_store.items():
        save_dataframe(d, outdir / f"window_trajectory_t{int(t0)}.csv")

    # clearance-only sweep (no passivation)
    clr_rows = []
    for cy in np.linspace(0.25, 0.60, 12):
        c = copy.deepcopy(cfg)
        c["interventions"] = [{"time": 8.0, "set": {"cY": float(cy)}}]
        res = run_simulation(c)
        m = summary_metrics(res)
        m["cY"] = float(cy)
        clr_rows.append(m)
    save_dataframe(pd.DataFrame(clr_rows), outdir / "clearance_sweep.csv")

    # 2D combination grid: clearance enhancement x passivation dose
    cYs = np.linspace(0.25, 0.60, 10)
    uPs = np.logspace(-3.3, -1.0, 10)
    grid = np.zeros((len(uPs), len(cYs)))
    for i, u in enumerate(uPs):
        for j, cy in enumerate(cYs):
            c = copy.deepcopy(cfg)
            c["interventions"] = [{"time": 8.0,
                                   "set": {"uP": float(u), "cY": float(cy)}}]
            res = run_simulation(c)
            grid[i, j] = summary_metrics(res)["mean_a_end"]
    grid_df = pd.DataFrame(grid, index=[f"{u:.5g}" for u in uPs],
                           columns=[f"{c:.3g}" for c in cYs])
    grid_df.index.name = "uP"
    grid_df.to_csv(outdir / "combination_grid.csv")
    return {"window": df, "grid": grid, "cYs": cYs, "uPs": uPs,
            "trajectories": traj_store}


# ---------------------------------------------------------------------------
# 5. global sensitivity (LHS + PRCC)
# ---------------------------------------------------------------------------
def exp_sensitivity(cfg: dict, outdir: Path, n: int = 400) -> dict:
    p = cfg["parameters"]
    names = ["lambda0", "beta", "cX", "kn", "ks", "ke", "cY", "kf", "cZ",
             "DY", "DX", "KY", "hY", "etaY", "etaN", "deltaA", "kon", "koff"]
    bounds = {k: (p[k] / 3.0, p[k] * 3.0) for k in names}
    # m and q are structural orders: narrow range
    bounds["m"] = (1.6, 2.4)
    bounds["q"] = (0.8, 1.2)
    names = names + ["m", "q"]
    res = global_sensitivity(cfg, names, bounds,
                             metrics=("Y_peak", "AUC_Y", "mean_a_end",
                                      "t_a_0.5", "frac_a_over_0.8_end"),
                             n_samples=n, seed=SEED)
    res["table"].to_csv(outdir / "prcc_table.csv", index=False)
    res["samples"].to_csv(outdir / "lhs_samples.csv", index=False)
    if res["ensemble"] is not None:
        ens = pd.DataFrame(res["ensemble"])
        ens.to_csv(outdir / "lhs_ensemble.csv", index=False)
    return res


# ---------------------------------------------------------------------------
# 6. identifiability
# ---------------------------------------------------------------------------
TRUE_PERTURB = {"kn": 0.011, "ks": 2.1, "ke": 0.045, "cY": 0.18,
                "etaY": 0.42, "KY": 0.09}


def exp_identifiability(cfg: dict, outdir: Path) -> dict:
    obs_times = np.arange(0.0, 40.01, 2.0)
    target = make_synthetic_target(cfg, TRUE_PERTURB, obs_times,
                                   columns=["mean_Y", "mean_Z", "mean_a"],
                                   noise_rel=0.03, seed=SEED)
    save_dataframe(target, outdir / "synthetic_target.csv")

    p = cfg["parameters"]
    fit_names = ["kn", "ks", "ke", "cY", "etaY", "KY"]
    bounds = {n: (p[n] / 4.0, p[n] * 4.0) for n in fit_names}
    fit = fit_parameters(cfg, fit_names, bounds, target, n_starts=3, seed=SEED)
    rec = {n: {"true": float(TRUE_PERTURB[n]),
               "fit": float(v),
               "rel_err": float(v / TRUE_PERTURB[n] - 1.0)}
           for n, v in zip(fit["names"], fit["x"])}
    save_json({"best": {n: float(v) for n, v in zip(fit["names"], fit["x"])},
               "sd": {n: float(v) for n, v in zip(fit["names"], fit["sd"])},
               "recovery": rec,
               "max_rel_err": float(max(abs(r["rel_err"]) for r in rec.values())),
               "cost": fit["cost"], "n_obs": fit["n_obs"]},
              outdir / "fit_result.json")
    ci = profile_likelihood_ci(cfg, fit, bounds, target)
    save_dataframe(ci, outdir / "profile_likelihood_ci.csv")
    corr = parameter_correlation(fit)
    corr.to_csv(outdir / "parameter_correlation.csv")

    # fitted trajectory + prediction band (parameter re-sampling)
    best_cfg = copy.deepcopy(cfg)
    for n, v in zip(fit["names"], fit["x"]):
        best_cfg["parameters"][n] = float(v)
    res = run_simulation(best_cfg)
    fit_df = timeseries_dataframe(res)
    save_dataframe(fit_df, outdir / "fit_trajectory.csv")

    rng = np.random.default_rng(SEED)
    band_Y, band_a = [], []
    for _ in range(12):
        c = copy.deepcopy(best_cfg)
        for n, v, s in zip(fit["names"], fit["x"], fit["sd"]):
            if np.isfinite(s) and s > 0:
                c["parameters"][n] = float(np.clip(
                    v + rng.normal(0, s), bounds[n][0], bounds[n][1]))
        r = run_simulation(c)
        d = timeseries_dataframe(r)
        band_Y.append(np.interp(obs_times, d["time"], d["mean_Y"]))
        band_a.append(np.interp(obs_times, d["time"], d["mean_a"]))
    band = pd.DataFrame({
        "time": obs_times,
        "Y_lo": np.quantile(np.array(band_Y), 0.05, axis=0),
        "Y_hi": np.quantile(np.array(band_Y), 0.95, axis=0),
        "a_lo": np.quantile(np.array(band_a), 0.05, axis=0),
        "a_hi": np.quantile(np.array(band_a), 0.95, axis=0),
    })
    save_dataframe(band, outdir / "fit_prediction_band.csv")
    return {"fit": fit, "ci": ci, "corr": corr, "target": target,
            "trajectory": fit_df, "band": band}


# ---------------------------------------------------------------------------
# 7. bifurcation and bistability
# ---------------------------------------------------------------------------
def critical_theta(p: Dict[str, float], cY: float | None = None) -> float:
    """Analytic critical passivation fraction (linear stability of the
    fibril-free state)."""
    p = dict(p)
    if cY is not None:
        p["cY"] = cY
    Xs = p["lambda0"] / p["cX"]
    return 1.0 - (p["cY"] + p["kf"]) * p["cZ"] / (
        Xs * (p["ks"] * p["kf"] + p["ke"] * (p["cY"] + p["kf"])))


def exp_bifurcation(cfg: dict, outdir: Path) -> dict:
    p = cfg["parameters"]
    KD = _KD(p)

    # --- one-parameter sweep in passivation fraction theta
    thetas = np.linspace(0.0, 0.97, 40)
    rows = []
    for th in thetas:
        P = P_of_theta(th, p)
        for s, lm in find_steady_states(p, P_fixed=P, n_starts=60):
            rows.append({"theta": th, "P": P, "X": s[0], "Y": s[1],
                         "Z": s[2], "a": s[5], "lmax": lm,
                         "stable": lm < 0})
    bif = pd.DataFrame(rows)
    save_dataframe(bif, outdir / "bifurcation_theta.csv")

    # --- critical curve in (cY, theta): numerically located saddle-node
    cYs = np.linspace(0.25, 0.55, 13)
    th_num, th_ana = [], []
    for cy in cYs:
        pc = dict(p)
        pc["cY"] = float(cy)
        th_ana.append(critical_theta(pc))
        # bisection on existence of a low-damage stable state
        lo, hi = 0.0, 0.995

        def exists(th):
            P = P_of_theta(th, pc)
            st = find_steady_states(pc, P_fixed=P, n_starts=20)
            return any((s[5] < 0.1) and (lm < 0) for s, lm in st)

        if exists(hi) and not exists(lo):
            for _ in range(14):
                mid = 0.5 * (lo + hi)
                if exists(mid):
                    hi = mid
                else:
                    lo = mid
            th_num.append(hi)
        else:
            th_num.append(np.nan)
    crit = pd.DataFrame({"cY": cYs, "theta_numeric": th_num,
                         "theta_analytic": th_ana})
    save_dataframe(crit, outdir / "critical_curve.csv")

    # --- bistability demonstration at high passivation (sustained infusion)
    traj = {}
    P_hi = P_of_theta(0.95, p)
    p_hi = dict(p)
    # sustained infusion keeps P at its steady state P_hi (dP/dt = cP(P_hi-P))
    p_hi["uP"] = p["cP"] * P_hi
    starts = {
        "healthy": np.array([p["lambda0"] / p["cX"] * 1.05, 0.0, 0.0, P_hi, 0.0, 0.0]),
        "diseased": np.array([1.2, 0.08, 0.15, P_hi, 0.95, 0.55]),
    }
    for name, y0 in starts.items():
        t, Y = simulate_node(cfg, p=p_hi, y0=y0, t_end=120, n=1200)
        traj[name] = pd.DataFrame({"time": t, "X": Y[:, 0], "Y": Y[:, 1],
                                   "Z": Y[:, 2], "P": Y[:, 3], "th": Y[:, 4],
                                   "a": Y[:, 5]})
        save_dataframe(traj[name], outdir / f"bistability_{name}.csv")

    # --- pseudo-potential (QSS reduction) at three passivation levels
    pots = {}
    a_grid = np.linspace(0.0, 1.0, 300)
    for label, th in [("none", 0.0), ("partial", 0.65), ("strong", 0.95)]:
        P = P_of_theta(th, p)
        f = np.array([qss_damage_rate(a, p, P) for a in a_grid])
        pots[label] = {"a": a_grid, "f": f,
                       "phi": pseudo_potential(a_grid, f), "theta": th}
        save_dataframe(pd.DataFrame({"a": a_grid, "f": f,
                                     "phi": pots[label]["phi"],
                                     "theta": th}),
                       outdir / f"potential_{label}.csv")
    return {"bifurcation": bif, "critical": crit, "traj": traj,
            "potentials": pots, "theta_c_analytic": critical_theta(p)}


# ---------------------------------------------------------------------------
# 8. spatial propagation snapshots (for Fig. 3)
# ---------------------------------------------------------------------------
def exp_spatial(cfg: dict, outdir: Path) -> dict:
    c = copy.deepcopy(cfg)
    c["graph"] = {"type": "connectome", "seed": 0}
    from .graph import build_graph
    g = build_graph(c["graph"])
    labels = g.labels
    c["seed_nodes"] = [i for i, lab in enumerate(labels)
                       if lab.endswith("entorhinal")]
    res = run_simulation(c)
    save_dataframe(timeseries_dataframe(res), outdir / "spatial_timeseries.csv")
    snap_times = [4.0, 10.0, 18.0, 28.0]
    snaps = {}
    for tp in snap_times:
        idx = int(np.argmin(np.abs(res.t - tp)))
        snaps[float(res.t[idx])] = {"a": res.a[idx], "Y": res.Y[idx],
                                    "Z": res.Z[idx], "th": res.th[idx]}
    # region x time matrix of damage
    mat = pd.DataFrame({"time": res.t})
    region_df = pd.DataFrame(res.a, columns=labels)
    region_df.insert(0, "time", res.t)
    region_df.to_csv(outdir / "spatial_region_time.csv", index=False)

    # onset time per region (first crossing of a = 0.3)
    onsets = []
    for j, lab in enumerate(labels):
        t_cross = threshold_time(res.t, res.a[:, j], 0.3)
        onsets.append({"region": lab, "onset": t_cross if t_cross else np.nan,
                       "a_end": float(res.a[-1, j])})
    onsets_df = pd.DataFrame(onsets).sort_values("onset")
    save_dataframe(onsets_df, outdir / "spatial_onsets.csv")

    # degree vs damage correlation
    strength = g.W.sum(axis=1)
    corr_df = pd.DataFrame({"region": labels, "strength": strength,
                            "a_end": res.a[-1], "Y_end": res.Y[-1]})
    corr_df.to_csv(outdir / "spatial_strength_damage.csv", index=False)
    r = float(np.corrcoef(np.log(strength + 1e-9), corr_df["a_end"])[0, 1])
    save_json({"pearson_r_strength_vs_damage": r}, outdir / "spatial_corr.json")
    return {"res": res, "labels": labels, "positions": g.positions,
            "snaps": snaps, "onsets": onsets_df, "strength": strength,
            "corr_r": r}


# ---------------------------------------------------------------------------
# master runner
# ---------------------------------------------------------------------------
def run_all(cfg_path: str | Path, outdir: str | Path,
            skip_sensitivity: bool = False, n_sens: int = 400) -> dict:
    outdir = ensure_dir(outdir)
    data_dir = ensure_dir(outdir / "data")
    cfg = load_config(cfg_path)
    save_json(cfg, data_dir / "used_config.json")

    results: Dict[str, object] = {}
    print("[1/8] baseline + topology ...", flush=True)
    results["baseline"] = exp_baseline(cfg, data_dir)
    print("[2/8] intervention strategies ...", flush=True)
    results["strategies"] = exp_interventions(cfg, data_dir)
    print("[3/8] dose-response ...", flush=True)
    results["dose"] = exp_dose_response(cfg, data_dir)
    print("[4/8] therapeutic window + combination ...", flush=True)
    results["therapy"] = exp_therapy(cfg, data_dir)
    if not skip_sensitivity:
        print(f"[5/8] global sensitivity (n={n_sens}) ...", flush=True)
        results["sensitivity"] = exp_sensitivity(cfg, data_dir, n=n_sens)
    else:
        print("[5/8] global sensitivity SKIPPED", flush=True)
    print("[6/8] identifiability ...", flush=True)
    results["identifiability"] = exp_identifiability(cfg, data_dir)
    print("[7/8] bifurcation ...", flush=True)
    results["bifurcation"] = exp_bifurcation(cfg, data_dir)
    print("[8/8] spatial propagation ...", flush=True)
    results["spatial"] = exp_spatial(cfg, data_dir)
    print("All experiments finished.", flush=True)
    return results
