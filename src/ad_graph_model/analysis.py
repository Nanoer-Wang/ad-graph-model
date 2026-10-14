"""Outcome metrics, time-series extraction and local sensitivity."""
from __future__ import annotations

import json
from typing import Dict

import numpy as np
import pandas as pd

from .model import aggregation_efficiency, langmuir_occupancy
from .simulate import SimulationResult, run_simulation


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def cumulative_trapezoid(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))])


def threshold_time(t: np.ndarray, y: np.ndarray, threshold: float):
    """First interpolated crossing time of ``threshold`` (None if never)."""
    hit = np.where(y >= threshold)[0]
    if len(hit) == 0:
        return None
    idx = int(hit[0])
    if idx == 0:
        return float(t[0])
    t0, t1 = float(t[idx - 1]), float(t[idx])
    y0, y1 = float(y[idx - 1]), float(y[idx])
    if abs(y1 - y0) < 1e-12:
        return float(t1)
    return float(t0 + (threshold - y0) * (t1 - t0) / (y1 - y0))


def node_positions(rows: int, cols: int) -> np.ndarray:
    coords = np.zeros((rows * cols, 2), dtype=float)
    for r in range(rows):
        for c in range(cols):
            coords[r * cols + c] = [c, r]
    return coords


def front_radius_series(res: SimulationResult, threshold: float,
                        seed_node: int) -> np.ndarray:
    coords = res.positions[:, :2]
    seed_xy = coords[int(seed_node)]
    dist = np.linalg.norm(coords - seed_xy, axis=1)
    radii = np.zeros(len(res.t), dtype=float)
    for k in range(len(res.t)):
        active = res.a[k] >= threshold
        radii[k] = float(dist[active].max()) if np.any(active) else 0.0
    return radii


def front_speed(t: np.ndarray, radius: np.ndarray) -> float:
    rmax = float(np.max(radius))
    if rmax <= 1e-12:
        return 0.0
    t10 = threshold_time(t, radius, 0.10 * rmax)
    t90 = threshold_time(t, radius, 0.90 * rmax)
    if t10 is None or t90 is None or t90 <= t10:
        return 0.0
    return float((0.80 * rmax) / (t90 - t10))


# --------------------------------------------------------------------------
# metrics
# --------------------------------------------------------------------------
def summary_metrics(res: SimulationResult) -> Dict[str, float]:
    mean_X = res.X.mean(axis=1)
    mean_Y = res.Y.mean(axis=1)
    mean_Z = res.Z.mean(axis=1)
    mean_a = res.a.mean(axis=1)
    mean_th = res.th.mean(axis=1)
    mean_P = res.P.mean(axis=1)

    idx_peak = int(np.argmax(mean_Y))
    auc_Y = cumulative_trapezoid(mean_Y, res.t)
    auc_a = cumulative_trapezoid(mean_a, res.t)
    auc_Z = cumulative_trapezoid(mean_Z, res.t)

    # aggregation efficiency: flux-weighted mean of (1 - theta) over nodes/time
    eff = aggregation_efficiency(res.th)
    mean_eff = eff.mean(axis=1)

    metrics = {
        "Y_peak": float(mean_Y[idx_peak]),
        "t_peak_Y": float(res.t[idx_peak]),
        "AUC_Y": float(auc_Y[-1]),
        "AUC_Z": float(auc_Z[-1]),
        "AUC_a": float(auc_a[-1]),
        "mean_X_end": float(mean_X[-1]),
        "mean_Y_end": float(mean_Y[-1]),
        "mean_Z_end": float(mean_Z[-1]),
        "mean_a_end": float(mean_a[-1]),
        "frac_a_over_0.5_end": float((res.a[-1] > 0.5).mean()),
        "frac_a_over_0.8_end": float((res.a[-1] > 0.8).mean()),
        "t_a_0.2": threshold_time(res.t, mean_a, 0.2) or float("nan"),
        "t_a_0.5": threshold_time(res.t, mean_a, 0.5) or float("nan"),
        "t_a_0.8": threshold_time(res.t, mean_a, 0.8) or float("nan"),
        "mean_passivation_end": float(mean_th[-1]),
        "mean_P_end": float(mean_P[-1]),
        "mean_eff_end": float(mean_eff[-1]),
        "eff_at_peak_Y": float(mean_eff[idx_peak]),
    }
    seed = int(res.config.get("seed_nodes", [0])[0])
    if res.kind == "sheet":
        r05 = front_radius_series(res, 0.5, seed)
        metrics["front_speed_0.5"] = front_speed(res.t, r05)
        metrics["front_radius_0.5_max"] = float(r05.max())
    return metrics


def timeseries_dataframe(res: SimulationResult) -> pd.DataFrame:
    return pd.DataFrame({
        "time": res.t,
        "mean_X": res.X.mean(axis=1),
        "mean_Y": res.Y.mean(axis=1),
        "mean_Z": res.Z.mean(axis=1),
        "mean_P": res.P.mean(axis=1),
        "mean_th": res.th.mean(axis=1),
        "mean_eff": aggregation_efficiency(res.th).mean(axis=1),
        "mean_a": res.a.mean(axis=1),
    })


def node_dataframe(res: SimulationResult, node: int) -> pd.DataFrame:
    return pd.DataFrame({
        "time": res.t,
        "X": res.X[:, node], "Y": res.Y[:, node], "Z": res.Z[:, node],
        "P": res.P[:, node], "th": res.th[:, node], "a": res.a[:, node],
    })


# --------------------------------------------------------------------------
# local sensitivity (kept for backwards compatibility / supplement)
# --------------------------------------------------------------------------
def _deepcopy_jsonable(obj):
    if isinstance(obj, dict):
        return {k: _deepcopy_jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_deepcopy_jsonable(v) for v in obj]
    return obj


def local_sensitivity(cfg: dict, metric_name: str = "Y_peak",
                      rel_step: float = 0.1) -> pd.DataFrame:
    base_res = run_simulation(cfg)
    base_metrics = summary_metrics(base_res)
    base_value = base_metrics[metric_name]

    rows = []
    for key, value in cfg.get("parameters", {}).items():
        if value == 0:
            continue
        cfg_up = _deepcopy_jsonable(cfg)
        cfg_dn = _deepcopy_jsonable(cfg)
        cfg_up["parameters"][key] = float(value) * (1.0 + rel_step)
        cfg_dn["parameters"][key] = max(float(value) * (1.0 - rel_step), 1e-12)

        up = summary_metrics(run_simulation(cfg_up))[metric_name]
        dn = summary_metrics(run_simulation(cfg_dn))[metric_name]
        sens = (up - dn) / (2.0 * rel_step * base_value + 1e-12)
        rows.append({"parameter": key, "baseline": value, "metric": metric_name,
                     "sensitivity": sens, "up": up, "down": dn})

    df = pd.DataFrame(rows)
    return df.sort_values("sensitivity", key=lambda s: np.abs(s),
                          ascending=False)
