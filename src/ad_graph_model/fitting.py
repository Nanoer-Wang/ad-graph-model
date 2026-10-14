"""Parameter identifiability: synthetic-target fitting, recovery and
profile-likelihood confidence intervals.

A synthetic surrogate dataset is generated from a reference parameter set
(``true`` values) sampled at observation times with heteroscedastic noise.
Selected parameters are then re-estimated by multi-start nonlinear least
squares.  Structural identifiability is summarised by the parameter
correlation matrix; practical identifiability by profile-likelihood
confidence intervals (chi-square threshold, 1 d.f.).
"""
from __future__ import annotations

import gc
import json
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import least_squares

from .simulate import run_simulation
from .analysis import timeseries_dataframe


def make_synthetic_target(cfg: dict, true_params: Dict[str, float],
                          obs_times: np.ndarray, columns: Sequence[str],
                          noise_rel: float = 0.05,
                          seed: int = 7) -> pd.DataFrame:
    """Simulate with ``true_params`` and add heteroscedastic noise."""
    c = json.loads(json.dumps(cfg))
    c["parameters"].update(true_params)
    res = run_simulation(c)
    df = timeseries_dataframe(res)
    rng = np.random.default_rng(seed)
    out = pd.DataFrame({"time": obs_times})
    for col in columns:
        vals = np.interp(obs_times, df["time"].to_numpy(), df[col].to_numpy())
        scale = np.maximum(np.abs(vals), 0.02)
        noisy = vals + rng.normal(0.0, noise_rel * scale)
        out[col] = np.clip(noisy, 0.0, None)
    return out


def _sim_columns(c: dict, names: Sequence[str], t_obs: np.ndarray,
                 cols: Sequence[str]) -> np.ndarray:
    """Run the model with parameters ``names = x`` and interpolate at t_obs."""
    res = run_simulation(c)
    df = timeseries_dataframe(res)
    return np.column_stack([
        np.interp(t_obs, df["time"].to_numpy(), df[col].to_numpy())
        for col in cols])


def fit_parameters(cfg: dict, param_names: Sequence[str],
                   bounds: Dict[str, Tuple[float, float]],
                   target_df: pd.DataFrame, n_starts: int = 3,
                   seed: int = 0, max_nfev: int = 30) -> dict:
    """Multi-start nonlinear least squares in **log parameter space**.

    Working in log space removes the scale disparity between parameters
    (e.g. ``k_n ~ 0.01`` vs ``k_s ~ 2``), which makes the Gauss-Newton
    linearisation, the Jacobian-based covariance and the warm-started
    inner optimisations of the profile likelihood far better conditioned.
    The returned ``x`` / ``sd`` / ``cov`` are mapped back to linear space
    with the delta method:  ``cov_lin ~ diag(x) cov_log diag(x)``.
    """
    names = list(param_names)
    lo_lin = np.array([bounds[n][0] for n in names])
    hi_lin = np.array([bounds[n][1] for n in names])
    lo_log = np.log(lo_lin)
    hi_log = np.log(hi_lin)
    cols = [c for c in target_df.columns if c != "time"]
    t_obs = target_df["time"].to_numpy(dtype=float)
    Y_obs = target_df[cols].to_numpy(dtype=float)
    scales = np.maximum(np.abs(Y_obs).max(axis=0), 0.05)

    # looser solver tolerances during fitting (final trajectories are
    # re-simulated at reference accuracy)
    fit_cfg = json.loads(json.dumps(cfg))
    fit_cfg.setdefault("simulation", {})["rtol"] = 1e-5
    fit_cfg["simulation"]["atol"] = 1e-7

    rng = np.random.default_rng(seed)
    x0s = []
    base_log = np.log(np.array([cfg["parameters"][n] for n in names]))
    x0s.append(base_log.copy())
    for _ in range(n_starts - 1):
        x0s.append(np.clip(base_log + rng.normal(0, 0.5, len(names)),
                           lo_log, hi_log))

    def residuals(xl):
        c = json.loads(json.dumps(fit_cfg))
        for n, v in zip(names, np.exp(xl)):
            c["parameters"][n] = float(v)
        sim = _sim_columns(c, names, t_obs, cols)
        gc.collect()      # release the integrator workspace between
        return ((sim - Y_obs) / scales).ravel()  # evaluations

    best = None
    for x0 in x0s:
        try:
            r = least_squares(residuals, x0=x0, bounds=(lo_log, hi_log),
                              max_nfev=max_nfev, xtol=1e-10, ftol=1e-10)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        raise RuntimeError("all fits failed")

    # Gauss-Newton covariance in log space, then delta method to linear
    J = best.jac
    dof = max(J.shape[0] - J.shape[1], 1)
    s2 = 2.0 * best.cost / dof
    try:
        cov_log = s2 * np.linalg.pinv(J.T @ J)
        sd_log = np.sqrt(np.clip(np.diag(cov_log), 0.0, None))
        x_lin = np.exp(best.x)
        cov_lin = (x_lin[:, None] * cov_log) * x_lin[None, :]
        sd = np.sqrt(np.clip(np.diag(cov_lin), 0.0, None))
    except Exception:
        cov_log = np.full((len(names), len(names)), np.nan)
        cov_lin = cov_log.copy()
        sd = np.full(len(names), np.nan)

    return {
        "names": names,
        "x": x_lin,
        "sd": sd,
        "cov": cov_lin,
        "cov_log": cov_log,
        "x_log": best.x.copy(),
        "cost": float(best.cost),
        "success": bool(best.success),
        "columns": list(cols),
        "n_obs": int(Y_obs.size),
    }


def profile_likelihood_ci(cfg: dict, fit: dict,
                          bounds: Dict[str, Tuple[float, float]],
                          target_df: pd.DataFrame, log_half_width: float = 0.8,
                          n_grid: int = 9,
                          inner_max_nfev: int = 6,
                          n_bisect: int = 8) -> pd.DataFrame:
    """True profile-likelihood 95% intervals (chi-square, 1 d.f.).

    For each parameter the SSE is scanned in **log space** while **all
    other parameters are re-optimised** (warm-started from the best fit,
    bounded, short runs).  A coarse grid locates the region below the
    chi-square threshold; the interval edges are then refined by
    bisection on the below/above transitions, so that intervals much
    narrower than the grid spacing are still resolved.
    """
    names = fit["names"]
    cols = [c for c in target_df.columns if c != "time"]
    t_obs = target_df["time"].to_numpy(dtype=float)
    Y_obs = target_df[cols].to_numpy(dtype=float)
    scales = np.maximum(np.abs(Y_obs).max(axis=0), 0.05)
    x_best = np.asarray(fit["x"], dtype=float)
    x_log = np.log(x_best)
    lo_lin = np.array([bounds[n][0] for n in names])
    hi_lin = np.array([bounds[n][1] for n in names])
    lo_log = np.log(lo_lin)
    hi_log = np.log(hi_lin)

    fit_cfg = json.loads(json.dumps(cfg))
    fit_cfg.setdefault("simulation", {})["rtol"] = 1e-5
    fit_cfg["simulation"]["atol"] = 1e-7

    def residuals(xl):
        c = json.loads(json.dumps(fit_cfg))
        for n, v in zip(names, np.exp(xl)):
            c["parameters"][n] = float(v)
        sim = _sim_columns(c, names, t_obs, cols)
        gc.collect()      # release the integrator workspace between
        return ((sim - Y_obs) / scales).ravel()  # evaluations

    def sse(xl):
        return float(np.sum(residuals(xl) ** 2))

    def profile_sse(i, g):
        """SSE with parameter i fixed at g and all others re-optimised."""
        x = x_log.copy()
        x[i] = g
        free = [j for j in range(len(names)) if j != i]

        def res_free(xf):
            xx = x.copy()
            xx[free] = xf
            return residuals(xx)

        r = least_squares(res_free, x0=x_log[free],
                          bounds=(lo_log[free], hi_log[free]),
                          max_nfev=inner_max_nfev,
                          xtol=1e-8, ftol=1e-8)
        x[free] = r.x
        return sse(x)

    sse_min = sse(x_log)
    thr = sse_min * (1.0 + stats.chi2.ppf(0.95, 1) / max(Y_obs.size, 1))
    rows = []
    for i, n in enumerate(names):
        grid = np.clip(x_log[i] + np.linspace(-log_half_width, log_half_width,
                                              n_grid), lo_log[i], hi_log[i])
        prof = np.array([profile_sse(i, g) for g in grid])
        below = prof <= thr
        if not below.any():
            ci_lo = ci_hi = float(x_best[i])
        else:
            # maximal runs of consecutive below-threshold grid points
            runs = []
            k = 0
            while k < len(below):
                if below[k]:
                    j = k
                    while j + 1 < len(below) and below[j + 1]:
                        j += 1
                    runs.append((k, j))
                    k = j + 1
                else:
                    k += 1
            # the run containing the best profile point
            k_best = int(np.argmin(prof))
            run = min(runs, key=lambda r: min(abs(k_best - r[0]),
                                              abs(k_best - r[1])))
            k0, k1 = run

            def refine_edge(k_left_below):
                """Bisect the below/above transition at grid index k.

                ``grid[k]`` is below threshold, ``grid[k + 1]`` is above;
                ``a`` converges to the below-threshold side.  Boundary
                indices are guarded (the grid may be clipped to the
                parameter bounds, producing repeated points).
                """
                kk = int(k_left_below)
                if kk < 0:
                    return grid[0]
                if kk + 1 >= len(grid):
                    return grid[-1]
                a, b = grid[kk], grid[kk + 1]
                for _ in range(n_bisect):
                    mid = 0.5 * (a + b)
                    if profile_sse(i, mid) <= thr:
                        a = mid
                    else:
                        b = mid
                return a  # a converges to the below-threshold side

            lo_edge = grid[k0] if k0 <= 0 else refine_edge(k0 - 1)
            hi_edge = grid[k1] if k1 >= len(grid) - 1 else refine_edge(k1)
            ci_lo = float(np.exp(lo_edge))
            ci_hi = float(np.exp(hi_edge))
        rows.append({"parameter": n, "best": float(x_best[i]),
                     "ci_lo": ci_lo, "ci_hi": ci_hi,
                     "sd_linear": float(fit["sd"][i])})
    return pd.DataFrame(rows)


def parameter_correlation(fit: dict) -> pd.DataFrame:
    """Correlation matrix from the linearised (Gauss-Newton) covariance."""
    cov = np.asarray(fit.get("cov"), dtype=float)
    sd = np.asarray(fit["sd"], dtype=float)
    if not np.all(np.isfinite(cov)) or not np.all(sd > 0):
        return pd.DataFrame(np.eye(len(sd)), index=fit["names"],
                            columns=fit["names"])
    denom = np.outer(sd, sd)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = np.where(denom > 0, cov / denom, 0.0)
    np.fill_diagonal(corr, 1.0)
    return pd.DataFrame(corr, index=fit["names"], columns=fit["names"])
