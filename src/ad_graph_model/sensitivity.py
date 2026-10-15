"""Global uncertainty and sensitivity analysis (Latin Hypercube + PRCC).

The PRCC (partial rank correlation coefficient) methodology follows the
standard systems-biology workflow (Marino et al., J Theor Biol 2008):

1. sample parameter vectors with Latin Hypercube Sampling,
2. simulate the model and record outcome metrics,
3. rank-transform inputs and outputs,
4. compute partial correlations of each input/output pair, controlling for
   all remaining inputs,
5. significance from the t-distribution with df = n - 2 - k.
"""
from __future__ import annotations

import json
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy import stats

from .simulate import run_simulation
from .analysis import summary_metrics


def latin_hypercube(n: int, d: int, seed: int = 0) -> np.ndarray:
    """LHS design in the unit hypercube [0,1]^d."""
    rng = np.random.default_rng(seed)
    cut = np.linspace(0, 1, n + 1)
    u = np.empty((n, d))
    for j in range(d):
        pts = (cut[:-1] + cut[1:]) / 2.0
        rng.shuffle(pts)
        u[:, j] = pts + rng.uniform(-0.5, 0.5, n) / n
    return np.clip(u, 0.0, 1.0)


def scale_samples(u: np.ndarray, bounds: Dict[str, Tuple[float, float]],
                  names: Sequence[str], log: bool = True) -> np.ndarray:
    """Map unit-cube samples to parameter values (log-uniform by default)."""
    out = np.empty_like(u)
    for j, name in enumerate(names):
        lo, hi = bounds[name]
        if log:
            out[:, j] = np.exp(np.log(lo) + u[:, j] * (np.log(hi) - np.log(lo)))
        else:
            out[:, j] = lo + u[:, j] * (hi - lo)
    return out


def prcc_matrix(X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Partial rank correlation coefficients and two-sided p-values.

    Parameters
    ----------
    X : (n, d) input samples
    Y : (n, m) output samples

    Returns
    -------
    prcc : (d, m) coefficients
    pval : (d, m) p-values
    """
    n, d = X.shape
    m = Y.shape[1]
    Xr = np.apply_along_axis(stats.rankdata, 0, X).astype(float)
    Yr = np.apply_along_axis(stats.rankdata, 0, Y).astype(float)

    prcc = np.zeros((d, m))
    pval = np.ones((d, m))
    df = n - 2 - (d - 1)
    if df <= 0:
        raise ValueError("n too small for PRCC with this many parameters")

    # residualize each input against all OTHER inputs
    R = np.empty((n, d), dtype=float)
    for j in range(d):
        others = np.delete(Xr, j, axis=1)
        A = np.column_stack([np.ones(n), others])
        coef, *_ = np.linalg.lstsq(A, Xr[:, j], rcond=None)
        R[:, j] = Xr[:, j] - A @ coef

    # For each (input j, output k) pair: residualize the output against all
    # inputs EXCEPT j, then correlate with R[:, j] (Marino et al. 2008).
    for j in range(d):
        others = np.delete(Xr, j, axis=1)
        A_y = np.column_stack([np.ones(n), others])
        for k in range(m):
            coef, *_ = np.linalg.lstsq(A_y, Yr[:, k], rcond=None)
            ry = Yr[:, k] - A_y @ coef
            num = float(np.dot(R[:, j], ry))
            den = float(np.linalg.norm(R[:, j]) * np.linalg.norm(ry))
            if den <= 0:
                continue
            r = num / den
            prcc[j, k] = r
            t = r * np.sqrt(df / max(1.0 - r * r, 1e-12))
            pval[j, k] = 2.0 * stats.t.sf(abs(t), df)
    return prcc, pval


def global_sensitivity(cfg: dict, param_names: Sequence[str],
                       bounds: Dict[str, Tuple[float, float]],
                       metrics: Sequence[str] = ("Y_peak", "AUC_Y",
                                                 "mean_a_end", "t_a_0.5",
                                                 "frac_a_over_0.8_end"),
                       n_samples: int = 512, seed: int = 1) -> dict:
    """LHS ensemble + PRCC for a list of outcome metrics."""
    names = list(param_names)
    u = latin_hypercube(n_samples, len(names), seed=seed)
    X = scale_samples(u, bounds, names, log=True)

    values = np.zeros((n_samples, len(metrics)))
    all_Y, all_a = [], []
    for i in range(n_samples):
        c = json.loads(json.dumps(cfg))
        for j, name in enumerate(names):
            c["parameters"][name] = float(X[i, j])
        res = run_simulation(c)
        m = summary_metrics(res)
        for k, met in enumerate(metrics):
            v = float(m.get(met, np.nan))
            if not np.isfinite(v):
                v = float(res.config["simulation"]["t_end"])
            values[i, k] = v
        from .analysis import timeseries_dataframe
        d = timeseries_dataframe(res)
        all_Y.append(d["mean_Y"].to_numpy())
        all_a.append(d["mean_a"].to_numpy())

    Y = values
    prcc, pval = prcc_matrix(X, Y)

    df = pd.DataFrame({"parameter": names})
    for k, met in enumerate(metrics):
        df[f"prcc_{met}"] = prcc[:, k]
        df[f"p_{met}"] = pval[:, k]
        df[f"sig_{met}"] = np.where(pval[:, k] < 0.05, "*", "")

    ens = None
    if all_Y:
        stackY = np.stack(all_Y)
        stack_a = np.stack(all_a)
        ens = {
            "time": np.asarray(res.t),
            "Y_q05": np.quantile(stackY, 0.05, axis=0),
            "Y_q50": np.quantile(stackY, 0.50, axis=0),
            "Y_q95": np.quantile(stackY, 0.95, axis=0),
            "a_q05": np.quantile(stack_a, 0.05, axis=0),
            "a_q50": np.quantile(stack_a, 0.50, axis=0),
            "a_q95": np.quantile(stack_a, 0.95, axis=0),
        }

    samples_df = pd.DataFrame(X, columns=names)
    for k, met in enumerate(metrics):
        samples_df[met] = Y[:, k]

    return {"table": df, "samples": samples_df, "ensemble": ens,
            "n_samples": n_samples, "metrics": list(metrics)}


def default_bounds(base: Dict[str, float], span: float = 2.0) -> Dict[str, Tuple[float, float]]:
    """Symmetric log bounds around baseline values."""
    return {k: (v / span, v * span) for k, v in base.items() if v > 0}
