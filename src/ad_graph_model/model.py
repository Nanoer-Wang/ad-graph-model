"""Extended amyloid-cascade model with fibril-interface passivation.

State variables per network node v:

    X_v  amyloid-beta monomer burden
    Y_v  toxic soluble oligomer burden
    Z_v  fibril / plaque burden
    P_v  passivator burden (interface-blocking agent)
    th_v fraction of fibril interfaces occupied by passivator (0 = free)
    a_v  neuronal damage in [0, 1]

Reference dynamics in model units (well-mixed node):

    dX/dt  = lam (1-a)(1+beta a) - cX X - kn X^m - ke X Z (1-th)
             + DX * Lap(X)
    dY/dt  = kn X^m + ks X Z^q (1-th) - (cY + kf) Y + DY * Lap(Y)
    dZ/dt  = kf Y + ke X Z (1-th) - cZ Z
    dP/dt  = u_P - cP P
    dth/dt = kon P (1-th) - koff th
    da/dt  = etaY Hill(Y; KY, hY) (1-a) + etaN * sum_u W_vu ReLU(a_u-a_v)
             - deltaA a

The passivation state th reduces both surface-mediated aggregation
processes (secondary nucleation and elongation) through the factor
(1 - th), i.e. the *interface activity*.  The quasi-steady-state
occupancy th* = P / (K_D + P), K_D = koff/kon, is a Langmuir isotherm.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np


@dataclass
class ParameterSchedule:
    base: Dict[str, float]
    interventions: List[Dict[str, Any]]

    def at(self, t: float) -> Dict[str, float]:
        params = dict(self.base)
        for item in self.interventions:
            if t >= float(item["time"]):
                for k, v in item.get("set", {}).items():
                    params[k] = float(v)
        return params


def hill(y, K: float, h: float):
    y = np.asarray(y, dtype=float)
    y_clip = np.clip(y, 0.0, None)
    return (y_clip ** h) / (K ** h + y_clip ** h + 1e-12)


def neighbor_damage(a: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Damage transmission from more damaged neighbours (a_u > a_v)."""
    diff = a[None, :] - a[:, None]
    relu = np.maximum(diff, 0.0)
    return (W * relu).sum(axis=1)


def diffusion(x: np.ndarray, W: np.ndarray) -> np.ndarray:
    """Row-normalized mixing, W x - x.

    Preserves constants and, for connected undirected raw weights, the
    strength-weighted sum. The unweighted sum is conserved only when W
    is doubly stochastic.
    """
    return W @ x - x


DEFAULT_PARAMETERS: Dict[str, float] = {
    # monomer production / clearance
    "lambda0": 0.22, "beta": 4.5, "cX": 0.28,
    # primary nucleation (cooperative, order m)
    "kn": 0.008, "m": 2.0,
    # surface-mediated aggregation on fibril interfaces
    "ks": 1.5, "q": 1.0, "ke": 0.03,
    # oligomer / fibril clearance and conversion
    "cY": 0.25, "kf": 0.06, "cZ": 0.114,
    # network transport
    "DX": 0.03, "DY": 0.06,
    # toxicity and damage
    "KY": 0.12, "hY": 4.0, "etaY": 0.30, "etaN": 0.15, "deltaA": 0.01,
    # passivator binding kinetics (K_D = koff/kon = 0.15)
    "kon": 2.0, "koff": 0.30, "cP": 0.02,
    # passivator infusion rate (interventions set "uP")
    "uP": 0.0,
}

DEFAULT_INITIAL: Dict[str, float] = {"X": 0.03, "Y": 0.0, "Z": 0.0,
                                     "P": 0.0, "th": 0.0, "a": 0.0}


def rhs(t: float, state: np.ndarray, n_nodes: int, W: np.ndarray,
        schedule: ParameterSchedule) -> np.ndarray:
    p = schedule.at(t)
    X = np.clip(state[0:n_nodes], 0.0, None)
    Y = np.clip(state[n_nodes:2 * n_nodes], 0.0, None)
    Z = np.clip(state[2 * n_nodes:3 * n_nodes], 0.0, None)
    P = np.clip(state[3 * n_nodes:4 * n_nodes], 0.0, None)
    th = np.clip(state[4 * n_nodes:5 * n_nodes], 0.0, 1.0)
    a = np.clip(state[5 * n_nodes:6 * n_nodes], 0.0, 1.0)

    activity = 1.0 - th                       # free (active) interface fraction
    prod = p["lambda0"] * (1.0 - a) * (1.0 + p["beta"] * a)

    dX = (prod - p["cX"] * X - p["kn"] * (X ** p["m"])
          - p["ke"] * X * Z * activity + p["DX"] * diffusion(X, W))
    dY = (p["kn"] * (X ** p["m"]) + p["ks"] * X * (Z ** p["q"]) * activity
          - (p["cY"] + p["kf"]) * Y + p["DY"] * diffusion(Y, W))
    dZ = p["kf"] * Y + p["ke"] * X * Z * activity - p["cZ"] * Z
    dP = p.get("uP", 0.0) - p["cP"] * P
    dth = p["kon"] * P * (1.0 - th) - p["koff"] * th
    tox = hill(Y, p["KY"], p["hY"])
    spread = neighbor_damage(a, W)
    da = p["etaY"] * tox * (1.0 - a) + p["etaN"] * spread - p["deltaA"] * a

    return np.concatenate([dX, dY, dZ, dP, dth, da], axis=0)


def aggregation_efficiency(th: np.ndarray) -> np.ndarray:
    """Interface activity = fraction of surface-mediated aggregation retained."""
    return 1.0 - np.clip(th, 0.0, 1.0)


def langmuir_occupancy(P: float, K_D: float) -> float:
    """Quasi-steady-state passivator occupancy of fibril interfaces."""
    return P / (K_D + P)
