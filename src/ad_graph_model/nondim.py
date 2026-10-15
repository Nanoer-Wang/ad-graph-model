"""Nondimensionalization and dimensionless control groups.

Scales (chosen once, used throughout the manuscript):

    X_s = lambda0 / cX      healthy monomer steady state
    Y_s = Z_s = X_s         mass-consistent concentration scale
    tau = cX * t            fast clearance time scale

Dimensionless groups reported by :func:`dimensionless_groups`.
"""
from __future__ import annotations

from typing import Dict

import numpy as np


def reference_scales(p: Dict[str, float]) -> Dict[str, float]:
    X_s = p["lambda0"] / p["cX"]
    return {"X_s": X_s, "Y_s": X_s, "Z_s": X_s, "tau": p["cX"]}


def dimensionless_groups(p: Dict[str, float]) -> Dict[str, float]:
    s = reference_scales(p)
    X_s, Y_s, Z_s = s["X_s"], s["Y_s"], s["Z_s"]
    cX = p["cX"]

    rho = p["kn"] * X_s ** (p["m"] - 1) / cX            # primary nucleation
    sigma = p["ks"] * X_s * Z_s ** (p["q"] - 1) / cX    # secondary nucleation
    eps = p["ke"] * Z_s / cX                            # elongation
    gamma = p["cY"] / cX                                # oligomer clearance
    phi = p["kf"] / cX                                  # oligomer->fibril
    zeta = p["cZ"] / cX                                 # fibril clearance
    kappa = p["DX"] / cX                                # monomer transport
    lam = p["DY"] / cX                                  # oligomer transport
    K = p["KY"] / Y_s                                   # toxicity threshold
    eY = p["etaY"] / cX                                 # toxicity->damage
    eN = p["etaN"] / cX                                 # nonlocal damage
    delta = p["deltaA"] / cX                            # damage repair
    KD = p["koff"] / p["kon"]                           # passivator K_D (in P units)
    # Cascade number: loop gain of the autocatalytic oligomer-fibril loop.
    #
    #   R = X_s [k_s k_f + k_e (cY + kf)] / [(cY + kf) cZ]
    #     = (secondary-nucleation amplification) + (elongation amplification)
    #
    # It is the composite that appears in the critical passivation fraction:
    # theta_c = 1 - 1/R.  Hence R > 1  <=>  the untreated fibril-poor state is
    # unstable, and R < 1  <=>  clearance dominates and no passivation is
    # needed.  (An earlier version used ks*Zs**(q-1) in the numerator, which
    # omits one factor of Xs and overestimates R by 27% at q = 1.)
    R = (X_s * (p["ks"] * p["kf"] + p["ke"] * (p["cY"] + p["kf"]))
         ) / ((p["cY"] + p["kf"]) * p["cZ"])
    return {
        "rho": rho, "sigma": sigma, "eps": eps, "gamma": gamma, "phi": phi,
        "zeta": zeta, "kappa": kappa, "lambda": lam, "K": K, "H": p["hY"],
        "eps_Y": eY, "eps_N": eN, "delta": delta, "K_D": KD,
        "R_cascade": R,
    }


GROUP_DESCRIPTIONS = {
    "rho": "primary-nucleation rate (monomer clearance units)",
    "sigma": "secondary-nucleation rate on fibril interfaces",
    "eps": "fibril elongation rate",
    "gamma": "oligomer clearance rate",
    "phi": "oligomer-to-fibril conversion rate",
    "zeta": "fibril clearance rate",
    "kappa": "monomer network-transport rate",
    "lambda": "oligomer network-transport rate",
    "K": "normalized toxicity threshold of the Hill response",
    "H": "Hill coefficient (toxicity ultrasensitivity)",
    "eps_Y": "toxicity-driven damage rate",
    "eps_N": "nonlocal (network) damage-transmission rate",
    "delta": "damage-repair rate",
    "K_D": "passivator dissociation constant (interface affinity)",
    "R_cascade": "cascade number: oligomer production vs removal at reference state",
}
