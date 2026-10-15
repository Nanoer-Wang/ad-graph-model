"""Steady-state, bifurcation and bistability analysis of the well-mixed node.

The network-averaged dynamics reduce (in the absence of spatial gradients)
to a single-node system.  Its steady states are found by multi-start Newton
iteration; stability follows from the eigenvalues of the Jacobian; saddle-node
folds are located by a scalar root find on det(J) along a branch.

A quasi-steady-state (QSS) reduction along the slow damage coordinate ``a``
is also provided: with (X, Y, Z) slaved to ``a``, the damage dynamics becomes
a scalar ODE whose roots are the fixed points and whose primitive is a
pseudo-potential whose wells are the stable states.
"""
from __future__ import annotations

import copy
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, fsolve

from .model import DEFAULT_PARAMETERS, hill, rhs

STATE_KEYS = ("X", "Y", "Z", "P", "th", "a")


# --------------------------------------------------------------------------
# single-node dynamics
# --------------------------------------------------------------------------
def _node_rhs(t, y, p):
    W = np.ones((1, 1))
    return rhs(t, y, 1, W, _ConstSchedule(p))


class _ConstSchedule:
    def __init__(self, p):
        self.p = dict(p)

    def at(self, t):
        return self.p


def node_f(y: np.ndarray, p: Dict[str, float]) -> np.ndarray:
    return _node_rhs(0.0, np.asarray(y, dtype=float), p)


def node_f4(y4: np.ndarray, p: Dict[str, float], P: float) -> np.ndarray:
    """Reduced node dynamics with passivator level P and occupancy th fixed
    at their quasi-steady values (P sustained, th = Langmuir)."""
    X, Y, Z, a = y4
    th = P / (p["koff"] / p["kon"] + P)
    act = 1.0 - th
    dX = (p["lambda0"] * (1 - a) * (1 + p["beta"] * a) - p["cX"] * X
          - p["kn"] * X ** p["m"] - p["ke"] * X * Z * act)
    dY = (p["kn"] * X ** p["m"] + p["ks"] * X * Z ** p["q"] * act
          - (p["cY"] + p["kf"]) * Y)
    dZ = p["kf"] * Y + p["ke"] * X * Z * act - p["cZ"] * Z
    tox = hill(np.array([Y]), p["KY"], p["hY"])[0]
    da = p["etaY"] * tox * (1.0 - a) - p["deltaA"] * a
    return np.array([dX, dY, dZ, da])


def node_jac4(y4: np.ndarray, p: Dict[str, float], P: float,
              h: float = 1e-7) -> np.ndarray:
    y4 = np.asarray(y4, dtype=float)
    n = len(y4)
    J = np.zeros((n, n))
    f0 = node_f4(y4, p, P)
    for j in range(n):
        dy = np.zeros(n)
        dy[j] = h * max(abs(y4[j]), 1.0)
        J[:, j] = (node_f4(y4 + dy, p, P) - f0) / dy[j]
    return J


def node_jac(y: np.ndarray, p: Dict[str, float], h: float = 1e-7) -> np.ndarray:
    y = np.asarray(y, dtype=float)
    n = len(y)
    J = np.zeros((n, n))
    f0 = node_f(y, p)
    for j in range(n):
        dy = np.zeros(n)
        dy[j] = h * max(abs(y[j]), 1.0)
        J[:, j] = (node_f(y + dy, p) - f0) / dy[j]
    return J


def _bounded_root(fun, y0, lo, hi, jac=None, n_restarts: int = 0):
    """Bounded least-squares root find (keeps variables in their physical box)."""
    from scipy.optimize import least_squares
    kw = {} if jac is None else {"jac": jac}
    r = least_squares(fun, y0, bounds=(lo, hi), xtol=1e-14,
                      ftol=1e-14, gtol=1e-14, max_nfev=2000, **kw)
    return r.x, float(np.max(np.abs(fun(r.x))))


def _healthy_start(p: Dict[str, float], P: float) -> np.ndarray:
    """Quasi-steady estimate of the low-damage (healthy) steady state.

    The exact healthy state is NOT (Xs, 0, 0, 0): primary nucleation
    continuously produces oligomer, so Y, Z > 0 and a small damage level
    balances toxicity against repair.  We estimate it from the slow-manifold
    scalings Y ~ kn X^m / (cY + kf), Z ~ kf Y / cZ, a ~ etaY H(Y) / (deltaA
    + etaY H(Y)) and X slightly below Xs.
    """
    th = P / (p["koff"] / p["kon"] + P)
    act = 1.0 - th
    X = p["lambda0"] / p["cX"]
    for _ in range(50):
        Y = p["kn"] * X ** p["m"] / (p["cY"] + p["kf"])
        Z = p["kf"] * Y / p["cZ"]
        Xn = (p["lambda0"] - p["kn"] * X ** p["m"]
              - p["ke"] * X * Z * act) / p["cX"]
        if abs(Xn - X) < 1e-12:
            X = Xn
            break
        X = 0.5 * X + 0.5 * Xn
    tox = hill(np.array([Y]), p["KY"], p["hY"])[0]
    a = p["etaY"] * tox / (p["deltaA"] + p["etaY"] * tox)
    return np.array([X, Y, Z, min(a, 1.0)])


def _slave_XYZ(p: Dict[str, float], P: float, a: float) -> np.ndarray:
    """(X, Y, Z) slaved to a fixed damage level ``a`` (healthy branch).

    Y and Z follow from the quasi-steady oligomer/fibril balances (fixed
    point, contractive near the healthy state); X solves the monotone
    monomer balance by bisection.
    """
    if abs(p["q"] - 1.0) < 1e-12:
        return fast_equilibrium_q1(a, p, P)
    th = P / (p["koff"] / p["kon"] + P)
    act = 1.0 - th
    prod = p["lambda0"] * (1.0 - a) * (1.0 + p["beta"] * a)

    def balance(X):
        Y, Z = 0.0, 0.0
        for _ in range(200):
            Yn = (p["kn"] * X ** p["m"]
                  + p["ks"] * X * Z ** p["q"] * act) / (p["cY"] + p["kf"])
            den = p["cZ"] - p["ke"] * X * act
            Zn = p["kf"] * Yn / den if den > 1e-12 else 0.0
            if abs(Yn - Y) < 1e-14 and abs(Zn - Z) < 1e-14:
                Y, Z = Yn, Zn
                break
            Y, Z = Yn, Zn
        return Y, Z

    def g(X):
        Y, Z = balance(X)
        return prod - p["cX"] * X - p["kn"] * X ** p["m"] - p["ke"] * X * Z * act

    lo, hi = 0.0, max(prod / p["cX"] * 1.5, 1e-6)
    glo = g(lo)
    if glo <= 0:
        X = lo
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if g(mid) > 0:
                lo = mid
            else:
                hi = mid
        X = 0.5 * (lo + hi)
    Y, Z = balance(X)
    return np.array([X, Y, Z])


def _healthy_branch_roots(p: Dict[str, float], P: float) -> List[np.ndarray]:
    """Low-damage steady states found by slaved reduction + bisection in a.

    The exact healthy state is not (Xs, 0, 0, 0) — primary nucleation keeps
    Y, Z positive — and plain least squares stalls near a ~ 0.01 because the
    damage equation is very flat there.  We therefore scan the damage
    coordinate with (X, Y, Z) slaved, bracket sign changes of da/dt, and
    refine by bisection; the resulting roots are polished in the full 4D
    system.
    """
    th = P / (p["koff"] / p["kon"] + P)
    act = 1.0 - th

    def da_of(a):
        X, Y, Z = _slave_XYZ(p, P, a)
        tox = hill(np.array([Y]), p["KY"], p["hY"])[0]
        return p["etaY"] * tox * (1.0 - a) - p["deltaA"] * a

    roots = []
    grid = np.linspace(0.0, 0.6, 241)
    vals = np.array([da_of(a) for a in grid])
    for i in range(len(grid) - 1):
        if vals[i] == 0.0:
            roots.append(grid[i])
            continue
        if vals[i] * vals[i + 1] < 0:
            lo, hi = grid[i], grid[i + 1]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if da_of(lo) * da_of(mid) <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(0.5 * (lo + hi))
    out = []
    for a in roots:
        X, Y, Z = _slave_XYZ(p, P, a)
        y0 = np.array([X, Y, Z, a])
        try:
            sol, resid = _bounded_root(lambda y: node_f4(y, p, P), y0,
                                       [0.0, 0.0, 0.0, 0.0],
                                       [np.inf, np.inf, np.inf, 1.0],
                                       jac=lambda y: node_jac4(y, p, P))
        except Exception:
            continue
        if resid <= 1e-9:
            out.append(sol)
    return out


def find_steady_states(p: Dict[str, float], P_fixed: float | None = None,
                       n_starts: int = 80, seed: int = 0
                       ) -> List[Tuple[np.ndarray, float]]:
    """Steady states of the well-mixed node.

    If ``P_fixed`` is given, the passivator level and interface occupancy are
    held at their sustained quasi-steady values and only (X, Y, Z, a) are
    solved (4-dimensional reduction).
    """
    rng = np.random.default_rng(seed)
    if P_fixed is None:
        lo = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
        hi = np.array([np.inf, np.inf, np.inf, np.inf, 1.0, 1.0])
        fun = lambda y: node_f(y, p)
        jac = lambda y: node_jac(y, p)
        starts = [np.zeros(6),
                  np.array([p["lambda0"] / p["cX"], 0.0, 0.0, 0.0, 0.0, 0.0]),
                  np.array([p["lambda0"] / p["cX"], 0.05, 0.05, 0.0, 0.0, 0.0]),
                  np.array([0.5, 0.10, 0.10, 0.0, 0.0, 0.30]),
                  np.array([0.3, 0.08, 0.20, 0.0, 0.0, 0.70])]
        for _ in range(n_starts):
            starts.append(rng.uniform([0, 0, 0, 0, 0, 0],
                                      [1.5, 1.0, 1.0, 2.0, 1.0, 1.0]))
    else:
        th = P_fixed / (p["koff"] / p["kon"] + P_fixed)
        lo = np.array([0.0, 0.0, 0.0, 0.0])
        hi = np.array([np.inf, np.inf, np.inf, 1.0])
        fun = lambda y: node_f4(y, p, P_fixed)
        jac = lambda y: node_jac4(y, p, P_fixed)
        starts = [np.zeros(4),
                  _healthy_start(p, P_fixed),
                  np.array([p["lambda0"] / p["cX"], 0.0, 0.0, 0.0]),
                  np.array([0.5, 0.10, 0.10, 0.30]),
                  np.array([0.3, 0.08, 0.20, 0.70])]
        # low-damage branch from the slaved reduction (robust near a ~ 0)
        starts.extend(_healthy_branch_roots(p, P_fixed))
        for _ in range(n_starts):
            starts.append(rng.uniform([0, 0, 0, 0], [1.5, 1.0, 1.0, 1.0]))

    found: List[Tuple[np.ndarray, float]] = []
    reduced: List[np.ndarray] = []
    for y0 in starts:
        try:
            sol, resid = _bounded_root(fun, y0, lo, hi, jac=jac)
        except Exception:
            continue
        if resid > 1e-9:
            continue
        if any(np.max(np.abs(sol - s)) < 1e-6 for s in reduced):
            continue
        eigs = np.linalg.eigvals(jac(sol))
        lmax = float(np.max(eigs.real))
        reduced.append(sol.copy())
        if P_fixed is not None:
            sol = np.array([sol[0], sol[1], sol[2], P_fixed, th, sol[3]])
        found.append((sol, lmax))
    return found

def continuation(cfg: dict, param_name: str, values: Sequence[float],
                 n_starts: int = 40) -> List[dict]:
    """Steady states along a one-parameter sweep (well-mixed node)."""
    out = []
    for v in values:
        p = dict(cfg["parameters"])
        p[param_name] = float(v)
        states = find_steady_states(p, n_starts=n_starts)
        out.append({"value": float(v), "states": states})
    return out


def find_fold(cfg: dict, param_name: str, branch: Sequence[Tuple[float, np.ndarray]],
              bracket: Tuple[float, float]) -> Optional[float]:
    """Refine a saddle-node fold on a branch by bisection on max Re(eig)."""
    def lmax(v: float) -> float:
        p = dict(cfg["parameters"])
        p[param_name] = float(v)
        y0 = _interp_branch(branch, v)
        sol, info, ier, _ = fsolve(node_f, y0, args=(p,), full_output=True)
        if ier != 1:
            return np.nan
        J = node_jac(sol, p)
        return float(np.max(np.linalg.eigvals(J).real))

    lo, hi = bracket
    flo, fhi = lmax(lo), lmax(hi)
    if not (np.isfinite(flo) and np.isfinite(fhi)):
        return None
    if flo * fhi > 0:
        return None
    try:
        return float(brentq(lmax, lo, hi, xtol=1e-8))
    except Exception:
        return None


def _interp_branch(branch, v):
    arr = np.array([[val, *y] for val, y in branch])
    order = np.argsort(arr[:, 0])
    arr = arr[order]
    return np.array([np.interp(v, arr[:, 0], arr[:, j]) for j in range(1, arr.shape[1])])


# --------------------------------------------------------------------------
# QSS reduction along the damage coordinate
# --------------------------------------------------------------------------
def fast_equilibrium_q1(a: float, p: Dict[str, float], P: float) -> np.ndarray:
    """Positive fast equilibrium for q=1, including monomer depletion.

    Eliminate Y and Z, then bracket the strictly decreasing monomer balance
    below its autocatalytic pole. Failed unconstrained root searches must
    never be interpreted as nonexistence of a slow manifold.
    """
    if abs(p["q"] - 1.0) > 1e-12:
        raise ValueError("Analytic fast elimination requires q=1")
    prod = p["lambda0"] * (1-a) * (1+p["beta"]*a)
    if prod <= 0:
        return np.zeros(3)
    act = 1.0 - P / (p["koff"] / p["kon"] + P)
    cy = p["cY"] + p["kf"]
    gain = act * (cy*p["ke"] + p["ks"]*p["kf"])
    pole = cy*p["cZ"] / gain if gain > 0 else np.inf

    def yz(x):
        z = p["kf"]*p["kn"]*x**p["m"] / (cy*p["cZ"]-gain*x)
        y = (p["cZ"]-p["ke"]*x*act)*z/p["kf"]
        return y, z

    def balance(x):
        _, z = yz(x)
        return prod-p["cX"]*x-p["kn"]*x**p["m"]-p["ke"]*x*z*act

    hi = min(prod/p["cX"], pole*(1-1e-12))
    x = brentq(balance, 0.0, hi, xtol=1e-14, rtol=1e-14)
    y, z = yz(x)
    return np.array([x, y, z])


def qss_damage_rate(a: float, p: Dict[str, float], P: float) -> float:
    """f(a) = da/dt with (X,Y,Z,th) slaved to a (th from sustained P)."""
    if abs(p["q"] - 1.0) < 1e-12:
        _, Y, _ = fast_equilibrium_q1(a, p, P)
        tox = hill(Y, p["KY"], p["hY"])
        return float(p["etaY"]*tox*(1-a)-p["deltaA"]*a)
    th = P / (p["koff"] / p["kon"] + P)
    act = 1.0 - th

    def fast(y):
        X, Y, Z = y
        dX = (p["lambda0"] * (1 - a) * (1 + p["beta"] * a) - p["cX"] * X
              - p["kn"] * X ** p["m"] - p["ke"] * X * Z * act)
        dY = (p["kn"] * X ** p["m"] + p["ks"] * X * Z ** p["q"] * act
              - (p["cY"] + p["kf"]) * Y)
        dZ = p["kf"] * Y + p["ke"] * X * Z * act - p["cZ"] * Z
        return np.array([dX, dY, dZ])

    lo = np.array([0.0, 0.0, 0.0])
    hi = np.array([np.inf, np.inf, np.inf])
    best, best_res = None, np.inf
    for y0 in (np.array([0.3, 0.0, 0.0]), np.array([0.8, 0.05, 0.05]),
               np.array([0.5, 0.2, 0.2]), np.array([0.6, 0.02, 0.3])):
        try:
            sol, resid = _bounded_root(fast, y0, lo, hi)
        except Exception:
            continue
        if resid < best_res:
            best, best_res = sol, resid
    if best is None or best_res > 1e-9:
        return np.nan
    X, Y, Z = best
    tox = hill(np.array([Y]), p["KY"], p["hY"])[0]
    return p["etaY"] * tox * (1.0 - a) - p["deltaA"] * a


def qss_branches(p: Dict[str, float], P: float = 0.0,
                 n_grid: int = 400) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Fixed points of the QSS damage ODE on a grid in a."""
    a_grid = np.linspace(0.0, 1.0, n_grid)
    f = np.array([qss_damage_rate(a, p, P) for a in a_grid])
    roots, kinds = [], []
    for i in range(len(a_grid) - 1):
        f0, f1 = f[i], f[i + 1]
        if not (np.isfinite(f0) and np.isfinite(f1)):
            continue
        if f0 == 0.0:
            roots.append(a_grid[i]); kinds.append("tangent" if abs(f1) < 1e-6 else "cross")
        elif f0 * f1 < 0:
            try:
                r = brentq(lambda x: qss_damage_rate(x, p, P), a_grid[i], a_grid[i + 1])
                roots.append(r)
                kinds.append("cross")
            except Exception:
                pass
    roots = np.array(roots)
    # stability: negative slope of f at the root => stable
    stable = np.zeros(len(roots), dtype=bool)
    for i, r in enumerate(roots):
        h = 1e-4
        fp = (qss_damage_rate(min(r + h, 1.0), p, P) - qss_damage_rate(max(r - h, 0.0), p, P)) / (2 * h)
        stable[i] = fp < 0
    return roots, np.array(stable), f


def pseudo_potential(a_grid: np.ndarray, f: np.ndarray) -> np.ndarray:
    """Phi(a) = -int f da (wells = stable fixed points).

    Nonfinite entries mark unresolved numerical evaluations, not proof of
    nonexistence. Integrate valid segments independently without bridging
    missing values. For q=1 the analytic elimination resolves the full range.
    """
    a_grid = np.asarray(a_grid, dtype=float)
    f = np.asarray(f, dtype=float)
    phi = np.full(f.shape, np.nan)
    valid = np.isfinite(f)
    if not valid.any():
        return phi
    change = np.flatnonzero(valid[1:] != valid[:-1]) + 1
    bounds = np.concatenate([[0], change, [len(f)]])
    for s, e in zip(bounds[:-1], bounds[1:]):
        if not valid[s] or e - s < 2:
            continue
        inc = 0.5 * (f[s + 1:e] + f[s:e - 1]) * np.diff(a_grid[s:e])
        phi[s:e] = -np.concatenate([[0.0], np.cumsum(inc)])
    return phi


# --------------------------------------------------------------------------
# full-node trajectories for bistability verification
# --------------------------------------------------------------------------
def simulate_node(cfg: dict, p: Optional[Dict[str, float]] = None,
                  y0: Optional[np.ndarray] = None, t_end: float = 60.0,
                  n: int = 1200) -> Tuple[np.ndarray, np.ndarray]:
    p = dict(p or cfg["parameters"])
    if y0 is None:
        y0 = np.array([p["lambda0"] / p["cX"] * 1.05, 0.0, 0.0, 0.0, 0.0, 0.0])
    t_eval = np.linspace(0.0, t_end, n)
    sol = solve_ivp(lambda t, y: _node_rhs(t, y, p), (0.0, t_end), y0,
                    t_eval=t_eval, method="LSODA", rtol=1e-8, atol=1e-10)
    return sol.t, sol.y.T
