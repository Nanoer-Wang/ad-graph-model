"""Simulation driver for the extended passivation model."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np
from scipy.integrate import solve_ivp

from .graph import GraphData, build_graph
from .model import DEFAULT_PARAMETERS, ParameterSchedule, rhs


@dataclass
class SimulationResult:
    t: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    Z: np.ndarray
    P: np.ndarray
    th: np.ndarray
    a: np.ndarray
    W: np.ndarray
    rows: int
    cols: int
    kind: str
    labels: list | None
    positions: np.ndarray
    config: Dict


def build_initial_state(cfg: Dict, n_nodes: int) -> np.ndarray:
    init = cfg.get("initial_conditions", {})
    X = np.full(n_nodes, float(init.get("X", 0.0)))
    Y = np.full(n_nodes, float(init.get("Y", 0.0)))
    Z = np.full(n_nodes, float(init.get("Z", 0.0)))
    P = np.full(n_nodes, float(init.get("P", 0.0)))
    th = np.full(n_nodes, float(init.get("th", 0.0)))
    a = np.full(n_nodes, float(init.get("a", 0.0)))

    seed_nodes = cfg.get("seed_nodes", [0])
    seed_values = cfg.get("seed_values", {})
    for idx in seed_nodes:
        idx = int(idx)
        if not (0 <= idx < n_nodes):
            raise IndexError(f"seed node {idx} out of range for n_nodes={n_nodes}")
        if "X" in seed_values:
            X[idx] = float(seed_values["X"])
        if "Y" in seed_values:
            Y[idx] = float(seed_values["Y"])
        if "Z" in seed_values:
            Z[idx] = float(seed_values["Z"])
        if "P" in seed_values:
            P[idx] = float(seed_values["P"])
        if "th" in seed_values:
            th[idx] = float(seed_values["th"])
        if "a" in seed_values:
            a[idx] = float(seed_values["a"])

    return np.concatenate([X, Y, Z, P, th, a], axis=0)


def run_simulation(cfg: Dict) -> SimulationResult:
    graph: GraphData = build_graph(cfg["graph"])
    n_nodes = graph.W.shape[0]
    y0 = build_initial_state(cfg, n_nodes)

    base = {k: float(v) for k, v in
            {**DEFAULT_PARAMETERS, **cfg.get("parameters", {})}.items()}
    schedule = ParameterSchedule(base=base,
                                 interventions=cfg.get("interventions", []))

    sim = cfg["simulation"]
    t_end = float(sim.get("t_end", 30.0))
    n_eval = int(sim.get("n_eval", 801))
    t_eval = np.linspace(0.0, t_end, n_eval)

    sol = solve_ivp(
        fun=lambda t, y: rhs(t, y, n_nodes=n_nodes, W=graph.W,
                             schedule=schedule),
        t_span=(0.0, t_end),
        y0=y0,
        t_eval=t_eval,
        method=str(sim.get("method", "LSODA")),
        rtol=float(sim.get("rtol", 1e-6)),
        atol=float(sim.get("atol", 1e-8)),
    )
    if not sol.success:
        raise RuntimeError(f"ODE solver failed: {sol.message}")

    X = np.clip(sol.y[0:n_nodes, :].T, 0.0, None)
    Y = np.clip(sol.y[n_nodes:2 * n_nodes, :].T, 0.0, None)
    Z = np.clip(sol.y[2 * n_nodes:3 * n_nodes, :].T, 0.0, None)
    P = np.clip(sol.y[3 * n_nodes:4 * n_nodes, :].T, 0.0, None)
    th = np.clip(sol.y[4 * n_nodes:5 * n_nodes, :].T, 0.0, 1.0)
    a = np.clip(sol.y[5 * n_nodes:6 * n_nodes, :].T, 0.0, 1.0)

    return SimulationResult(t=sol.t, X=X, Y=Y, Z=Z, P=P, th=th, a=a,
                            W=graph.W, rows=graph.rows, cols=graph.cols,
                            kind=graph.kind, labels=graph.labels,
                            positions=graph.positions, config=cfg)
