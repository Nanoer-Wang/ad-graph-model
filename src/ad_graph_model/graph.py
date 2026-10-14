"""Graph/network construction for the amyloid cascade model.

Supports three network types used in the manuscript:

1. ``sheet``       : 2D cortical sheet abstraction (regular grid with
                     distance-decaying weights).
2. ``smallworld``  : Watts-Strogatz small-world network (topology control).
3. ``connectome``  : Surrogate structural connectome with 68 Desikan-like
                     regions, distance-dependent weights and hub structure,
                     calibrated to macro-scale connectome statistics
                     (see graph_stats).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import numpy as np


@dataclass
class GraphData:
    W: np.ndarray                 # weighted adjacency (row-normalized)
    positions: np.ndarray         # (n, 2) or (n, 3) coordinates for plotting
    rows: int
    cols: int
    kind: str
    labels: list | None = None    # region names for connectome
    stats: Dict | None = None     # connectome statistics


# --------------------------------------------------------------------------
# 2D sheet
# --------------------------------------------------------------------------
def build_grid_graph(rows: int, cols: int, diagonal: bool = True,
                     weight_mode: str = "distance") -> np.ndarray:
    n = rows * cols
    W = np.zeros((n, n), dtype=float)

    def idx(r: int, c: int) -> int:
        return r * cols + c

    for r in range(rows):
        for c in range(cols):
            i = idx(r, c)
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    if not diagonal and abs(dr) + abs(dc) != 1:
                        continue
                    rr, cc = r + dr, c + dc
                    if 0 <= rr < rows and 0 <= cc < cols:
                        dist = (dr * dr + dc * dc) ** 0.5
                        w = 1.0 if weight_mode == "binary" else 1.0 / dist
                        W[i, idx(rr, cc)] = w
    return W


# --------------------------------------------------------------------------
# Small-world network (Watts-Strogatz with geometric weighting)
# --------------------------------------------------------------------------
def build_smallworld(n: int, k: int, p: float, seed: int = 0) -> np.ndarray:
    rng = np.random.default_rng(seed)
    W = np.zeros((n, n), dtype=float)
    half = k // 2
    for i in range(n):
        for j in range(1, half + 1):
            w = 1.0 / j  # distance-decaying ring weights
            W[i, (i + j) % n] = w
            W[i, (i - j) % n] = w
    # Watts-Strogatz rewiring
    for i in range(n):
        for j in range(1, half + 1):
            if rng.random() < p:
                target = (i + j) % n
                if W[i, target] > 0:
                    W[i, target] = 0.0
                    new = rng.integers(0, n)
                    if new != i:
                        W[i, new] = 1.0 / max(abs(new - i), 1)
    W = np.maximum(W, W.T)
    return W


# --------------------------------------------------------------------------
# Surrogate structural connectome (68 Desikan-like regions)
# --------------------------------------------------------------------------
# Approximate 3D coordinates (x: posterior->anterior, y: inferior->superior,
# z: left->right) for the 34 Desikan-Killiany regions per hemisphere.
_REGION_COORDS: Dict[str, tuple] = {
    "bankssts": (-0.55, 0.05, 0.0),
    "caudalanteriorcingulate": (-0.15, 0.55, 0.0),
    "caudalmiddlefrontal": (-0.55, 0.75, 0.0),
    "cuneus": (0.55, 0.75, 0.0),
    "entorhinal": (-0.45, -0.45, 0.0),
    "frontalpole": (0.85, 0.75, 0.0),
    "fusiform": (0.30, -0.75, 0.0),
    "inferiorparietal": (0.15, 0.65, 0.0),
    "inferiortemporal": (0.35, -0.55, 0.0),
    "insula": (-0.35, 0.10, 0.0),
    "isthmuscingulate": (-0.02, 0.42, 0.0),
    "lateraloccipital": (0.60, -0.05, 0.0),
    "lateralorbitofrontal": (0.55, -0.30, 0.0),
    "lingual": (0.55, -0.80, 0.0),
    "medialorbitofrontal": (0.60, -0.35, 0.0),
    "middletemporal": (0.40, -0.45, 0.0),
    "paracentral": (-0.05, 0.85, 0.0),
    "parahippocampal": (-0.35, -0.55, 0.0),
    "parsopercularis": (-0.65, 0.45, 0.0),
    "parsorbitalis": (-0.60, -0.15, 0.0),
    "parstriangularis": (-0.60, 0.35, 0.0),
    "pericalcarine": (0.70, 0.10, 0.0),
    "postcentral": (-0.25, 0.60, 0.0),
    "posteriorcingulate": (-0.10, 0.35, 0.0),
    "precentral": (-0.45, 0.75, 0.0),
    "precuneus": (0.05, 0.55, 0.0),
    "rostralanteriorcingulate": (0.30, 0.55, 0.0),
    "rostralmiddlefrontal": (-0.25, 0.85, 0.0),
    "superiorfrontal": (-0.10, 0.90, 0.0),
    "superiorparietal": (0.20, 0.80, 0.0),
    "superiortemporal": (0.05, -0.35, 0.0),
    "supramarginal": (0.05, 0.50, 0.0),
    "temporalpole": (0.80, -0.50, 0.0),
    "transversetemporal": (-0.50, -0.20, 0.0),
}

# Regions treated as high-degree hubs (association/limbic core).
_HUB_REGIONS = {"precuneus", "posteriorcingulate", "isthmuscingulate",
                "rostralanteriorcingulate", "superiorfrontal", "middletemporal",
                "inferiorparietal", "lateraloccipital"}

# Earliest-affected regions in sporadic AD (Braak I-II transentorhinal).
_EARLY_REGIONS = {"entorhinal", "parahippocampal", "transversetemporal",
                  "bankssts"}


def _desikan_labels() -> list:
    names = sorted(_REGION_COORDS.keys())
    labels = [f"L-{n}" for n in names] + [f"R-{n}" for n in names]
    return labels


def _desikan_positions() -> np.ndarray:
    names = sorted(_REGION_COORDS.keys())
    coords = []
    for n in names:
        x, y, _ = _REGION_COORDS[n]
        coords.append((x, y, -0.25))
    for n in names:
        x, y, _ = _REGION_COORDS[n]
        coords.append((x, y, 0.25))
    return np.asarray(coords, dtype=float)


def build_connectome(n_scale: float = 1.0, d0: float = 0.45, seed: int = 0,
                     cross_hemisphere: float = 0.45, hub_boost: float = 2.5,
                     d_max: float = 0.60, w_min: float = 0.02,
                     k_nn: int = 3, p_rewire: float = 0.10) -> tuple:
    """Surrogate structural connectome with realistic macro-scale statistics.

    Construction (exponential-distance rule with spatial rewiring, cf.
    Roberts et al. 2016; Vázquez-Rodríguez et al. 2019):

    1. Pairwise weight  w_ij = exp(-d_ij/d0) * cross_hemisphere * hub_boost
       * lognormal noise, restricted to d_ij <= d_max and w_ij >= w_min,
       which yields a sparse graph with realistic density (~0.2) rather
       than a fully connected noise matrix.
    2. Watts-Strogatz-style rewiring with probability ``p_rewire`` replaces
       a fraction of short geometric edges by random long-range edges
       (association pathways), producing small-world behaviour: high
       clustering with a short characteristic path length.
    3. k-nearest-neighbour backstop: every region keeps at least ``k_nn``
       neighbours so that no region is isolated.

    The returned matrix is un-normalized; ``build_graph`` applies row
    normalization (conservative transport) afterwards.
    """
    rng = np.random.default_rng(seed)
    labels = _desikan_labels()
    pos = _desikan_positions()
    n = len(labels)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = float(np.linalg.norm(pos[i] - pos[j]))
            D[i, j] = d
            D[j, i] = d
    is_hub = np.array([lab[2:] in _HUB_REGIONS for lab in labels])
    cross = np.array([[labels[i][0] != labels[j][0] for j in range(n)]
                      for i in range(n)])

    def _base(i: int, j: int) -> float:
        b = np.exp(-D[i, j] / d0)
        if cross[i, j]:
            b *= cross_hemisphere
        if is_hub[i] or is_hub[j]:
            b *= hub_boost
        return float(b)

    W = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] > d_max:
                continue
            w = _base(i, j) * float(np.exp(rng.normal(0.0, 0.45)))
            if w < w_min:
                continue
            W[i, j] = w
            W[j, i] = w
    # WS-style rewiring -> long-range shortcuts (association fibers)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n) if W[i, j] > 0]
    for (i, j) in edges:
        if rng.random() >= p_rewire:
            continue
        if (W[i] > 0).sum() <= k_nn or (W[j] > 0).sum() <= k_nn:
            continue
        new = int(rng.integers(0, n))
        if new in (i, j) or W[i, new] > 0:
            continue
        w = W[i, j]
        W[i, j] = 0.0
        W[j, i] = 0.0
        W[i, new] = w
        W[new, i] = w
    # kNN backstop: guarantee a minimum degree of k_nn for every region
    for i in range(n):
        if (W[i] > 0).sum() >= k_nn:
            continue
        for j in np.argsort(D[i]):
            if j == i or W[i, j] > 0:
                continue
            w = max(_base(i, j), w_min)
            W[i, j] = w
            W[j, i] = w
            if (W[i] > 0).sum() >= k_nn:
                break
    return W, pos, labels


def graph_stats(W: np.ndarray) -> Dict:
    """Binary-adjacency graph statistics.

    Clustering coefficient (Watts-Strogatz):
        C_i = 2 T_i / [k_i (k_i - 1)],  T_i = (A^3)_ii / 2
    Path length: mean geodesic distance in hops (Floyd-Warshall).
    """
    n = W.shape[0]
    A = (W > 0).astype(float)
    np.fill_diagonal(A, 0.0)
    deg = A.sum(axis=1)
    A2 = A @ A
    tri = np.diag(A2 @ A) / 2.0
    with np.errstate(divide="ignore", invalid="ignore"):
        cc = np.where(deg > 1, 2.0 * tri / (deg * (deg - 1.0)), 0.0)
    dist = np.where(A > 0, 1.0, np.inf)
    np.fill_diagonal(dist, 0.0)
    for kk in range(n):
        dist = np.minimum(dist, dist[:, kk][:, None] + dist[kk, :][None, :])
    finite = np.isfinite(dist) & (dist > 0)
    L = float(dist[finite].mean()) if finite.any() else float("nan")
    return {
        "n_nodes": int(n),
        "mean_degree": float(deg.mean()),
        "density": float(A.sum() / (n * (n - 1))),
        "mean_clustering": float(np.nanmean(cc)),
        "path_length": L,
    }


# --------------------------------------------------------------------------
# Utilities
# --------------------------------------------------------------------------
def row_normalize(W: np.ndarray) -> np.ndarray:
    row_sum = W.sum(axis=1, keepdims=True)
    out = W.copy()
    mask = row_sum[:, 0] > 0
    out[mask] = out[mask] / row_sum[mask]
    return out


def build_graph(graph_cfg: dict) -> GraphData:
    gtype = graph_cfg.get("type", "sheet")
    if gtype == "grid":  # backward-compatible alias
        gtype = "sheet"
    if gtype == "sheet":
        rows = int(graph_cfg.get("rows", 20))
        cols = int(graph_cfg.get("cols", 20))
        diagonal = bool(graph_cfg.get("diagonal", True))
        weight_mode = str(graph_cfg.get("weight_mode", "distance"))
        W = build_grid_graph(rows, cols, diagonal, weight_mode)
        positions = np.array([[c, r] for r in range(rows) for c in range(cols)],
                             dtype=float)
        if graph_cfg.get("normalize_weights", True):
            W = row_normalize(W)
        return GraphData(W=W, positions=positions, rows=rows, cols=cols,
                         kind="sheet", labels=None, stats=None)
    if gtype == "smallworld":
        n = int(graph_cfg.get("n", 68))
        k = int(graph_cfg.get("k", 6))
        p = float(graph_cfg.get("p", 0.12))
        W = build_smallworld(n, k, p, seed=int(graph_cfg.get("seed", 0)))
        theta = np.linspace(0, 2 * np.pi, n, endpoint=False)
        positions = np.stack([np.cos(theta), np.sin(theta)], axis=1)
        if graph_cfg.get("normalize_weights", True):
            W = row_normalize(W)
        stats = graph_stats(W)
        return GraphData(W=W, positions=positions, rows=1, cols=n,
                         kind="smallworld", labels=None, stats=stats)
    if gtype == "connectome":
        W, pos, labels = build_connectome(
            d0=float(graph_cfg.get("d0", 0.45)),
            seed=int(graph_cfg.get("seed", 0)),
            cross_hemisphere=float(graph_cfg.get("cross_hemisphere", 0.45)),
        )
        if graph_cfg.get("normalize_weights", True):
            W = row_normalize(W)
        stats = graph_stats(W)
        return GraphData(W=W, positions=pos, rows=1, cols=len(labels),
                         kind="connectome", labels=labels, stats=stats)
    raise ValueError(f"Unsupported graph type: {gtype}")
