"""Rebuild the network-level spatial statistics used by Section 4.2 / Fig. 3.

Two problems in the original pipeline are fixed here.

1.  The hub test used ``strength = W.sum(axis=1)`` on the **row-normalised**
    adjacency, so every node had strength identically 1 (to floating-point
    precision).  The reported correlation was therefore computed on numerical
    noise and carries no information.  Here the hubness measures are computed
    from the *un-normalised* connectome: weighted strength, degree, and
    eigenvector centrality.

2.  The manuscript quotes "geodesic distance from the entorhinal seed predicts
    onset (r = 0.60, p = 7e-8)" but no code ever computed it.  Here the
    shortest-path hop distance from the seed is computed by BFS on the binary
    adjacency and correlated with the regional onset times.

Outputs (results_v2/data/):
    spatial_hub_metrics.csv   region, degree, strength_raw, eigenvector, hop
    network_distance_stats.json  Pearson/Spearman r and p for every pairing
"""

from __future__ import annotations

import json
import sys
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "src"))

from ad_graph_model.graph import build_connectome, graph_stats  # noqa: E402

DATA = ROOT / "results_v2" / "data"


def hop_distances(A: np.ndarray, seeds: list[int]) -> np.ndarray:
    """Unweighted shortest-path distance (in hops) from the nearest seed."""
    n = A.shape[0]
    dist = np.full(n, np.inf)
    q: deque[int] = deque()
    for s in seeds:
        dist[s] = 0.0
        q.append(s)
    while q:
        u = q.popleft()
        for v in np.nonzero(A[u])[0]:
            if dist[v] == np.inf:
                dist[v] = dist[u] + 1.0
                q.append(int(v))
    return dist


def eigenvector_centrality(W: np.ndarray, iters: int = 500) -> np.ndarray:
    """Power iteration on the symmetric weighted adjacency."""
    n = W.shape[0]
    v = np.ones(n) / n
    for _ in range(iters):
        w = W @ v
        nrm = np.linalg.norm(w)
        if nrm == 0:
            break
        v_new = w / nrm
        if np.linalg.norm(v_new - v) < 1e-12:
            v = v_new
            break
        v = v_new
    return v / v.sum()


def partial_corr(x, y, z):
    """Pearson correlation between x and y with the linear effect of z removed."""
    x, y, z = map(lambda v: np.asarray(v, dtype=float), (x, y, z))
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(z)
    x, y, z = x[ok], y[ok], z[ok]
    Z = np.column_stack([np.ones_like(z), z])
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    n = len(x)
    dof = n - 3
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    p = float(2 * stats.t.sf(abs(t), dof))
    return {"n": n, "partial_r": r, "partial_p": p}


def partial_corr_multi(y, x, zs):
    """Partial correlation of x and y controlling for several covariates."""
    y, x = np.asarray(y, float), np.asarray(x, float)
    Z = np.column_stack([np.asarray(z, float) for z in zs])
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(Z).all(axis=1)
    x, y, Z = x[ok], y[ok], Z[ok]
    Z = np.column_stack([np.ones(len(Z)), Z])
    rx = x - Z @ np.linalg.lstsq(Z, x, rcond=None)[0]
    ry = y - Z @ np.linalg.lstsq(Z, y, rcond=None)[0]
    r = float(np.corrcoef(rx, ry)[0, 1])
    dof = len(x) - Z.shape[1] - 1
    t = r * np.sqrt(dof / max(1e-12, 1 - r ** 2))
    return {"n": int(len(x)), "partial_r": r,
            "partial_p": float(2 * stats.t.sf(abs(t), dof))}


def corr_report(x: np.ndarray, y: np.ndarray) -> dict:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[ok], np.asarray(y)[ok]
    pr, pp = stats.pearsonr(x, y)
    sr, sp = stats.spearmanr(x, y)
    return {"n": int(ok.sum()), "pearson_r": float(pr), "pearson_p": float(pp),
            "spearman_r": float(sr), "spearman_p": float(sp)}


def main() -> None:
    W_raw, pos, labels = build_connectome(seed=0)
    A = (W_raw > 0).astype(float)
    np.fill_diagonal(A, 0.0)

    degree = A.sum(axis=1)
    strength_raw = W_raw.sum(axis=1)
    eig = eigenvector_centrality(W_raw)
    seed_idx = [i for i, lab in enumerate(labels) if lab.endswith("entorhinal")]
    hop = hop_distances(A, seed_idx)

    onsets = pd.read_csv(DATA / "spatial_onsets.csv")
    tab = pd.DataFrame({"region": labels, "degree": degree,
                        "strength_raw": strength_raw,
                        "eigenvector": eig, "hop": hop})
    tab = tab.merge(onsets[["region", "onset", "a_end"]], on="region",
                    how="left")
    tab.to_csv(DATA / "spatial_hub_metrics.csv", index=False)

    out = {"seed_regions": [labels[i] for i in seed_idx],
           "n_regions": int(len(labels)),
           "graph_stats_unweighted": {k: float(v) for k, v in
                                      graph_stats(W_raw).items()}}
    for pred, col in [("hop", "hop"), ("degree", "degree"),
                      ("strength_raw", "strength_raw"),
                      ("eigenvector", "eigenvector")]:
        out[f"{pred}_vs_onset"] = corr_report(tab[col], tab["onset"])
        out[f"{pred}_vs_a_end"] = corr_report(tab[col], tab["a_end"])
    # log strength (the manuscript reports the log-scaled hub test)
    out["log_strength_vs_a_end"] = corr_report(np.log(tab["strength_raw"]),
                                               tab["a_end"])
    # partial correlations: hubness vs onset / damage, controlling for hop
    # distance from the seed (the manuscript's "once distance is accounted
    # for" claim in Section 5.2).
    for name, col in [("degree", "degree"), ("strength_raw", "strength_raw"),
                      ("eigenvector", "eigenvector")]:
        out[f"{name}_vs_onset_given_hop"] = partial_corr(
            tab[col], tab["onset"], tab["hop"])
        out[f"{name}_vs_a_end_given_hop"] = partial_corr(
            tab[col], tab["a_end"], tab["hop"])
    out["hop_vs_onset_given_strength"] = partial_corr(
        tab["hop"], tab["onset"], np.log(tab["strength_raw"]))
    # Euclidean (physical) distance from the seed: separates "hubs happen to
    # lie farther from the entorhinal cortex" from a genuine topological
    # hub effect.
    euclid = np.min(np.linalg.norm(pos - pos[seed_idx][:, None, :], axis=2),
                    axis=0)
    tab["euclid"] = euclid
    out["euclid_vs_onset"] = corr_report(euclid, tab["onset"])
    out["hop_vs_onset_given_euclid"] = partial_corr(
        tab["hop"], tab["onset"], euclid)
    out["strength_vs_onset_given_euclid"] = partial_corr(
        tab["strength_raw"], tab["onset"], euclid)
    out["strength_vs_onset_given_hop_euclid"] = partial_corr_multi(
        tab["strength_raw"], tab["onset"], [tab["hop"], euclid])
    out["degree_vs_onset_given_hop_euclid"] = partial_corr_multi(
        tab["degree"], tab["onset"], [tab["hop"], euclid])

    with open(DATA / "network_distance_stats.json", "w") as fh:
        json.dump(out, fh, indent=1)

    print(json.dumps(out, indent=1))
    print("\nwrote", DATA / "spatial_hub_metrics.csv")
    print("wrote", DATA / "network_distance_stats.json")


if __name__ == "__main__":
    main()
