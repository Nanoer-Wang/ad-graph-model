from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .simulate import SimulationResult


plt.rcParams.update(
    {
        "font.size": 11,
        "axes.linewidth": 1.0,
        "figure.dpi": 200,
        "savefig.dpi": 300,
    }
)


def save_node_trajectories(res: SimulationResult, node: int, outpath: str | Path) -> None:
    outpath = Path(outpath)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(res.t, res.X[:, node], label="X (Aβ monomer)", linewidth=2)
    ax.plot(res.t, res.Y[:, node], label="Y (toxic oligomer)", linewidth=2)
    ax.plot(res.t, res.Z[:, node], label="Z (plaque burden)", linewidth=2)
    ax.plot(res.t, res.a[:, node], label="a (degree)", linewidth=1)
    ax.set_xlabel("Time")
    ax.set_ylabel("State value")

    ax.legend(frameon=False)
    fig.text(
        0.01, 0.98, "A",
        transform=fig.transFigure,
        fontsize=20,
        fontweight="bold",
        va="top",
        ha="left",
        color="black"
    )
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def save_global_means(res: SimulationResult, outpath: str | Path) -> None:
    outpath = Path(outpath)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(res.t, res.X.mean(axis=1), label="Mean X", linewidth=2)
    ax.plot(res.t, res.Y.mean(axis=1), label="Mean Y", linewidth=2)
    ax.plot(res.t, res.Z.mean(axis=1), label="Mean Z", linewidth=2)
    ax.plot(res.t, res.a.mean(axis=1), label="Mean a", linewidth=2)
    ax.set_xlabel("Time")
    ax.set_ylabel("Global mean")

    ax.legend(frameon=False)
    fig.text(
        0.01, 0.98, "B",
        transform=fig.transFigure,
        fontsize=20,
        fontweight="bold",
        va="top",
        ha="left",
        color="black"
    )
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)

def _cumulative_trapezoid(y: np.ndarray, t: np.ndarray) -> np.ndarray:
    return np.concatenate([
        [0.0],
        np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(t))
    ])


def save_cumulative_burden(res: SimulationResult, outpath: str | Path) -> None:
    outpath = Path(outpath)

    mean_Y = res.Y.mean(axis=1)
    mean_Z = res.Z.mean(axis=1)
    mean_a = res.a.mean(axis=1)

    auc_Y = _cumulative_trapezoid(mean_Y, res.t)
    auc_Z = _cumulative_trapezoid(mean_Z, res.t)
    auc_a = _cumulative_trapezoid(mean_a, res.t)

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(res.t, auc_Y, label="AUC(Y)", linewidth=2)
    ax.plot(res.t, auc_Z, label="AUC(Z)", linewidth=2)
    ax.plot(res.t, auc_a, label="AUC(a)", linewidth=2)

    ax.set_xlabel("Time")
    ax.set_ylabel("Cumulative burden")

    ax.legend(frameon=False)
    fig.text(
        0.01, 0.98, "C",
        transform=fig.transFigure,
        fontsize=20,
        fontweight="bold",
        va="top",
        ha="left",
        color="black"
    )
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)


def _threshold_time(t: np.ndarray, y: np.ndarray, threshold: float) -> float | None:
    above = np.where(y >= threshold)[0]

    if len(above) == 0:
        return None

    idx = above[0]

    if idx == 0:
        return float(t[0])

    t0, t1 = t[idx - 1], t[idx]
    y0, y1 = y[idx - 1], y[idx]

    return float(t0 + (threshold - y0) * (t1 - t0) / (y1 - y0))


def save_threshold_times(res: SimulationResult, outpath: str | Path) -> None:
    outpath = Path(outpath)

    mean_a = res.a.mean(axis=1)
    thresholds = [0.2, 0.5, 0.8]

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.plot(res.t, mean_a, color="red", linewidth=2)

    for threshold in thresholds:
        t_cross = _threshold_time(res.t, mean_a, threshold)

        ax.axhline(
            threshold,
            color="gray",
            linestyle="--",
            linewidth=1,
            alpha=0.7
        )

        if t_cross is not None:
            ax.axvline(
                t_cross,
                color="gray",
                linestyle="--",
                linewidth=1,
                alpha=0.7
            )

            ax.scatter(
                [t_cross],
                [threshold],
                color="red",
                s=30,
                zorder=3
            )

            ax.text(
                t_cross,
                threshold,
                f"{t_cross:.2f}",
                fontsize=9,
                ha="left",
                va="bottom"
            )

    ax.set_xlabel("Time")
    ax.set_ylabel("Mean damage a")

    ax.set_ylim(0, 1.0)
    fig.text(
        0.01, 0.98, "D",
        transform=fig.transFigure,
        fontsize=20,
        fontweight="bold",
        va="top",
        ha="left",
        color="black"
    )
    fig.tight_layout()
    fig.savefig(outpath, bbox_inches="tight")
    plt.close(fig)





def _reshape(values: np.ndarray, rows: int, cols: int) -> np.ndarray:
    return values.reshape(rows, cols)


def save_heatmaps(res: SimulationResult, variable: str, time_points: Sequence[float], outdir: str | Path) -> None:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    mapping = {"X": res.X, "Y": res.Y, "Z": res.Z, "a": res.a}
    if variable not in mapping:
        raise ValueError(f"Unknown variable: {variable}")
    arr = mapping[variable]
    vmin = float(arr.min())
    vmax = float(arr.max())

    for tp in time_points:
        idx = int(np.argmin(np.abs(res.t - tp)))
        grid = _reshape(arr[idx], res.rows, res.cols)
        fig, ax = plt.subplots(figsize=(6, 5))
        im = ax.imshow(grid, origin="lower", vmin=vmin, vmax=vmax)
        ax.set_title(f"{variable} at t={res.t[idx]:.2f}")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        fig.tight_layout()
        fig.savefig(outdir / f"heatmap_{variable}_{idx:04d}.png")
        plt.close(fig)


def save_sensitivity_bar(df: pd.DataFrame, outpath: str | Path) -> None:
    outpath = Path(outpath)
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.bar(df["parameter"], df["sensitivity"])
    ax.set_ylabel("Normalized local sensitivity")
    ax.set_title("Parameter sensitivity ranking")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)


def save_fit_comparison(target_df: pd.DataFrame, sim_df: pd.DataFrame, outpath: str | Path) -> None:
    outpath = Path(outpath)
    fig, ax = plt.subplots(figsize=(8, 5))
    for col in target_df.columns:
        if col == "time":
            continue
        sim_col = col.replace("_target", "")
        ax.plot(target_df["time"], target_df[col], linestyle="--", linewidth=2, label=f"{col}")
        if sim_col in sim_df.columns:
            ax.plot(sim_df["time"], sim_df[sim_col], linewidth=2, label=f"sim_{sim_col}")
    ax.set_xlabel("Time")
    ax.set_ylabel("Value")
    ax.set_title("Target vs fitted simulation")
    ax.legend(frameon=False, ncol=2)
    fig.tight_layout()
    fig.savefig(outpath)
    plt.close(fig)





