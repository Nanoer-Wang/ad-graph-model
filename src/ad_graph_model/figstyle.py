"""Publication figure style and shared plotting helpers.

Style targets scientific journals: colour-blind-safe palette (Okabe-Ito),
spines removed, light grid, bold panel labels, vector PDF + 600 dpi PNG.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

# ---------------------------------------------------------------------------
# palette (Okabe-Ito colour-blind safe)
# ---------------------------------------------------------------------------
C = {
    "X": "#0072B2",      # monomer - blue
    "Y": "#D55E00",      # toxic oligomer - vermillion
    "Z": "#009E73",      # fibril - green
    "a": "#C0272D",      # damage - deep red
    "P": "#E69F00",      # passivator - orange
    "th": "#56B4E9",     # passivation occupancy - sky
    "gray": "#7F7F7F",
    "dark": "#222222",
    "light": "#E8E8E8",
    "purple": "#9467BD",
}

STRATEGY_COLORS = {
    "Baseline": "#4D4D4D",
    "Passivation": "#D55E00",
    "Clearance": "#0072B2",
    "Nucleation": "#009E73",
    "Propagation": "#9467BD",
    "Combination": "#C0272D",
}

# perceptually-uniform heatmaps
CMAP_DAMAGE = "YlOrRd"
CMAP_Y = "inferno"
CMAP_Z = "viridis"
CMAP_TH = "cividis"

PANEL_KW = dict(fontsize=15, fontweight="bold")


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.5,
        "axes.titlesize": 11,
        "axes.labelsize": 11,
        "axes.linewidth": 0.9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "grid.linewidth": 0.6,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "figure.dpi": 110,
        "savefig.dpi": 600,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "mathtext.fontset": "dejavusans",
    })


def panel_label(ax, label: str, dx: float = -0.13, dy: float = 1.07) -> None:
    ax.text(dx, dy, label, transform=ax.transAxes, va="top", ha="left",
            **PANEL_KW)


def save_fig(fig, outbase: Path, formats=("png", "pdf")) -> list:
    outbase = Path(outbase)
    outbase.parent.mkdir(parents=True, exist_ok=True)
    paths = []
    for fmt in formats:
        p = outbase.with_suffix(f".{fmt}")
        fig.savefig(p, format=fmt)
        paths.append(p)
    plt.close(fig)
    return paths


def damage_cmap():
    return LinearSegmentedColormap.from_list(
        "damage", ["#F7F7F7", "#FBD8C0", "#E8845A", "#C0272D", "#6E0E14"])


def arrow(ax, xy_from, xy_to, color=None, lw=1.6, style="-|>", ms=11,
          connectionstyle="arc3,rad=0.0", ls="-"):
    ax.annotate("", xy=xy_to, xytext=xy_from,
                arrowprops=dict(arrowstyle=style, lw=lw, color=color or C["dark"],
                                mutation_scale=ms, linestyle=ls,
                                connectionstyle=connectionstyle))


def sig_stars(p: float) -> str:
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."
