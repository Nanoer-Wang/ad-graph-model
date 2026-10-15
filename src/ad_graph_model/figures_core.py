"""Publication figures 1-4: schematic, baseline dynamics, spatial
propagation, bifurcation/bistability.

Every function returns the list of written file paths.  Data are taken from
the experiment outputs under ``data_dir`` (produced by ``experiments.py``).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import json

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Patch
from matplotlib.lines import Line2D

from .figstyle import (C, CMAP_DAMAGE, CMAP_TH, CMAP_Y, CMAP_Z, STRATEGY_COLORS,
                       damage_cmap, panel_label, save_fig, setup_style)
from .model import hill, langmuir_occupancy

setup_style()


# ===========================================================================
# Figure 1 - model schematic
# ===========================================================================
def fig1_schematic(outdir: Path, baseline_csv: Path,
                   passivation_csv: Path) -> list:
    outdir = Path(outdir)
    fig = plt.figure(figsize=(13.2, 6.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[2.15, 1.0, 1.0],
                          height_ratios=[1.0, 1.0], hspace=0.32, wspace=0.22)
    axA = fig.add_subplot(gs[:, 0])
    axB = fig.add_subplot(gs[0, 1:])
    axC = fig.add_subplot(gs[1, 1:])

    # ---------------- panel A: molecular module --------------------------
    axA.set_xlim(0, 10.4)
    axA.set_ylim(0, 6.6)
    axA.axis("off")
    axA.text(0.0, 6.45, "A", fontsize=17, fontweight="bold", va="top")

    def node(x, y, label, color, r=0.46, fs=11):
        axA.add_patch(Circle((x, y), r, facecolor=color, edgecolor="white",
                             lw=1.6, zorder=3, alpha=0.92))
        axA.text(x, y, label, ha="center", va="center", fontsize=fs,
                 color="white", fontweight="bold", zorder=4)

    def rarrow(p0, p1, color=None, lw=1.8, rad=0.0, ls="-", ms=13):
        axA.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>",
                                      mutation_scale=ms, lw=lw,
                                      color=color or C["dark"], linestyle=ls,
                                      connectionstyle=f"arc3,rad={rad}",
                                      zorder=2, shrinkA=2, shrinkB=2))

    def inhibit(p0, p1):
        axA.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>",
                                      mutation_scale=12, lw=1.5,
                                      color="#B02020", linestyle=(0, (4, 2)),
                                      zorder=2, shrinkA=2, shrinkB=2))

    node(1.35, 5.15, "X", C["X"])
    node(1.35, 2.55, "Y", C["Y"])
    axA.add_patch(FancyBboxPatch((4.15, 2.1), 2.0, 0.9,
                                 boxstyle="round,pad=0.08,rounding_size=0.18",
                                 facecolor=C["Z"], edgecolor="white", lw=1.6,
                                 zorder=3, alpha=0.92))
    axA.text(5.15, 2.55, "Z", ha="center", va="center", fontsize=12,
             color="white", fontweight="bold", zorder=4)
    node(8.35, 5.15, "P", C["P"])
    node(8.35, 2.55, "a", C["a"], r=0.42)

    axA.text(1.35, 5.85, "monomer", ha="center", fontsize=9.5, color=C["dark"])
    axA.text(1.35, 1.85, "toxic oligomer", ha="center", fontsize=9.5, color=C["dark"])
    axA.text(5.15, 1.72, "fibril / plaque", ha="center", fontsize=9.5, color=C["dark"])
    axA.text(8.35, 5.85, "passivator", ha="center", fontsize=9.5, color=C["dark"])
    axA.text(8.35, 1.85, "neuronal damage", ha="center", fontsize=9.5, color=C["dark"])

    # primary nucleation
    rarrow((1.35, 4.68), (1.35, 3.05), color=C["X"], lw=2.0)
    axA.text(0.28, 3.9, "primary nucleation\n$k_n X^m$", fontsize=9,
             color=C["X"], ha="left")
    # secondary nucleation (blocked by theta)
    rarrow((4.15, 2.9), (1.85, 2.9), color=C["Y"], lw=2.0, rad=-0.25)
    axA.text(3.0, 3.42, "secondary nucleation\n$k_s X Z^q (1-\\theta)$",
             fontsize=9, color=C["Y"], ha="center")
    inhibit((6.35, 3.6), (4.6, 3.15))
    axA.text(6.5, 3.75, "interface passivation $\\theta$\nblocks both surface\n"
             "reactions", fontsize=8.6, color="#B02020", ha="left")
    # elongation
    rarrow((1.9, 4.85), (4.3, 3.15), color=C["Z"], lw=2.0, rad=0.22)
    axA.text(2.35, 4.42, "elongation\n$k_e X Z (1-\\theta)$", fontsize=9,
             color=C["Z"], ha="center")
    inhibit((6.35, 3.35), (4.75, 3.05))
    # sequestration
    rarrow((1.85, 2.3), (4.1, 2.3), color=C["Z"], lw=1.8, rad=0.18)
    axA.text(2.95, 1.62, "sequestration $k_f Y$", fontsize=9, color=C["Z"],
             ha="center")
    # passivator binding
    rarrow((8.0, 4.75), (6.3, 3.05), color=C["P"], lw=2.0, rad=-0.15)
    axA.text(7.6, 4.05, "binding $k_{on}P(1-\\theta)-k_{off}\\theta$\n"
             "$\\theta^* = P/(K_D+P)$,  $K_D=k_{off}/k_{on}$",
             fontsize=8.6, color=C["P"], ha="left")
    # toxicity -> damage
    rarrow((1.9, 2.55), (7.9, 2.55), color=C["a"], lw=2.2, rad=-0.08)
    axA.text(4.9, 2.98, "Hill toxicity  $\\eta_Y\\, \\mathcal{H}(Y)\\,(1-a)$",
             fontsize=9.5, color=C["a"], ha="center")
    # feedback a -> X
    rarrow((8.35, 2.1), (1.6, 4.75), color=C["a"], lw=1.6, rad=0.32, ls=(0, (5, 3)))
    axA.text(4.4, 0.62, "failed-clearance feedback:  $\\lambda_0(1-a)(1+\\beta a)$",
             fontsize=9, color=C["a"], ha="center")
    # repair
    rarrow((8.35, 2.05), (9.7, 2.05), color=C["gray"], lw=1.4)
    axA.text(9.72, 2.05, "repair\n$\\delta_a a$", fontsize=8.6, color=C["gray"],
             va="center")
    # network spread
    rarrow((8.9, 2.9), (9.9, 2.9), color=C["purple"], lw=1.6)
    axA.text(9.9, 3.25, "network spread\n$\\eta_N \\sum_u W_{vu}"
             "\\mathrm{ReLU}(a_u-a_v)$", fontsize=8.2, color=C["purple"], ha="right")
    axA.text(0.0, 0.18, "Clearance: $c_X X$, $(c_Y+k_f)Y$, $c_Z Z$;   "
             "transport: $D_X, D_Y$ along network edges",
             fontsize=9, color=C["gray"], ha="left")

    # ---------------- panel B: network module ----------------------------
    rng = np.random.default_rng(4)
    n = 14
    ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
    pos = np.stack([np.cos(ang), np.sin(ang)], axis=1)
    axB.set_xlim(-1.55, 1.75)
    axB.set_ylim(-1.72, 1.72)
    axB.axis("off")
    axB.text(-1.5, 1.68, "B", fontsize=17, fontweight="bold", va="top")
    axB.set_title("Brain-network propagation", fontsize=11, pad=2)
    edges = [(i, (i + 1) % n) for i in range(n)] + [(0, 7), (3, 10), (5, 12)]
    for i, j in edges:
        axB.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]],
                 color="#BBBBBB", lw=1.0, zorder=1)
    dmg = np.zeros(n)
    dmg[0] = 0.9
    dmg[[1, 13, 7]] = 0.45
    dmg[[2, 12, 6]] = 0.2
    sizes = 120 + 900 * dmg
    axB.scatter(pos[:, 0], pos[:, 1], s=sizes, c=dmg, cmap=damage_cmap(),
                vmin=0, vmax=1, edgecolor="white", lw=1.2, zorder=3)
    # passivated node shield
    axB.add_patch(Circle(pos[3], 0.34, facecolor="none", edgecolor=C["P"],
                         lw=2.4, zorder=4))
    axB.annotate("passivated\ninterface", xy=pos[3], xytext=(0.78, -1.18),
                 fontsize=8.6, color=C["P"], ha="center",
                 arrowprops=dict(arrowstyle="->", color=C["P"], lw=1.2))
    axB.annotate("damage front", xy=pos[6], xytext=(-1.45, 1.02),
                 fontsize=8.6, color=C["a"], ha="center",
                 arrowprops=dict(arrowstyle="->", color=C["a"], lw=1.2))
    axB.text(0, -1.52, "node colour and size = neuronal damage $a$;  "
             "orange ring = passivated interface",
             fontsize=8.4, ha="center", color=C["gray"])

    # ---------------- panel C: intervention timeline ---------------------
    base = pd.read_csv(baseline_csv)
    pas = pd.read_csv(passivation_csv)
    axC.set_title("Intervention timeline and readouts", fontsize=11, pad=2)
    axC.plot(base["time"], base["mean_a"], color=C["a"], lw=2.4,
             label="untreated")
    axC.plot(pas["time"], pas["mean_a"], color=C["P"], lw=2.4,
             label="interface passivation")
    axC.axvline(8.0, color=C["dark"], lw=1.2, ls="--")
    axC.text(8.4, 0.94, "intervention\n($t=8$)", fontsize=8.6, va="top")
    axC.annotate("", xy=(34, 0.93), xytext=(34, 0.36),
                 arrowprops=dict(arrowstyle="<->", color=C["dark"], lw=1.1))
    axC.text(33.4, 0.62, "therapeutic\nbenefit", fontsize=8.6, ha="right",
             color=C["dark"])
    axC.set_xlabel("Time (model units)")
    axC.set_ylabel("Mean neuronal damage $a$")
    axC.set_xlim(0, 40)
    axC.set_ylim(0, 1.02)
    axC.legend(loc="upper left", fontsize=8.6)
    panel_label(axC, "C")

    paths = save_fig(fig, outdir / "fig1_model_schematic")
    return paths


# ===========================================================================
# Figure 2 - baseline dynamics
# ===========================================================================
def fig2_baseline(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    ts = pd.read_csv(data_dir / "baseline_timeseries.csv")
    seed = pd.read_csv(data_dir / "baseline_seed_node.csv")
    metrics = pd.read_json(data_dir / "baseline_metrics.json", typ="series")

    fig, axes = plt.subplots(2, 2, figsize=(11.4, 8.2))
    axA, axB, axC, axD = axes.ravel()

    # A: seed-region trajectories
    axA.plot(seed["time"], seed["Y"], color=C["Y"], lw=2.2, label="$Y$ oligomer")
    axA.plot(seed["time"], seed["Z"], color=C["Z"], lw=2.2, label="$Z$ fibril")
    axA.plot(seed["time"], seed["X"], color=C["X"], lw=1.8, label="$X$ monomer")
    axA.set_xlabel("Time")
    axA.set_ylabel("Burden (seed region)")
    axA.set_xlim(0, 40)
    axA2 = axA.twinx()
    axA2.plot(seed["time"], seed["a"], color=C["a"], lw=2.4, label="$a$ damage")
    axA2.set_ylabel("Damage $a$", color=C["a"])
    axA2.tick_params(axis="y", colors=C["a"])
    axA2.set_ylim(0, 1.05)
    axA2.spines["right"].set_visible(True)
    axA2.grid(False)
    h1, l1 = axA.get_legend_handles_labels()
    h2, l2 = axA2.get_legend_handles_labels()
    # legends are placed below the axes so that they never overlap the curves
    axA.legend(h1 + h2, l1 + l2, loc="upper center",
               bbox_to_anchor=(0.5, -0.20), ncol=4, fontsize=8.6)
    axA.set_title("Seed-region (entorhinal-like) dynamics", fontsize=10.5)
    panel_label(axA, "A")

    # B: global means with phase shading
    axB.plot(ts["time"], ts["mean_X"], color=C["X"], lw=2.0, label="$X$")
    axB.plot(ts["time"], ts["mean_Y"], color=C["Y"], lw=2.2, label="$Y$")
    axB.plot(ts["time"], ts["mean_Z"], color=C["Z"], lw=2.0, label="$Z$")
    axB.plot(ts["time"], ts["mean_a"], color=C["a"], lw=2.4, label="$a$")
    for x0, x1, lab in [(0, 14, "silent"), (14, 24, "amplification"),
                        (24, 40, "clinical")]:
        axB.axvspan(x0, x1, color="#F2F2F2" if lab == "silent" else
                    ("#FBEEE6" if lab == "amplification" else "#F9E3E1"),
                    zorder=0)
    # phase labels sit in a dedicated strip above the plotting area
    for x0, x1, lab in [(0, 14, "silent"), (14, 24, "amplification"),
                        (24, 40, "clinical")]:
        axB.text((x0 + x1) / 2, 1.16, lab, ha="center", va="bottom",
                 fontsize=8.4, color=C["gray"])
    axB.set_xlabel("Time")
    axB.set_ylabel("Network mean")
    axB.set_xlim(0, 40)
    axB.set_ylim(0, 1.28)
    axB.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axB.legend(loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=4,
               fontsize=8.6)
    axB.set_title("Network-averaged cascade", fontsize=10.5)
    panel_label(axB, "B")

    # C: phase portrait Y vs Z
    axC.plot(seed["Z"], seed["Y"], color=C["dark"], lw=1.6, alpha=0.85)
    for tp in [0, 8, 16, 24, 32, 39]:
        idx = int(np.argmin(np.abs(seed["time"] - tp)))
        axC.scatter(seed["Z"][idx], seed["Y"][idx], s=42, zorder=3,
                    color=plt.cm.viridis(tp / 40), edgecolor="white", lw=0.8)
        axC.annotate(f"t={tp}", (seed["Z"][idx], seed["Y"][idx]),
                     textcoords="offset points", xytext=(6, 4), fontsize=8)
    axC.set_xlabel("Fibril burden $Z$")
    axC.set_ylabel("Toxic oligomer burden $Y$")
    axC.set_title("Phase portrait: nucleation lag, then runaway", fontsize=10.5)
    panel_label(axC, "C")

    # D: failed-clearance feedback hump
    a = np.linspace(0, 1, 400)
    lam = 0.22
    beta = 4.5
    prod = lam * (1 - a) * (1 + beta * a)
    axD.plot(a, prod, color=C["a"], lw=2.4,
             label=r"production $\lambda_0(1-a)(1+\beta a)$")
    for Xs, lab in [(0.79, "healthy"), (1.35, "damaged")]:
        axD.axhline(0.28 * Xs, ls="--", lw=1.4, color=C["X"],
                    label=f"clearance $c_XX$ ({lab})")
    axD.axvline((beta - 1) / (2 * beta), color=C["gray"], ls=":", lw=1.2)
    axD.text((beta - 1) / (2 * beta) + 0.015, 0.06,
             "max production\n$a^*=(\\beta-1)/2\\beta$", fontsize=8.4,
             color=C["gray"])
    axD.set_xlabel("Damage $a$")
    axD.set_ylabel("Monomer production / clearance rate")
    axD.set_xlim(0, 1)
    axD.set_ylim(0, 0.62)
    axD.legend(fontsize=8.4, loc="upper right")
    axD.set_title("Failed-clearance feedback loop", fontsize=10.5)
    panel_label(axD, "D")

    fig.tight_layout()
    return save_fig(fig, outdir / "fig2_baseline_dynamics")


# ===========================================================================
# Figure 3 - spatial propagation on the surrogate connectome
# ===========================================================================
def fig3_spatial(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    region = pd.read_csv(data_dir / "spatial_region_time.csv")
    ts = pd.read_csv(data_dir / "spatial_timeseries.csv")
    onsets = pd.read_csv(data_dir / "spatial_onsets.csv")
    strength = pd.read_csv(data_dir / "spatial_strength_damage.csv")
    labels = [c for c in region.columns if c != "time"]
    # positions rebuilt deterministically (same construction as graph.py)
    from .graph import _desikan_labels, _desikan_positions, build_graph
    pos = _desikan_positions()
    px = pos[:, 0] + 0.5 * pos[:, 2]
    py = pos[:, 1]
    t = region["time"].to_numpy()
    A = region[labels].to_numpy().T          # regions x time
    Y = None

    # actual connectome edges (deterministic rebuild, seed 0): the 120
    # strongest weighted connections, drawn faintly behind the nodes
    g = build_graph({"type": "connectome", "seed": 0})
    W = g.W
    iu = np.triu_indices(len(labels), 1)
    w_flat = W[iu]
    top = np.argsort(w_flat)[::-1][:120]
    edge_list = [(iu[0][k], iu[1][k]) for k in top]

    fig = plt.figure(figsize=(13.4, 9.8))
    gs = fig.add_gridspec(3, 4, height_ratios=[1.0, 1.05, 0.9],
                          hspace=0.52, wspace=0.34)
    snap_axes = [fig.add_subplot(gs[0, i]) for i in range(4)]
    axE = fig.add_subplot(gs[1, :2])
    axF = fig.add_subplot(gs[1, 2])
    axG = fig.add_subplot(gs[1, 3])
    axH = fig.add_subplot(gs[2, :])

    # ---- A-D: connectome snapshots
    snap_times = [4.0, 10.0, 18.0, 28.0]
    dmg_cmap = damage_cmap()
    for ax, tp, lab in zip(snap_axes, snap_times, "ABCD"):
        idx = int(np.argmin(np.abs(t - tp)))
        av = A[:, idx]
        for (i, j) in edge_list:
            ax.plot([px[i], px[j]], [py[i], py[j]],
                    color="#C9C9C9", lw=0.5, alpha=0.55, zorder=1)
        ax.scatter(px, py, s=26 + 300 * av, c=av, cmap=dmg_cmap, vmin=0,
                   vmax=1, edgecolor="white", lw=0.5, zorder=3)
        ax.set_xlim(-1.35, 1.85)
        ax.set_ylim(-1.15, 1.25)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"t = {t[idx]:.0f}", fontsize=10)
        panel_label(ax, lab, dx=-0.02, dy=1.13)
    snap_axes[0].set_ylabel("68-region surrogate connectome", fontsize=9)

    # ---- E: region x time heatmap (sorted by onset)
    order = onsets["region"].tolist()
    mat = region[order].to_numpy().T
    im = axE.imshow(mat, aspect="auto", cmap=dmg_cmap, vmin=0, vmax=1,
                    extent=[t[0], t[-1], len(order) - 0.5, -0.5],
                    interpolation="nearest")
    # 68 region names cannot be legibly placed on a single heat-map axis, so
    # every fourth region is labelled and the rest are ticked.
    axE.set_yticks(np.arange(len(order)))
    axE.set_yticklabels([r if (i % 4 == 0) else "" for i, r in enumerate(order)],
                        fontsize=7.0)
    axE.tick_params(axis="y", length=2.5, pad=1.5)
    for s in ("top", "right"):
        axE.spines[s].set_visible(False)
    axE.set_xlabel("Time")
    axE.set_ylabel("region (ordered by onset; every 4th labelled)", fontsize=8.5)
    axE.set_title("Regional damage propagation (ordered by onset)", fontsize=10)
    cb = fig.colorbar(im, ax=axE, fraction=0.03, pad=0.01)
    cb.set_label("damage $a$", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    panel_label(axE, "E", dx=-0.02, dy=1.12)

    # ---- F: topology robustness
    for name, fname in [("2D sheet", "topology_sheet_timeseries.csv"),
                        ("small-world", "topology_smallworld_timeseries.csv"),
                        ("connectome", "topology_connectome_timeseries.csv")]:
        d = pd.read_csv(data_dir / fname)
        axF.plot(d["time"], d["mean_a"], lw=2.0,
                 color={"2D sheet": C["X"], "small-world": C["Z"],
                        "connectome": C["a"]}[name], label=name)
    axF.set_xlabel("Time")
    axF.set_ylabel("Mean damage $a$")
    axF.set_xlim(0, 40)
    axF.set_ylim(0, 1.02)
    axF.legend(fontsize=8, loc="lower right")
    axF.set_title("Topology robustness", fontsize=10)
    panel_label(axF, "F")

    # ---- G: hub test on the un-normalised connectome
    # The row-normalised adjacency has unit row sums, so node strength read
    # from it is identically 1 for every region.  Hubness is therefore taken
    # from the raw connectome (weighted strength before row normalisation).
    try:
        hub = pd.read_csv(data_dir / "spatial_hub_metrics.csv")
        stats_json = json.loads(
            (data_dir / "network_distance_stats.json").read_text())
    except (FileNotFoundError, ValueError):
        hub = pd.read_csv(data_dir / "spatial_strength_damage.csv")
        hub = hub.rename(columns={"strength": "strength_raw"})
        stats_json = None

    axG.scatter(hub["strength_raw"], hub["a_end"], s=30, color=C["gray"],
                alpha=0.85, edgecolor="white", lw=0.5, zorder=3,
                label="terminal damage")
    axG.set_xscale("log")
    axG.set_xlabel("Node strength (un-normalised, log)")
    axG.set_ylabel("Damage at $t=40$", color=C["gray"])
    axG.tick_params(axis="y", colors=C["gray"])
    axG.set_ylim(0.90, 1.0)

    axG2 = axG.twinx()
    axG2.scatter(hub["strength_raw"], hub["onset"], s=30, color=C["Z"],
                 alpha=0.9, marker="^", edgecolor="white", lw=0.5, zorder=4,
                 label="onset time")
    axG2.set_ylabel("Onset time ($a=0.3$)", color=C["Z"])
    axG2.tick_params(axis="y", colors=C["Z"])
    axG2.grid(False)
    axG2.spines["right"].set_visible(True)

    if stats_json is not None:
        rd = stats_json["strength_raw_vs_a_end"]
        ro = stats_json["strength_raw_vs_onset"]
        sub = (f"endpoint $r$={rd['pearson_r']:+.2f}, $p$={rd['pearson_p']:.2f}\n"
               f"onset $r$={ro['pearson_r']:+.2f}, $p$={ro['pearson_p']:.3f}")
    else:  # pragma: no cover - legacy fallback
        r = float(np.corrcoef(hub["strength_raw"], hub["a_end"])[0, 1])
        sub = f"endpoint $r$={r:+.2f}"
    axG.set_title("Hub test: no hub vulnerability", fontsize=10)
    axG.text(0.03, 0.03, sub, transform=axG.transAxes, fontsize=7.6,
             va="bottom", ha="left", color=C["dark"],
             bbox=dict(fc="white", ec="#CCCCCC", lw=0.6, alpha=0.9, pad=2.5))
    h1, l1 = axG.get_legend_handles_labels()
    h2, l2 = axG2.get_legend_handles_labels()
    axG.legend(h1 + h2, l1 + l2, loc="upper center",
               bbox_to_anchor=(0.5, -0.24), ncol=2, fontsize=7.6)
    panel_label(axG, "G")

    # ---- H: onset ranking
    top = onsets.dropna(subset=["onset"]).head(14).iloc[::-1]
    axH.barh(top["region"], top["onset"], color=C["a"], alpha=0.85)
    axH.set_xlabel("Time to reach $a=0.3$")
    axH.set_title("Earliest-affected regions (Braak-like ordering)", fontsize=10)
    axH.tick_params(axis="y", labelsize=8)
    panel_label(axH, "H")

    return save_fig(fig, outdir / "fig3_spatial_propagation")


# ===========================================================================
# Figure 4 - bifurcation and bistability
# ===========================================================================
def _f_prime_at(d: pd.DataFrame, r: float) -> float:
    """Central-difference slope of f(a) at a root, from the grid data."""
    a = d["a"].to_numpy()
    f = d["f"].to_numpy()
    i = int(np.argmin(np.abs(a - r)))
    lo = max(i - 1, 0)
    hi = min(i + 1, len(a) - 1)
    if hi == lo or not (np.isfinite(f[lo]) and np.isfinite(f[hi])):
        return 1.0
    return float((f[hi] - f[lo]) / (a[hi] - a[lo]))


def fig4_bifurcation(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    bif = pd.read_csv(data_dir / "bifurcation_theta.csv")
    crit = pd.read_csv(data_dir / "critical_curve.csv")
    traj_h = pd.read_csv(data_dir / "bistability_healthy.csv")
    traj_d = pd.read_csv(data_dir / "bistability_diseased.csv")
    pot = {lab: pd.read_csv(data_dir / f"potential_{lab}.csv")
           for lab in ("none", "partial", "strong")}

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    axA, axB, axC, axD = axes.ravel()

    # ---- A: bifurcation in theta
    for stable, ls, lab in [(True, "-", "stable"), (False, "--", "unstable")]:
        sub = bif[bif["stable"] == stable].sort_values("theta")
        axA.plot(sub["theta"], sub["a"], ls, color=C["dark"], lw=2.0,
                 label=f"{lab} steady states")
    axA.scatter(bif.loc[bif["stable"], "theta"],
                bif.loc[bif["stable"], "a"], s=14, color=C["dark"], zorder=3)
    th_c = float(crit["theta_analytic"].iloc[0])
    axA.axvline(th_c, color=C["P"], lw=1.6, ls=":")
    # numeric saddle-node fold at the reference clearance (full system)
    num_ref = crit["theta_numeric"].dropna()
    th_on = float(num_ref.iloc[0]) if len(num_ref) else float("nan")
    if np.isfinite(th_on):
        axA.axvline(th_on, color=C["a"], lw=1.6, ls="--")
        axA.annotate("", xy=(th_c, 0.955), xytext=(th_on, 0.955),
                     arrowprops=dict(arrowstyle="<->", color=C["dark"],
                                     lw=1.2))
        axA.text(0.5 * (th_c + th_on), 0.975,
                 f"feedback penalty {th_on - th_c:.2f}", fontsize=8.2,
                 color=C["dark"], ha="center", va="bottom")
    # thresholds are reported in the legend rather than as rotated text, so
    # that they can never collide with the branches when the data change.
    hA, lA = axA.get_legend_handles_labels()
    hA += [Line2D([], [], color=C["P"], ls=":", lw=1.6),
           Line2D([], [], color=C["a"], ls="--", lw=1.6)]
    lA += [f"subsystem $\\theta_c$ = {th_c:.2f}",
           f"full-system fold = {th_on:.2f}"]
    axA.set_xlabel("Passivated interface fraction $\\theta$")
    axA.set_ylabel("Steady-state damage $a^*$")
    axA.set_xlim(0, 1.0)
    axA.set_ylim(0, 1.08)
    axA.legend(hA, lA, fontsize=8.2, loc="upper right")
    axA.set_title("Steady states vs interface passivation", fontsize=10.5)
    panel_label(axA, "A")

    # ---- B: critical boundary in (cY, theta)
    axB.plot(crit["cY"], crit["theta_analytic"], color=C["X"], lw=2.2,
             label="analytic $\\theta_c(c_Y)$ (subsystem, $a=0$)")
    num = crit.dropna(subset=["theta_numeric"]).sort_values("cY")
    if len(num) >= 2:
        xs = num["cY"].to_numpy()
        ys = num["theta_numeric"].to_numpy()
        axB.plot(xs, ys, color=C["a"], lw=1.6, ls="--", zorder=3)
        axB.scatter(xs, ys, s=46, color=C["a"], zorder=4,
                    label="numeric saddle-node (full system)")
        # full-system arrest region: above the numeric boundary
        grid = np.linspace(crit["cY"].min(), crit["cY"].max(), 200)
        num_i = np.interp(grid, xs, ys)
        ana_i = np.interp(grid, crit["cY"], crit["theta_analytic"])
        axB.fill_between(grid, num_i, 1.0, color="#EAF3EA", zorder=0)
        # feedback gap: subsystem arrested but full system still runaway
        axB.fill_between(grid, ana_i, num_i, color="#FDF0E3", zorder=0)
        hB, lB = axB.get_legend_handles_labels()
        hB += [Patch(facecolor="#EAF3EA", edgecolor="#CCCCCC", lw=0.5),
               Patch(facecolor="#FDF0E3", edgecolor="#CCCCCC", lw=0.5),
               Patch(facecolor="#FBE9E7", edgecolor="#CCCCCC", lw=0.5)]
        lB += ["cascade arrested (full system)",
               "feedback gap: subsystem only",
               "runaway (both)"]
        axB.legend(hB, lB, fontsize=8.0, loc="lower left")
    else:
        axB.scatter(num["cY"], num["theta_numeric"], s=46, color=C["a"],
                    zorder=3, label="numeric saddle-node")
        axB.fill_between(crit["cY"], crit["theta_analytic"], 1.0,
                         color="#EAF3EA", zorder=0)
        hB, lB = axB.get_legend_handles_labels()
        hB += [Patch(facecolor="#EAF3EA", edgecolor="#CCCCCC", lw=0.5)]
        lB += ["cascade arrested (full system)"]
        axB.legend(hB, lB, fontsize=8.0, loc="lower left")
    axB.set_xlabel("Oligomer clearance rate $c_Y$")
    axB.set_ylabel("Critical passivation $\\theta_c$")
    axB.set_xlim(0.25, 0.55)
    axB.set_ylim(0, 1.0)
    axB.set_title("Therapeutic boundary: clearance vs passivation",
                  fontsize=10.5)
    panel_label(axB, "B")

    # ---- C: bistability trajectories
    axC.plot(traj_h["time"], traj_h["a"], color=C["X"], lw=2.4,
             label="healthy initial state")
    axC.plot(traj_d["time"], traj_d["a"], color=C["a"], lw=2.4,
             label="diseased initial state")
    axC.set_xlabel("Time")
    axC.set_ylabel("Damage $a$")
    axC.set_xlim(0, 120)
    axC.set_ylim(0, 1.0)
    axC.legend(fontsize=8.6, loc="center right")
    axC.set_title("Bistability at strong passivation ($\\theta=0.95$)",
                  fontsize=10.5)
    panel_label(axC, "C")

    # ---- D: pseudo-potential (defined only on the slow manifold)
    gap_span = {}
    for lab, col in [("none", C["a"]), ("partial", C["P"]), ("strong", C["X"])]:
        d = pot[lab]
        a_v = d["a"].to_numpy()
        phi_v = d["phi"].to_numpy()
        f_v = d["f"].to_numpy()
        th_v = float(d["theta"].iloc[0])
        axD.plot(a_v, phi_v, color=col, lw=2.2,
                 label=f"$\\theta$ = {th_v:.2f}")
        # shade the runaway gap (no positive QSS equilibrium)
        nan = ~np.isfinite(f_v)
        idx = np.flatnonzero(nan)
        if len(idx):
            runs = []
            start = prev = idx[0]
            for i in idx[1:]:
                if i != prev + 1:
                    runs.append((start, prev))
                    start = i
                prev = i
            runs.append((start, prev))
            for s, e in runs:
                if e - s >= 2:
                    axD.axvspan(a_v[s], a_v[e], color="#BBBBBB", alpha=0.18,
                                lw=0)
                    gap_span[lab] = (a_v[s], a_v[e])
        # stable fixed points of the slow damage ODE (sign changes of f)
        roots = []
        for i in range(len(f_v) - 1):
            if np.isfinite(f_v[i]) and np.isfinite(f_v[i + 1]) \
                    and f_v[i] * f_v[i + 1] < 0:
                roots.append(0.5 * (a_v[i] + a_v[i + 1]))
        stable_roots = [r for r in roots
                        if (qss_slope := _f_prime_at(d, r)) < 0]
        axD.scatter(stable_roots,
                    [np.interp(r, a_v, phi_v) for r in stable_roots],
                    s=42, color=col, zorder=4, edgecolor="white", lw=0.8)
    if gap_span:
        axD.text(0.5, 0.02, "grey bands: runaway gap\n(no slow manifold)",
                 transform=axD.transAxes, fontsize=8, color="#666666",
                 ha="center")
    axD.set_xlabel("Damage coordinate $a$")
    axD.set_ylabel("Pseudo-potential $\\Phi(a)$")
    axD.set_xlim(0, 1)
    axD.legend(fontsize=8.6, loc="upper right", frameon=True, framealpha=0.95,
               facecolor="white", edgecolor="none")
    axD.set_title("Energy landscape of the slow damage variable",
                  fontsize=10.5)
    panel_label(axD, "D")

    fig.tight_layout()
    return save_fig(fig, outdir / "fig4_bifurcation_bistability")
