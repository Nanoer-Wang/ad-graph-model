"""Publication figures 5-8: dose-response, therapy design, global
sensitivity, parameter identifiability."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit

from .figstyle import (C, STRATEGY_COLORS, damage_cmap, panel_label, save_fig,
                       setup_style, sig_stars)

setup_style()


def _hill4(x, bottom, top, ec50, hill):
    return bottom + (top - bottom) / (1.0 + (ec50 / np.maximum(x, 1e-12)) ** hill)


def _fit_ec50(x, y):
    """Inhibitory 4-parameter logistic fit (data-driven bounds, multi-start).

    y = bottom + (top - bottom) / (1 + (ec50/x)^hill),  hill > 0
    (bottom: response at zero dose, top: saturating response).
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    hi = float(np.nanmax(y)) * 1.3
    best = None
    try:
        for ec50_0 in np.geomspace(x.min(), x.max(), 5):
            for h0 in (0.8, 2.0, 5.0):
                popt, _ = curve_fit(_hill4, x, y,
                                    p0=[float(np.nanmax(y)),
                                        float(np.nanmin(y)),
                                        float(ec50_0), h0],
                                    bounds=([0, 0, x.min() / 100.0, 0.2],
                                            [hi, hi, x.max() * 100.0, 20.0]),
                                    maxfev=40000)
                resid = y - _hill4(x, *popt)
                ss = float(np.sum(resid ** 2))
                if best is None or ss < best[0]:
                    best = (ss, popt)
        return best[1]
    except Exception:
        return None


# ===========================================================================
# Figure 5 - passivation dose-response
# ===========================================================================
def fig5_dose(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    dose = pd.read_csv(data_dir / "dose_response.csv")
    kin = pd.read_csv(data_dir / "dose_kinetics.csv")
    cfg_p = {"KD": 0.15}

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    axA, axB, axC, axD = axes.ravel()

    # ---- A: Langmuir isotherm + binding kinetics inset
    P = np.logspace(-3, 1.2, 400)
    th = P / (cfg_p["KD"] + P)
    axA.plot(P, th, color=C["th"], lw=2.6)
    axA.axhline(0.5, color=C["gray"], ls=":", lw=1.2)
    axA.axvline(cfg_p["KD"], color=C["gray"], ls=":", lw=1.2)
    axA.scatter([cfg_p["KD"]], [0.5], s=55, color=C["dark"], zorder=3)
    axA.annotate("$K_D$ = 0.15", (cfg_p["KD"], 0.5),
                 textcoords="offset points", xytext=(10, -16), fontsize=9)
    axA.set_xscale("log")
    axA.set_xlabel("Sustained passivator level $P$")
    axA.set_ylabel("Interface occupancy $\\theta$")
    axA.set_xlim(1e-3, 10 ** 1.2)
    axA.set_ylim(0, 1.02)
    axA.set_title("Langmuir binding of passivator to fibril interfaces",
                  fontsize=10.5)
    panel_label(axA, "A")
    axins = axA.inset_axes([0.52, 0.12, 0.45, 0.42])
    for u, col in zip([0.003, 0.006, 0.012, 0.025],
                      ["#F4C7A3", "#E89B6A", "#D55E00", "#8C3B05"]):
        d = kin[np.isclose(kin["uP"], u)]
        axins.plot(d["time"], d["mean_th"], lw=1.8, color=col,
                   label=f"$u_P$={u}")
    axins.set_xlim(0, 40)
    axins.set_ylim(0, 1.0)
    axins.set_title("binding kinetics $\\theta(t)$", fontsize=8)
    axins.tick_params(labelsize=7)
    axins.legend(fontsize=6.4, frameon=False, loc="lower right")

    # ---- B: mechanism - oligomer trajectories + efficiency
    for u, col in zip([0.003, 0.006, 0.012, 0.025],
                      ["#F4C7A3", "#E89B6A", "#D55E00", "#8C3B05"]):
        d = kin[np.isclose(kin["uP"], u)]
        axB.plot(d["time"], d["mean_Y"], lw=2.2, color=col,
                 label=f"$u_P$ = {u}")
    axB.set_xlabel("Time")
    axB.set_ylabel("Mean oligomer burden $Y$")
    axB.set_xlim(0, 40)
    axB.legend(fontsize=8.4, title="passivation dose", title_fontsize=8.4)
    axB.set_title("Suppression of the toxic-oligomer wave", fontsize=10.5)
    panel_label(axB, "B")
    axb2 = axB.twinx()
    d0 = kin[np.isclose(kin["uP"], 0.012)]
    axb2.plot(d0["time"], 1 - d0["mean_th"], color=C["th"], lw=1.6, ls="--")
    axb2.set_ylabel("interface activity $1-\\theta$", color=C["th"], fontsize=9)
    axb2.tick_params(axis="y", colors=C["th"], labelsize=8)
    axb2.set_ylim(0, 1.05)
    axb2.grid(False)
    axb2.spines["right"].set_visible(True)

    # ---- C: dose-response of toxicity metrics
    for col_y, lab, mk, key in [("Y_peak", "$Y_{peak}$", "o", "Y"),
                                ("AUC_Y", "AUC($Y$)", "s", "Z")]:
        fit = _fit_ec50(dose["uP"].to_numpy(), dose[col_y].to_numpy())
        if fit is not None:
            xs = np.logspace(np.log10(dose["uP"].min()),
                             np.log10(dose["uP"].max()), 200)
            axC.plot(xs, _hill4(xs, *fit), color=C[key],
                     lw=1.4, alpha=0.6)
            axC.scatter([fit[2]], [_hill4(fit[2], *fit)], marker="v", s=60,
                        color=C[key], edgecolor="white", lw=0.6, zorder=4)
            axC.annotate(f"EC$_{{50}}$ = {fit[2]:.4g}",
                         (fit[2], _hill4(fit[2], *fit)),
                         textcoords="offset points", xytext=(8, 10), fontsize=8.4,
                         color=C[key])
        axC.scatter(dose["uP"], dose[col_y], marker=mk, s=42,
                    color=C[key],
                    label=lab, zorder=3, edgecolor="white", lw=0.5)
    axC.set_xscale("log")
    axC.set_xlabel("Passivation dose $u_P$ (infusion rate)")
    axC.set_ylabel("Toxicity metric")
    axC.legend(fontsize=8.6)
    axC.set_title("Toxicity dose-response (4-parameter logistic fit)",
                  fontsize=10.5)
    panel_label(axC, "C")

    # ---- D: damage dose-response at two horizons (the central caution)
    # A finite observation horizon conflates genuine arrest with critical
    # slowing down: the t = 40 snapshot shows a sharp sigmoidal transition
    # whose EC50/Hill are artefacts; the t = 200 durable outcome is graded.
    doseL = pd.read_csv(data_dir / "dose_response_long.csv")
    th = doseL["theta_ss"].to_numpy()
    a40 = doseL["mean_a_end40"].to_numpy()
    a200 = doseL["mean_a_end"].to_numpy()
    order = np.argsort(th)

    # t = 40 horizon: 4-parameter logistic fit (spurious sharpness)
    fit40 = _fit_ec50(th, a40)
    xs = np.linspace(0.10, 1.0, 300)
    if fit40 is not None:
        axD.plot(xs, _hill4(xs, *fit40), color=C["gray"], lw=1.5, alpha=0.75,
                 zorder=2)
        axD.scatter([fit40[2]], [_hill4(fit40[2], *fit40)], marker="v", s=70,
                    color=C["dark"], zorder=4)
        axD.annotate(f"logistic fit ($t$ = 40):\nEC$_{{50}}$ = {fit40[2]:.2f}, "
                     f"Hill = {fit40[3]:.0f}",
                     xy=(fit40[2], _hill4(fit40[2], *fit40)),
                     xytext=(0.30, 0.42), textcoords=axD.transAxes,
                     fontsize=8.4, color=C["dark"],
                     arrowprops=dict(arrowstyle="->", lw=1.0,
                                     color=C["gray"]))
    axD.plot(th[order], a40[order], "o--", ms=4.5, lw=1.2, color=C["gray"],
             alpha=0.9, label="$t$ = 40 snapshot", zorder=3)

    # t = 200 horizon: graded durable decline (monotone interpolation)
    from scipy.interpolate import PchipInterpolator
    pch = PchipInterpolator(th[order], a200[order])
    axD.plot(xs, pch(xs), color=C["a"], lw=2.4, zorder=3,
             label="$t$ = 200 durable")
    axD.plot(th, a200, "o", ms=4.5, color=C["a"], zorder=4,
             markeredgecolor="white", markeredgewidth=0.5)

    # saddle-node fold of the full system (Fig. 4), read from the continuation
    crit = pd.read_csv(data_dir / "critical_curve.csv")
    num = crit.dropna(subset=["theta_numeric"])
    th_sn = float(num["theta_numeric"].iloc[0]) if len(num) else float("nan")
    if np.isfinite(th_sn):
        axD.axvline(th_sn, color=C["purple"], ls=":", lw=1.6, zorder=2)
        axD.annotate(f"full-system fold\n$\\theta_{{SN}}$ = {th_sn:.2f}",
                     (th_sn, 0.06),
                     textcoords="offset points", xytext=(6, 0), fontsize=8,
                     color=C["purple"])

    # the unravelling of an intermediate dose
    i83 = int(np.argmin(np.abs(th - 0.827)))
    axD.annotate("", xy=(th[i83], a40[i83]), xytext=(th[i83], a200[i83]),
                 arrowprops=dict(arrowstyle="<->", lw=1.4, color=C["dark"]))
    axD.annotate(f"{a40[i83]:.2f} at $t$ = 40\n$\\rightarrow$ "
                 f"{a200[i83]:.2f} at $t$ = 200",
                 (th[i83], 0.5 * (a40[i83] + a200[i83])),
                 textcoords="offset points", xytext=(-112, -4), fontsize=8.2,
                 color=C["dark"])

    # durable midpoint: the theta* at which the monotone durable curve reaches
    # the midpoint of its own range (root of pch(th) - mid, by bisection).
    lo_t, hi_t = float(th[order][0]), float(th[order][-1])
    y_lo, y_hi = float(pch(hi_t)), float(pch(lo_t))
    mid = 0.5 * (y_lo + y_hi)
    f = lambda x: float(pch(x)) - mid          # noqa: E731
    a_, b_ = lo_t, hi_t
    for _ in range(80):
        m_ = 0.5 * (a_ + b_)
        if f(a_) * f(m_) <= 0:
            b_ = m_
        else:
            a_ = m_
    th_mid = 0.5 * (a_ + b_)
    axD.scatter([th_mid], [mid], marker="D", s=42, color=C["dark"], zorder=5)
    axD.annotate(f"durable midpoint\n$\\theta^*$ $\\approx$ {th_mid:.2f}",
                 (th_mid, mid),
                 textcoords="offset points", xytext=(-18, 16), fontsize=8.2)

    axD.set_xlabel("Sustained interface occupancy $\\theta^*$")
    axD.set_ylabel("Terminal mean damage $\\langle a\\rangle_{end}$")
    axD.set_xlim(0.10, 1.0)
    axD.set_ylim(-0.03, 1.03)
    axD.legend(fontsize=8.6, loc="upper right")
    axD.set_title("Damage dose-response: horizon dependence", fontsize=10.5)
    panel_label(axD, "D")

    fig.tight_layout()
    return save_fig(fig, outdir / "fig5_passivation_dose_response")


# ===========================================================================
# Figure 6 - therapeutic window and combination therapy
# ===========================================================================
def fig6_therapy(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    win = pd.read_csv(data_dir / "treatment_window_long.csv")
    grid = pd.read_csv(data_dir / "combination_grid_long.csv", index_col=0)
    clr = pd.read_csv(data_dir / "clearance_sweep_long.csv")
    doseL = pd.read_csv(data_dir / "dose_response_long.csv")
    sm_long = pd.read_json(data_dir / "strategy_metrics_long.json")
    base_dmg = float(sm_long["Baseline"]["a200"])   # durable baseline
    cYs = np.array([float(c) for c in grid.columns])
    uPs = np.array([float(u) for u in grid.index])
    G = grid.to_numpy()
    cP, KD = 0.02, 0.15

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    axA, axB, axC, axD = axes.ravel()

    # ---- A: onset sweep trajectories on the durable horizon (passivation)
    onsets_show = [0.0, 8.0, 16.0, 24.0]
    cmap = plt.cm.YlOrRd(np.linspace(0.25, 0.9, len(onsets_show)))
    for onset, col in zip(onsets_show, cmap):
        d = pd.read_csv(data_dir / f"window_trajectory_long_t{int(onset)}.csv")
        axA.plot(d["time"], d["mean_a"], lw=2.2, color=col,
                 label=f"onset $t$ = {int(onset)}")
    dc = pd.read_csv(data_dir / "window_trajectory_long_comb_t0.csv")
    axA.plot(dc["time"], dc["mean_a"], lw=2.6, color=C["dark"],
             label="combination ($t$ = 0)")
    axA.axvline(40, color=C["gray"], ls=":", lw=1.4)
    axA.annotate("$t$ = 40\nsnapshot", (40, 0.965), fontsize=8,
                 color=C["gray"], ha="center")
    # durable single-agent plateau, read from the onset scan
    plat = float(win.loc[win["strategy"] == "Passivation",
                         "mean_a_end"].mean())
    axA.axhline(plat, color=C["gray"], ls="--", lw=1.0, alpha=0.7)
    axA.annotate(f"single-agent durable plateau $\\approx$ {plat:.2f}",
                 (0.985, plat + 0.045), xycoords=("axes fraction", "data"),
                 fontsize=8, color=C["gray"], ha="right")
    traj0 = pd.read_csv(data_dir / "window_trajectory_long_t0.csv")
    a40_pass = float(np.interp(40.0, traj0["time"], traj0["mean_a"]))
    axA.annotate("", xy=(46, plat), xytext=(46, a40_pass),
                 arrowprops=dict(arrowstyle="<->", lw=1.3, color=C["gray"]))
    axA.annotate("apparent protection\nunravels", (52, 0.5 * (a40_pass + plat)),
                 fontsize=8, color=C["gray"], ha="left")
    axA.set_xlabel("Time")
    axA.set_ylabel("Mean damage $a$")
    axA.set_xlim(0, 200)
    axA.set_ylim(0, 1.05)
    axA.legend(fontsize=8.0, loc="upper center",
               bbox_to_anchor=(0.5, -0.17), ncol=3,
               title="passivation onset", title_fontsize=8.0)
    axA.set_title("Onset delay on the durable horizon ($t$ = 200)",
                  fontsize=10.5)
    panel_label(axA, "A")

    # ---- B: onset curves for three strategies at two horizons
    horizon_handles = [
        Line2D([], [], color=C["dark"], lw=2.2, label="durable ($t$ = 200)"),
        Line2D([], [], color=C["dark"], lw=1.3, ls="--",
               label="$t$ = 40 snapshot"),
    ]
    for sname, col in [("Passivation", STRATEGY_COLORS["Passivation"]),
                       ("Clearance", STRATEGY_COLORS["Clearance"]),
                       ("Combination", STRATEGY_COLORS["Combination"])]:
        sub = win[win["strategy"] == sname].sort_values("onset")
        axB.plot(sub["onset"], sub["mean_a_end"], "-o", ms=3.5, lw=2.2,
                 color=col, label=sname)
        axB.plot(sub["onset"], sub["mean_a_end40"], "--", lw=1.3, color=col,
                 alpha=0.85)
    axB.axhline(base_dmg, color=C["gray"], ls=":", lw=1.2)
    axB.annotate(f"untreated ({base_dmg:.2f})", (1.0, base_dmg + 0.03),
                 fontsize=8, color=C["gray"])
    comb = win[win["strategy"] == "Combination"].sort_values("onset")
    ARREST_THRESHOLD = 0.14          # durable-arrest criterion used in Sec. 4.5
    ok = comb[comb["mean_a_end"] <= ARREST_THRESHOLD]
    if len(ok):
        t_last = float(ok["onset"].max())
        axB.axvspan(ok["onset"].min(), t_last, color=C["Z"],
                    alpha=0.10, lw=0)
        axB.annotate(f"combination durably arrests\n"
                     f"($\\langle a\\rangle_{{end}} \\leq {ARREST_THRESHOLD:.2f}$) "
                     f"up to onset $t$ $\\approx$ {t_last:.0f}",
                     (11, 0.30), fontsize=8.2, color="#00756A", ha="center")
    axB.annotate("single agents flatten:\ncritical slowing only",
                 (23, 0.85), fontsize=8.2, color=C["gray"], ha="center")
    axB.set_xlabel("Intervention onset time")
    axB.set_ylabel("Terminal mean damage $\\langle a\\rangle_{end}$")
    axB.set_ylim(-0.03, 1.05)
    leg1 = axB.legend(fontsize=8.6, loc="lower left")
    axB.add_artist(leg1)
    axB.legend(handles=horizon_handles, fontsize=8, loc="upper right")
    axB.set_title("Therapeutic window at two horizons", fontsize=10.5)
    panel_label(axB, "B")

    # ---- C: durable combination design landscape
    im = axC.pcolormesh(cYs, uPs, G, cmap="YlOrRd", shading="auto",
                        vmin=0, vmax=1)
    cs = axC.contour(cYs, uPs, G, levels=[0.15], colors=["#1A1A1A"],
                     linewidths=1.8)
    axC.clabel(cs, fmt="%.2f", fontsize=8.5)
    axC.set_yscale("log")
    axC.set_xlabel("Oligomer clearance $c_Y$")
    axC.set_ylabel("Passivation dose $u_P$")
    cb = fig.colorbar(im, ax=axC, fraction=0.046, pad=0.03)
    cb.set_label("durable damage ($t$ = 200)", fontsize=8)

    def _uP_to_th(u):
        P = np.asarray(u, dtype=float) / cP
        return P / (KD + P)

    def _th_to_uP(t):
        t = np.clip(np.asarray(t, dtype=float), 0.0, 0.999)
        P = KD * t / (1.0 - t)
        return P * cP

    sec = axC.secondary_yaxis("right", functions=(_uP_to_th, _th_to_uP))
    sec.set_yticks([0.2, 0.4, 0.6, 0.8, 0.95])
    sec.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "0.95"])
    sec.set_ylabel("sustained occupancy $\\theta^*$", fontsize=9)

    axC.scatter([0.50], [0.012], marker="*", s=220, color=C["dark"], zorder=5,
                edgecolor="white", lw=0.7)
    axC.annotate("reference\nintervention", (0.50, 0.012),
                 textcoords="offset points", xytext=(12, -2), fontsize=8)
    axC.annotate("durable-arrest ridge\n(0.15 contour)", (0.40, 0.042),
                 fontsize=8.2, color="white", ha="center",
                 bbox=dict(fc="black", alpha=0.35, pad=2.0, lw=0))
    axC.set_title("Combination design landscape (durable, $t$ = 200)",
                  fontsize=10.5)
    panel_label(axC, "C")

    # ---- D: Bliss synergy analysis on durable outcomes
    def benefit_pass(u):
        i = int(np.argmin(np.abs(doseL["uP"].to_numpy() - u)))
        return 1.0 - float(doseL["mean_a_end"].iloc[i]) / base_dmg

    def benefit_clr(cy):
        i = int(np.argmin(np.abs(clr["cY"].to_numpy() - cy)))
        return 1.0 - float(clr["mean_a_end"].iloc[i]) / base_dmg

    syn = []
    for i, u in enumerate(uPs):
        for j, cy in enumerate(cYs):
            bp = benefit_pass(u)
            bc = benefit_clr(cy)
            bliss = bp + bc - bp * bc
            obs = 1.0 - G[i, j] / base_dmg
            syn.append({"uP": u, "cY": cy, "bliss": bliss, "observed": obs,
                        "excess": obs - bliss})
    syn = pd.DataFrame(syn)
    axD.scatter(syn["bliss"], syn["observed"], s=34, color=C["X"], alpha=0.8,
                edgecolor="white", lw=0.5, zorder=3)
    lim = [-0.05, 1.08]
    axD.plot(lim, lim, ls="--", color=C["gray"], lw=1.2)
    axD.set_xlim(lim)
    axD.set_ylim(lim)
    axD.set_xlabel("Bliss-predicted benefit")
    axD.set_ylabel("Observed combination benefit")
    r_bliss = float(np.corrcoef(syn["bliss"], syn["observed"])[0, 1])
    frac_above = float((syn["excess"] > 1e-12).mean())
    mean_exc = float(syn["excess"].mean())
    # reference operating point
    i_ref = int(np.argmin(np.abs(uPs - 0.012)))
    j_ref = int(np.argmin(np.abs(cYs - 0.50)))
    ref = syn.iloc[i_ref * len(cYs) + j_ref]
    axD.scatter([ref["bliss"]], [ref["observed"]], marker="*", s=200,
                color=C["dark"], zorder=5, edgecolor="white", lw=0.7)
    axD.annotate(f"reference:\nBliss {ref['bliss']:.2f} $\\rightarrow$ "
                 f"observed {ref['observed']:.2f}",
                 (ref["bliss"], ref["observed"]),
                 textcoords="offset points", xytext=(-104, 6), fontsize=8.2,
                 color=C["dark"])
    axD.set_title(f"Bliss test on durable outcomes "
                  f"({int(round(frac_above * 100))}% above line; "
                  f"mean excess +{mean_exc:.2f})", fontsize=10.5)
    panel_label(axD, "D")
    syn.to_csv(outdir / "fig6_synergy_long.csv", index=False)

    fig.tight_layout()
    return save_fig(fig, outdir / "fig6_therapy_design")


# ===========================================================================
# Figure 7 - global sensitivity (LHS + PRCC)
# ===========================================================================
def fig7_sensitivity(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    prcc = pd.read_csv(data_dir / "prcc_table.csv")
    ens = pd.read_csv(data_dir / "lhs_ensemble.csv")
    samples = pd.read_csv(data_dir / "lhs_samples.csv")
    base = pd.read_csv(data_dir / "baseline_timeseries.csv")

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.6))
    axA, axB, axC, axD = axes.ravel()

    metric = "mean_a_end"
    sub = prcc[["parameter", f"prcc_{metric}", f"p_{metric}"]].copy()
    sub["abs"] = sub[f"prcc_{metric}"].abs()
    sub = sub.sort_values("abs", ascending=True).tail(12)
    colors = [C["a"] if v > 0 else C["X"] for v in sub[f"prcc_{metric}"]]
    axA.barh(sub["parameter"], sub[f"prcc_{metric}"], color=colors, alpha=0.9)
    for y, (_, r) in enumerate(sub.iterrows()):
        star = sig_stars(r[f"p_{metric}"])
        if star != "n.s.":
            axA.text(r[f"prcc_{metric}"] + (0.02 if r[f"prcc_{metric}"] > 0 else -0.02),
                     y, star, va="center",
                     ha="left" if r[f"prcc_{metric}"] > 0 else "right",
                     fontsize=9, color=C["dark"])
    axA.axvline(0, color=C["dark"], lw=0.8)
    axA.set_xlabel(f"PRCC ({metric})")
    axA.set_title("Global sensitivity ranking (LHS, n=400)", fontsize=10.5)
    axA.set_xlim(-1.15, 1.15)
    panel_label(axA, "A")

    # ---- B: PRCC heatmap across metrics
    metrics = ["Y_peak", "AUC_Y", "mean_a_end", "t_a_0.5",
               "frac_a_over_0.8_end"]
    heat = prcc.set_index("parameter")[[f"prcc_{m}" for m in metrics]]
    heat.columns = metrics
    heat = heat.loc[heat.abs().max(axis=1).sort_values(ascending=False).index]
    vmax = float(np.nanmax(np.abs(heat.to_numpy())))
    im = axB.imshow(heat.to_numpy(), cmap="coolwarm", vmin=-vmax, vmax=vmax,
                    aspect="auto")
    axB.set_xticks(range(len(metrics)))
    axB.set_xticklabels(metrics, rotation=30, ha="right", fontsize=8.4)
    axB.set_yticks(range(len(heat.index)))
    axB.set_yticklabels(heat.index, fontsize=8.4)
    for i in range(heat.shape[0]):
        for j in range(heat.shape[1]):
            v = heat.to_numpy()[i, j]
            if abs(v) > 0.35:
                axB.text(j, i, sig_stars(float(prcc.set_index("parameter")
                                               .loc[heat.index[i], f"p_{metrics[j]}"])),
                         ha="center", va="center", fontsize=6.6, color="white")
    cb = fig.colorbar(im, ax=axB, fraction=0.046, pad=0.03)
    cb.set_label("PRCC", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    axB.set_title("PRCC across outcome metrics", fontsize=10.5)
    axB.grid(False)
    panel_label(axB, "B")

    # ---- C: ensemble uncertainty bands (twin axes: damage | oligomer)
    axC.fill_between(ens["time"], ens["a_q05"], ens["a_q95"], color=C["a"],
                     alpha=0.20, label="90% ensemble ($a$)")
    axC.plot(ens["time"], ens["a_q50"], color=C["a"], lw=2.0,
             label="median $a$")
    axC.plot(base["time"], base["mean_a"], color=C["dark"], lw=1.4, ls="--",
             label="baseline $a$")
    axC.set_xlabel("Time")
    axC.set_ylabel("Damage $a$", color=C["a"])
    axC.set_xlim(0, 40)
    axC.set_ylim(0, 1.02)
    axC.tick_params(axis="y", colors=C["a"])
    axc2 = axC.twinx()
    axc2.fill_between(ens["time"], ens["Y_q05"], ens["Y_q95"], color=C["Y"],
                      alpha=0.20)
    axc2.plot(ens["time"], ens["Y_q50"], color=C["Y"], lw=2.0,
              label="median $Y$")
    axc2.plot(base["time"], base["mean_Y"], color=C["Y"], lw=1.2, ls=":",
              label="baseline $Y$")
    axc2.set_yscale("log")
    axc2.set_ylim(1e-4, 1e2)
    axc2.set_ylabel("Oligomer burden $Y$ (log)", color=C["Y"])
    axc2.tick_params(axis="y", colors=C["Y"], labelsize=8.5)
    axc2.grid(False)
    axc2.spines["right"].set_visible(True)
    h1, l1 = axC.get_legend_handles_labels()
    h2, l2 = axc2.get_legend_handles_labels()
    axC.legend(h1 + h2, l1 + l2, fontsize=7.8, loc="center left")
    axC.set_title("Parameter-uncertainty ensemble", fontsize=10.5)
    panel_label(axC, "C")

    # ---- D: top driver scatter
    top_param = heat.index[0]
    axD.scatter(samples[top_param], samples["mean_a_end"], s=16, color=C["purple"],
                alpha=0.55, edgecolor="none")
    x = samples[top_param].to_numpy()
    y = samples["mean_a_end"].to_numpy()
    b, a_ = np.polyfit(np.log(x), y, 1)
    xs = np.logspace(np.log10(x.min()), np.log10(x.max()), 100)
    axD.plot(xs, a_ + b * np.log(xs), color=C["dark"], lw=1.8)
    r = float(np.corrcoef(np.log(x), y)[0, 1])
    axD.set_xscale("log")
    axD.set_xlabel(f"{top_param} (log scale)")
    axD.set_ylabel("Final damage $\\langle a\\rangle_{end}$")
    axD.set_title(f"Dominant driver: {top_param} (r = {r:.2f})", fontsize=10.5)
    panel_label(axD, "D")

    fig.tight_layout()
    return save_fig(fig, outdir / "fig7_global_sensitivity")


# ===========================================================================
# Figure 8 - parameter identifiability
# ===========================================================================
def fig8_identifiability(data_dir: Path, outdir: Path) -> list:
    data_dir = Path(data_dir)
    target = pd.read_csv(data_dir / "synthetic_target.csv")
    traj = pd.read_csv(data_dir / "fit_trajectory.csv")
    band = pd.read_csv(data_dir / "fit_prediction_band.csv")
    ci = pd.read_csv(data_dir / "local_gaussian_intervals.csv")
    corr = pd.read_csv(data_dir / "parameter_correlation.csv", index_col=0)
    fit_json = pd.read_json(data_dir / "fit_result.json")

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.4))
    axA, axB, axC, axD = axes.ravel()

    # ---- A: synthetic data + fitted curve + band
    axA.scatter(target["time"], target["mean_Y"], s=34, color=C["Y"],
                label="synthetic data $Y$", zorder=3, edgecolor="white",
                lw=0.5)
    axA.scatter(target["time"], target["mean_a"], s=34, color=C["a"],
                label="synthetic data $a$", zorder=3, edgecolor="white",
                lw=0.5)
    axA.plot(traj["time"], traj["mean_Y"], color=C["Y"], lw=2.0,
             label="fitted $Y$")
    axA.plot(traj["time"], traj["mean_a"], color=C["a"], lw=2.0,
             label="fitted $a$")
    axA.fill_between(band["time"], band["Y_lo"], band["Y_hi"], color=C["Y"],
                     alpha=0.15)
    axA.fill_between(band["time"], band["a_lo"], band["a_hi"], color=C["a"],
                     alpha=0.15)
    axA.set_xlabel("Time")
    axA.set_ylabel("Network mean")
    axA.set_xlim(0, 40)
    axA.set_ylim(0, float(max(target["mean_a"].max(), target["mean_Y"].max())) * 1.15)
    axA.legend(fontsize=8.2, loc="upper left")
    axA.set_title("Synthetic benchmark and fitted trajectory", fontsize=10.5)
    panel_label(axA, "A")

    # ---- B: parameter recovery
    true = {"kn": 0.011, "ks": 2.1, "ke": 0.045, "cY": 0.18, "etaY": 0.42,
            "KY": 0.09}
    rec = {n: float(v) for n, v in fit_json["best"].items()}
    names = list(true.keys())
    tv = [true[n] for n in names]
    rv = [rec[n] for n in names]
    lim = [min(min(tv), min(rv)) * 0.6, max(max(tv), max(rv)) * 1.7]
    axB.plot(lim, lim, ls="--", color=C["gray"], lw=1.2)
    axB.scatter(tv, rv, s=80, color=C["X"], zorder=3, edgecolor="white", lw=0.8)
    for n, x0, y0 in zip(names, tv, rv):
        axB.annotate(n, (x0, y0), textcoords="offset points", xytext=(7, 4),
                     fontsize=8.6)
    axB.set_xscale("log")
    axB.set_yscale("log")
    axB.set_xlim(lim)
    axB.set_ylim(lim)
    axB.set_xlabel("True value")
    axB.set_ylabel("Recovered value")
    axB.set_title("Parameter recovery (multi-start least squares)", fontsize=10.5)
    panel_label(axB, "B")

    # ---- C: profile-likelihood CIs
    ci = ci.iloc[::-1]
    y = np.arange(len(ci))
    for i, (_, r) in enumerate(ci.iterrows()):
        axB_ = None
        axC.plot([r["ci_lo"], r["ci_hi"]], [i, i], color=C["X"], lw=2.6,
                 solid_capstyle="round")
        axC.scatter([r["best"]], [i], s=55, color=C["dark"], zorder=3)
        t = true.get(r["parameter"], np.nan)
        if np.isfinite(t):
            axC.scatter([t], [i], marker="|", s=110, color=C["a"], zorder=4,
                        linewidths=2.0)
    axC.set_yticks(y)
    axC.set_yticklabels(ci["parameter"], fontsize=9)
    axC.set_xscale("log")
    axC.set_xlabel("Parameter value (log scale)")
    axC.set_title("Local Gaussian uncertainty (descriptive)", fontsize=10.5)
    handles = [Line2D([], [], color=C["X"], lw=2.6, label="local interval"),
               Line2D([], [], marker="o", color="none", markerfacecolor=C["dark"],
                      markersize=6, label="best fit"),
               Line2D([], [], marker="|", color=C["a"], markersize=10,
                      markeredgewidth=2, linestyle="none", label="true value")]
    axC.legend(handles=handles, fontsize=8, loc="lower right")
    panel_label(axC, "C")

    # ---- D: parameter correlation
    im = axD.imshow(corr.to_numpy(), cmap="coolwarm", vmin=-1, vmax=1)
    axD.set_xticks(range(len(corr)))
    axD.set_xticklabels(corr.columns, rotation=35, ha="right", fontsize=8.4)
    axD.set_yticks(range(len(corr)))
    axD.set_yticklabels(corr.index, fontsize=8.4)
    for i in range(len(corr)):
        for j in range(len(corr)):
            v = corr.to_numpy()[i, j]
            if abs(v) > 0.3:
                axD.text(j, i, f"{v:.2f}", ha="center", va="center",
                         fontsize=7.4,
                         color="white" if abs(v) > 0.65 else C["dark"])
    cb = fig.colorbar(im, ax=axD, fraction=0.046, pad=0.03)
    cb.set_label("correlation", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    axD.grid(False)
    axD.set_title("Local parameter correlations",
                  fontsize=10.5)
    panel_label(axD, "D")

    fig.tight_layout()
    return save_fig(fig, outdir / "fig8_identifiability")
