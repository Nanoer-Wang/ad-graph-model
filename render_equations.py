"""Render manuscript display equations to high-resolution PNGs.

Matplotlib mathtext is used so that equations appear as crisp images inside
the Word manuscript (no LaTeX toolchain required).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.rcParams["mathtext.fontset"] = "stix"

EQUATIONS = {
    "eq_monomer": r"$\dfrac{dX_v}{dt}=\lambda_0(1-a_v)(1+\beta a_v)-c_X X_v-k_n X_v^{m}-k_e X_v Z_v(1-\theta_v)+D_X\sum_{u}L_{vu}X_u$",
    "eq_oligomer": r"$\dfrac{dY_v}{dt}=k_n X_v^{m}+k_s X_v Z_v^{q}(1-\theta_v)-(c_Y+k_f)Y_v+D_Y\sum_{u}L_{vu}Y_u$",
    "eq_fibril": r"$\dfrac{dZ_v}{dt}=k_f Y_v+k_e X_v Z_v(1-\theta_v)-c_Z Z_v$",
    "eq_passivator": r"$\dfrac{dP_v}{dt}=u_P-c_P P_v$",
    "eq_occupancy": r"$\dfrac{d\theta_v}{dt}=k_{on}P_v(1-\theta_v)-k_{off}\theta_v$",
    "eq_damage": r"$\dfrac{da_v}{dt}=\eta_Y\,\mathcal{H}(Y_v;K_Y,h_Y)(1-a_v)+\eta_N\sum_{u}W_{vu}\,\mathrm{ReLU}(a_u-a_v)-\delta_a a_v$",
    "eq_hill": r"$\mathcal{H}(Y;K_Y,h_Y)=\dfrac{Y^{h_Y}}{K_Y^{h_Y}+Y^{h_Y}}$",
    "eq_langmuir": r"$\theta^{*}=\dfrac{P}{K_D+P},\qquad K_D=\dfrac{k_{off}}{k_{on}}$",
    "eq_efficiency": r"$\eta_{\mathrm{agg}}(t)=1-\langle\theta\rangle(t)$",
    "eq_cascade_number": (
        r"$\mathcal{R}=\dfrac{X_s\left[k_sk_f+k_e(c_Y+k_f)\right]}{(c_Y+k_f)c_Z}$"
        r"=\dfrac{k_sX_sk_f}{(c_Y+k_f)c_Z}+\dfrac{k_eX_s}{c_Z}$"
    ),
    "eq_theta_c_from_R": r"$\theta_c=1-\dfrac{1}{\mathcal{R}}$",
    "eq_critical": r"$(c_Y+k_f)\left[c_Z-k_eX_s(1-\theta)\right]=k_sX_sk_f(1-\theta)$",
    "eq_theta_c": r"$\theta_c=1-\dfrac{(c_Y+k_f)c_Z}{X_s\left[k_sk_f+k_e(c_Y+k_f)\right]}$",
    "eq_graph_laplacian": r"$L_{vu}=W_{vu}\ (u\neq v),\quad L_{vv}=-1,\qquad \sum_{u}W_{vu}=1$",
    "eq_prcc": r"$\mathrm{PRCC}(p,m)=\mathrm{corr}\left[r(p),r(m)\mid r(\mathbf{p}_{-p})\right]$",
    "eq_sse": r"$\mathrm{SSE}(\vartheta)=\sum_{i}\left(\dfrac{y_i^{\mathrm{sim}}(\vartheta)-y_i^{\mathrm{obs}}}{\sigma_i}\right)^2$",
    "eq_bliss": r"$E_{\mathrm{Bliss}}=E_P+E_C-E_PE_C$",
    "eq_hysteresis": r"$\Delta a_{\mathrm{thr}}=\theta_{c,2}-\theta_{c,1}>0$",
}


REVISED_EQUATIONS = {
    'eq_monomer': r'$\frac{dX_v}{dt}=\lambda_0(1-a_v)(1+\beta a_v)-c_XX_v-k_nX_v^m$'+'\n'+r'$\qquad{}-k_eX_vZ_v(1-\theta_v)+D_X\sum_uL_{vu}X_u$',
    'eq_oligomer': r'$\frac{dY_v}{dt}=k_nX_v^m+k_sX_vZ_v^q(1-\theta_v)$'+'\n'+r'$\qquad{}-(c_Y+k_f)Y_v+D_Y\sum_uL_{vu}Y_u$',
    'eq_fibril': r'$\frac{dZ_v}{dt}=k_fY_v+k_eX_vZ_v(1-\theta_v)-c_ZZ_v$',
    'eq_passivator': r'$\frac{dP_v}{dt}=u_P-c_PP_v$',
    'eq_occupancy': r'$\frac{d\theta_v}{dt}=k_{\mathrm{on}}P_v(1-\theta_v)-k_{\mathrm{off}}\theta_v$',
    'eq_damage': r'$\frac{da_v}{dt}=\eta_Y\mathcal{H}(Y_v;K_Y,h_Y)(1-a_v)-\delta_a a_v$'+'\n'+r'$\qquad{}+\eta_N\sum_uW_{vu}[a_u-a_v]_+$',
    'eq_hill': r'$\mathcal{H}(Y;K_Y,h_Y)=\frac{Y^{h_Y}}{K_Y^{h_Y}+Y^{h_Y}},\qquad [z]_+=\max(z,0)$',
    'eq_efficiency': r'$\eta_{\mathrm{agg}}(t)=1-\frac{1}{N}\sum_v\theta_v(t)$',
    'eq_langmuir': r'$\theta^*=\frac{P}{K_D+P},\qquad K_D=\frac{k_{\mathrm{off}}}{k_{\mathrm{on}}}$',
    'eq_linear': r'$\dot{Y}=k_sX_s(1-\theta)Z-(c_Y+k_f)Y$'+'\n'+r'$\dot{Z}=k_fY+[k_eX_s(1-\theta)-c_Z]Z$',
    'eq_cascade_number': r'$\mathcal{R}=\frac{X_s[k_sk_f+k_e(c_Y+k_f)]}{(c_Y+k_f)c_Z}$',
    'eq_theta_c': r'$\theta_c=1-\frac{1}{\mathcal{R}}=1-\frac{(c_Y+k_f)c_Z}{X_s[k_sk_f+k_e(c_Y+k_f)]}$',
    'eq_fast': r'$Z(X)=\frac{k_fk_nX^m}{Cc_Z-bX(Ck_e+k_sk_f)}$'+'\n'+r'$Y(X)=\frac{(c_Z-k_eXb)Z(X)}{k_f}$',
    'eq_potential': r'$f(a;\theta)=\eta_Y\mathcal{H}(Y(a;\theta))(1-a)-\delta_a a$'+'\n'+r'$\Phi(a;\theta)=-\int_0^a f(s;\theta)\,ds$',
    'eq_bliss': r'$E_{\mathrm{Bliss}}=E_P+E_C-E_PE_C$',
}


def render_all(outdir: Path, fontsize: float = 11.5, dpi: int = 900) -> list:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, tex in REVISED_EQUATIONS.items():
        tex = tex.replace(r'\frac', r'\dfrac')
        fig = plt.figure(figsize=(0.1, 0.1))
        t = fig.text(0, 0, tex, fontsize=fontsize, color="#1A1A1A")
        fig.canvas.draw()
        bbox = t.get_window_extent()
        plt.close(fig)
        w = bbox.width / fig.dpi + 0.15
        h = bbox.height / fig.dpi + 0.12
        fig = plt.figure(figsize=(w, h))
        fig.text(0.5, 0.5, tex, fontsize=fontsize, color="#1A1A1A",
                 ha="center", va="center")
        p = outdir / f"{name}.png"
        fig.savefig(p, dpi=dpi, transparent=True, bbox_inches="tight",
                    pad_inches=0.03)
        plt.close(fig)
        paths.append(p)
    return paths


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="results_v2/equations")
    args = parser.parse_args()
    paths = render_all(args.outdir)
    print(f"Rendered {len(paths)} equations to {args.outdir}")


if __name__ == "__main__":
    main()
