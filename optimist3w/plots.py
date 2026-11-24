"""Figures : conception, budget d'erreur, ajustement, vérification."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

BLUE, ORANGE, GREEN, GREY = "#2563eb", "#ea580c", "#059669", "#6b7280"
CYC = [BLUE, ORANGE, GREEN, "#7c3aed", "#db2777", "#ca8a04", "#0891b2"]
plt.rcParams.update({"figure.dpi": 140, "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False,
                     "axes.spines.right": False, "font.size": 9})


def design_figure(res, path):
    c_ = res.cands
    w = np.array([c.w for c in c_]) * 1e6
    fe = np.array([c.feasible for c in c_])
    vis = np.array([c.visibility for c in c_])
    sx, sz = np.array([c.S_kx for c in c_]), np.array([c.S_kz for c in c_])
    ux = 100 * np.array([c.single["sigma_kx"] for c in c_])
    uz = 100 * np.array([c.single["sigma_kz"] for c in c_])
    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    ax[0].semilogx(w, 100 * vis, "o-", color=GREEN)
    ax[0].set_ylabel("visibilité R_film/θ_ruban (%)")
    ax[1].semilogx(w, sz, "o-", color=ORANGE, label="kz")
    ax[1].semilogx(w, sx, "o-", color=BLUE, label="kx")
    ax[1].set_ylabel("sensibilité −∂ln θ/∂ln k")
    ax[1].legend(frameon=False)
    ax[2].loglog(w[fe], ux[fe], "o-", color=BLUE, label="σ(kx), ruban seul")
    ax[2].loglog(w[fe], uz[fe], "o-", color=ORANGE, label="σ(kz), ruban seul")
    ax[2].axhline(100 * res.best["sigma_kx"], color=BLUE, ls="--", lw=1, label="σ(kx), combinaison")
    ax[2].axhline(100 * res.best["sigma_kz"], color=ORANGE, ls="--", lw=1, label="σ(kz), combinaison")
    ax[2].set_ylabel("incertitude a posteriori (%)")
    ax[2].legend(frameon=False, fontsize=7)
    for a in ax:
        a.set_xlabel("largeur 2b (µm)")
        if (~fe).any():
            a.axvspan(w[~fe].min() * 0.9, res.cfg.design.w_min * 1e6, color=GREY, alpha=0.15, lw=0)
        for c in res.best_cands:
            a.axvline(c.w * 1e6, color=GREEN, ls="--", lw=1)
    t = res.cfg.stack.layers[res.cfg.target]
    fig.suptitle(f"Conception — cible {t.name} ; gris : non réalisable ; tirets verts : rubans retenus", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def budget_figure(bud, path, top=6):
    fig, axes = plt.subplots(2, 1, figsize=(9, 4))
    for a, key in zip(axes, ("kx", "kz")):
        b = bud[key]
        items = sorted(b["parts"].items(), key=lambda t: -t[1])
        main = items[:top] + [("autres", sum(v for _, v in items[top:]))]
        left = 0.0
        for j, (n, v) in enumerate(main):
            fr = 100 * v / b["variance"]
            a.barh(0, fr, left=left, color=CYC[j % len(CYC)], edgecolor="white", label=f"{n} ({fr:.1f} %)")
            left += fr
        a.set_xlim(0, 100)
        a.set_yticks([])
        a.grid(False)
        a.set_title(f"σ({key}) = {100 * np.sqrt(b['variance']):.3g} % : origine de la variance", loc="left", fontsize=9)
        a.legend(fontsize=6.5, frameon=False, bbox_to_anchor=(1.01, 1), loc="upper left")
    axes[1].set_xlabel("part de la variance (%)")
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def fit_figure(data, fr, k, path):
    n = len(data)
    fig, ax = plt.subplots(2, n, figsize=(4.6 * n, 6), squeeze=False, sharex="col",
                           gridspec_kw={"height_ratios": [2.2, 1]})
    for j, (d, m) in enumerate(zip(data, fr.models)):
        a = ax[0, j]
        a.errorbar(d.freqs, d.T.real, d.sigma, fmt="o", ms=3.5, color=BLUE, label="mesure, en phase")
        a.errorbar(d.freqs, -d.T.imag, d.sigma, fmt="s", ms=3.5, color=ORANGE, label="mesure, −quadrature")
        a.semilogx(d.freqs, m.real, color=BLUE, label="modèle ajusté")
        a.semilogx(d.freqs, -m.imag, color=ORANGE)
        a.set_title(f"{d.name} (2b = {d.w * 1e6:.3g} µm)", fontsize=9)
        a.set_ylabel("ΔT à 2ω (K)")
        if j == 0:
            a.legend(fontsize=7, frameon=False)
        r = ax[1, j]
        r.semilogx(d.freqs, (d.T.real - m.real) / d.sigma, "o-", ms=3, color=BLUE)
        r.semilogx(d.freqs, (d.T.imag - m.imag) / d.sigma, "s-", ms=3, color=ORANGE)
        r.axhspan(-2, 2, color=GREEN, alpha=0.1, lw=0)
        r.set_ylabel("résidu / σ")
        r.set_xlabel("fréquence électrique f (Hz)")
    fig.suptitle(f"Ajustement : kx = {k['kx']:.4g} ± {100 * k['sigma_kx_rel']:.3g} %, kz = {k['kz']:.4g} ± "
                 f"{100 * k['sigma_kz_rel']:.3g} % (W/m/K), corrélation {k['corr']:+.2f}, χ² réduit = {fr.chi2_red:.2f}",
                 fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)


def verify_figure(rows, path):
    fig, ax = plt.subplots(figsize=(7, 3.8))
    for j, r in enumerate(rows):
        ax.semilogx(r["f"], r["ecart_sigma"], "o-", color=CYC[j % len(CYC)], label=f"2b = {r['w'] * 1e6:.3g} µm")
    ax.axhspan(-1, 1, color=GREEN, alpha=0.12, lw=0)
    ax.set_xlabel("fréquence électrique f (Hz)")
    ax.set_ylabel("|T_EF − T_analytique| / σ_mesure")
    ax.set_title("Vérification : écart entre modèles rapporté au bruit (vert : < 1σ)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
