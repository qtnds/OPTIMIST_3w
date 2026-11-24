"""Fichiers de mesure (format laboratoire), conversion tension -> température, simulation, méthode de la pente.

Un fichier CSV par ruban. En-tête de métadonnées (lignes « # clé = valeur ») :
    width_m      largeur 2b mesurée au MEB
    L_taps_m     distance entre prises de tension
    R_taps_ohm   résistance entre prises (quatre fils, à T0)
    alpha_1_K    TCR étalonné
puis les colonnes : f_Hz, X3_V, Y3_V, sX3_V, sY3_V, V1_V, I_A
(X3, Y3 : composantes à 3ω ; V1 : amplitude à 1ω entre prises ; I : courant efficace ; s : écarts-types).

Conversion (Phase F du protocole) :
    ΔT = 2 (X3 + i Y3) / (α V1),   p = I² R_prises / L_prises (W/m),   σ_T = 2 σ_V3 / (α V1).
"""

from __future__ import annotations

import glob
import json
from pathlib import Path

import numpy as np

from .config import Config
from .identify import Dataset, apply, build_params, model

COLS = ["f_Hz", "X3_V", "Y3_V", "sX3_V", "sY3_V", "V1_V", "I_A"]


def read(path: str) -> Dataset:
    meta = {}
    rows = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "=" in line:
                    k, v = line[1:].split("=", 1)
                    meta[k.strip()] = float(v)
                continue
            if line.startswith("f_Hz"):
                continue
            rows.append([float(v) for v in line.split(",")])
    a = np.array(rows)
    f, X, Y, sX, sY, V1, I = a.T
    alpha = meta["alpha_1_K"]
    T = 2 * (X + 1j * Y) / (alpha * V1)
    sig = 2 * np.sqrt(0.5 * (sX ** 2 + sY ** 2)) / (alpha * V1)
    p = float(np.mean(I ** 2) * meta["R_taps_ohm"] / meta["L_taps_m"])
    return Dataset(Path(path).stem, meta["width_m"], f, T, sig, p)


def read_dir(pattern: str) -> list[Dataset]:
    files = sorted(glob.glob(pattern if any(c in pattern for c in "*?") else str(Path(pattern) / "*.csv")))
    if not files:
        raise FileNotFoundError(f"aucun fichier de mesure pour {pattern}")
    return [read(f) for f in files]


def simulate(cfg: Config, cands, out: Path, offsets=(-0.2, 0.15), seed=1, use_fem=False, max_nodes=None) -> dict:
    """Mesures synthétiques au format laboratoire sur les rubans retenus, avec une vérité connue :
    kx, kz décalés de `offsets`, nuisances, largeurs réelles et TCR tirés selon leurs a priori.
    use_fem : la « réalité » est le modèle EF complet (ruban épais, extrapolé en maillage) au lieu du modèle
    analytique ; c'est le test réaliste de la chaîne, l'ajustement utilisant le modèle analytique."""
    rng = np.random.default_rng(seed)
    params = build_params(cfg)
    G = len(params)
    x = np.zeros(G)
    for k, p in enumerate(params):
        if p.target:
            x[k] = np.log(1 + (offsets[0] if p.where[0] == "kx" else offsets[1]))
        else:
            x[k] = rng.normal(0, p.sigma if p.kind == "log" else 1.0)
    out.mkdir(parents=True, exist_ok=True)
    truth = {"x": x.tolist(), "names": [p.name for p in params], "heaters": []}
    hz = cfg.heater
    for c in cands:
        dw, ds = rng.normal(0, 1.0), rng.normal(0, hz.sigma_scale)
        Rt = hz.rho_e / (c.w * hz.h) * c.L / 2                      # résistance entre prises (largeur nominale)
        I = c.I
        V1 = np.full(len(c.freqs), np.sqrt(2) * I * Rt)
        p_lin = I ** 2 * Rt / (c.L / 2)
        if use_fem:
            from .fem import MAX_NODES, richardson
            st, hz_t = apply(cfg, params, x)
            wt = c.w + dw * hz.w_unc
            T = np.exp(ds) * richardson(st, hz_t, wt, c.freqs, p_lin, 3.0, max_nodes or MAX_NODES)[2]
        else:
            T = model(cfg, params, x, c.w, c.freqs, p_lin, dw=dw, ds=ds)
        V3 = 0.5 * V1 * hz.tcr * T
        sV = np.sqrt((cfg.rel_noise * np.abs(V3)) ** 2 + cfg.v_noise ** 2)
        X = V3.real + rng.normal(0, sV)
        Y = V3.imag + rng.normal(0, sV)
        tag = f"{c.w * 1e6:.4g}".replace(".", "p")
        name = out / f"ruban_{tag}um.csv"
        with open(name, "w", encoding="utf-8") as fh:
            fh.write(f"# width_m = {c.w:.6e}\n# L_taps_m = {c.L / 2:.6e}\n# R_taps_ohm = {Rt:.6e}\n"
                     f"# alpha_1_K = {hz.tcr:.6e}\n" + ",".join(COLS) + "\n")
            for row in zip(c.freqs, X, Y, sV, sV, V1, np.full(len(c.freqs), I)):
                fh.write(",".join(f"{v:.8e}" for v in row) + "\n")
        truth["heaters"].append({"file": name.name, "dw_sigma": float(dw), "ln_scale": float(ds)})
    L = cfg.stack.layers[cfg.target]
    st, _ = apply(cfg, params, x)
    truth["kx_vrai"] = st.layers[cfg.target].kx
    truth["kz_vrai"] = st.layers[cfg.target].kz
    truth["kx_nominal"], truth["kz_nominal"] = L.kx, L.kz
    (out / "verite.json").write_text(json.dumps(truth, indent=2, ensure_ascii=False), encoding="utf-8")
    return truth


def slope_method(ds: Dataset, f_max: float | None = None) -> dict:
    """Méthode de la pente (Cahill) : X = −p/(2πk)·ln(2ω) + cte en régime linéaire (basses fréquences)."""
    m = np.ones(len(ds.freqs), bool) if f_max is None else ds.freqs <= f_max
    x = np.log(2 * 2 * np.pi * ds.freqs[m])
    y = ds.T.real[m]
    A = np.column_stack([x, np.ones_like(x)])
    coef, res, *_ = np.linalg.lstsq(A, y, rcond=None)
    cov = np.linalg.inv(A.T @ A) * (np.mean(ds.sigma[m] ** 2))
    k = -ds.p / (2 * np.pi * coef[0])
    return {"k": float(k), "sigma_k_rel": float(np.sqrt(cov[0, 0]) / abs(coef[0])),
            "k_quadrature": float(-ds.p / (4 * np.mean(ds.T.imag[m]))), "n_points": int(m.sum())}