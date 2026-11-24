"""Problème inverse : paramètres, sensibilités, conception optimale (Fisher) et ajustement (MAP).

Paramètres θ (vecteur d'écarts au nominal) :
  - inconnues     : ln kx et ln kz de la couche cible (sans a priori) ;
  - nuisances     : ln k des autres couches, ln ρc (cible, substrat), R d'interface, R ruban/couche,
                    ln R_back (contact arrière) ;
  - par ruban     : largeur réelle (en unités de w_unc) et ln du facteur d'échelle (étalonnage TCR, L).
Les paramètres « abs » sont normalisés par leur incertitude a priori (bon conditionnement).

Mesure d'un ruban : T(f) = e^s · response(empilement(θ), largeur(θ), f, p), composantes X et Y.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .analytic import penetration, response
from .config import Config, Stack

NO_PRIOR = 1e3


@dataclass
class Param:
    name: str
    kind: str        # 'log' | 'abs'
    where: tuple
    sigma: float
    target: bool = False


def build_params(cfg: Config) -> list[Param]:
    st, t = cfg.stack, cfg.target
    n = len(st.layers)
    tn = st.layers[t].name
    ps = [Param(f"kx {tn}", "log", ("kx", t), NO_PRIOR, True), Param(f"kz {tn}", "log", ("kz", t), NO_PRIOR, True),
          Param(f"ρc {tn}", "log", ("C", t), cfg.prior_C)]
    for i, l in enumerate(st.layers):
        if i != t:
            ps.append(Param(f"k {l.name}#{i}", "log", ("k", i), l.sigma_k))
    if t != n - 1:
        ps.append(Param(f"ρc {st.layers[-1].name}", "log", ("C", n - 1), cfg.prior_C))
    for j in range(n - 1):
        ps.append(Param(f"R {st.layers[j].name}/{st.layers[j + 1].name}", "abs", ("rint", j), st.sigma_r_int[j]))
    ps.append(Param("R ruban/couche", "abs", ("rh",), cfg.heater.sigma_r_h))
    if st.bottom == "contact":
        ps.append(Param("R contact arrière", "log", ("rback",), st.sigma_ln_rback))
    return ps


def apply(cfg: Config, params: list[Param], x: np.ndarray):
    """Empilement et ruban modifiés par les écarts x (paramètres globaux)."""
    st = cfg.stack.copy()
    from dataclasses import replace
    hz = replace(cfg.heater)
    for p, v in zip(params, x):
        if v == 0.0:
            continue
        kind = p.where[0]
        if kind in ("kx", "kz", "k", "C"):
            L = st.layers[p.where[1]]
            f = np.exp(v)
            if kind in ("kx", "k"):
                L.kx *= f
            if kind in ("kz", "k"):
                L.kz *= f
            if kind == "C":
                L.rho *= f
        elif kind == "rint":
            st.r_int[p.where[1]] += v * p.sigma
        elif kind == "rh":
            hz.r_h += v * p.sigma
        elif kind == "rback":
            st.r_back *= np.exp(v)
    return st, hz


def model(cfg, params, x, w, freqs, p, dw=0.0, ds=0.0) -> np.ndarray:
    st, hz = apply(cfg, params, x)
    return np.exp(ds) * response(st, hz, 0.5 * (w + dw * cfg.heater.w_unc), freqs, p)


def jacobian(cfg, params, w, freqs, p, h=1e-3):
    """T0 (nf,) et dérivées complexes (nf, G+2) : globaux, largeur (par w_unc), échelle (ln)."""
    G = len(params)
    T0 = model(cfg, params, np.zeros(G), w, freqs, p)
    cols = []
    for k in range(G):
        e = np.zeros(G)
        e[k] = h
        cols.append((model(cfg, params, e, w, freqs, p) - model(cfg, params, -e, w, freqs, p)) / (2 * h))
    hw = 0.1
    cols.append((model(cfg, params, np.zeros(G), w, freqs, p, dw=hw)
                 - model(cfg, params, np.zeros(G), w, freqs, p, dw=-hw)) / (2 * hw))
    cols.append(T0.copy())
    return T0, np.column_stack(cols)


def precision(params, cfg, m: int) -> np.ndarray:
    pr = [1.0 / p.sigma ** 2 if p.kind == "log" else 1.0 for p in params]
    return np.array(pr + [1.0, 1.0 / cfg.heater.sigma_scale ** 2] * m)


def posterior(Js: list[np.ndarray], params, cfg):
    G, m = len(params), len(Js)
    N = G + 2 * m
    Jn = np.zeros((N, N))
    for k, J in enumerate(Js):
        idx = list(range(G)) + [G + 2 * k, G + 2 * k + 1]
        Jn[np.ix_(idx, idx)] += J
    P = precision(params, cfg, m)
    return np.linalg.inv(Jn + np.diag(P)), Jn, P


def target_stats(cov, params) -> dict:
    ix, iz = [k for k, p in enumerate(params) if p.target]
    sx, sz = np.sqrt(cov[ix, ix]), np.sqrt(cov[iz, iz])
    return {"sigma_kx": float(sx), "sigma_kz": float(sz), "corr": float(cov[ix, iz] / (sx * sz))}


def budget(Js, names, params, cfg) -> dict:
    """Décomposition exacte de Var(kx), Var(kz) entre bruit et chaque a priori."""
    cov, Jn, P = posterior(Js, params, cfg)
    ix, iz = [k for k, p in enumerate(params) if p.target]
    out = {}
    for lab, i in (("kx", ix), ("kz", iz)):
        parts = {"bruit de mesure": float((cov @ Jn @ cov)[i, i])}
        for j, nm in enumerate(names):
            parts[nm] = parts.get(nm, 0.0) + float(cov[i, j] ** 2 * P[j])
        out[lab] = {"variance": float(cov[i, i]), "parts": parts}
    return out


# ------------------------------------------------------------------- conception
@dataclass
class Candidate:
    w: float
    feasible: bool
    L: float
    freqs: np.ndarray
    P: float
    I: float
    J_peak: float
    theta_dc: float
    V1: float
    limit: str
    T: np.ndarray
    sigma: np.ndarray
    J: np.ndarray
    visibility: float
    S_kx: float
    S_kz: float
    single: dict = field(default_factory=dict)


def window(cfg: Config, w: float):
    d, sub = cfg.design, cfg.stack.layers[-1]
    D = sub.kx / sub.C
    f_hi = min(d.f_max, d.f_lockin / 3.0)
    lam = float(penetration(D, d.f_min))
    L = min(max(d.L_min, 40 * w, 4.7 * lam), d.L_max)
    f_lo = d.f_min if L >= 4.7 * lam else D / (4 * np.pi * (L / 4.7) ** 2)
    return L, np.geomspace(f_lo, f_hi, d.n_freq)


def evaluate(cfg: Config, params, w: float) -> Candidate:
    hz, d = cfg.heater, cfg.design
    L, freqs = window(cfg, w)
    T1, D1 = jacobian(cfg, params, w, freqs, 1.0)
    st = cfg.stack
    th_dc = float(response(st, hz, w / 2, [0.0], 1.0)[0].real) if st.bottom in ("contact", "isothermal") else np.inf
    area = w * hz.h
    Rp = hz.rho_e / area
    lim = {"P visée": d.P_target, "électromigration": (d.j_max * area / np.sqrt(2)) ** 2 * Rp}
    if np.isfinite(th_dc):
        lim["ΔT_DC"] = d.dT_dc_max / th_dc
        lim["linéarité"] = d.lin_max / (hz.tcr * th_dc)
    limit = min(lim, key=lim.get)
    P = lim[limit]
    I = np.sqrt(P / Rp)
    V1 = np.sqrt(2) * I * Rp * L / 2
    T, D = P * T1, P * D1
    sig = np.sqrt((cfg.rel_noise * np.abs(T)) ** 2 + (2 * cfg.v_noise / (hz.tcr * V1)) ** 2)
    W = np.vstack([D.real, D.imag]) / np.concatenate([sig, sig])[:, None]
    k = len(freqs) // 2
    ix, iz = [j for j, p in enumerate(params) if p.target]
    vis = np.nan
    if cfg.target < len(st.layers) - 1:
        ctrl = response(st.without(cfg.target), hz, w / 2, freqs[k:k + 1], 1.0)[0]
        vis = float((T1[k] - ctrl).real / T1[k].real)
    return Candidate(w, w >= d.w_min * (1 - 1e-9), L, freqs, P, I, np.sqrt(2) * I / area, P * th_dc, V1, limit, T, sig,
                     W.T @ W, vis, float(-D1[k, ix].real / T1[k].real), float(-D1[k, iz].real / T1[k].real))


@dataclass
class DesignResult:
    cfg: Config
    params: list
    cands: list
    best_by_n: dict
    budget: dict

    @property
    def best(self):
        return self.best_by_n[max(self.best_by_n)]

    @property
    def best_cands(self):
        return [self.cands[i] for i in self.best["idx"]]


def score(u, objective):
    return {"kx": u["sigma_kx"], "kz": u["sigma_kz"]}.get(objective, max(u["sigma_kx"], u["sigma_kz"]))


def local_names(ws):
    return sum(([f"largeur ruban {w * 1e6:.3g} µm", f"étalonnage ruban {w * 1e6:.3g} µm"] for w in ws), [])


def design(cfg: Config) -> DesignResult:
    params = build_params(cfg)
    d = cfg.design
    L = cfg.stack.layers[cfg.target]
    w_narrow = 0.2 * L.d * np.sqrt(L.kx / L.kz)       # B = 0.1 (Dames §6.1)
    widths = np.unique(np.concatenate([np.geomspace(min(d.w_min, 0.5 * w_narrow), d.w_max, d.n_widths), [d.w_min]]))
    cands = [evaluate(cfg, params, w) for w in widths]
    for c in cands:
        c.single = target_stats(posterior([c.J], params, cfg)[0], params)
    feas = [i for i, c in enumerate(cands) if c.feasible]
    best = {}
    for n in range(1, d.n_heaters + 1):
        for combo in itertools.combinations(feas, n):
            u = target_stats(posterior([cands[i].J for i in combo], params, cfg)[0], params)
            s = score(u, d.objective)
            if n not in best or s < best[n]["score"]:
                best[n] = dict(idx=list(combo), score=s, **u)
    res = DesignResult(cfg, params, cands, best, {})
    bc = res.best_cands
    res.budget = budget([c.J for c in bc], [p.name for p in params] + local_names([c.w for c in bc]), params, cfg)
    return res


def model_bias(cfg: Config, params, x: np.ndarray, data, deltas) -> np.ndarray:
    """Biais linéarisé des paramètres si les données contiennent un écart de modèle δ (complexe, K) :
    biais = Cov · Σ_rubans Dᵀ W δ, avec la covariance a posteriori de l'ajustement."""
    st, hz = apply(cfg, params, x[:len(params)])
    from dataclasses import replace as _r
    c2 = _r(cfg, stack=st, heater=hz)
    G = len(params)
    Js, g = [], np.zeros(G + 2 * len(data))
    for j, (d, dl) in enumerate(zip(data, deltas)):
        _, D = jacobian(c2, params, d.w, d.freqs, d.p)
        W = np.vstack([D.real, D.imag]) / np.concatenate([d.sigma, d.sigma])[:, None]
        r = np.concatenate([dl.real, dl.imag]) / np.concatenate([d.sigma, d.sigma])
        Js.append(W.T @ W)
        idx = list(range(G)) + [G + 2 * j, G + 2 * j + 1]
        g[idx] += W.T @ r
    cov = posterior(Js, params, cfg)[0]
    return cov @ g


# -------------------------------------------------------------------- ajustement
@dataclass
class Dataset:
    name: str
    w: float            # largeur mesurée au MEB (m)
    freqs: np.ndarray
    T: np.ndarray       # ΔT complexe mesuré (K)
    sigma: np.ndarray   # incertitude par composante (K)
    p: float            # W/m


@dataclass
class FitResult:
    names: list
    x: np.ndarray
    sigma: np.ndarray
    cov: np.ndarray
    chi2_red: float
    models: list
    params: list
    restart_agree: bool

    def k_values(self, cfg):
        L = cfg.stack.layers[cfg.target]
        ix, iz = [k for k, p in enumerate(self.params) if p.target]
        return {"kx": float(L.kx * np.exp(self.x[ix])), "kz": float(L.kz * np.exp(self.x[iz])),
                "sigma_kx_rel": float(self.sigma[ix]), "sigma_kz_rel": float(self.sigma[iz]),
                "corr": float(self.cov[ix, iz] / (self.sigma[ix] * self.sigma[iz])),
                "kx_identifiable": bool(self.sigma[ix] < 1.0), "kz_identifiable": bool(self.sigma[iz] < 1.0)}


def fit(cfg: Config, data: list[Dataset], verbose=True, corrections: dict | None = None) -> FitResult:
    """corrections : {nom du fichier: rapport complexe T_EF/T_analytique par fréquence} (sortie de verify)."""
    params = build_params(cfg)
    G, m = len(params), len(data)
    prior_sd = np.array([np.inf if p.target else (p.sigma if p.kind == "log" else 1.0) for p in params]
                        + [1.0, cfg.heater.sigma_scale] * m)
    names = [p.name for p in params] + local_names([d.w for d in data])
    corr = [np.asarray(corrections[d.name]) if corrections and d.name in corrections else 1.0 for d in data]

    def predict(x):
        return [c * model(cfg, params, x[:G], d.w, d.freqs, d.p, dw=x[G + 2 * j], ds=x[G + 2 * j + 1])
                for j, (d, c) in enumerate(zip(data, corr))]

    def residuals(x, free):
        full = base.copy()
        full[free] = x
        r = []
        for T, d in zip(predict(full), data):
            r += [(T.real - d.T.real) / d.sigma, (T.imag - d.T.imag) / d.sigma]
        pri = np.isfinite(prior_sd) & free
        r.append(full[pri] / prior_sd[pri])
        return np.concatenate(r)

    # bornes : inconnues à un facteur 100 du nominal, nuisances à ±5 σ de leur a priori
    span = np.where(np.isfinite(prior_sd), 5.0 * prior_sd, np.log(100.0))
    lo, hi = -span, span

    def solve(x0, free):
        nonlocal base
        base = x0.copy()
        x_in = np.clip(x0[free], lo[free] * 0.999, hi[free] * 0.999)
        sol = least_squares(residuals, x_in, args=(free,), bounds=(lo[free], hi[free]), method="trf",
                            diff_step=1e-3, xtol=1e-8, ftol=1e-8)
        out = x0.copy()
        out[free] = sol.x
        return out, sol

    ix, iz = [k for k, p in enumerate(params) if p.target]
    base = np.zeros(G + 2 * m)
    # étape 1 : niveau global de la cible (kx et kz liés)
    x = np.zeros(G + 2 * m)
    f1 = np.zeros(G + 2 * m, bool); f1[[ix, iz]] = True
    sol1 = least_squares(lambda a: residuals(np.array([a[0], a[0]]), f1), [0.0], bounds=([-4.6], [4.6]),
                         method="trf", diff_step=1e-3)
    x[[ix, iz]] = sol1.x[0]
    # étape 2 : kx et kz libres
    x, _ = solve(x, f1)
    # étape 3 : ajustement global avec a priori (MAP)
    allf = np.ones(G + 2 * m, bool)
    x, sol = solve(x, allf)
    # relance depuis le nominal pour contrôler l'unicité du minimum
    x2, sol2 = solve(np.zeros(G + 2 * m), allf)
    agree = abs(sol2.cost - sol.cost) <= 1e-3 * max(sol.cost, 1.0) and np.allclose(x2[[ix, iz]], x[[ix, iz]], atol=1e-2)
    if sol2.cost < sol.cost:
        x, sol = x2, sol2
    Jm = sol.jac
    cov = np.linalg.pinv(Jm.T @ Jm)
    n_meas = sum(2 * len(d.freqs) for d in data)
    chi2 = float(np.sum(residuals(x, allf)[:n_meas] ** 2) / max(n_meas - 2, 1))
    if verbose:
        print(f"  ajustement : {G + 2 * m} paramètres, {n_meas} mesures, χ² réduit = {chi2:.3g}, "
              f"relance depuis le nominal {'concordante' if agree else 'DIFFÉRENTE (minimum local ?)'}")
    return FitResult(names, x, np.sqrt(np.diag(cov)), cov, chi2, predict(x), params, agree)
