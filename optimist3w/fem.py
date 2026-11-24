"""Éléments finis 2D (coupe transverse x, z) : vérification des hypothèses du modèle analytique.

Le modèle analytique suppose un ruban à flux uniforme, sans épaisseur, corrigé de sa masse et de son
contact de façon globale (éq. 58). Ici le ruban est maillé avec son épaisseur et sa conductivité réelles,
la source Joule est répartie dans sa section, et sa température est moyennée sur la section, comme la lit
la mesure 3ω. L'écart entre les deux modèles, comparé au bruit de mesure, dit si l'ajustement analytique
est suffisant.

Équation résolue (régime harmonique, θ̃ à 2ω) : (K + i·4πf·M) θ̃ = F, éléments P1, masse condensée,
demi-section (symétrie x = 0), bords adiabatiques, face arrière en contact R_back ou isotherme.
Les résistances d'interface et de contact sont représentées par des couches fictives de 1 nm de
conductivité d/R (équivalence exacte en flux 1D, capacité négligeable).
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .analytic import penetration
from .config import Heater, Layer, Stack

_D_FICT = 1e-9


def _graded(a, b, h0, r):
    if b - a <= 1e-12:
        return np.array([])
    pts, x, h = [], a, h0
    while x + h < b - 1e-3 * h:
        x += h
        pts.append(x)
        h *= r
    pts.append(b)
    return np.array(pts)


def _layers_with_resistances(stack: Stack, heater: Heater) -> list[Layer]:
    out = []
    if heater.r_h > 0:
        k = _D_FICT / heater.r_h
        out.append(Layer("R_h", _D_FICT, k, k, 1.0, 1.0, 0.0))
    for i, L in enumerate(stack.layers):
        out.append(L)
        if i < len(stack.r_int) and stack.r_int[i] > 0:
            k = _D_FICT / stack.r_int[i]
            out.append(Layer("R_int", _D_FICT, k, k, 1.0, 1.0, 0.0))
    return out


MAX_NODES = 80_000        # plafond par défaut : ~0,3 Go de mémoire pour la factorisation LU


def _axes(stack: Stack, heater: Heater, w: float, freqs, refine: float):
    layers = _layers_with_resistances(stack, heater)
    b, h = w / 2, heater.h
    H = sum(l.d for l in layers)
    sub = stack.layers[-1]
    X = max(3 * H, 40 * b, 5 * float(penetration(sub.kx / sub.C, min(np.atleast_1d(freqs)))))
    if stack.bottom == "contact":
        X = max(X, 5 * np.sqrt(sub.kx * H * stack.r_back))
    r = 1.0 + 0.24 / refine
    nb = max(int(10 * refine), 4)
    xs = np.concatenate([np.linspace(0, b, nb + 1), _graded(b, X, b / nb, r)])
    hz = min(b / 5, H) / refine
    zs, z = list(np.linspace(-h, 0, max(int(3 * refine), 2) + 1)), 0.0
    for i, L in enumerate(layers):
        if L.name in ("R_h", "R_int"):                          # couche fictive : flux 1D, 2 mailles suffisent
            zs += list(np.linspace(z, z + L.d, 3)[1:])
        elif L.d <= 8 * hz or i < len(layers) - 1:
            n = int(np.clip(max(np.ceil(L.d / hz), 3 * refine), 2, 80))
            zs += list(np.linspace(z, z + L.d, n + 1)[1:])
        else:
            zs += list(_graded(z, z + L.d, min(hz, L.d / 8), r))
        z += L.d
    return layers, xs, np.array(zs)


def mesh_size(stack: Stack, heater: Heater, w: float, freqs, refine: float) -> int:
    """Nombre de nœuds (majorant) du maillage, calculé sans rien assembler."""
    _, xs, zs = _axes(stack, heater, w, freqs, refine)
    return len(xs) * len(zs)


def safe_refine(stack: Stack, heater: Heater, w: float, freqs, refine: float, max_nodes: int = MAX_NODES) -> float:
    """Plus grande finesse ≤ refine dont le maillage tient sous max_nodes nœuds."""
    r = refine
    while r > 0.5 and mesh_size(stack, heater, w, freqs, r) > max_nodes:
        r *= 0.85
    return r


def richardson(stack: Stack, heater: Heater, w: float, freqs, p: float = 1.0, refine: float = 3.0,
               max_nodes: int = MAX_NODES):
    """Calculs à r et 2r puis extrapolation d'ordre 1 (2·T(2r) − T(r)). La finesse 2r est plafonnée pour
    que le maillage reste sous max_nodes nœuds (protection mémoire) ; r suit (r = 2r/2).
    Renvoie (T(r), T(2r), T_extrapolée, r_effectif)."""
    r2 = safe_refine(stack, heater, w, freqs, 2 * refine, max_nodes)
    r1 = r2 / 2
    if r2 < 2 * refine - 1e-9:
        print(f"  [EF] 2b = {w * 1e6:.3g} µm : finesse réduite {refine:g} -> {r1:.2g} "
              f"(plafond {max_nodes} nœuds, --max_nodes pour changer)", flush=True)
    t1 = heater_response(stack, heater, w, freqs, p, r1, max_nodes=None)
    t2 = heater_response(stack, heater, w, freqs, p, r2, max_nodes=None)
    return t1, t2, 2 * t2 - t1, r1


def heater_response(stack: Stack, heater: Heater, w: float, freqs, p: float = 1.0, refine: float = 2.0,
                    max_nodes: int | None = MAX_NODES) -> np.ndarray:
    if max_nodes:
        refine = safe_refine(stack, heater, w, freqs, refine, max_nodes)
    layers, xs, zs = _axes(stack, heater, w, freqs, refine)
    b = w / 2
    bounds = np.cumsum([0.0] + [l.d for l in layers])
    zc, xc = 0.5 * (zs[:-1] + zs[1:]), 0.5 * (xs[:-1] + xs[1:])
    nx = len(xs)
    lay_z = np.where(zc < 0, -1, np.searchsorted(bounds, zc, side="right") - 1)
    J, I = np.meshgrid(np.arange(len(zs) - 1), np.arange(nx - 1), indexing="ij")
    keep = ~((lay_z[J] == -1) & (xc[I] > b))                  # pas de mailles d'air à côté du ruban
    J, I = J[keep], I[keep]
    a, bb, c, d = J * nx + I, J * nx + I + 1, (J + 1) * nx + I, (J + 1) * nx + I + 1
    tris = np.concatenate([np.column_stack([a, bb, d]), np.column_stack([a, d, c])])
    rg = np.concatenate([lay_z[J], lay_z[J]])
    used = np.unique(tris)
    remap = -np.ones(nx * len(zs), int)
    remap[used] = np.arange(len(used))
    XX, ZZ = np.meshgrid(xs, zs)
    nodes = np.column_stack([XX.ravel(), ZZ.ravel()])[used]
    tris = remap[tris]
    N = len(nodes)

    P = nodes[tris]
    x_, z_ = P[..., 0], P[..., 1]
    det = (x_[:, 1] - x_[:, 0]) * (z_[:, 2] - z_[:, 0]) - (x_[:, 2] - x_[:, 0]) * (z_[:, 1] - z_[:, 0])
    A = 0.5 * np.abs(det)
    gx = np.stack([z_[:, 1] - z_[:, 2], z_[:, 2] - z_[:, 0], z_[:, 0] - z_[:, 1]], 1) / det[:, None]
    gz = np.stack([x_[:, 2] - x_[:, 1], x_[:, 0] - x_[:, 2], x_[:, 1] - x_[:, 0]], 1) / det[:, None]
    del P, x_, z_
    props = np.array([(heater.k, heater.k, heater.C)] + [(l.kx, l.kz, l.C) for l in layers])
    kx, kz, C = props[rg + 1].T
    Ke = A[:, None, None] * (kx[:, None, None] * gx[:, :, None] * gx[:, None, :]
                             + kz[:, None, None] * gz[:, :, None] * gz[:, None, :])
    rows, cols = np.repeat(tris, 3, 1).ravel(), np.tile(tris, (1, 3)).ravel()
    K = sp.csr_matrix((Ke.ravel(), (rows, cols)), shape=(N, N))
    del Ke, rows, cols, gx, gz
    m = np.zeros(N)
    np.add.at(m, tris.ravel(), np.repeat(C * A / 3, 3))
    he = rg == -1
    q = 0.5 * p / A[he].sum()                                  # demi-section : p/2 sur l'aire du demi-ruban
    F = np.zeros(N)
    np.add.at(F, tris[he].ravel(), np.repeat(q * A[he] / 3, 3))
    wH = np.zeros(N)
    np.add.at(wH, tris[he].ravel(), np.repeat(A[he] / 3, 3))
    wH /= A[he].sum()

    bot = remap[(len(zs) - 1) * nx + np.arange(nx)]
    free = np.arange(N)
    if stack.bottom == "contact":
        Lb = np.diff(xs)
        Kb = Lb[:, None, None] / (6 * stack.r_back) * np.array([[2.0, 1.0], [1.0, 2.0]])
        e = np.column_stack([bot[:-1], bot[1:]])
        K = K + sp.csr_matrix((Kb.ravel(), (np.repeat(e, 2, 1).ravel(), np.tile(e, (1, 2)).ravel())), shape=(N, N))
    elif stack.bottom == "isothermal":
        free = np.setdiff1d(free, bot)
    Kf = K.tocsr()[free][:, free].tocsc()
    mf = m[free]
    out = []
    for f in np.atleast_1d(freqs):
        # matrice symétrique : l'ordonnancement MMD(AᵀA+A) réduit le remplissage LU (~ −25 % mémoire, ×2 vitesse)
        lu = spla.splu((Kf + sp.diags(1j * 4 * np.pi * f * mf)).tocsc(), permc_spec="MMD_AT_PLUS_A")
        th = np.zeros(N, complex)
        th[free] = lu.solve(F[free].astype(complex))
        del lu
        out.append(wH @ th)
    return np.array(out)