"""Modèle analytique multicouche : moteur rapide de la conception et de l'ajustement.

Pour chaque nombre d'onde η (transformée en cosinus selon x), chaque couche obéit à
    kz T'' = (kx η² + i·2ω·ρc) T     =>     B = √((kx/kz) η² + i·2ω ρc / kz),
et l'impédance de surface Z = T/φ se transmet du bas vers le haut (Borca-Tasciuc 2001) :
    Z_haut = (Z_bas + tanh(Bd)/(kz B)) / (1 + kz B Z_bas tanh(Bd)),   Z -> Z + R à chaque interface.
Face arrière : contact Z_bas = R_back, isotherme 0, adiabatique Z = 1/(kz B tanh), semi-infini 1/(kz B).

Ruban de demi-largeur b, flux uniforme p/2b, température moyenne (chauffe et mesure sur 2b) :
    T0 = (p/π) ∫ Z(η) sin²(ηb)/(ηb)² dη.
Puis contact R_h et masse du ruban (Acquaroli, éq. 58) :
    T_h = (T0 + R_h p/2b) / (1 + i·2ω·ρc_h·h·(R_h + T0·2b/p)).
Les fréquences sont électriques (f) ; l'échauffement oscille à 2f (pulsation 4πf). f = 0 donne le régime continu.
"""

from __future__ import annotations

import numpy as np

from .config import Heater, Stack

_trapz = getattr(np, "trapezoid", None) or np.trapz


def _tanh(x):
    out = np.ones_like(x, dtype=complex)
    m = x.real < 20.0
    out[m] = np.tanh(x[m])
    return out


def impedance(stack: Stack, lam: np.ndarray, w_th: np.ndarray) -> np.ndarray:
    lam, w = lam[None, :], np.atleast_1d(w_th)[:, None]
    n = len(stack.layers)
    Z = None
    for i in range(n - 1, -1, -1):
        L = stack.layers[i]
        B = np.sqrt((L.kx / L.kz) * lam ** 2 + 1j * w * L.C / L.kz)
        kB = L.kz * B
        if i == n - 1:
            if stack.bottom == "semi-infinite":
                Z = 1.0 / kB
                continue
            t = _tanh(B * L.d)
            if stack.bottom == "adiabatic":
                Z = 1.0 / (kB * t)
                continue
            Zb = stack.r_back if stack.bottom == "contact" else 0.0
        else:
            Zb = Z + stack.r_int[i]
            t = _tanh(B * L.d)
        Z = (Zb + t / kB) / (1.0 + kB * Zb * t)
    return Z


def response(stack: Stack, heater: Heater, b: float, freqs, p: float = 1.0, with_heater: bool = True) -> np.ndarray:
    """Oscillation complexe de température du ruban (K) pour une puissance p (W/m) ; Im < 0 (retard)."""
    f = np.atleast_1d(np.asarray(freqs, dtype=float))
    w_th = 4.0 * np.pi * f
    d_min = min(l.d for l in stack.layers)
    u1 = np.concatenate([np.logspace(-9, 0, 300, endpoint=False), np.linspace(1.0, 200.0, 1600)])
    u_max = 10.0 ** np.ceil(np.log10(max(2e3, 1e3 * b / d_min)))
    u2 = np.logspace(np.log10(200.0), np.log10(u_max), 500)
    I1 = _trapz(impedance(stack, u1 / b, w_th) * (np.sin(u1) ** 2 / u1 ** 2)[None], u1, axis=1)
    I2 = _trapz(impedance(stack, u2 / b, w_th) * (0.5 / u2 ** 2)[None], u2, axis=1)
    top = stack.layers[0]
    tail = b / (4.0 * np.sqrt(top.kx * top.kz) * u_max ** 2)
    T0 = p / (np.pi * b) * (I1 + I2 + tail)
    if not with_heater:
        return T0
    q = p / (2.0 * b)
    return (T0 + heater.r_h * q) / (1.0 + 1j * w_th * heater.C * heater.h * (heater.r_h + T0 / q))


def penetration(D: float, f) -> np.ndarray:
    """λ = √(D/2ω) = √(D/4πf)."""
    return np.sqrt(D / (4.0 * np.pi * np.asarray(f, dtype=float)))
