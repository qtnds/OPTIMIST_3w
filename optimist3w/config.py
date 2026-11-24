"""Description de l'expérience : empilement, interfaces, ruban, contraintes de conception, bruit.

Tout est lu dans un fichier JSON (voir examples/hemt.json). Les couches sont données du haut vers le
bas ; le ruban repose sur la couche 0 ; la dernière couche est le substrat. Chaque grandeur mal connue
porte son incertitude a priori, qui sert à la conception (budget d'erreur) et à l'ajustement.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace

from .units import parse_length


@dataclass
class Layer:
    name: str
    d: float            # épaisseur (m)
    kx: float           # W/m/K dans le plan
    kz: float           # W/m/K hors plan
    rho: float          # kg/m3
    cp: float           # J/kg/K
    sigma_k: float      # incertitude relative a priori sur k (si la couche n'est pas la cible)

    @property
    def C(self) -> float:
        return self.rho * self.cp


@dataclass
class Stack:
    layers: list[Layer]
    r_int: list[float]          # n-1 résistances d'interface (m²K/W), entre couches i et i+1
    sigma_r_int: list[float]    # incertitudes absolues
    bottom: str                 # contact | isothermal | adiabatic | semi-infinite
    r_back: float               # m²K/W (contact)
    sigma_ln_rback: float

    def copy(self) -> "Stack":
        return replace(self, layers=[replace(l) for l in self.layers], r_int=list(self.r_int))

    def without(self, t: int) -> "Stack":
        """Témoin de la méthode différentielle : même empilement sans la couche t (et son interface du dessous)."""
        s = self.copy()
        del s.layers[t]
        if s.r_int:
            del s.r_int[min(t, len(s.r_int) - 1)]
        return s


@dataclass
class Heater:
    metal: str
    h: float            # épaisseur du métal (m)
    k: float            # W/m/K (vérification EF)
    C: float            # J/m3/K
    rho_e: float        # Ohm.m
    tcr: float          # 1/K
    r_h: float          # résistance ruban/couche 0 (m²K/W)
    sigma_r_h: float
    sigma_scale: float  # incertitude relative d'étalonnage (TCR, longueur)
    w_unc: float        # incertitude sur la largeur réalisée (m)


@dataclass
class Design:
    w_min: float
    w_max: float
    n_widths: int
    n_heaters: int
    f_min: float
    f_max: float
    f_lockin: float
    n_freq: int
    P_target: float     # W/m
    j_max: float        # A/m² crête
    lin_max: float      # α·θ_DC max
    dT_dc_max: float    # K
    L_min: float        # longueur totale du ruban (m)
    L_max: float
    objective: str = "max"


@dataclass
class Config:
    stack: Stack
    heater: Heater
    target: int
    design: Design
    v_noise: float
    rel_noise: float
    prior_C: float
    meta: dict = field(default_factory=dict)


METALS = {
    "Pt": dict(k=50.0, rho=21450.0, cp=133.0, rho_e=2.0e-7, tcr=2.5e-3),
    "Au": dict(k=200.0, rho=19300.0, cp=129.0, rho_e=3.5e-8, tcr=3.0e-3),
    "Al": dict(k=180.0, rho=2700.0, cp=897.0, rho_e=3.5e-8, tcr=3.6e-3),
}


def load_config(path: str, target=None, w_min=None, w_unc=None) -> Config:
    with open(path, encoding="utf-8") as fh:
        c = json.load(fh)
    layers = [Layer(l["name"], parse_length(l["d"]), float(l.get("kx", l.get("k"))), float(l.get("kz", l.get("k"))),
                    float(l["rho"]), float(l["cp"]), float(l.get("sigma_k", 0.1))) for l in c["layers"]]
    n = len(layers)
    itf = c.get("interfaces", [{}] * (n - 1))
    if len(itf) != n - 1:
        raise ValueError(f"'interfaces' doit contenir {n - 1} entrées")
    b = c.get("bottom", {})
    stack = Stack(layers, [float(i.get("R", 0.0)) for i in itf], [float(i.get("sigma", 1e-8)) for i in itf],
                  b.get("type", "contact"), float(b.get("R", 1e-5)), float(b.get("sigma_ln", 1.0)))
    h = c["heater"]
    m = dict(METALS[h.get("metal", "Pt")])
    m.update({k: h[k] for k in ("k", "rho", "cp", "rho_e", "tcr") if k in h})
    heater = Heater(h.get("metal", "Pt"), parse_length(h.get("h", "100nm")), m["k"], m["rho"] * m["cp"], m["rho_e"],
                    m["tcr"], float(h.get("R_contact", 2e-8)), float(h.get("sigma_R_contact", 1e-8)),
                    float(h.get("sigma_scale", 0.02)), parse_length(w_unc or h.get("w_unc", "0.1um"), "um"))
    d = c.get("design", {})
    design = Design(parse_length(w_min or d.get("w_min", "1um"), "um"), parse_length(d.get("w_max", "100um"), "um"),
                    int(d.get("n_widths", 20)), int(d.get("n_heaters", 3)), float(d.get("f_min", 20.0)),
                    float(d.get("f_max", 3e4)), float(d.get("f_lockin", 1.02e5)), int(d.get("n_freq", 12)),
                    float(d.get("P_target", 10.0)), float(d.get("j_max", 1e10)), float(d.get("lin_max", 0.01)),
                    float(d.get("dT_dc_max", 10.0)), parse_length(d.get("L_min", "0.5mm"), "mm"),
                    parse_length(d.get("L_max", "5mm"), "mm"), d.get("objective", "max"))
    t = int(target if target is not None else c.get("target", 0)) % n
    nz = c.get("noise", {})
    return Config(stack, heater, t, design, float(nz.get("v_noise", 1e-8)), float(nz.get("rel_noise", 1e-3)),
                  float(c.get("prior_C", 0.05)), {"path": path})
