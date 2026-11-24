"""Tests : python -m pytest -q (environ 1 min)."""

import json
from dataclasses import replace

import numpy as np
import pytest

from optimist3w import data as dt
from optimist3w import identify as idf
from optimist3w.analytic import response
from optimist3w.config import Heater, Layer, Stack, load_config
from optimist3w.fem import heater_response

HZ = Heater("Pt", 100e-9, 50.0, 2.85e6, 2e-7, 2.5e-3, 0.0, 1e-9, 0.02, 0.1e-6)


def _stack(layers, r=None, bottom="semi-infinite"):
    n = len(layers)
    return Stack(layers, r or [0.0] * (n - 1), [1e-9] * (n - 1), bottom, 1e-5, 1.0)


def test_cahill_linear_regime():
    """Substrat seul : pente −p/(2πk) en ln(2ω) et quadrature −p/(4k) (Cahill 1990)."""
    k = 1.38
    st = _stack([Layer("q", 1e-2, k, k, 1.6e6, 1.0, 0.0)])
    T = response(st, HZ, 10e-6, [0.5, 2.0], with_heater=False)
    slope = (T[1].real - T[0].real) / np.log(4.0)
    assert abs(slope / (-1 / (2 * np.pi * k)) - 1) < 0.01
    assert abs(T[0].imag / (-1 / (4 * k)) - 1) < 0.02


def test_film_and_interface_are_series_resistances():
    """Ruban large : un film isolant ajoute d/(2bk), une interface R/(2b) (éq. 47 d'Acquaroli)."""
    sub = Layer("s", 1e-2, 1000.0, 1000.0, 2e6, 1.0, 0.0)
    film = Layer("f", 100e-9, 1.0, 1.0, 2e6, 1.0, 0.0)
    b = 50e-6
    bare = response(_stack([sub]), HZ, b, [100.0], with_heater=False)
    with_film = response(_stack([film, sub]), HZ, b, [100.0], with_heater=False)
    with_r = response(_stack([film, sub], [1e-8]), HZ, b, [100.0], with_heater=False)
    assert abs((with_film - bare)[0].real / (100e-9 / (2 * b)) - 1) < 0.01
    assert abs((with_r - with_film)[0].real / (1e-8 / (2 * b)) - 1) < 0.02


def test_fem_matches_analytic():
    """EF (ruban réel, résistances par couches fictives) et analytique à mieux que 2 %."""
    cfg = load_config("examples/hemt.json")
    f = np.array([50.0, 5000.0])
    for w in (2e-6, 20e-6):
        a = response(cfg.stack, cfg.heater, w / 2, f)
        n = heater_response(cfg.stack, cfg.heater, w, f, 1.0, 3.0)
        assert np.all(np.abs(n - a) / np.abs(a) < 0.02)


def test_budget_is_exact():
    cfg = load_config("examples/hemt.json")
    params = idf.build_params(cfg)
    c = idf.evaluate(cfg, params, 20e-6)
    b = idf.budget([c.J], [p.name for p in params] + idf.local_names([c.w]), params, cfg)
    for key in ("kx", "kz"):
        assert abs(sum(b[key]["parts"].values()) / b[key]["variance"] - 1) < 1e-6


def test_slope_method_on_reference(tmp_path):
    cfg = load_config("examples/silice_reference.json")
    params = idf.build_params(cfg)
    cand = idf.evaluate(cfg, params, 10e-6)
    truth = dt.simulate(cfg, [cand], tmp_path, (0.0, 0.0), seed=3)
    r = dt.slope_method(dt.read_dir(str(tmp_path))[0], f_max=60.0)
    scale = np.exp(truth["heaters"][0]["ln_scale"])
    assert abs(r["k"] * scale / 1.38 - 1) < 0.02


def test_fit_recovers_truth(tmp_path):
    """Chaîne complète sur un film isolant de 200 nm sur Si : simulation puis ajustement de kx, kz."""
    conf = {"target": 0,
            "layers": [{"name": "film", "d": "200nm", "kx": 1.5, "kz": 1.0, "rho": 2200, "cp": 740},
                       {"name": "Si", "d": "500um", "k": 148, "rho": 2330, "cp": 712, "sigma_k": 0.03}],
            "interfaces": [{"R": 0, "sigma": 1e-10}], "bottom": {"type": "contact", "R": 1e-5, "sigma_ln": 1.0},
            "heater": {"metal": "Pt", "h": "50nm", "R_contact": 0, "sigma_R_contact": 1e-10, "w_unc": "5nm"},
            "design": {"w_min": "40nm", "n_freq": 6, "f_min": 100, "f_max": 10000}}
    p = tmp_path / "cfg.json"
    p.write_text(json.dumps(conf))
    cfg = load_config(str(p))
    params = idf.build_params(cfg)
    cands = [idf.evaluate(cfg, params, w) for w in (60e-9, 20e-6)]
    truth = dt.simulate(cfg, cands, tmp_path / "mes", (-0.2, 0.15), seed=2)
    fr = idf.fit(cfg, dt.read_dir(str(tmp_path / "mes")), verbose=False)
    k = fr.k_values(cfg)
    assert k["kz_identifiable"] and k["kx_identifiable"]
    assert abs(np.log(k["kz"] / truth["kz_vrai"])) < 3 * k["sigma_kz_rel"]
    assert abs(np.log(k["kx"] / truth["kx_vrai"])) < 3 * k["sigma_kx_rel"]
