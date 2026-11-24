"""OPTIMIST_3w — mise en œuvre du protocole de mesure 3ω (conception, ajustement, vérification).

  python main.py design   --config examples/hemt.json                 # phase A2 : géométrie des rubans
  python main.py simulate --config examples/hemt.json --out mesures   # mesures synthétiques (test de la chaîne)
  python main.py fit      --config examples/hemt.json --data mesures  # phase G : identification de kx, kz
  python main.py verify   --config examples/hemt.json --data mesures  # phase G6 : contrôle par éléments finis
  python main.py slope    --data reference.csv --f_max 200            # phase F3 : validation du banc (pente)
"""

import argparse
import json
import os

# Une seule tâche de calcul par bibliothèque : évite de saturer tous les cœurs (et la machine).
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "2")
from pathlib import Path

import numpy as np

from optimist3w import data as dt
from optimist3w import identify as idf
from optimist3w import plots
from optimist3w.analytic import response
from optimist3w.config import load_config
from optimist3w.units import parse_length, parse_list


def _pct(v):
    return f"{100 * v:.3g} %" if v < 10 else "> 1000 % (non identifiable)"


def _cfg(a):
    return load_config(a.config, a.target, a.w_min, a.w_unc)


def cmd_design(a):
    cfg = _cfg(a)
    res = idf.design(cfg)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t = cfg.stack.layers[cfg.target]
    print(f"Cible : {t.name} (couche {cfg.target}) ; lithographie 2b ≥ {cfg.design.w_min * 1e6:.3g} µm "
          f"(± {cfg.heater.w_unc * 1e6:.3g} µm)\n")
    print(f"{'2b (µm)':>9} {'réal.':>5} {'visib.':>7} {'S_kx':>7} {'S_kz':>7} {'σkx seul':>10} {'σkz seul':>10}  limite")
    for c in res.cands:
        print(f"{c.w * 1e6:9.3g} {'oui' if c.feasible else 'non':>5} {100 * c.visibility:6.1f}% {c.S_kx:7.3f} "
              f"{c.S_kz:7.3f} {min(100 * c.single['sigma_kx'], 999):9.3g}% {min(100 * c.single['sigma_kz'], 999):9.3g}%  {c.limit}")
    print()
    for n, b in res.best_by_n.items():
        ws = ", ".join(f"{res.cands[i].w * 1e6:.3g}" for i in b["idx"])
        print(f"{n} ruban(s) : 2b = [{ws}] µm -> σ(kx) = {_pct(b['sigma_kx'])}, σ(kz) = {_pct(b['sigma_kz'])}")
    print("\nGéométrie retenue :")
    geo = []
    for c in res.best_cands:
        g = {"2b_um": c.w * 1e6, "L_totale_mm": c.L * 1e3, "prises_a_+-L/4_mm": c.L / 4 * 1e3,
             "f_Hz": [float(c.freqs[0]), float(c.freqs[-1]), len(c.freqs)], "I_eff_mA": c.I * 1e3,
             "P_W_par_m": c.P, "J_MA_cm2": c.J_peak * 1e-10, "theta_DC_K": c.theta_dc, "V1_V": c.V1,
             "limite": c.limit}
        geo.append(g)
        print("  " + ", ".join(f"{k} {v:.3g}" if isinstance(v, float) else f"{k} {v}" for k, v in g.items()))
    print("\nBudget d'erreur :")
    for key in ("kx", "kz"):
        b = res.budget[key]
        top = sorted(b["parts"].items(), key=lambda x: -x[1])[:4]
        print(f"  σ({key}) = {_pct(np.sqrt(b['variance']))} : "
              + ", ".join(f"{n} {100 * v / b['variance']:.0f} %" for n, v in top))
    json.dump({"rubans": geo, "best_by_n": res.best_by_n, "budget": res.budget},
              open(out / "conception.json", "w"), indent=2, ensure_ascii=False, default=float)
    plots.design_figure(res, out / "conception.png")
    plots.budget_figure(res.budget, out / "budget.png")
    print(f"\n-> {out / 'conception.json'}, conception.png, budget.png")
    return res


def cmd_simulate(a):
    cfg = _cfg(a)
    if a.widths:
        params = idf.build_params(cfg)
        cands = [idf.evaluate(cfg, params, parse_length(w, "um")) for w in parse_list(a.widths)]
    else:
        print("Conception préalable pour choisir les rubans...")
        cands = idf.design(cfg).best_cands
    truth = dt.simulate(cfg, cands, Path(a.out), tuple(float(v) for v in parse_list(a.truth)), a.seed, a.fem, a.max_nodes)
    print(f"{len(cands)} fichiers de mesure écrits dans {a.out} ; vérité : kx = {truth['kx_vrai']:.4g}, "
          f"kz = {truth['kz_vrai']:.4g} W/m/K (voir verite.json)")


def cmd_fit(a):
    cfg = _cfg(a)
    data = dt.read_dir(a.data)
    print(f"{len(data)} rubans lus : " + ", ".join(f"{d.name} (2b = {d.w * 1e6:.3g} µm)" for d in data))
    corr = None
    if a.correction:
        corr = {k: np.array([complex(*v) for v in vals]) for k, vals in json.load(open(a.correction)).items()}
        print("Correction de modèle appliquée :", a.correction)
    fr = idf.fit(cfg, data, corrections=corr)
    k = fr.k_values(cfg)
    for lab in ("kx", "kz"):
        if k[f"{lab}_identifiable"]:
            print(f"{lab} = {k[lab]:.4g} W/m/K ± {100 * k[f'sigma_{lab}_rel']:.3g} %")
        else:
            print(f"{lab} : NON IDENTIFIABLE avec ces rubans (incertitude > 100 %, valeur sans signification)")
    print(f"corrélation kx-kz : {k['corr']:+.2f}")
    truth = Path(a.data) / "verite.json"
    if truth.exists():
        t = json.load(open(truth))
        print(f"vérité (simulation) : kx = {t['kx_vrai']:.4g}, kz = {t['kz_vrai']:.4g} W/m/K")
        for lab in ("kx", "kz"):
            if k[f"{lab}_identifiable"]:
                print(f"  écart {lab} : {100 * (k[lab] / t[lab + '_vrai'] - 1):+.2f} % "
                      f"= {np.log(k[lab] / t[lab + '_vrai']) / k[f'sigma_{lab}_rel']:+.2f} σ")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    res = {**k, "chi2_reduit": fr.chi2_red, "relance_concordante": fr.restart_agree,
           "parametres": [{"nom": n, "ecart": float(x), "sigma": float(s)} for n, x, s in zip(fr.names, fr.x, fr.sigma)]}
    json.dump(res, open(out / "ajustement.json", "w"), indent=2, ensure_ascii=False)
    plots.fit_figure(data, fr, k, out / "ajustement.png")
    print(f"-> {out / 'ajustement.json'}, ajustement.png")


def cmd_verify(a):
    """Phase G6 : le modèle d'ajustement (analytique) est-il assez fidèle ? Comparaison au modèle EF complet
    aux valeurs ajustées, extrapolé en maillage (Richardson, convergence d'ordre 1), puis propagation de
    l'écart au résultat : biais(kx, kz) comparé à leur incertitude."""
    from optimist3w.fem import richardson
    cfg = _cfg(a)
    data = dt.read_dir(a.data)
    params = idf.build_params(cfg)
    G = len(params)
    x = np.zeros(G + 2 * len(data))
    fj = Path(a.out) / "ajustement.json"
    sig = None
    if fj.exists():
        js = json.load(open(fj))
        x = np.array([p["ecart"] for p in js["parametres"]])
        sig = {"kx": js["sigma_kx_rel"], "kz": js["sigma_kz_rel"]}
        print("Valeurs ajustées lues dans", fj)
    st, hz = idf.apply(cfg, params, x[:G])
    rows, deltas, corr = [], [], {}
    for j, d in enumerate(data):
        w = d.w + x[G + 2 * j] * cfg.heater.w_unc
        s = np.exp(x[G + 2 * j + 1])
        ta = s * response(st, hz, w / 2, d.freqs, d.p)
        t1, t2, tf, _ = richardson(st, hz, w, d.freqs, d.p, a.refine, a.max_nodes)   # Richardson, ordre 1
        t1, t2, tf = s * t1, s * t2, s * tf
        deltas.append(tf - ta)
        corr[d.name] = [[v.real, v.imag] for v in tf / ta]
        rows.append({"w": d.w, "f": d.freqs, "ecart_sigma": np.abs(tf - ta) / d.sigma})
        print(f"  {d.name} : écart EF/analytique {100 * np.max(np.abs(tf - ta) / np.abs(ta)):.2f} % "
              f"(= {np.max(np.abs(tf - ta) / d.sigma):.1f} σ de bruit ; incertitude EF {100 * np.max(np.abs(t2 - t1) / np.abs(ta)):.2f} %)")
    bias = idf.model_bias(cfg, params, x, data, deltas)
    ix, iz = [k for k, p in enumerate(params) if p.target]
    print("Effet sur le résultat (biais linéarisé dû au modèle analytique) :")
    for lab, i in (("kx", ix), ("kz", iz)):
        ref = sig[lab] if sig else np.nan
        verdict = ("négligeable" if abs(bias[i]) < 0.3 * ref else "à corriger : relancer fit --correction") \
            if np.isfinite(ref) and ref < 1 else "sans objet (non identifiable)"
        print(f"  {lab} : biais {100 * bias[i]:+.2f} % (incertitude {100 * ref:.3g} %) -> {verdict}")
    out = Path(a.out)
    json.dump(corr, open(out / "correction.json", "w"), indent=1)
    plots.verify_figure(rows, out / "verification.png")
    print(f"-> {out / 'correction.json'} (rapport T_EF/T_analytique), verification.png")


def cmd_slope(a):
    import glob
    for f in [Path(p) for p in sorted(glob.glob(a.data))]:
        r = dt.slope_method(dt.read(str(f)), a.f_max)
        print(f"{f.name} : k (pente) = {r['k']:.4g} W/m/K ± {100 * r['sigma_k_rel']:.2g} %, "
              f"k (quadrature) = {r['k_quadrature']:.4g} W/m/K, {r['n_points']} points")


def main(argv=None):
    p = argparse.ArgumentParser(description="OPTIMIST_3w : protocole 3ω, de la conception à l'identification de kx, kz")
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("design", "simulate", "fit", "verify", "slope"):
        s = sub.add_parser(name)
        if name != "slope":
            s.add_argument("--config", required=True, help="fichier JSON de l'expérience (examples/hemt.json)")
            s.add_argument("--target", type=int, default=None, help="couche à caractériser (0 = haut)")
            s.add_argument("--w_min", default=None, help="largeur minimale réalisable (remplace le JSON)")
            s.add_argument("--w_unc", default=None, help="incertitude sur la largeur (remplace le JSON)")
            s.add_argument("--out", default="resultats", help="dossier de sortie")
        if name in ("fit", "verify", "slope"):
            s.add_argument("--data", required=True, help="dossier ou motif des fichiers CSV de mesure")
        if name == "simulate":
            s.add_argument("--widths", default=None, help="largeurs imposées, ex. [2um,20um] (sinon : conception)")
            s.add_argument("--truth", default="[-0.2,0.15]", help="écarts relatifs vrais de kx, kz")
            s.add_argument("--seed", type=int, default=1)
            s.add_argument("--fem", action="store_true", help="réalité simulée par éléments finis (test réaliste)")
        if name == "verify":
            s.add_argument("--refine", type=float, default=3.0, help="finesse EF (calcul à refine et 2×refine)")
        if name in ("verify", "simulate"):
            s.add_argument("--max_nodes", type=int, default=80000,
                           help="plafond de nœuds du maillage EF (~0,3 Go pour 80000) ; finesse réduite au-delà")
        if name == "fit":
            s.add_argument("--correction", default=None, help="correction.json produit par verify")
        if name == "slope":
            s.add_argument("--f_max", type=float, default=None, help="borne haute du régime linéaire (Hz)")
    a = p.parse_args(argv)
    {"design": cmd_design, "simulate": cmd_simulate, "fit": cmd_fit, "verify": cmd_verify, "slope": cmd_slope}[a.cmd](a)


if __name__ == "__main__":
    main()