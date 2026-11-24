# OPTIMIST_3w

Mise en œuvre du protocole de mesure 3ω (voir `protocole_mesure_3w.md`) : conception des rubans, mise en forme des mesures, identification de kx et kz d'une couche par ajustement, et vérification du modèle par éléments finis.

## Installation

```bash
cd OPTIMIST_3w
python3 -m venv .venv
source .venv/bin/activate          # Windows : .venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirement.txt
python -m pytest -q                # environ 10 s
```

## Utilisation, phase par phase

| Phase du protocole | Commande | Sortie |
|---|---|---|
| A2 Conception | `python main.py design --config examples/hemt.json` | largeurs, longueur, prises, fréquences, courant ; incertitudes prédites, budget d'erreur |
| (test de la chaîne) | `python main.py simulate --config examples/hemt.json --out mesures [--fem]` | fichiers CSV au format laboratoire, vérité connue dans `verite.json` |
| F3 Validation du banc | `python main.py slope --data "ref/*.csv" --f_max 60` | k par la pente (et par la quadrature) sur l'échantillon de référence |
| G1–G5 Identification | `python main.py fit --config examples/hemt.json --data mesures` | kx, kz, incertitudes, corrélation, χ² réduit ; figure données/modèle/résidus |
| G6 Vérification | `python main.py verify --config examples/hemt.json --data mesures` | écart entre modèle analytique et modèle EF, biais induit sur kx, kz, `correction.json` |
| G6 Correction (si demandée) | `python main.py fit ... --correction resultats/correction.json` | ajustement corrigé de l'écart de modèle |

Options communes : `--target` (couche visée, 0 = haut), `--w_min` et `--w_unc` (lithographie, remplacent le JSON), `--out` (dossier de résultats, `resultats` par défaut).

Exemple complet :

```bash
python main.py design   --config examples/hemt.json --w_min 20nm --w_unc 5nm
python main.py simulate --config examples/hemt.json --w_min 20nm --w_unc 5nm --fem --out mesures
python main.py fit      --config examples/hemt.json --w_min 20nm --w_unc 5nm --data mesures
python main.py verify   --config examples/hemt.json --w_min 20nm --w_unc 5nm --data mesures
```

## Fichiers d'entrée

**Configuration (JSON)** : `examples/hemt.json` et `examples/silice_reference.json`.

- `layers` : couches du haut vers le bas, avec `d`, `kx` et `kz` (ou `k`), `rho`, `cp`, et `sigma_k` (incertitude relative a priori si la couche n'est pas la cible).
- `interfaces` : n−1 résistances `R` (m²K/W) avec leur incertitude `sigma`.
- `bottom` : `contact` (`R`, `sigma_ln`), `isothermal`, `adiabatic` ou `semi-infinite`.
- `heater` : métal, épaisseur `h`, contact `R_contact` et `sigma_R_contact`, `sigma_scale` (étalonnage du TCR), `w_unc` (incertitude sur la largeur).
- `design` : contraintes de conception ; `noise` : `v_noise` (V) et `rel_noise`.

**Mesures (un CSV par ruban)** :

```
# width_m = 1.0e-06        largeur mesurée au MEB
# L_taps_m = 9.8e-04       distance entre prises de tension
# R_taps_ohm = 1.9e+03     résistance entre prises à T0 (quatre fils)
# alpha_1_K = 2.5e-03      TCR étalonné
f_Hz,X3_V,Y3_V,sX3_V,sY3_V,V1_V,I_A
```

X3 et Y3 sont les composantes à 3ω, sX3 et sY3 leurs écarts-types sur les répétitions. V1 est l'amplitude à 1ω entre prises, I le courant efficace. La conversion est ΔT = 2(X3 + iY3)/(αV1) et p = I²R/L.

## Ce que fait le code

**Modèle d'ajustement (`analytic.py`)** : la solution multicouche anisotrope de Borca-Tasciuc, sous forme d'impédance,

  Z_haut = (Z_bas + tanh(Bd)/(kzB)) / (1 + kzB·Z_bas·tanh(Bd)),  avec B = √((kx/kz)η² + i2ωρc/kz).

Une résistance d'interface s'ajoute en série (Z → Z + R). La température moyenne du ruban vaut T0 = (p/π)∫Z(η) sin²(ηb)/(ηb)² dη. Le contact et la masse du ruban sont ajoutés par l'équation 58 d'Acquaroli. Un balayage complet prend quelques millisecondes.

**Conception (`identify.design`)** : pour chaque largeur candidate, le code détermine la fenêtre de fréquences (L ≥ 4,7 λ_S, 3f ≤ f_lockin) et le courant (électromigration, linéarité α·θ_DC, ΔT_DC, puissance visée). Il calcule ensuite les sensibilités et la matrice de Fisher. Pour chaque combinaison de 1 à N rubans, il évalue la borne de Cramér-Rao, puis décompose exactement la variance de kx et de kz entre le bruit et chaque paramètre mal connu (budget d'erreur).

**Ajustement (`identify.fit`)** : moindres carrés non linéaires avec a priori (maximum a posteriori), sur tous les rubans, en phase et en quadrature. Paramètres :
- inconnues : kx, kz ;
- nuisances : autres couches, interfaces, contact, ρc, contact arrière ;
- propres à chaque ruban : largeur réelle et étalonnage.

L'ajustement se fait par étapes : kx et kz liés, puis kx et kz libres, puis tous les paramètres. Il est relancé depuis le nominal pour contrôler l'unicité du minimum. Les inconnues sont bornées à un facteur 100 du nominal, les nuisances à ±5 σ. Une inconnue dont l'incertitude dépasse 100 % est déclarée non identifiable.

**Vérification (`fem.py`, `verify`)** : modèle éléments finis 2D, avec le ruban épais et conducteur et les résistances d'interface représentées par des couches fictives de 1 nm. Il est calculé à deux finesses de maillage et extrapolé (Richardson, convergence d'ordre 1 constatée). L'écart avec le modèle analytique est propagé au résultat : biais ≈ Cov·JᵀW·δ. Si le biais dépasse 0,3 σ, on relance `fit --correction`.

## Ce que montrent les exemples

Résultats obtenus avec les valeurs de `examples/hemt.json`, qui sont des ordres de grandeur.

- **SiN, photolithographie (2b ≥ 1 µm)** : kz à environ 32 %, kx non identifiable. Les deux tiers de la variance de kz viennent du contact ruban/SiN, inséparable de la résistance du film avec une seule épaisseur : il faut une méthode différentielle ou une série d'épaisseurs.
- **SiN avec résistances d'interface bien connues et lithographie électronique** : kz à environ 7 %, kx à environ 50 %. Le modèle analytique biaise kz d'environ −4 %, ce que `verify` détecte et que `fit --correction` retire.
- **Référence en silice** : la méthode de la pente retrouve k à l'erreur d'étalonnage du TCR près.

## Organisation

```
main.py                 commandes design, simulate, fit, verify, slope
optimist3w/config.py    lecture de la configuration (empilement, interfaces, ruban, contraintes, a priori)
optimist3w/analytic.py  modèle multicouche (ajustement et conception)
optimist3w/identify.py  paramètres, Fisher, conception, budget, ajustement, biais de modèle
optimist3w/data.py      CSV de mesure, conversion V -> ΔT, simulation, méthode de la pente
optimist3w/fem.py       éléments finis 2D (vérification)
optimist3w/plots.py     figures
examples/               hemt.json, silice_reference.json
tests/                  limites de Cahill et de résistance en série, EF/analytique, budget, pente, chaîne complète
```

## Mémoire et temps de calcul

Seules `simulate --fem` et `verify` utilisent les éléments finis. Le maillage est plafonné à `--max_nodes` nœuds (80 000 par défaut, soit environ 0,3 Go). Au-delà, typiquement pour les rubans de largeur inférieure à 0,5 µm, la finesse est réduite automatiquement et un message `[EF] ... finesse réduite` s'affiche. Le calcul utilise 2 cœurs au plus.

Sur la chaîne d'exemple (lithographie 20 nm), le pic mémoire est d'environ 0,3 Go et chaque commande dure moins de 30 s. Sur une machine très limitée, baisser le plafond, par exemple `--max_nodes 40000`.

## Limites

- Modèle 2D : ruban infiniment long. Les plots et les prises sont traités par le critère de longueur, sans calcul 3D.
- Couches planes et homogènes ; propriétés thermiques indépendantes de T ; pertes de surface négligées (justifié par l'éq. 60 d'Acquaroli).
- Les valeurs des matériaux dans les exemples sont des ordres de grandeur, et les conclusions dépendent fortement des incertitudes a priori, en particulier sur les résistances d'interface.