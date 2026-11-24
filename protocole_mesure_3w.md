# Protocole expérimental — identification de kx et kz par la méthode 3ω

Ce protocole décrit la chaîne complète qui mène de la conception d'un transducteur 3ω à l'identification des conductivités thermiques dans le plan (kx) et hors plan (kz) d'une couche d'un empilement. Il est rédigé pour l'empilement HEMT SiN 50 nm / GaN 5 nm / AlGaN 20 nm / GaN 300 nm / AlN 60 nm / SiC 300 µm, mais s'applique à tout empilement de couches minces.

Il s'appuie sur la littérature (Cahill 1990 ; Borca-Tasciuc, Kumar et Chen 2001 ; Dames 2013 ; Zhang 2018 ; Acquaroli) et sur le code OPTIMIST_3w, utilisé pour la conception (`--design`) et pour la vérification des hypothèses.

---

## 0. Principe de l'identification

Un ruban métallique déposé sur l'échantillon sert à la fois de source de chaleur et de thermomètre. On injecte un courant I(t) = I₀ cos ωt entre les plots externes et on mesure la tension entre les prises internes. L'effet Joule chauffe le ruban à 2ω, sa résistance oscille à 2ω et la tension contient une composante à 3ω proportionnelle à l'oscillation de température :

$$\Delta T_{2\omega}=\frac{2\,V_{3\omega}}{\alpha\,V_{1\omega}} .$$

La conductivité n'est pas mesurée directement : elle est identifiée en ajustant un modèle de conduction de la chaleur sur la réponse mesurée ΔT₂ω(f), composantes en phase et en quadrature, de tous les rubans à la fois. La méthode de la pente de Cahill et la méthode différentielle sont des cas particuliers de cet ajustement, avec un modèle simplifié. Pour séparer kx et kz, il faut au moins deux rubans de largeurs différentes. Un ruban large (B = (b/d)·√(kz/kx) > 30) est sensible à kz seul ; un ruban étroit (B < 0,1) est sensible à kx et kz.

### Choix du modèle

Les trois approches suivantes résolvent la même physique. Elles diffèrent par leur coût de calcul et par les configurations qu'elles traitent.

| Modèle | Coût d'une évaluation | Ce qu'il traite | Rôle |
|---|---|---|---|
| Analytique multicouche (Borca-Tasciuc, formalisme matriciel, correction du ruban) | ~ms par balayage en fréquence | couches planes infinies, anisotropie, résistances de Kapitza, masse et contact du ruban | moteur de l'ajustement |
| Éléments finis (OPTIMIST_3w) | 0,1 à 1 s par balayage | géométrie quelconque : épaisseur du ruban, redistribution de la chaleur dans le métal, plots en 3D | conception de l'expérience, vérification des hypothèses |
| Volumes finis | comparable aux éléments finis | idem ; résistances d'interface simples à introduire | alternative équivalente aux éléments finis |

Le modèle rapide sert à ajuster. Le modèle complet sert à concevoir l'expérience, puis à vérifier que les simplifications du modèle rapide ne biaisent pas le résultat au-delà de l'incertitude de mesure.

---

## Phase A — Avant la fabrication

### A1. Caractérisation indépendante de l'empilement

Ces grandeurs sont les entrées du modèle et les a priori des paramètres de nuisance.

- Épaisseur de chaque couche : ellipsométrie, réflectivité X, coupe FIB-SEM ou TEM. Une erreur relative sur l'épaisseur d se reporte directement sur kz (R = d/kz).
- Rugosité de surface : AFM.
- Masse volumique et capacité thermique de chaque couche : littérature.
- Conductivité du substrat : mesure 3ω sur substrat nu avec le même masque si possible.
- Conductivité des couches qui ne sont pas la cible mais pèsent sur le signal : mesure séparée (TDTR, échantillon dédié) ou valeur de littérature avec une incertitude réaliste. Pour la cible GaN 300 nm, la conductivité du SiN supérieur représente 65 % de la variance de kz.

### A2. Conception de la géométrie

Lancer l'optimisation avec les contraintes réelles du procédé :

```bash
python main.py --config examples/hemt_gan_sic.json --design --target=<indice> \
               --w_min=<résolution litho> --w_unc=<incertitude sur la largeur> \
               --j_max=<A/m²> --f_lockin=<Hz> --fit_demo
```

Relever :

- les largeurs 2b retenues ;
- la longueur totale L et la position des prises de tension ;
- la fenêtre de fréquences ;
- le courant ou la puissance linéique ;
- les incertitudes prédites σ(kx), σ(kz) et le budget d'erreur.

Si l'incertitude prédite est insuffisante, agir avant de fabriquer sur le premier poste du budget : changer de lithographie, mesurer indépendamment une couche, ajouter un témoin.

### A3. Dessin du masque

- Rubans en configuration quatre plots : plots de courant I+ et I− aux extrémités, prises de tension V+ et V− à ±L/4 de la longueur totale.
- Longueur totale L ≥ 4,7 λ_S à la fréquence minimale (Dames, tableau 3), avec λ_S = √(D_S / 4πf) dans le substrat.
- Prises de tension fines et perpendiculaires au ruban ; plots éloignés de la zone mesurée.
- Plusieurs exemplaires de chaque largeur.
- Espacement entre rubans d'au moins 5 λ_S à la fréquence minimale, pour éviter le couplage thermique.
- Motifs de contrôle : lignes de mesure de largeur au MEB, marches de mesure d'épaisseur, croix d'alignement.
- Masque identique sur l'échantillon et sur les témoins de la méthode différentielle (substrat nu, empilement sans la couche visée, série d'épaisseurs).

### A4. Isolation électrique

Le ruban doit être isolé électriquement de l'échantillon. Sur le HEMT, le SiN joue ce rôle et le gaz 2D n'est pas contacté. Après une étape de délaminage, redéposer un diélectrique connu (Al₂O₃ ALD par exemple) et le déposer aussi sur le témoin.

---

## Phase B — Fabrication et métrologie

### B1. Préparation de surface

Nettoyage aux solvants (acétone, isopropanol), plasma O₂ léger. Sur échantillon délaminé : contrôle de propreté à la spectroscopie Auger et de rugosité à l'AFM.

### B2. Réalisation des rubans par lift-off

1. Enduction d'une bicouche de résine (profil en casquette).
2. Insolation et développement.
3. Évaporation par canon à électrons de Ti (accroche, 5 à 10 nm) puis Pt, épaisseur contrôlée par microbalance à quartz.
4. Lift-off.
5. Recuit sous N₂ pour stabiliser la résistance et le TCR.

### B3. Métrologie de ce qui a été fabriqué

Ces valeurs remplacent les valeurs nominales dans le modèle.

| Grandeur | Méthode | Usage |
|---|---|---|
| Largeur 2b | MEB en plusieurs points (moyenne, écart-type) | paramètre et a priori de largeur, souvent premier poste d'erreur sur kx |
| Épaisseur du métal | profilomètre ou AFM | masse thermique du ruban |
| Résistance R₀ entre prises | quatre fils | puissance linéique |
| Isolement ruban/échantillon | mesure à haute impédance | doit dépasser le GΩ |

---

## Phase C — Étalonnages

### C1. Coefficient de température (TCR) de chaque ruban

1. Échantillon sous vide sur platine régulée en température.
2. Mesure de R(T) en quatre fils à faible courant, sur 5 à 7 paliers (par exemple 290 à 330 K), stabilisation à chaque palier.
3. Ajustement linéaire : α = (1/R₀)·dR/dT.
4. L'incertitude sur α devient l'a priori du facteur d'échelle de ce ruban : une erreur sur α se reporte entièrement sur ΔT.

### C2. Échauffement propre

Tracer R en fonction de I² en courant continu pour vérifier l'ordre de grandeur de θ_DC et le bon fonctionnement du ruban.

---

## Phase D — Montage de mesure

### D1. Montage thermique

- Échantillon sous vide (pas de convection).
- Porte-échantillon régulé, T₀ stable à ±10 mK.
- Contact thermique en face arrière reproductible.

### D2. Montage électrique

| Élément | Réalisation |
|---|---|
| Source | sortie sinusoïdale de la détection synchrone, en série avec une résistance grande devant R₀ (courant quasi imposé) |
| Mesure du courant | résistance série connue à faible TCR, lue à 1ω |
| Mesure de tension | entre V+ et V−, entrée différentielle |
| Suppression du fondamental | pont de Wheatstone ou montage différentiel avec résistance de référence à faible TCR (V₃ω ≈ V₁ω / 1000) |
| Lectures | X₃ et Y₃ à 3f ; V₁ω à f |

### D3. Mesure du fond parasite

Remplacer le ruban par une résistance de référence à très faible TCR et mesurer le signal à 3ω résiduel (distorsion de la source et de l'électronique). C'est le plancher de mesure : le soustraire ou vérifier qu'il est négligeable.

---

## Phase E — Acquisition

### E1. Courant de mesure

Courant issu de la conception, en vérifiant α·θ_DC < 1 % et la densité de courant sous la limite d'électromigration.

### E2. Vérification de linéarité

À fréquence fixe, mesurer V₃ω pour 3 ou 4 courants. En échelle log-log, la pente doit valoir 3 (V₃ω ∝ I³). Un écart signale une non-linéarité, une distorsion ou une mauvaise compensation du fondamental.

### E3. Balayage en fréquence de chaque ruban

De f_min à f_max, sur une échelle logarithmique (8 à 15 points), avec 3f ≤ f_lockin. À chaque fréquence :

1. constante de temps de la détection synchrone ≥ 10 périodes ;
2. attente de 5 constantes de temps et d'au moins quelques périodes thermiques (2 à 3 périodes à 1 kHz d'après l'étude transitoire, davantage aux basses fréquences) ;
3. enregistrement de X₃, Y₃, V₁ω et I ;
4. N répétitions (environ 10).

Parcourir aussi les fréquences en sens inverse pour détecter une dérive.

### E4. Témoins et référence

Même protocole sur :

- l'échantillon de référence de conductivité connue (silice fondue par exemple) ;
- le substrat nu ;
- le témoin sans la couche visée (méthode différentielle).

### E5. Suivi de dérive

Relever T₀ et R₀ avant et après chaque balayage.

---

## Phase F — Mise en forme des données

### F1. Conversion en température et en puissance

$$\Delta T_{2\omega}=\frac{2\,(X_3+iY_3)}{\alpha\,V_{1\omega}},\qquad p=\frac{I_\text{eff}^2\,R_\text{prises}}{L_\text{prises}}\quad(\text{W/m}).$$

Utiliser la même convention (amplitude ou valeur efficace) pour V₃ω et V₁ω. Vérifier le signe sur l'échantillon de référence : la quadrature doit être négative (la température est en retard sur la puissance).

### F2. Incertitude par point

Écart-type des répétitions pour chaque fréquence et chaque composante. Il sert de poids dans l'ajustement.

### F3. Validation du banc

Sur l'échantillon de référence, la méthode de la pente (pente de la composante en phase en fonction de ln ω, égale à −p/(2πk)) doit redonner la conductivité de littérature à quelques pour cent. Ne pas poursuivre tant que ce n'est pas le cas.

---

## Phase G — Identification (problème inverse)

### G1. Modèle direct

Construit avec les dimensions mesurées (b, h, L, épaisseurs). Les paramètres sont classés en trois groupes :

| Groupe | Paramètres | Traitement |
|---|---|---|
| Inconnues | kx et kz de la couche visée | sans a priori |
| Nuisances | k des autres couches, résistances d'interface, contact ruban/couche, ρc, facteur d'échelle de chaque ruban, largeur réelle de chaque ruban | a priori gaussien (valeur, incertitude) |
| Fixés | paramètres sans influence sur le signal | valeur a priori, non ajustés |

### G2. Vérification d'identifiabilité

Avant tout ajustement, calculer la sensibilité du signal à chaque paramètre autour des valeurs nominales (`--design`). Un paramètre sans influence sur le signal doit être fixé : laissé libre, il rend l'ajustement instable.

### G3. Ajustement par étapes (valeurs de départ)

1. Conductivité du substrat : pente basse fréquence sur substrat nu ou sur le ruban large de l'échantillon.
2. kz de la couche, interfaces voisines incluses : méthode différentielle avec le ruban large,
   $$k_f=\frac{d_f}{2}\Big[\Big(\frac{b\,T_{2\omega}}{p}\Big)_\text{éch.}-\Big(\frac{b\,T_{2\omega}}{p}\Big)_\text{témoin}\Big]^{-1}_\text{moy}.$$
3. kx : première estimation par le ruban étroit, kz fixé à la valeur précédente.

### G4. Ajustement global

Tous les rubans, toutes les fréquences, composantes en phase et en quadrature simultanément ; moindres carrés non linéaires avec a priori (Levenberg-Marquardt), paramètres en logarithme :

$$\min_{\theta}\;\sum_{\text{rubans}}\sum_f\frac{\big|T_\text{modèle}(f;\theta)-T_\text{mesure}(f)\big|^2}{\sigma_f^2}\;+\;\sum_\text{nuisances}\Big(\frac{p-p_0}{\sigma_p}\Big)^2 .$$

Partir des valeurs de G3 ; relancer depuis plusieurs points de départ pour vérifier l'unicité du minimum.

### G5. Incertitudes

- Covariance au minimum : Cov = (JᵀWJ + Λ)⁻¹, avec J la matrice des dérivées du modèle, W les poids et Λ les précisions a priori. On en tire σ(kx), σ(kz) et leur corrélation.
- Ajouter les erreurs systématiques non ajustées, en premier lieu l'épaisseur d.
- Option : Monte-Carlo, en retirant les nuisances selon leurs a priori et en refaisant l'ajustement.

### G6. Validation du résultat

- Résidus aléatoires, sans tendance en fréquence ; χ² réduit proche de 1.
- Cohérence entre rubans : kz identique quand chaque largeur est ajustée séparément.
- Cas limite : la pente basse fréquence redonne la conductivité du substrat.
- Contrôle par éléments finis : recalculer la réponse aux valeurs ajustées avec le modèle complet (ruban épais, 3D si nécessaire). L'écart avec le modèle d'ajustement doit rester sous l'incertitude de mesure ; sinon, corriger ou changer de modèle.

### G7. Compte rendu

- kx et kz avec incertitudes et corrélation.
- Hypothèses : ce qu'inclut kz « effectif » (interfaces), modèle utilisé.
- Valeurs a priori des nuisances et leur origine.
- Figures de l'ajustement et des résidus.

---

## Points décisifs pour l'empilement HEMT

- **A1** : la conductivité du SiN supérieur, mesurée indépendamment, conditionne l'identification du GaN (65 % de la variance de kz).
- **G2** : avant fabrication, l'analyse d'identifiabilité montre qu'en l'état le GaN n'est identifiable qu'en mesurant la couche qui le masque, ou après délaminage.
- **B3** : la largeur réelle des rubans étroits est le premier poste d'erreur sur kx ; elle doit être mesurée au MEB.

---

## Annexe — Fiche de saisie

### Échantillon

| Champ | Valeur |
|---|---|
| Identifiant | |
| Empilement (couches, épaisseurs mesurées) | |
| Étape de délaminage | |
| Témoins associés | |

### Rubans

| Ruban | 2b nominal | 2b MEB (moy. ± σ) | h métal | L prises | R₀ (Ω) | α (1/K) ± σ | Isolement |
|---|---|---|---|---|---|---|---|
| | | | | | | | |

### Conditions de mesure

| Champ | Valeur |
|---|---|
| T₀ consigne / stabilité | |
| Pression | |
| Courant I_eff | |
| Fond parasite à 3ω | |
| Pente log V₃ω / log I (linéarité) | |

### Balayage (un tableau par ruban)

| f (Hz) | X₃ (V) | Y₃ (V) | V₁ω (V) | I (A) | N rép. | σ(X₃) | σ(Y₃) | ΔT phase (K) | ΔT quad. (K) |
|---|---|---|---|---|---|---|---|---|---|
| | | | | | | | | | |

### Résultat

| Grandeur | Valeur | Incertitude | Corrélation kx–kz | χ² réduit |
|---|---|---|---|---|
| kx (W/m/K) | | | | |
| kz (W/m/K) | | | | |

---

## Références

- D. G. Cahill, *Rev. Sci. Instrum.* 61, 802 (1990).
- T. Borca-Tasciuc, A. R. Kumar, G. Chen, *Rev. Sci. Instrum.* 72, 2139 (2001).
- C. Dames, *Annual Review of Heat Transfer* 16, 7 (2013).
- D. Zhang, mémoire de master, Utah State University (2018).
- L. N. Acquaroli, *3-omega method for thermal properties of thin film multilayers*.
