# `optim_adrs.py` — Optimisation linéaire pour ADRS 1D (maillages fixes & adaptatifs)

Ce projet implémente l’optimisation d’une **ADRS 1D** (Advection–Diffusion–Réaction avec Source) en exploitant la **linéarité** de l’équation d’état **par rapport aux contrôles**.  
On compare deux pipelines :
- **Référence** : maillage **uniforme fixe** (assez fin).  
- **Adaptatif** : chaque solution élémentaire est calculée sur **son maillage adapté**, puis **interpolée** vers une **grille commune** pour assembler les intégrales.

---

## 🔧 Problème

Sur \(x\in[0,1]\), on approche le stationnaire (par pseudo-temps explicite) de :
\[
u_t + V\,u_x - K\,u_{xx} + \lambda\,u \;=\; \sum_{i=1}^{N} x_i\,g_i(x),
\]
où \(x_i\) sont les **contrôles** (poids des sources) et \(g_i\) des **gaussiennes** centrées \(c_i\).

**Fonctionnelle à minimiser** :
\[
J(x) \;=\; \int_0^1 \big(u(x)-\text{Target}(x)\big)^2\,dx.
\]

---

## 🧠 Linéarité → base élémentaire \(\{T_i\}\)

En notant :
- \(T_0\) la solution pour contrôle nul \(x=0\),
- \(T_i\) la solution pour le contrôle unité \(e_i\) (1 sur i, 0 ailleurs),

on approxime (très bien dans ce cadre linéaire) :
\[
u(x) \;\approx\; T_0 \;+\; \sum_{i=1}^N x_i\,T_i.
\]

D’où le **système normal de moindres carrés** :
\[
A_{ij} = \langle T_i, T_j \rangle,\qquad
B_i = \langle \text{Target}-T_0,\; T_i \rangle,\qquad
A\,x^\star = B.
\]

---

## 🧩 Problème pratique : maillages différents (adaptation)

Si chaque \(T_i\) est calculé sur **son propre maillage adapté** \(x^{(i)}\), les produits
\(\langle T_i,T_j\rangle\) et \(\langle \text{Target}-T_0, T_i\rangle\) **ne sont pas** directement évaluables.

**Solution** :  
1) Construire une **grille de référence commune** \(x^{\text{ref}}\) **suffisamment fine** (pas plus petit que les pas locaux / facteur de sécurité).  
2) **Interpôler** toutes les fonctions \(T_i\), \(T_0\), \(\text{Target}\) vers \(x^{\text{ref}}\).  
3) Évaluer les **intégrales** avec une quadrature simple (trapèzes).  
4) **(Optionnel)** : ajouter un **test de Cauchy** sur la quadrature (raffiner \(x^{\text{ref}}\) jusqu’à stabilisation à \(10^{-8}\) par ex.), pour garantir que l’erreur d’intégration est **négligeable** vis-à-vis de l’erreur de modèle.

> Règle pratique : la quadrature (grille commune) doit être **bien plus fine** que les maillages des \(T_i\), de sorte que l’erreur d’intégration soit **1–2 ordres de grandeur** en dessous de l’erreur finale tolérée (cf. poly de cours).

---

## 🛠️ Ce que fait `optim_adrs.py`

- **Résolution stationnaire** (pseudo-temps explicite, stabilisation advection par viscosité numérique).  
- **Référence fixe** : calcule \(T_0\), \(\text{Target}\), \(\{T_i\}\) sur **maillage uniforme fin**.  
- **Adaptatif** : pour chaque \(i\),  
  1) construit un **maillage adapté** (métrique \(\sqrt{|u_{xx}|}\) + équidistribution),  
  2) résout \(T_i\) dessus,  
  3) **interpôle** toutes les fonctions sur une **grille commune** pour assembler \(A,B\).  
- **Résout** \(A\,x=B\) → obtient \(x^\star\) (fixe) et \(x^\star_{\text{adapt}}\) (adapté).  
- **Compare** les vecteurs optimaux, les reconstructions \(u_0 + \sum x_i T_i\) et l’erreur \(L^2\) vs `Target`.  
- **Trace** la **surface \(J(x_1,x_2)\)** (les autres \(x_k\) figés à l’optimum fixe).

---

## ▶️ Utilisation

### Prérequis
```bash
pip install numpy matplotlib
```

### Lancer
```bash
python optim_adrs.py
```

### Paramètres utiles (en tête du fichier)
- `N_ctrl` (nombre de contrôles), `centers` (positions), `alpha` (largeur gaussienne),
- `NX_ref` (taille maillage uniforme de référence),
- `NX_adapt` (taille visée des maillages adaptés),
- `CFL`, `Tmax`, `eps_rel`, `NTmax` (intégration pseudo-temps),
- facteurs de construction de la **grille commune** (`refine_factor`, `NXmax`).

---

## 📈 Sorties / Figures

- **Carte \(J(x_1,x_2)\)** (contours), avec marqueurs des optima **fixe** vs **adapté**.  
- **Barres** des composantes du contrôle optimal \(x^\star\) : **fixe** vs **adapté**.  
- **Comparaison** des reconstructions \(u^\star\) (fixe/adapté) à la cible `Target`.  
- Log des **erreurs \(L^2\)** finales.

---

## ✅ Validation recommandée

1) **Rafraîchir la cible** : `Target` = solution associée à un vecteur « vérité » `xcible` (fourni dans le code).  
2) **Monter** `NX_ref` jusqu’à stabiliser \(\|u^\star - \text{Target}\|_{L^2}\) (critère Cauchy).  
3) **Comparer** `xopt_fix` vs `xopt_adapt` :
   - Si la **grille commune** est assez fine et l’**adaptation** bien posée, les deux doivent être **proches** (differences mineures liées à l’interpolation et au pas de temps).  
4) **Tester** la sensibilité à `refine_factor` (grille commune) pour garantir que l’erreur d’intégration est **négligeable**.

---

## 📝 Notes / Astuces

- **Stabilisation advection** : viscosité numérique \(x_\nu = K + \tfrac{1}{2}|V|\,\Delta x\).  
- **CFL** conservateur : \(\Delta t \le \text{CFL}\cdot \min(\Delta x/|V|,\ \Delta x^2/(2K))\).  
- **Adaptation** : la métrique \(\sqrt{|u_{xx}|}\) concentre les nœuds là où la solution est **courbe**.  
- **Quadrature** : si la **grille commune** est trop grossière, \(A\) et \(B\) sont **biaisés** → optimum faux. Toujours vérifier la **stabilité** des intégrales (test de Cauchy simple).

---

## 📌 Résumé

- On **exploite la linéarité** en contrôles : \(u \approx T_0+\sum x_iT_i\) et on résout \(A\,x=B\).  
- Avec **maillages adaptés**, on **interpôle** tout vers une **grille commune** fine pour **intégrer** \(A_{ij}\), \(B_i\) sans biais.  
- On **compare** l’optimum obtenu avec celui d’un **maillage fixe de référence** (validation).  
- On **visualise** la **surface \(J(x_1,x_2)\)** et les **écarts** entre solutions optimales.

