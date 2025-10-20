# Intégration Riemann vs « Lebesgue-genre » + Adaptation a posteriori

## 🎯 Objectif
Comparer trois approches pour approximer l’intégrale de `f` sur `[0,1]` avec une précision cible `≤ 1e-3` (référence fournie `I_ref ≈ 6.94`).

Fonction test :
```
f(x) = a*x^2 + b*x + c*sin(4*pi*x) + 10*exp(-100*(x-0.5)^2)
a = 0.5, b = 10, c = 3,  Left = 0, Right = 1
```

## 🧩 Méthodes

### 1) Riemann (pas uniforme en x)
- Maillage uniforme : `x_i = Left + i*(Right-Left)/N`.
- Milieux des intervalles (midpoint rule), pas `dx` constant.
- Contrôle d’arrêt par **suite de Cauchy** : on double `N` jusqu’à ce que `|I_{2N} - I_N| < tol/2` (et, si disponible, `|I_{2N} - I_ref| < tol`).

### 2) « Lebesgue-genre » (pas uniforme en y)
- Échantillonnage très fin en `x` (grand `Nx_fine`) pour obtenir le nuage `(x, f(x))`.
- Discrétisation **uniforme en y** en `Ny` bacs.
- Pour chaque bac `[y_k, y_{k+1})`, on approxime la **mesure en x** (compte d’échantillons dans le bac × `dx_fine`).
- Intégrale approchée : `sum( y_center_k * measure_x_k )`.
- On augmente `Ny` jusqu’au critère de **Cauchy** (et `|I - I_ref| < tol` si la référence est connue).

### 3) Adaptation par « métrique » (contrôle de l’erreur d’interpolation)
- Idée classique P1 : erreur locale ~ `|f''(x)| * h(x)^2`.  
  Pour équilibrer l’erreur, imposer `sqrt(|f''(x)|) * h(x) ≈ cste`.
- Définir la **métrique** `M(x) = sqrt(|f''(x)| + eps)` (eps petit pour éviter 0).
- Placer les nœuds `x_k` tels que la **longueur métrique** soit uniforme :  
  `∫_0^{x_k} M(s) ds = k/N * ∫_0^1 M(s) ds`.
- Intégrer ensuite par midpoint sur ce **maillage adapté**.  
→ En pratique, cette méthode atteint la tolérance avec **moins de points** quand `f` a des zones « raides » (ici la bosse gaussienne).

---

## 📈 Observations attendues
- **Riemann (x-uniforme)** : simple et robuste ; `N` peut devoir être élevé à cause de la bosse gaussienne.
- **Lebesgue-genre (y-uniforme)** : exploite la distribution des **valeurs** de `f`, mais dépend d’un pré-échantillonnage fin en `x`.
- **Adaptation métrique** : maillage concentré autour des fortes courbures (`|f''|` grand) → **meilleur compromis précision / nombre de points**.

---

## 🛠️ Bonnes pratiques
- Fixer une **tolérance** (`tol = 1e-3`) et utiliser la **suite de Cauchy** pour décider de l’arrêt.
- Si `I_ref` est inconnu, conserver uniquement le **critère Cauchy**.
- Si `f''` n’est pas disponible analytiquement, utiliser une **approximation numérique** (différences finies).
- Journaliser `N`, `Ny`, la tolérance et la valeur trouvée pour la **reproductibilité**.

