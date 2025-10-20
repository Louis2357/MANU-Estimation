# ADRS instationnaire — Maillages uniformes & adaptatifs

Ce dépôt illustre l’**équation d’Advection–Diffusion–Réaction avec Source** (ADRS) en 1D, avec **solution exacte instationnaire fabriquée** et étude d’erreur sur **maillages uniformes** et **adaptatifs**.

## 📐 Problème
On résout sur :

$$ x\in[0,1] $$

$$ t\in[0,T] $$

$$
u_t + V\,u_x - K\,u_{xx} + \lambda\,u = f(x,t).
$$

On choisit une solution exacte **séparable** :

$$
u_{\mathrm{ex}}(t,x) = u(t)\,v(x),\quad u(t)=\sin(4\pi t).
$$

Alors

$$
u_t = 4\pi\cos(4\pi t)\,v(x),\quad
u_x = u(t)\,v_x(x),\quad
u_{xx} = u(t)\,v_{xx}(x),
$$

et le terme source

$$
f(x,t)= u_t + V\,u_x - K\,u_{xx} + \lambda\,u.
$$

> Comme \(u(t)\) oscille, le **résidu ne tend pas vers 0** (problème instationnaire).

---

## 📁 Scripts
- **`adrs_insta.py`** — maillages **uniformes**.
  - Trace l’**erreur \(L^2\)** à \(T/2\) et à \(T\) pour **plusieurs \(N_X\)** (RK3 par défaut).
  - Trace l’**erreur au milieu du domaine** au cours du temps pour **RK1..RK4**.
- **`adrs_insta_multiple_mesh_adap.py`** — maillages **adaptatifs** (métriques).
  - Source **dépendante du temps** (formule ci-dessus).
  - **Snapshots** de la solution à plusieurs instants (\(t/T=0.25,0.5,0.75,1.0\)).
  - **Critère d’arrêt mixte** (erreur \(L^2\) et taille du maillage).
  - **Métriques** : stationnaire (`final`), moyenne en temps (`avg`), maximum temporel (`max`).

---

## ⚙️ Discrétisation & stabilité
- **Espace** : différences finies **centrées**, CL **Dirichlet 0**. Stabilisation advection : viscosité numérique \(x_\nu=K+\tfrac{1}{2}|V|\Delta x\).
- **Temps** : **Runge–Kutta** (RK1..RK4) explicites.
- **CFL conservateur** :
  
$$ 
\Delta t \le \mathrm{CFL}\cdot \min\!\Big(\tfrac{\Delta x}{|V|},\,\tfrac{\Delta x^2}{2K}\Big).
$$

---

## 🧪 Mesures d’erreur
- **Norme \(L^2\)** à un instant \(t\) : \(\displaystyle \|e(t)\|_{L^2}\approx\sqrt{\sum_i (u_i-u_{\mathrm{ex},i})^2\,\Delta x}\) (trapèzes).
- **Erreur au milieu** : \(|u(x=L/2,t)-u_{\mathrm{ex}}(L/2,t)|\).

---

## 🔧 Métrologie pour l’adaptation (script 2)
Pour une interpolation P1, l’erreur locale \(\sim |u_{xx}|\,h^2\). On impose
\[\sqrt{|u_{xx}(t,x)|}\,h(x)\approx \text{cste},\]
soit une **métrique** \(M(t,x)=\sqrt{|u_{xx}(t,x)|+\varepsilon}\) et une **équidistribution** de \(\int M\,dx\).

**Modes de métrique :**
- `final` : \(M(x)=\sqrt{|u_{xx}(T,x)|}\) (stationnaire).
- `avg` : moyenne en temps \(\frac{1}{N_t}\sum_j \sqrt{|u_{xx}(t_j,x)|}\).
- `max` : enveloppe \(\max_j \sqrt{|u_{xx}(t_j,x)|}\) (intersection temporelle).

**Critère d’arrêt mixte (adaptation)** : poursuivre tant que **les deux** conditions ne sont pas **simultanément** vraies :  
1) \(\|e(T)\|_{L^2}\le \text{tol}\) et 2) \(N\le N_{\max}\).

---

## 📈 Attendus
- **Uniforme** : l’erreur \(L^2\) décroit avec \(N_X\) ; RK d’ordre plus élevé → meilleure précision temporelle (erreur milieu plus faible).
- **Adaptatif** : `final` place les nœuds là où la solution **finale** est courbe ; `avg` et `max` protègent les pics **pendant l’évolution**. Le **critère mixte** évite les arrêts trop précoces/tardifs.

---

## 📝 Astuces
- Vérifier la **CFL** après changement de paramètres.
- Si \(u_{xx}\) analytique indisponible, utiliser **différences finies**.
- Pour comparer `final/avg/max`, garder la **même tolérance** et une progression identique de `N`.
- Ajouter `plt.savefig("figure_...png", dpi=200)` pour exporter les figures.
