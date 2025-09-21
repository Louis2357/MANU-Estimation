# Séance 1 – Méthode d’Euler explicite pour une EDO

## 📘 Problème étudié
On considère l’équation différentielle ordinaire suivante :

$$
u'(t) = -\lambda\,u(t),\qquad u(0)=u_0=1,\qquad \lambda=1
$$

La solution exacte est connue :

\[
u_{\text{exact}}(t) = u_0 e^{-\lambda t}
\]

---

## 🔹 Schéma d’Euler explicite
Le schéma d’Euler explicite s’écrit :

\[
u_{n+1} = u_n + \Delta t \cdot f(t_n, u_n)
\]

où \( f(t,u) = -\lambda u \).  
Ainsi :

\[
u_{n+1} = u_n - \Delta t \, \lambda u_n
\]

ou encore :

\[
u_{n+1} = (1 - \lambda \Delta t)\, u_n
\]

---

## 📊 Objectifs du code
Le script `Euler_ODE_Errors.py` permet de :

1. Résoudre l’équation différentielle par Euler explicite pour un pas de temps fixé (\(\Delta t = 1s\)) sur un intervalle de 1 minute.  
2. Comparer la solution numérique avec la solution exacte.  
3. Tracer :
   - la solution exacte et la solution numérique,  
   - l’erreur en fonction du temps,  
   - les erreurs \(L^2\) de la solution et de sa dérivée en fonction du pas de temps (\(\Delta t \in [1, 0.001]\)).

---

## ⚙️ Données utilisées
- Durée de la simulation : \(T = 60 \, s\)  
- Pas de temps pour comparaison : \(\Delta t = 1s\)  
- Étude de convergence : 20 valeurs de \(\Delta t\) décroissantes de 1 à 0.001  
- Paramètres :
  - \(u_0 = 1\)  
  - \(\lambda = 1\)

---

## 📂 Organisation
