import numpy as np
import matplotlib.pyplot as plt

# Paramètres du problème
lam = 1.0
u0 = 1.0
T = 60.0  # 1 minute = 60 secondes

# ------------------------------------------------------
# Fonction pour Euler explicite
def euler_explicit(dt, T, lam, u0):
    N = int(T/dt)
    t = np.linspace(0, N*dt, N+1)
    u = np.zeros(N+1)
    u[0] = u0
    for n in range(N):
        u[n+1] = (1 - dt*lam)*u[n]
    return t, u

# Solution exacte
def u_exact(t, lam=1.0, u0=1.0):
    return u0*np.exp(-lam*t)

# ------------------------------------------------------
# 1) Comparaison pour dt=1s
dt = 1.0
t, u_num = euler_explicit(dt, T, lam, u0)
u_ex = u_exact(t, lam, u0)
err = np.abs(u_num - u_ex)

plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.plot(t, u_ex, 'k-', label="Solution exacte")
plt.plot(t, u_num, 'ro--', label="Euler explicite")
plt.xlabel("Temps (s)")
plt.ylabel("u(t)")
plt.title("Solution exacte vs. numérique (dt=1s)")
plt.legend()

plt.subplot(1,2,2)
plt.plot(t, err, 'b-')
plt.xlabel("Temps (s)")
plt.ylabel("Erreur absolue")
plt.title("Erreur en fonction du temps (dt=1s)")

plt.tight_layout()
plt.savefig("figure1.png")
plt.show()

# ------------------------------------------------------
# 2) Étude de convergence : erreur L2 en fonction du pas
dts = np.logspace(0, -3, 20)   # 20 valeurs de dt de 1 à 0.001
err_L2 = []
err_deriv_L2 = []

for dt in dts:
    t, u_num = euler_explicit(dt, T, lam, u0)
    u_ex = u_exact(t, lam, u0)

    # Erreur L2 sur la fonction
    e = u_num - u_ex
    L2 = np.sqrt(np.sum(e**2) * dt / T)
    err_L2.append(L2)

    # Approx dérivée numérique par différences finies
    du_num = np.diff(u_num)/dt
    du_ex = -lam * u_ex[:-1]
    e_der = du_num - du_ex
    L2_der = np.sqrt(np.sum(e_der**2) * dt / T)
    err_deriv_L2.append(L2_der)

# ------------------------------------------------------
# Tracé des erreurs L2
plt.figure(figsize=(12,5))

plt.subplot(1,2,1)
plt.loglog(dts, err_L2, 'o-')
plt.xlabel("Pas de temps Δt")
plt.ylabel("Erreur L2 sur u")
plt.title("Erreur L2 de la solution")

plt.subplot(1,2,2)
plt.loglog(dts, err_deriv_L2, 'o-r')
plt.xlabel("Pas de temps Δt")
plt.ylabel("Erreur L2 sur u'")
plt.title("Erreur L2 de la dérivée")

plt.tight_layout()
plt.savefig("figure2.png")
plt.show()

