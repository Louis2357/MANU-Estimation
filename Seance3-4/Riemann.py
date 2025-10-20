import numpy as np
import matplotlib.pyplot as plt

def integrale_riemann(f, a, b, N, method="midpoint"):
    """
    Intégrale de Riemann de f sur [a,b] avec N sous-intervalles (pas uniforme en x).
    method ∈ {"left", "right", "midpoint", "trapz"}.
    """
    if N <= 0:
        raise ValueError("N doit être > 0.")
    dx = (b - a) / N

    if method == "left":
        x = a + dx * np.arange(0, N)          # points à gauche
        return np.sum(f(x)) * dx, x, dx
    elif method == "right":
        x = a + dx * np.arange(1, N+1)        # points à droite
        return np.sum(f(x)) * dx, x, dx
    elif method == "midpoint":
        x = a + dx * (np.arange(0, N) + 0.5)  # milieux
        return np.sum(f(x)) * dx, x, dx
    elif method == "trapz":
        # règle des trapèzes (N sous-intervalles => N+1 points)
        x = np.linspace(a, b, N+1)
        fx = f(x)
        return np.trapz(fx, x), x, dx
    else:
        raise ValueError("method doit être 'left', 'right', 'midpoint' ou 'trapz'.")

# ------------------ Exemple d'utilisation ------------------
if __name__ == "__main__":
    # Données de l’exercice
    acoef, bcoef, ccoef = 0.5, 10.0, 3.0
    A, B = 0.0, 1.0
    def f_ex(x):
        return (acoef*x**2
                + bcoef*x
                + ccoef*np.sin(4*np.pi*x)
                + 10.0*np.exp(-100.0*(x-0.5)**2))

    N = 200          # nb de sous-intervalles
    method = "midpoint"  # "left" | "right" | "midpoint" | "trapz"

    I, xsamp, dx = integrale_riemann(f_ex, A, B, N=N, method=method)
    print(f"Intégrale ({method}, N={N}) ≈ {I:.6f}")

    # ----- TRACÉ -----
    xplot = np.linspace(A, B, 2000)
    yplot = f_ex(xplot)

    plt.figure(figsize=(9,4))
    plt.plot(xplot, yplot, label="f(x)", color="red")

    # afficher les points d’échantillonnage (sauf pour trapz où xsamp sont les nœuds)
    if method in {"left", "right", "midpoint"}:
        plt.scatter(xsamp, f_ex(xsamp), s=12, alpha=0.7, label=f"échantillons ({method})")

    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.title(f"f(x) sur [{A},{B}] — Riemann {method}, N={N},  ∫≈ {I:.6f}")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()
