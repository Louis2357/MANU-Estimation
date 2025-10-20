import numpy as np

def integrale_lebesgue(f, L, R, N, n_points_y=40000):
    """
    Approximation numérique de l'intégrale de Lebesgue de f sur [L, R].
    N : nombre de niveaux pour discrétiser l'image de f (en y).
    n_points_y : nombre de points pour estimer les mesures (grille en x).
    """

    # Évaluer f sur un grand nombre de points pour estimer les mesures
    # Plus ce nombre est grand, meilleure est l'approximation
    x = np.linspace(L, R, n_points_y)
    y = f(x)

    # Niveaux de discrétisation de l'image
    y_min, y_max = np.min(y), np.max(y)
    levels = np.linspace(y_min, y_max, N + 1)   # N+1 niveaux (y0, y1, ..., yN)

    # Calculer la mesure de chaque "tranche" { x | y_i ≤ f(x) < y_{i+1} }
    integral = 0.0
    longueur_cellule = (R - L) / (n_points_y - 1)   # maille en x
    for i in range(N):
        mask = (y >= levels[i]) & (y < levels[i+1])
        # Mesure de Lebesgue = longueur ≈ (#points) * Δx
        measure = np.sum(mask) * longueur_cellule
        # Contribution ≈ (valeur représentative de la tranche) * mesure
        integral += levels[i] * measure

    return integral
