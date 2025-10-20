import numpy as np

def lebesgue_adapt(fct, fct_xx, L, R, itermax_lebesgue=10, tol=1e-3):
    """
    Intégration 'façon Lebesgue' avec pas local h(x) piloté par une métrique issue de f''(x).
    On rafraîchit la grille adaptative et on compare l'intégrale obtenue à la première passe.

    Paramètres
    ----------
    fct : callable
        f(x)
    fct_xx : callable
        f''(x) (donnée ici ; pas de point-fixe car dérivée seconde fournie)
        # pas de point-fixe car dérivée seconde donnée, comment se passer de fct_xx,
        # même pb qu'avec une EDP (note de la version d'origine)
    L, R : float
        Bords du domaine d'intégration
    itermax_lebesgue : int
        Nombre max de passes (pour tester la stabilité de l'intégrale)
    tol : float
        Seuil d'arrêt sur l'écart |IL[0] - IL[n]| (< tol)

    Retour
    ------
    IL  : np.ndarray
        Valeurs de l'intégrale à chaque passe
    nptL: np.ndarray
        Nombre de sous-intervalles utilisés à chaque passe
    eps : np.ndarray
        Valeur d'epsilon utilisée (constante ici) par passe
    """
    # bornes de pas local
    hmin = (R - L) / 100.0
    hmax = (R - L) / 3.0
    epsilon0 = 0.01  # valeur initiale indicative (non utilisée ensuite)

    nptL = np.zeros(itermax_lebesgue, dtype=int)
    eps  = np.zeros(itermax_lebesgue)
    IL   = np.zeros(itermax_lebesgue)
    
    fpp = np.zeros_like(u)
    fpp[1:-1] = (u[2:] - 2.0*u[1:-1] + u[:-2]) / (h*h)
    # bords : on copie la valeur intérieure (ou on peut mettre des schémas décentrés d'ordre 1)
    fpp[0]  = fpp[1]
    fpp[-1] = fpp[-2]

    for npt in range(itermax_lebesgue):
        # dans le code d'origine, epsilon est remis à 0.9 à chaque passe
        epsilon = 0.9

        x = L
        u = fct(x)

        while x < R:
            #uxx = fct_xx(x)
            fpp = np.zeros_like(u)
            fpp[1:-1] = (u[2:] - 2.0*u[1:-1] + u[:-2]) / (h*h)
            # bords : on copie la valeur intérieure (ou on peut mettre des schémas décentrés d'ordre 1)
            fpp[0]  = fpp[1]
            fpp[-1] = fpp[-2]
            
            uxx = fpp

            # métrique basée sur f'' :
            # metric = min( max(|f''|/epsilon, 1/hmax^2), 1/hmin^2 )
            metric = min(max(abs(uxx) / epsilon, 1.0 / (hmax**2)), 1.0 / (hmin**2))

            # pas local : h_loc = min( sqrt(1/metric), R - x )
            hloc = min(np.sqrt(1.0 / metric), R - x)
            # print(hloc**2 * metric)  # debug comme dans l'original

            u0 = u
            x  = x + hloc
            u  = fct(x)

            # intégration locale (trapèze)
            IL[npt] += hloc * (u + u0) / 2.0
            nptL[npt] += 1
            eps[npt] = epsilon

        # contrôle de stabilité par rapport à la première passe
        if npt > 0:
            error = abs(IL[0] - IL[npt])
            if error < tol:
                print("Approximate Lebesgues integral  epsilon, npt, error:",
                      epsilon, nptL[npt], error)
                break

        print("Approximate Lebesgues integral  epsilon, npt, error:",
              epsilon, nptL[npt], (abs(IL[0] - IL[npt]) if npt > 0 else np.nan))

    return IL[:npt+1], nptL[:npt+1], eps[:npt+1]


# ---------------- EXEMPLE D'UTILISATION ----------------
if __name__ == "__main__":
    # Exemple : f(x) = a x^2 + b x + c sin(4πx) + 10 exp(-100 (x-0.5)^2)
    a, b, c = 0.5, 10.0, 3.0
    def fct(x):
        return (a*x**2
                + b*x
                + c*np.sin(4.0*np.pi*x)
                + 10.0*np.exp(-100.0*(x - 0.5)**2))

    def fct_xx(x):
        r = x - 0.5
        return (2.0*a
                - 16.0*(np.pi**2)*c*np.sin(4.0*np.pi*x)
                + (400000.0*r**2 - 2000.0)*np.exp(-100.0*r**2))

    L, R = 0.0, 1.0
    IL, nptL, eps = lebesgue_adapt(fct, fct_xx, L, R, itermax_lebesgue=10, tol=1e-3)
    print("IL (par passe)  :", IL)
    print("nptL (par passe):", nptL)
    print("eps (par passe) :", eps)

    
