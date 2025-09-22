import numpy as np
import matplotlib.pyplot as plt

# -------------------------
# Paramètres physiques
# -------------------------
L = 1.0
v = 1.0
nu = 0.01
lam = 1.0

# -------------------------
# Solution exacte et source
# u_ex(s) = exp(-10 (s - L/2)^2)
# f(s) = v u' - nu u'' + lam u
# -------------------------
def u_exact(s):
    r = s - 0.5*L
    return np.exp(-10.0*r**2)

def f_source(s):
    r = s - 0.5*L
    u = np.exp(-10.0*r**2)
    up = -20.0*r*u
    upp = (-20.0 + 400.0*r**2)*u
    return v*up - nu*upp + lam*u

# -------------------------
# Discrétisation spatiale
# -------------------------
def build_grid(N):
    # N points, h = L/(N-1), s_i = i*h, i=0..N-1
    h = L/(N-1)
    s = np.linspace(0.0, L, N)
    return s, h

# -------------------------
# Schémas en espace
# -------------------------
def upwind_first_derivative(u, h, v):
    """ u_s avec schéma amont (upwind) au centre du domaine.
        Bord gauche/droit seront traités via CL. """
    N = len(u)
    du = np.zeros_like(u)
    if v >= 0:
        # amont gauche
        # du[i] = (u[i]-u[i-1])/h pour i=1..N-1 ; du[0] fixé par CL (non utilisé)
        du[1:] = (u[1:] - u[:-1]) / h
        du[0] = du[1]
    else:
        # amont droit
        # du[i] = (u[i+1]-u[i])/h pour i=0..N-2 ; du[N-1] fixé par CL (non utilisé)
        du[:-1] = (u[1:] - u[:-1]) / h
        du[-1] = du[-2]
    return du

def laplacian_central(u, h):
    N = len(u)
    lap = np.zeros_like(u)
    lap[1:-1] = (u[0:-2] - 2.0*u[1:-1] + u[2:]) / (h*h)
    # bords traités via CL
    lap[0]  = lap[1]
    lap[-1] = lap[-2]
    return lap

# -------------------------
# Conditions aux limites
# Gauche: Dirichlet u(0)=u_l (= u_ex(0) ici)
# Droite: Neumann u_s(L)=0  -> (u_N - u_{N-1})/h = 0 -> u_N = u_{N-1}
# -------------------------
def apply_bc(u, h, u_l):
    u[0] = u_l           # Dirichlet gauche
    u[-1] = u[-2]        # Neumann homogène à droite

# -------------------------
# Pas de temps CFL
# dt <= 1 / ( |v|/h + 2nu/h^2 + lam )
# -------------------------
def cfl_dt(h, safety=0.9):
    denom = abs(v)/h + 2.0*nu/(h*h) + lam
    return safety/denom

# -------------------------
# Marche en temps vers stationnaire
# On intègre: u_t = -v u_s + nu u_ss - lam u + f
# -------------------------
def solve_stationary(N, tol=1e-6, itmax=2_000_000, report=False):
    s, h = build_grid(N)
    u = np.zeros(N)
    f = f_source(s)

    # CI cohérente (option): u^0 = u_ex
    # u = u_exact(s).copy()

    # CL gauche = valeur exacte (cohérence)
    u_l = u_exact(0.0)

    dt = cfl_dt(h, safety=0.9)
    res_hist = []

    # initial res
    u_old = u.copy()
    apply_bc(u, h, u_l)
    du = upwind_first_derivative(u, h, v)
    lap = laplacian_central(u, h)
    rhs = -v*du + nu*lap - lam*u + f
    res0 = np.sqrt(h*np.sum(rhs**2))
    if res0 == 0: res0 = 1.0
    res = res0
    res_hist.append(res/res0)

    n=0
    while n < itmax and (res/res0) > tol:
        n += 1
        # un pas d'Euler explicite
        du  = upwind_first_derivative(u, h, v)
        lap = laplacian_central(u, h)
        rhs = -v*du + nu*lap - lam*u + f
        u   = u + dt*rhs

        apply_bc(u, h, u_l)

        # convergence L2 du résidu (ou ||u^{n+1}-u^n||)
        res = np.sqrt(h*np.sum(rhs**2))
        res_hist.append(res/res0)

        # (option) break si la variation solution est petite
        # diff = np.sqrt(h*np.sum((u - u_old)**2))
        # u_old = u.copy()

    if report:
        print(f"N={N}, steps={n}, res/res0={res/res0:.3e}, dt={dt:.3e}, h={h:.3e}")
    return s, u, res_hist, h

# -------------------------
# Normes L2 et H1-semi (dérivée)
# -------------------------
def l2_norm(vec, h):
    return np.sqrt(h*np.sum(vec**2))

def h1_seminorm(u, uex, h):
    # ||(u-uex)'||_L2
    du   = np.zeros_like(u)
    duex = np.zeros_like(u)
    # centrée au centre du domaine
    du[1:-1]   = (u[2:]   - u[:-2])/(2*h)
    duex[1:-1] = (uex[2:] - uex[:-2])/(2*h)
    du[0], du[-1] = du[1], du[-2]
    duex[0], duex[-1] = duex[1], duex[-2]
    return l2_norm(du - duex, h)

# -------------------------
# (A) Une exécution avec N=100 + figures demandées
# -------------------------
N = 100
s, u, res_hist, h = solve_stationary(N, tol=1e-8, report=True)
uex = u_exact(s)

# Figures: solution, convergence, (u-uex)
plt.figure(figsize=(15,4))

plt.subplot(1,3,1)
plt.plot(s, uex, '-', label='u_exact')
plt.plot(s, u , '--', label='u_h')
plt.xlabel('s'); plt.ylabel('u'); plt.title('Solution')
plt.legend(); plt.grid(True)

plt.subplot(1,3,2)
it = np.arange(len(res_hist))
plt.semilogy(it, res_hist, '-')
plt.xlabel('itération'); plt.ylabel('||R||/||R0||')
plt.title('Convergence vers le stationnaire'); plt.grid(True)

plt.subplot(1,3,3)
plt.plot(s, u-uex, '-')
plt.xlabel('s'); plt.ylabel('u_h - u_exact')
plt.title('Erreur ponctuelle'); plt.grid(True)

plt.tight_layout()
plt.show()

# -------------------------
# (B) Étude de convergence sur 5 maillages
# N = 3, 6, 12, 24, 48  (par ex.)
# -------------------------
Ns = [3, 6, 12, 24, 48]
hs, errL2, errH1 = [], [], []
for Ntest in Ns:
    s, u, res_hist, h = solve_stationary(Ntest, tol=1e-8)
    uex = u_exact(s)
    hs.append(h)
    errL2.append(l2_norm(u - uex, h))
    errH1.append(h1_seminorm(u, uex, h))

# Tracé erreurs vs h
plt.figure(figsize=(12,4))

plt.subplot(1,2,1)
plt.loglog(hs, errL2, 'o-')
plt.xlabel('h'); plt.ylabel('||u_h - u_exact||_{L2}')
plt.title('Erreur L2 vs h'); plt.grid(True, which='both')

plt.subplot(1,2,2)
plt.loglog(hs, errH1, 'o-')
plt.xlabel('h'); plt.ylabel('||u_h - u_exact||_{H1-semi}')
plt.title('Erreur H1(semi) vs h'); plt.grid(True, which='both')

plt.tight_layout()
plt.show()
