import numpy as np
import matplotlib.pyplot as plt

# Paramètres
Lx, Ly = 1.0, 1.0
Nx, Ny = 50, 50
dx, dy = Lx/Nx, Ly/Ny
x = np.linspace(0, Lx, Nx+1)
y = np.linspace(0, Ly, Ny+1)
X, Y = np.meshgrid(x, y, indexing="ij")

Tmax = 0.1
dt = 1e-4
Nt = int(Tmax/dt)

# Paramètres physiques
v1, v2 = 1.0, 0.5
nu = 0.01
lam = 1.0

# Source gaussienne
Tc, k = 1.0, 50.0
sc1, sc2 = 0.5, 0.5
def f_source(t, X, Y):
    d2 = (X-sc1)**2 + (Y-sc2)**2
    return Tc*np.exp(-k*d2)

# Condition initiale
u = np.zeros((Nx+1, Ny+1))

# Fonction solution exacte (exemple fictif ici, tu peux définir une vraie solution test)
def u_exact(t, X, Y):
    return np.exp(-lam*t) * f_source(0, X, Y)

# Boucle en temps
for n in range(Nt):
    # gradients
    u_x = (np.roll(u, -1, axis=0) - np.roll(u, 1, axis=0)) / (2*dx)
    u_y = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2*dy)

    # laplacien
    u_xx = (np.roll(u, -1, axis=0) - 2*u + np.roll(u, 1, axis=0)) / dx**2
    u_yy = (np.roll(u, -1, axis=1) - 2*u + np.roll(u, 1, axis=1)) / dy**2
    lap = u_xx + u_yy

    # PDE : Euler explicite
    rhs = -v1*u_x - v2*u_y + nu*lap - lam*u + f_source(n*dt, X, Y)
    u = u + dt*rhs

    # Conditions de Dirichlet sur bords entrants
    u[0,:]  = 0 if v1>0 else u[0,:]
    u[-1,:] = 0 if v1<0 else u[-1,:]
    u[:,0]  = 0 if v2>0 else u[:,0]
    u[:,-1] = 0 if v2<0 else u[:,-1]

# Comparaison avec une "solution exacte"
u_ex = u_exact(Tmax, X, Y)

# Erreur L2
err_u = np.sqrt(np.sum((u-u_ex)**2)*dx*dy)
err_grad = np.sqrt(np.sum(((np.gradient(u)[0]-np.gradient(u_ex)[0])**2 + 
                           (np.gradient(u)[1]-np.gradient(u_ex)[1])**2)*dx*dy))

print("Erreur L2 sur u =", err_u)
print("Erreur L2 sur grad u =", err_grad)

# Figures
plt.figure(figsize=(15,5))
plt.subplot(1,3,1)
plt.pcolormesh(X,Y,u, shading='auto')
plt.colorbar()
plt.title("Solution numérique")

plt.subplot(1,3,2)
plt.pcolormesh(X,Y, np.abs(u-u_ex), shading='auto')
plt.colorbar()
plt.title("Erreur |u-uex|")

plt.subplot(1,3,3)
grad_diff = np.sqrt((np.gradient(u)[0]-np.gradient(u_ex)[0])**2 +
                    (np.gradient(u)[1]-np.gradient(u_ex)[1])**2)
plt.pcolormesh(X,Y,grad_diff, shading='auto')
plt.colorbar()
plt.title("Erreur sur ||grad u||")

plt.tight_layout()
plt.show()

