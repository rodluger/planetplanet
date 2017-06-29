import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import rebound
import numpy as np
import matplotlib;
import matplotlib.pyplot as plt
G = 11466.9811868

def run_sim(e):
    sim = rebound.Simulation()
    sim.G = G
    sim.add(m=25785.694672633395)
    sim.add(primary=sim.particles[0],m=1.29, a=257.629687, e = 0.01, inc = np.pi / 2, Omega = 0, omega = 1.32)
    sim.add(primary=sim.particles[0],m=1.70, a=352.86464, e = e, inc = np.pi / 2, Omega = 0, omega = 1.89)
    sim.integrate(10.)
    return sim.particles[1].x

N=400
x_exact = np.zeros((N))
e_grid = np.linspace(0.01, 0.2, N)
for i,e in enumerate(e_grid):
    x_exact[i] = run_sim(e)

def run_sim_var(e):
    sim = rebound.Simulation()
    sim.G = G
    sim.add(m=25785.694672633395)
    sim.add(primary=sim.particles[0],m=1.29, a=257.629687, e = 0.01, inc = np.pi / 2, Omega = 0, omega = 1.32)
    sim.add(primary=sim.particles[0],m=1.70, a=352.86464, e = e, inc = np.pi / 2, Omega = 0, omega = 1.89)
    var_de = sim.add_variation()
    var_dde = sim.add_variation(order=2, first_order=var_de)
    var_de.vary(2, "e")
    var_dde.vary(2, "e")

    sim.integrate(10.)
    return sim.particles[1].x, var_de.particles[1].x, var_dde.particles[1].x

e_0 = 0.1
x, dxde, ddxdde = run_sim_var(e_0)
x_1st_order = np.zeros(N)
x_2nd_order = np.zeros(N)
for i,e in enumerate(e_grid):
    x_1st_order[i] = x + (e-e_0)*dxde
    x_2nd_order[i] = x + (e-e_0)*dxde + 0.5*(e-e_0)*(e-e_0)*ddxdde

fig = plt.figure(figsize=(6,4))
ax = plt.subplot(111)
ax.set_xlim(e_grid[0],e_grid[-1])
ax.set_xlabel("initial eccentricity of the outer planet")
ax.set_ylabel("$x$ position of inner planet after 10 orbits")
ax.plot(e_grid, x_exact, "-", color="black", lw=2)
ax.plot(e_grid, x_1st_order, "--", color="green")
ax.plot(e_grid, x_2nd_order, ":", color="blue")
ax.plot(e_0, x, "ro",ms=10);

plt.show()