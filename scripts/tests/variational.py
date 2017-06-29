import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
G = 11466.9811868 # REARTH^3 / MEARTH / day^2

def Integrate(system, t0 = 0, tf = 10, timestep = 0.001):
  '''
  
  '''
  
  # Add all the particles
  sim = rebound.Simulation()  
  sim.G = G
  sim.dt = timestep
  sim.add(m = system.A._m)
  for planet in system.bodies[1:]:
    sim.add(primary = sim.particles[0], 
            m = planet._m, 
            a = planet.a,
            e = planet.ecc, 
            inc = planet._inc, 
            omega = planet._w,
            Omega = planet._Omega, 
            M = 2. * np.pi / planet.per * ((t0 - planet.t0) % planet.per))
  
  # Vary the eccentricity of `b`
  var_de = sim.add_variation()
  var_dde = sim.add_variation(order = 2, first_order = var_de)
  var_de.vary(1, "e")
  var_dde.vary(1, "e")
  e0 = sim.particles[1].e
  
  # Integrate to tf
  sim.integrate(tf - t0)
    
  # DEBUG
  i = 1
  j = 2
  
  # The relative (x, y) positions of the two planets and their time and eccentricity derivatives
  x0 = sim.particles[i].x - sim.particles[j].x
  vx0 = sim.particles[i].vx - sim.particles[j].vx
  ax0 = sim.particles[i].ax - sim.particles[j].ax
  
  dxde = var_de.particles[i].x - var_de.particles[j].x
  dvxde = var_de.particles[i].vx - var_de.particles[j].vx
  daxde = var_de.particles[i].ax - var_de.particles[j].ax
  
  ddxdde = var_dde.particles[i].x - var_dde.particles[j].x
  ddvxdde = var_dde.particles[i].vx - var_dde.particles[j].vx
  ddaxdde = var_dde.particles[i].ax - var_dde.particles[j].ax

  y0 = sim.particles[i].y - sim.particles[j].y
  vy0 = sim.particles[i].vy - sim.particles[j].vy
  ay0 = sim.particles[i].ay - sim.particles[j].ay
  
  dyde = var_de.particles[i].y - var_de.particles[j].y
  dvyde = var_de.particles[i].vy - var_de.particles[j].vy
  dayde = var_de.particles[i].ay - var_de.particles[j].ay
  
  ddydde = var_dde.particles[i].y - var_dde.particles[j].y
  ddvydde = var_dde.particles[i].vy - var_dde.particles[j].vy
  ddaydde = var_dde.particles[i].ay - var_dde.particles[j].ay
  
  # The positions, velocities, and accelerations at a given eccentricity
  def _x(e):
    return x0 + (e - e0) * dxde + 0.5 * (e - e0) ** 2 * ddxdde
  def _vx(e):
    return vx0 + (e - e0) * dvxde + 0.5 * (e - e0) ** 2 * ddvxdde
  def _ax(e):
    return ax0 + (e - e0) * daxde + 0.5 * (e - e0) ** 2 * ddaxdde
  def _y(e):
    return y0 + (e - e0) * dyde + 0.5 * (e - e0) ** 2 * ddydde
  def _vy(e):
    return vy0 + (e - e0) * dvyde + 0.5 * (e - e0) ** 2 * ddvydde
  def _ay(e):
    return ay0 + (e - e0) * dayde + 0.5 * (e - e0) ** 2 * ddaydde
  
  # The positions as a function of time and eccentricity
  def x(t, e):
    return _x(e) + _vx(e) * t + 0.5 * _ax(e) * t ** 2
  def y(t, e):
    return _y(e) + _vy(e) * t + 0.5 * _ay(e) * t ** 2
  
  # The separation of the two particles as a function of time and eccentricity
  def r(t, e):
    return np.sqrt(x(e, t) ** 2 + y(e, t) ** 2)
  
  return r

# Initialize
np.random.seed(45)
trappist1 = Trappist1(sample = True, nbody = True)



# Simulation start and end (ends at an occultation)
t0 = 7702.677442913073 
tf = 7854.02 + 0.006

# Compute the variation
r = Integrate(trappist1, t0 = t0, tf = tf)

# Compute the distance of closest approach on a grid of time and initial eccentricity
logemin, logemax, ne = np.log10(trappist1.b.ecc) - 2, np.log10(trappist1.b.ecc) + 2, 100
tmin, tmax, nt = -0.25, 0.25, 100
t, loge = np.meshgrid(np.linspace(tmin, tmax, nt), np.linspace(logemin, logemax, ne))
distance = r(t, 10 ** loge) / (trappist1.b.r + trappist1.c.r)

# Plot the heatmap
pl.imshow(distance, origin = 'lower', extent = (tmin, tmax, logemin, logemax), aspect = 'auto')

# Occultation contour
pl.contour(distance, [0, 1], origin = 'lower', lw = 1, extent = (tmin, tmax, logemin, logemax), aspect = 'auto', colors = 'w', linestyles = '--')

# Label the true point
pl.plot(0, np.log10(trappist1.b.ecc), 'wo')
pl.show()
