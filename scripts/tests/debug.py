import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Trappist1
import rebound
import numpy as np
import matplotlib
import matplotlib.pyplot as pl
from tqdm import tqdm
G = 11466.9811868 # REARTH^3 / MEARTH / day^2
seed = 45

def BruteForce(vary = 'm', t = 1):
  '''
  
  '''
  
  # Set up the grid
  if vary == 'm':
    default = 1.
    grid = np.linspace(0.1, 10, 100)
  elif vary == 'e':
    default = 0.001
    grid = np.linspace(0.0001, 0.01, 100)
  else:
    raise ValueError("Invalid value for `vary`.")
  x = np.zeros_like(grid)
  
  # Brute force
  for i in tqdm(range(100)):
    
    # Get system params
    system = Trappist1(sample = True, seed = seed)
    
    # Set the mass or eccentricity of b
    if vary == 'm':
      system.b.m = grid[i]
    elif vary == 'e':
      system.b.ecc = grid[i]
    else:
      raise ValueError("Invalid value for `vary`.")
  
    # Set up the simulation
    sim = rebound.Simulation()  
    sim.G = G
    sim.dt = 0.0001
    sim.add(m = system.A._m)
    for planet in system.bodies[1:]:
      sim.add(primary = sim.particles[0], 
              m = planet._m, 
              a = planet.a,
              e = planet.ecc, 
              inc = planet._inc, 
              omega = planet._w,
              Omega = planet._Omega, 
              M = planet.M(0))
  
    # Integrate and save the x position of `c`
    sim.integrate(t)
    x[i] = sim.particles[2].x
  
  # Plot
  pl.plot(grid, x, 'g-', label = 'Brute force')

def Variational(vary = 'm', t = 1):
  '''
  
  '''

  # Get system params
  system = Trappist1(sample = True, seed = seed)
  
  # Set up the grid
  if vary == 'm':
    default = 1.
    grid = np.linspace(0.1, 10, 100)
    system.b.m = default
  elif vary == 'e':
    default = 0.001
    grid = np.linspace(0.0001, 0.01, 100)
    system.b.ecc = default
  else:
    raise ValueError("Invalid value for `vary`.")

  # Set up the simulation
  sim = rebound.Simulation()  
  sim.G = G
  sim.dt = 0.0001
  sim.add(m = system.A._m)
  for planet in system.bodies[1:]:
    sim.add(primary = sim.particles[0], 
            m = planet._m, 
            a = planet.a,
            e = planet.ecc, 
            inc = planet._inc, 
            omega = planet._w,
            Omega = planet._Omega, 
            M = planet.M(0))

  # Vary the mass or eccentricity of `b`
  var = sim.add_variation()
  var.vary(1, vary)

  # Integrate and compute the x position of `c`
  # and its derivative w/ respect to the mass or eccentricity
  sim.integrate(t)
  x0 = sim.particles[2].x
  dxdp = var.particles[2].x
  
  # Taylor expand the position of `c`
  x = x0 + (grid - default) * dxdp
  
  # Plot
  pl.plot(grid, x, 'b-', label = 'Variational')
  pl.plot(default, x0, 'ro')

BruteForce('m')
Variational('m')
pl.legend()
pl.show()
