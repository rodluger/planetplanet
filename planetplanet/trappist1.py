#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`ppo.py` - Python interface to C
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .ppo import Star, Planet, System
import numpy as np
np.random.seed(1234)
MSUN = 1.988416e30
RSUN = 6.957e8
G = 6.67428e-11
MEARTH = 5.9722e24
REARTH = 6.3781e6
DAYSEC = 86400.
AUM = 1.49598e11

def Trappist1(nl = 11, polyeps1 = 1e-8, polyeps2 = 1e-15, ttvs = True, uncertainty = True):
  '''
  
  '''
  
  # Account for the uncertainty?
  if not uncertainty:
    N = lambda mu, sigma: mu
  else: 
    N = lambda mu, sigma: mu + sigma * np.random.randn()
    
  # Instantiate the star
  mstar = N(0.0802, 0.0073)
  rstar = N(0.117, 0.0036)
  star = Star('star', m = mstar, r = rstar)
  
  # Parameters from Gillon et al. (2017) and Luger et al. (2017)
  # Mass for `h` is currently unconstrained, so basing it loosely on 
  # the mass distribution for `d`, which has a similar radius.
  planets = [None for i in range(7)]
  names = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
  periods = [(1.51087081, 0.60e-6), (2.4218233, 0.17e-5), (4.049610, 0.63e-4), (6.099615, 0.11e-4), 
             (9.206690, 0.15e-4), (12.35294, 0.12e-3), (18.767, 0.004)]
  transits = [(7322.51736, 0.00010), (7282.80728, 0.00019), (7670.14165, 0.00035), (7660.37859, 0.00038),
              (7671.39767, 0.00023), (7665.34937, 0.00021), (7662.55284, 0.00037)]
  masses = [(0.85, 0.72), (1.38, 0.61), (0.41, 0.27), (0.62, 0.58), (0.68, 0.18), (1.34, 0.88), (0.4, 0.3)]
  inclinations = [(89.65, 0.245), (89.67, 0.17), (89.75, 0.16), (89.86, 0.11), (89.680, 0.034),
                  (89.710, 0.025), (89.80, 0.075)]
  depths = [(0.7266, 0.0088), (0.687, 0.010), (0.367, 0.017), (0.519, 0.026), (0.673, 0.023), 
            (0.782, 0.027), (0.752, 0.032)]

  # Instantiate the planets
  for i in range(7):
  
    # Period and time of transit
    per = N(*periods[i])
    trn0 = N(*transits[i])
  
    # Positive mass
    m = 0
    while m <= 0:
      m = N(*masses[i])
  
    # Inclination in range [0, 90]
    inc = N(*inclinations[i])
    if inc > 90:
      inc = 180 - inc
  
    # Semi-major axis in AU from Kepler's law
    a = ((per * DAYSEC) ** 2 * G * (mstar * MSUN + m * MEARTH) / (4 * np.pi ** 2)) ** (1. / 3.) / AUM
  
    # Radius from Rp / Rstar
    mu = np.sqrt(depths[i][0] / 100)
    sig = 0.5 * depths[i][1] / 100 / mu
    RpRs = N(mu, sig)
    r = RpRs * rstar * RSUN / REARTH
  
    # Instantiate!
    planets[i] = Planet(names[i], m = m, per = per, inc = inc, a = a, r = r, trn0 = trn0, nl = nl)

  # Return the system
  return System(star, *planets, polyeps1 = polyeps1, polyeps2 = polyeps2, ttvs = ttvs)