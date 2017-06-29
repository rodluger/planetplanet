#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`ppo.py` - Python interface to C
----------------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from .ppo import Star, Planet, System
import numpy as np
import matplotlib.pyplot as pl
import os
from tqdm import tqdm
MSUN = 1.988416e30
LSUN = 3.846e26
RSUN = 6.957e8
G = 6.67428e-11
MEARTH = 5.9722e24
REARTH = 6.3781e6
DAYSEC = 86400.
AUM = 1.49598e11
SBOLTZ = 5.670367e-8

__all__ = ['Trappist1']

def Trappist1(sample = True, airless = True, distance = 12, seed = None, **kwargs):
  '''
  
  '''
  
  # Randomizer seed
  if seed is not None:
    np.random.seed(seed)
  
  # Account for the uncertainty?
  if not sample:
    N = lambda mu, sigma: mu
  else: 
    N = lambda mu, sigma: mu + sigma * np.random.randn()
    
  # Instantiate the star; radius from Burgasser & Mamajek (2017)
  mstar = N(0.0802, 0.0073)
  rstar = N(0.121, 0.003)
  teff = (N(0.000524, 0.000034) * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
  star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k', **kwargs)
  
  # Parameters from Gillon et al. (2017) and Luger et al. (2017)
  # Mass for `h` is currently unconstrained, so basing it loosely on 
  # the mass distribution for `d`, which has a similar radius.
  planets = [None for i in range(7)]
  names = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
  
  periods = [(1.51087081, 0.60e-6), 
             (2.4218233, 0.17e-5), 
             (4.049610, 0.63e-4), 
             (6.099615, 0.11e-4), 
             (9.206690, 0.15e-4), 
             (12.35294, 0.12e-3), 
             (18.767, 0.004)]
  
  transits = [(7322.51736, 0.00010), 
              (7282.80728, 0.00019), 
              (7670.14165, 0.00035), 
              (7660.37859, 0.00038),
              (7671.39767, 0.00023), 
              (7665.34937, 0.00021), 
              (7662.55284, 0.00037)]
  
  masses = [(0.85, 0.72), 
            (1.38, 0.61), 
            (0.41, 0.27), 
            (0.62, 0.58), 
            (0.68, 0.18), 
            (1.34, 0.88), 
            (0.4, 0.3)]
            
  inclinations = [(89.65, 0.245), 
                  (89.67, 0.17), 
                  (89.75, 0.16), 
                  (89.86, 0.11), 
                  (89.680, 0.034),
                  (89.710, 0.025), 
                  (89.80, 0.075)]
                  
  depths = [(0.7266, 0.0088), 
            (0.687, 0.010), 
            (0.367, 0.017), 
            (0.519, 0.026), 
            (0.673, 0.023), 
            (0.782, 0.027), 
            (0.752, 0.032)]
  
  # These are eyeballed from Supplementary Figure 6 in Luger et al. (2017).
  # These are likely quite biased and model-specific. Need to re-think them.
  eccentricities = [(0.0005, 0.0001), 
                    (0.004, 0.001), 
                    (0.0004, 0.0003), 
                    (0.007, 0.0005), 
                    (0.009, 0.001), 
                    (0.004, 0.001), 
                    (0.003, 0.001)]
  
  # These we're just going to fix for now. We have no prior constraints on them.
  albedos = [(0.3, 0), (0.3, 0), (0.3, 0), (0.3, 0), (0.3, 0), (0.3, 0), (0.3, 0)]
  tnights = [(40., 0), (40., 0), (40., 0), (40., 0), (40., 0), (40., 0), (40., 0)]
  
  # Colors for plotting
  colors = ['firebrick', 'coral', 'gold', 'mediumseagreen', 'turquoise', 'cornflowerblue', 'midnightblue']
  
  # Instantiate the planets
  for i in range(7):
  
    # Period and time of transit
    per = N(*periods[i])
    t0 = N(*transits[i])
  
    # Positive mass
    m = 0
    while m <= 0:
      m = N(*masses[i])
  
    # Inclination in range [0, 90]
    inc = N(*inclinations[i])
    if inc > 90:
      inc = 180 - inc
    
    # Longitude of ascending node in degrees
    # A standard deviation of 0.3 is what Eric Agol got
    # in his Monte Carlo runs
    if i == 0:
      Omega = 0
    else:
      Omega = 0.3 * np.random.randn()
    
    # Longitude of pericenter (uniform over [0-360 deg])
    if sample:
      w = 360. * np.random.rand()
    else:
      w = 0.
    
    # Eccentricity
    ecc = 1
    while (ecc < 0) or (ecc >= 1):
      ecc = N(*eccentricities[i])

    # Radius from Rp / Rstar
    mu = np.sqrt(depths[i][0] / 100)
    sig = 0.5 * depths[i][1] / 100 / mu
    RpRs = N(mu, sig)
    r = RpRs * rstar * RSUN / REARTH
  
    # Albedo, night side temperature, effective temperature
    albedo = N(*albedos[i])
    tnight = N(*tnights[i])
  
    # Instantiate!
    planets[i] = Planet(names[i], m = m, per = per, inc = inc, r = r, t0 = t0, 
                        Omega = Omega, w = w, ecc = ecc, color = colors[i], 
                        tnight = tnight, albedo = albedo, 
                        airless = airless, **kwargs)

  # Return the system
  system = System(star, distance = distance, *planets, **kwargs)
  system._reset()
  return system