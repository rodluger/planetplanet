#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
jwst_example.py
---------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.constants import *
from planetplanet.photo import Star, Planet, System
import matplotlib.pyplot as pl
import numpy as np

def retrograde_bc():
  '''
  A retrograde `c` occults `b` for 180 minutes.
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = False, oversample = 10, airless = True, seed = 1234)

  # Get the next occultation
  t = system.next_occultation(2000, system.b, occultor = system.c)
  system.b.limbdark = [0]
  time = np.arange(t - 0.15, t + 0.15, 5. / 1440.)

  # Compute and plot the light curve
  system.compute(time)
  system.plot_lightcurve(15.)

  # Observe it (one exposure)
  system.observe(stack = 1)
  pl.show()

# Instantiate the star
mstar = 0.0802
rstar = 0.121
teff = (0.000524 * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k')
  
# Instantiate `b`
RpRs = np.sqrt(0.7266 / 100)
r = RpRs * rstar * RSUN / REARTH    
b = Planet('b', m = 0.85, per = 1.51087081, inc = 89.65, r = r, t0 = 0, 
           Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., albedo = 0.3, 
           airless = True)

# Instantiate `c`
RpRs = np.sqrt(0.687 / 100)
r = RpRs * rstar * RSUN / REARTH
c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67, r = r, t0 = 0, 
           Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., albedo = 0.3, 
           airless = True)

# Instantiate the system
system = System(star, b, c, distance = 12)

#
time = np.arange(252.75, 253.50, 10 * MINUTE)

# Compute and plot the light curve
system.compute(time)
system.plot_lightcurve(15.)

# Observe it (one exposure)
system.observe(stack = 1)
pl.show()