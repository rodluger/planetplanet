#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit.py
----------

A simple transit light curve.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def u1(lam):
  '''
  A really silly linear limb darkening law with a linear
  wavelength dependence.
  
  '''
  
  lam = np.atleast_1d(lam)
  result = 0.5 * (1 - (lam - 5) / 10) + 0.5
  result[lam < 5] = 1.
  result[lam > 15] = 0.5
  return result
  
# Instantiate the star
star = Star('A', m = 0.1, r = 0.1, nz = 11, color = 'k', limbdark = [u1])

# Planet b
b = Planet('b', m = 1, per = 2, inc = 90.4, r = 2., t0 = 0, 
           nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

# Compute the light curve, no optimization
system = System(star, b, batmanopt = False)
time = np.arange(-0.025, 0.025, 1 * SECOND)
system.compute(time)
flux1 = system.A.flux[:,0] / system.A.flux[0,0]

# Compute the light curve w/ batman optimization
system = System(star, b, batmanopt = True)
time = np.arange(-0.025, 0.025, 1 * SECOND)
system.compute(time)
flux2 = system.A.flux[:,0] / system.A.flux[0,0]

# Plot it
pl.plot(system.A.time, flux1, label = 'Standard')
pl.plot(system.A.time, flux2, label = 'Batman')
pl.legend()
pl.show()