#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
flower.py
---------

Computes a mutual transit among four planets with longitudes
of ascending node at right angles to each other.

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet.photo import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

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
star = Star('A', m = 0.1, r = 0.1, nz = 21, color = 'k', limbdark = [u1])

# Planet b
b = Planet('b', m = 1, per = 3, inc = 89.6, r = 5., t0 = 0, 
           nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

# Planet c
c = Planet('c', m = 1, per = 3 + 1e-5, inc = 89.6, r = 5., t0 = 0, 
           nz = 11, Omega = 90, w = 0., ecc = 0., phasecurve = False, color = 'b')

# Planet c
d = Planet('d', m = 1, per = 3 + 2e-5, inc = 89.6, r = 5., t0 = 0, 
           nz = 11, Omega = 180, w = 0., ecc = 0., phasecurve = False, color = 'b')
           
# Planet c
e = Planet('e', m = 1, per = 3 + 3e-5, inc = 89.6, r = 5., t0 = 0, 
           nz = 11, Omega = 270, w = 0., ecc = 0., phasecurve = False, color = 'b')

# System
system = System(star, b, c, d, e)

# Get the occultation light curves
time = np.linspace(-0.02, 0.02, 100)
system.compute(time)
system.plot_occultation('A', 0.) #, gifname = 'flower')
pl.show()