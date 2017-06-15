#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
mutual_transit.py
-----------------

Computes and plots a hypothetical mutual transit event, where two large 
planets transit the star and occult each other simultaneously.

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
star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [u1])

# Planet b
b = Planet('b', m = 1, per = 3, inc = 89.8, r = 3., trn0 = 0, 
           nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

# Planet c
c = Planet('c', m = 1, per = 30, inc = 90., r = 3., trn0 = 0., 
           nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'b')

# System
system = System(star, b, c)

# Get the occultation light curves
time = np.linspace(-0.06, 0.06, 1000)
system.compute(time)
system.plot_occultation('A', 0.) #, gifname = 'mutual')
pl.show()