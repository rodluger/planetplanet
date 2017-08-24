#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
mutual_transit.py |github|
--------------------------

Computes and plots a hypothetical mutual transit event, where two large 
planets transit the star and occult each other simultaneously.

  .. plot::
     :align: center
     
     from scripts import mutual_transit
     import matplotlib.pyplot as pl
     mutual_transit.plot()
     pl.show()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/mutual_transit.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet import Planet, Star, System
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

def plot():
  '''
  
  '''
  
  # Instantiate the star
  star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [u1])

  # Planet b
  b = Planet('b', m = 1, per = 3, inc = 89.8, r = 3., t0 = 0., 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'r')

  # Planet c
  c = Planet('c', m = 1, per = 30, inc = 90., r = 3., t0 = 0., 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False, color = 'b')

  # System
  system = System(star, b, c)

  # Get the occultation light curves
  time = np.linspace(-0.6, 0.6, 10000)
  system.compute(time)
  fig, axlc, axxz, axim = system.plot_occultation('A', -0.05)

  return fig, axlc, axxz, axim
  
if __name__ == '__main__':
  plot()
  pl.show()