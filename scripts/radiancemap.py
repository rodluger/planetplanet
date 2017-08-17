#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
radiancemap.py |github|
-----------------------

Illustrates how to specify different surface radiance maps.

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/radiancemap.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, unicode_literals
from planetplanet.photo import Trappist1
from planetplanet.constants import *
from planetplanet.photo import RadiativeEquilibriumMap
import matplotlib.pyplot as pl
import numpy as np
from numba import cfunc
import timeit, builtins

@cfunc("float64(float64, float64)")
def BandedCloudsMap(lam, z):
  '''
  
  '''
  
  # Planet c
  albedo = 0.3
  tnight = 40.
  irrad = 2.27 * SEARTH
  
  # Compute the temperature
  temp = ((irrad * np.cos(5 * z) ** 2 * (1 - albedo)) / SBOLTZ) ** 0.25
  
  # Compute the radiance
  a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
  b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
  return a / (np.exp(b) - 1.)

def compute(radiancemap = RadiativeEquilibriumMap, Phi = 0, Lambda = 0):
  '''
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = True, phasecurve = True, airless = True, nbody = True, seed = 999)

  # Give `c` a large latitudinal offset in its hotspot just for fun
  system.c.Phi = Phi
  system.c.Lambda = Lambda
  system.c.nz = 99
  
  # Give `c` a custom radiance map
  system.c.radiancemap = radiancemap

  # Compute an occultation by `b`
  system.quiet = True
  time = np.linspace(9552.9364, 9552.9564, 500)
  system.compute(time)
  
  return system

def plot(**kwargs):
  '''
  
  '''
  
  # Compute
  system = compute(**kwargs)
  
  # Plot the occultation
  fig, axlc, axxz, axim = system.plot_occultation('c', 9552.95, nz = 99, draw_ellipses = False, draw_terminator = False)
  return fig, axlc, axxz, axim

if __name__ == '__main__':
  plot(radiancemap = RadiativeEquilibriumMap, Lambda = -30)
  plot(radiancemap = BandedCloudsMap)
  pl.show()