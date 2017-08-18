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
from planetplanet.photo import RadiativeEquilibriumMap, BandedCloudsMap, UniformMap
import matplotlib.pyplot as pl
import numpy as np
from numba import cfunc
import timeit, builtins

def compute(radiancemap = None, Phi = 0, Lambda = 0):
  '''
  
  '''
  
  # Instantiate the Trappist-1 system
  system = Trappist1(sample = True, phasecurve = True, airless = True, nbody = True, seed = 999)

  # Give `c` a large latitudinal offset in its hotspot just for fun
  system.c.Phi = Phi
  system.c.Lambda = Lambda
  system.c.nz = 99
  system.b.r /= 2
  
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
  plot(radiancemap = None, Lambda = -30)
  plot(radiancemap = RadiativeEquilibriumMap(), Lambda = -30)
  plot(radiancemap = BandedCloudsMap())
  pl.show()