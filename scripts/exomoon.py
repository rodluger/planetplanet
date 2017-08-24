#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
exomoon.py |github|
-------------------

A silly example of a star-planet-moon system. The
orbital parameters are unrealistic, but this showcases
how the code can handle moons.

  .. plot::
     :align: center
     
     from scripts import exomoon
     import matplotlib.pyplot as pl
     exomoon.plot()
     pl.show()

  .. warning:: The ability to model moons is still in beta and hasn't been thoroughly tested. \
               Please `let us know <https://github.com/rodluger/planetplanet/issues>`_ if you find any bugs!
              

  .. role:: raw-html(raw)
     :format: html
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/exomoon.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet import Planet, Star, Moon, System
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def plot():
  '''
  
  '''
  
  # Instantiate the star
  star = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', limbdark = [0.4, 0.26])

  # Planet b
  b = Planet('b', m = 50., per = 1, inc = 90., r = 1., t0 = 0, 
             nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = True, color = 'r')
           
  # Moon
  m = Moon('bI', host = 'b', m = 0., per = 0.1, inc = 90., r = 0.5, t0 = 0.5, 
           nz = 11, Omega = 10, w = 0., ecc = 0., phasecurve = True, color = 'b')

  # Compute the light curve
  system = System(star, b, m, nbody = True, integrator = 'ias15', timestep = MINUTE)
  time = np.arange(-1.1, 1.1, 0.1 * MINUTE)
  system.compute(time, lambda2 = 100)

  # Plot
  system.plot_lightcurve(wavelength = 100)

if __name__ == '__main__':
  plot()
  pl.show()