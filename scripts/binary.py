#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
binary.py |github|
------------------

A binary star system.

  .. role:: raw-html(raw)
     :format: html
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/binary.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet import Planet, Star, Moon, System
from planetplanet.constants import *
from planetplanet.photo.maps import UniformMap
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''

    '''

    pass

def plot():
    '''

    '''

    # Instantiate the primary
    A = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'k', 
             limbdark = [0.4, 0.26])
    
    # Instantiate the secondary
    B = Star('B', m = 0.1, r = 0.1, nz = 31, t0 = 0,
             Omega = 0, w = 0, ecc = 0, inc = 90,
             per = 3, color = 'k', 
             limbdark = [0.4, 0.26])

    # Instantiate a planet
    b = Planet('b', m = 1, r = 1, nz = 31, t0 = 0,
              Omega = 0, w = 0, ecc = 0, inc = 90,
              per = 1, color = 'k', radiancemap = UniformMap())

    # Compute the light curve
    system = System(A, B, b, nbody = True, integrator = 'ias15', timestep = MINUTE)
    time = np.arange(0, 10., 0.1 * MINUTE)
    system.compute(time, lambda2 = 100)
    
    system.plot_lightcurve(wavelength = 100, interactive = True)

if __name__ == '__main__':

  plot()