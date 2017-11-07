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
    A = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'r', 
             limbdark = [0.4, 0.26])
    
    # Instantiate the secondary
    B = Star('B', m = 0.05, r = 0.1, nz = 31, t0 = -0.05,
             Omega = 0, w = 0, ecc = 0, inc = 87.2,
             per = 2, color = 'g', 
             limbdark = [0.4, 0.26])

    # Instantiate a planet
    b = Planet('b', host = 'A', m = 300, r = 2, nz = 31, t0 = 0.,
              Omega = 0, w = 0, ecc = 0, inc = 90,
              per = 0.2, color = 'b', radiancemap = UniformMap())
    
    # Instantiate a moon
    bI = Moon('bI', host = 'b', m = 0.1, r = 1, nz = 31, t0 = 0,
              Omega = 0, w = 0, ecc = 0, inc = 45,
              per = 0.05, color = 'b', radiancemap = UniformMap())
    
    # Compute the light curve
    system = System(A, B, b, bI, nbody = True, integrator = 'ias15', 
                    timestep = MINUTE)
    time = np.arange(-1., 10., 0.1 * MINUTE)
    system.compute(time, lambda2 = 100)
    system.plot_occultation('A', 0.)
    pl.show()

if __name__ == '__main__':

  plot()