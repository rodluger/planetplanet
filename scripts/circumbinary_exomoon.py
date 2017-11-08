#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
circumbinary_exomoon.py |github|
--------------------------------

A (crazy) example of the code's ability to handle nested multi-body systems.
Here we have a moon orbiting a planet orbiting one star in a binary system.
All bodies transit the primary at `t = 0`. This system isn't stable (and
the moon isn't even bound to the planet), but it highlights what the code can
do.

  .. plot::
     :align: center
     
     from scripts import circumbinary_exomoon
     circumbinary_exomoon._test()

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

    plot()

def plot():
    '''

    '''

    # Instantiate the primary (an M dwarf)
    A = Star('A', m = 0.1, r = 0.1, nz = 31, color = 'r', 
             limbdark = [0.4, 0.26])
    
    # Instantiate the secondary (a brown dwarf)
    B = Star('B', m = 0.05, r = 0.05, nz = 31, t0 = -0.05,
             Omega = 0, w = 0, ecc = 0, inc = 88.,
             per = 2, color = 'g', 
             limbdark = [0.4, 0.26])

    # Instantiate a planet (a hot Jupiter orbiting the primary)
    b = Planet('b', host = 'A', m = 300, r = 2, nz = 1, t0 = 0.,
              Omega = 0, w = 0, ecc = 0, inc = 90,
              per = 0.2, color = 'b', radiancemap = UniformMap())
    
    # Instantiate its moon (a puffy Mars-sized thing)
    bI = Moon('bI', host = 'b', m = 0.1, r = 1, nz = 1, t0 = 0,
              Omega = 0, w = 0, ecc = 0, inc = 45,
              per = 0.05, color = 'dodgerblue', radiancemap = UniformMap())
    
    # Compute the light curve
    system = System(A, B, b, bI, nbody = True, integrator = 'ias15', 
                    timestep = MINUTE)
    time = np.arange(-1., 10., 0.1 * MINUTE)
    system.compute(time, lambda2 = 100)
    system.plot_occultation('A', 0.)
    pl.show()

if __name__ == '__main__':

  plot()