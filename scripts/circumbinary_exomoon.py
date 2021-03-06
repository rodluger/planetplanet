#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
circumbinary_exomoon.py |github|
--------------------------------

A (crazy) example of the code's ability to handle nested multi-body systems.
Here we have a moon orbiting a planet orbiting two stars in a binary system.
All bodies transit the primary close to `t = 0`. This system is probably not
stable, but it highlights what the code can do!

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

    # Instantiate the primary
    A = Star('A', m = 0.1, r = 0.12, nz = 31, color = 'r',
             teff = 2500,
             limbdark = [0.4, 0.26],
             x0 = 22.5, y0 = 0, z0 = 10,
             vx0 = 0, vy0 = -100, vz0 = 2500,
             cartesian = True)

    # Instantiate the secondary
    B = Star('B', m = 0.1, r = 0.08, nz = 31, teff = 3500,
             color = 'g', limbdark = [0.9, 0.26],
             x0 = -22.5, y0 = 20, z0 = -10,
             vx0 = 0, vy0 = 100, vz0 = -2500,
             cartesian = True)

    # Instantiate a planet (a hot Jupiter orbiting the binary)
    planet = Planet('planet', m = 300, r = 2, nz = 1, host = None,
                    color = 'b', radiancemap = UniformMap(),
                    x0 = -20, y0 = 5, z0 = -500,
                    vx0 = 1000, vy0 = 0, vz0 = 300,
                    cartesian = True)

    # Instantiate its moon (a puffy Mars-sized thing)
    moon = Moon('moon', host = 'planet', m = 0.1, r = 1, nz = 1, t0 = 0.04,
                Omega = 0, w = 0, ecc = 0, inc = 85,
                per = 0.1, color = 'dodgerblue', radiancemap = UniformMap())

    # Compute the light curve
    system = System(A, B, planet, moon, nbody = True, integrator = 'ias15',
                    timestep = MINUTE)
    time = np.arange(-0.02, 0.03, 0.1 * MINUTE)
    system.compute(time, lambda1 = 0.5, lambda2 = 0.9)
    system.plot_occultation('A', time = 0., wavelength = 0.75,
                            full_lightcurve = True)
    pl.show()

if __name__ == '__main__':

  plot()
