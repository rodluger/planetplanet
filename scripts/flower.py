#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
flower.py |github|
------------------

Computes a mutual transit among four planets with longitudes
of ascending node at right angles to each other. The premise
is silly, but this showcases the ability of the code to handle
mutual transits.

  .. plot::
     :align: center

     from scripts import flower
     flower._test()

  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/flower.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet import Planet, Star, System
import matplotlib.pyplot as pl
import numpy as np
np.random.seed(1234)

def _test():
    '''

    '''

    plot()

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
    star = Star('A', m = 0.1, r = 0.1, nz = 21, color = 'k', limbdark = [u1])

    # Planet b
    b = Planet('b', m = 1, per = 3, inc = 89.6, r = 5., t0 = 0,
               nz = 11, Omega = 0, w = 0., ecc = 0., phasecurve = False,
               color = 'r')

    # Planet c
    c = Planet('c', m = 1, per = 3 + 1e-5, inc = 89.6, r = 5., t0 = 0,
               nz = 11, Omega = 90, w = 0., ecc = 0., phasecurve = False,
               color = 'b')

    # Planet c
    d = Planet('d', m = 1, per = 3 + 2e-5, inc = 89.6, r = 5., t0 = 0,
               nz = 11, Omega = 180, w = 0., ecc = 0., phasecurve = False,
               color = 'g')

    # Planet c
    e = Planet('e', m = 1, per = 3 + 3e-5, inc = 89.6, r = 5., t0 = 0,
               nz = 11, Omega = 270, w = 0., ecc = 0., phasecurve = False,
               color = 'y')

    # System
    system = System(star, b, c, d, e)

    # Get the occultation light curves
    time = np.linspace(-0.02, 0.02, 100)
    system.compute(time)
    system.plot_occultation('A', 0.) #, gifname = 'flower')
    pl.show()

if __name__ == '__main__':

    plot()
