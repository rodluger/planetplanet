#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
transit.py |github|
-------------------

A simple transit light curve. Here we compare it to one generated with :py:obj:`batman`.

  .. plot::
     :align: center
     
     from scripts import transit
     transit._test()
     
  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/transit.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import os, sys
sys.path.insert(1, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from planetplanet import Planet, Star, System
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''
    
    '''
    
    plot()

def plot():
    '''
    
    '''
    
    # Instantiate the star
    star = Star('star', m = 0.1, r = 0.1, nz = 31, color = 'k', 
                limbdark = [0.4, 0.26])

    # Planet b
    planet = Planet('planet', m = 1, per = 0.5, inc = 90.4, r = 2., t0 = 0, 
                    nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

    # Compute the light curve, no optimization
    system = System(star, planet, batmanopt = False)
    time = np.arange(-0.025, 0.025, 0.1 * MINUTE)
    system.compute(time, lambda1 = 0.5, lambda2 = 2.)
    flux1 = system.star.flux[:,0] / system.star.flux[0,0]
    
    # Compute the light curve w/ batman optimization
    system = System(star, planet, batmanopt = True)
    system.compute(time, lambda1 = 0.5, lambda2 = 2.)
    flux2 = system.star.flux[:,0] / system.star.flux[0,0]

    # Plot it
    fig = pl.figure()
    pl.plot(system.star.time, flux1, label = 'planetplanet')
    pl.plot(system.star.time, flux2, '--', label = 'batman')
    pl.legend()

    # Plot it interactively
    _, ax, _, _ = system.plot_occultation('star', 0, spectral = False)
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    pl.show()

if __name__ == '__main__':

    plot()