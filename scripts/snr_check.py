#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
snr_check.py |github|
---------------------

Computes the SNRs on the secondary eclipses of all TRAPPIST-1 planets
in the 15 micron JWST MIRI filter. This is used for benchmarking the 
code. The benchmarked values are

  .. code-block:: python

    b: SNR = 4.250
    c: SNR = 2.752
    d: SNR = 0.931
    e: SNR = 0.863
    f: SNR = 0.652
    g: SNR = 0.502
    h: SNR = 0.117

  .. role:: raw-html(raw)
     :format: html
  
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/snr_check.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                unicode_literals
from planetplanet.constants import *
from planetplanet.photo import LimbDarkenedMap
from planetplanet import Trappist1, Star, Planet, System, jwst
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''
    
    '''
    
    snrs = compute()
    
def compute():
    '''
    Simulate an observation of a secondary eclipse of each of
    the TRAPPIST-1 planets with JWST MIRI at 15 microns.
    
    .. plot::
         :align: center
         
         from scripts import jwst_example
         import matplotlib.pyplot as pl
         jwst_example.Triple_bc()
         pl.show()

    '''
    
    snrs = np.zeros(7)
    
    # Loop over the seven planets
    for planet in range(7):
    
        # Instantiate the system, no random draw
        system = Trappist1(distance = 12, oversample = 30, nbody = False, 
                           sample = False, quiet = True)
        
        # Get the body object
        body = system.bodies[planet + 1]
        
        # Set the t0 and radiance map for this planet
        # and force zero eccentricity.
        body.ecc = 0.
        body.t0 = 0.
        body.radiancemap = LimbDarkenedMap()
        
        # Fudge the inclinations of all other planets so 
        # that they don't transit or eclipse. This ensures
        # we get a clean secondary eclipse light curve.
        for b in system.bodies[1:]:
            if b.name != body.name:
                b.inc = 10.
        
        # Compute the light curve in the vicinity of secondary eclipse
        tecl = body.per / 2
        time = np.arange(tecl - 120 * MINUTE, tecl + 120 * MINUTE, 5 * MINUTE)
        system.compute(time, lambda1 = 10, lambda2 = 20, R = 100)

        # Let's re-center the time array for a prettier x axis
        system.A.time_hr -= tecl
        system.A.time -= tecl
    
        # Observe it (one exposure)
        system.observe(stack = 1, filter = 'f1500w')
        pl.close()
        
        snrs[planet] = system.filter.lightcurve.event_SNRs[0]
        
    return snrs

if __name__ == '__main__':
    
    snrs = compute()
    for planet, snr in zip(['b', 'c', 'd', 'e', 'f', 'g', 'h'], snrs):
        print("%s: SNR = %.3f" % (planet, snr))