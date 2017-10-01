#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
snr_check.py |github|
---------------------

Computes the SNRs on the secondary eclipses of all TRAPPIST-1 planets
in the 15 micron JWST MIRI filter. This is used for benchmarking the 
code. The benchmarked values are

  .. code-block:: python
    
        planetplanet    analytic     ratio
    b:     4.250          4.488      0.947
    c:     2.752          2.950      0.933
    d:     0.931          0.978      0.953
    e:     0.863          0.900      0.958
    f:     0.652          0.691      0.942
    g:     0.502          0.539      0.931
    h:     0.117          0.122      0.958

Note that the `planetplanet` estimates are lower than the `analytic` estimates
by a few percent because the latter neglects thermal/background noise.

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
    
    planetplanet_snr()
    analytic_snr()
    
def planetplanet_snr():
    '''
    Simulate an observation of a secondary eclipse of each of
    the TRAPPIST-1 planets with JWST MIRI at 15 microns and
    compute the SNR.
    
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

def analytic_snr():
    '''
    Compute the SNR analytically with back-of-the envelope
    calculations. This is a sanity check on the values we get
    with planetplanet.
    
    '''
    
    # Constants
    PLANCK = 6.62606983e-27
    GRAV = 6.67384e-8
    BOLTZMANN = 1.3806488e-16
    STEFANBOLTZMANN = 5.670373e-5
    MHYDROGEN = 1.00794*1.66053886e-24
    YEAR = 3.15569e7
    C = 29979245800e0
    AU = 1.49598e13
    PC = 3.08568025e18
    MSUN = 1.98892e33
    RSUN = 6.955e10
    MEARTH = 5.9742e27
    REARTH = 6.3781e8
    MMOON = 7.342e25
    RMOON = 1737.1e5
    MJUPITER = 1.8986e30
    RJUPITER = 6.9911e9
    AJUPITER = 4.84143144246472090 # AU
    KILOGRAM = 1e3
    METER = 1e2

    # area of JWST in cm^2
    area = 25.0 * 1e4  

    # Star properties
    teff = 2513.
    rs = 0.121

    # Eclipse duration in minutes
    tec = np.array([36.4,42.37,49.13,57.21,62.6,68.4,76.7])

    # Filter throughput
    eps = 0.3

    # Planet temperatures, zero albedo
    tp = np.array([400.,342.,288.,251.,219.,199.,172.7])

    # Planet radii, Earth radius
    rp = np.array([1.086,1.056,0.772,0.918,1.045,1.127,0.755])

    # Distance in parsecs
    d = 12.0

    # Filter wavelength, microns
    lam = 15.0

    # Filter width, microns
    dlam = 3.0

    # Compute the signal
    nu = C / (1e-4 * lam)
    dnu = dlam / lam * nu
    x = PLANCK * nu / (BOLTZMANN * tp)
    photons = 2. / (lam*1e-4) ** 2 / (np.exp(x) - 1.0) * np.pi * \
              (rp * REARTH / d / PC) ** 2 * dnu * area * eps * tec * 60.

    # Compute the noise
    x = PLANCK * nu / (BOLTZMANN * teff)
    photon_star = 2. / (lam * 1e-4) ** 2 / (np.exp(x) - 1.0) * np.pi * \
                 (rs * RSUN / d / PC) ** 2 * dnu * area * eps * tec * 60.

    # Compute the SNR
    return photons / np.sqrt(photon_star)

if __name__ == '__main__':
    
    snrs1 = planetplanet_snr()
    snrs2 = analytic_snr()
    print("    planetplanet    analytic     ratio")
    for planet, snr1, snr2 in zip(['b', 'c', 'd', 'e', 'f', 'g', 'h'], snrs1, snrs2):
        print("%s:     %.3f          %.3f      %.3f" % (planet, snr1, snr2, snr1/snr2))