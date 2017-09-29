#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
wasp47.py |github|
------------------

Estimates the signal-to-noise on 15 stacked
secondary eclipses for the innermost planet (b)
in the WASP-47 system with JWST.

  .. plot::
     :align: center

     from scripts import wasp47
     wasp47._test()

  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/wasp47.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet.constants import *
from planetplanet import Star, Planet, System, jwst
from planetplanet import Wasp47
import matplotlib.pyplot as pl
import numpy as np
import astropy.units as u

def _test():
    '''

    '''

    fig1, _, fig2, _ = compute()
    pl.close(fig2)
    pl.show()

def _old_compute():
    '''

    '''

    # number of occultations observed
    Nocc = 1

    # WASP-47 stellar and observation parameters
    nout = 4.0                  # obseved out of transit durations [tint]
    lammin = 0.4                # min wavelength [um]
    lammax = 30.0               # max wavelength [um]
    Tstar = 5576.               # stellar temperature [K]
    Rs = 1.150                  # stellar radius [Solar Radii]
    d = 200.                    # distance to system [pc]

    # Additional params for WASP-47b
    Rp = 13.1                   # planet radius [Earth Radii]
    r = 0.0520                  # semi-major axis [AU]
    A = 0.0                     # Planet albedo
    e = 1.0                     # Planet emissivity
    i = 89.02                   # Orbital inclination [degrees]
    P = 4.15912                 # Orbital period [days]

    # Additional params for WASP-47d
    Rp = 3.71                   # planet radius [Earth Radii]
    r = 0.088                   # semi-major axis [AU]
    A = 0.0                     # Planet albedo
    e = 1.0                     # Planet emissivity
    i = 89.22                   # Orbital inclination [degrees]
    P = 9.0304                  # Orbital period [days]

    # Additional params for WASP-47e
    Rp = 1.87                   # planet radius [Earth Radii]
    r = 0.0173                  # semi-major axis [AU]
    A = 0.0                     # Planet albedo
    e = 1.0                     # Planet emissivity
    i = 86.2                    # Orbital inclination [degrees]
    P = 0.78961                 # Orbital period [days]

    # Convert some units to km
    Rs_km = Rs * u.solRad.in_units(u.km)
    Rp_km = Rp * u.earthRad.in_units(u.km)
    r_km = r * u.AU.in_units(u.km)
    P_mins = P * 24. * 60.

    # Planet temperature [K]
    Tplan = Tstar * ((1.-A)/e)**0.25 * (0.5*Rs_km/r_km)**0.5

    # transit duration
    tdur_mins = (P_mins / np.pi) * np.arcsin(np.sqrt((Rp_km + Rs_km) ** 2
                                   - r_km / (Rs_km * np.cos(i))) / r_km)

    # integration time [seconds]
    tdur = Nocc * tdur_mins * 60.0

    print("Planet Temperature : %.1f K" %Tplan)
    print("Transit Duration : %.2f mins" %(tdur_mins))

    # Create custom Kepler filter
    #k_band = jwst.create_tophat_filter(0.45, 0.8, dlam=0.1, Tput=0.5, name="Kepler")

    # Load Spitzer filters
    filters = jwst.get_spitzer_filter_wheel()
    atel_spitzer = np.pi * (0.85 / 2.) ** 2

    # Estimate for JWST
    fig1, ax1 = jwst.estimate_eclipse_snr(tint = tdur, nout = nout,
                                          lammin = lammin, lammax = lammax,
                                          Tstar = Tstar, Tplan = Tplan,
                                          Rs = Rs, Rp = Rp, d = d,
                                          thermal = True, filters = "MIRI"
                                         )

    # Estimate for an OST-like telescope, but using the JWST filters
    fig2, ax2 = jwst.estimate_eclipse_snr(tint = tdur, nout = nout,
                                          lammin = lammin, lammax = lammax,
                                          Tstar = Tstar, Tplan = Tplan,
                                          Rs = Rs, Rp = Rp, d = d,
                                          atel = 144., thermal = False)

    return fig1, ax1, fig2, ax2

def compute(seed = None):
    '''
    
    '''
    
    # Draw a system
    system = Wasp47(sample = True, seed = seed)

    # Re-instantiate and compute
    system = Wasp47(sample = True, seed = seed, oversample = 1)
    time = np.arange(0., 250., 5 * MINUTE)
    system.compute(time, lambda1 = 5, lambda2 = 10, R = 100)
    
    # HACK: Delete stellar occultations
    system.A.occultor *= 0
    system.A.occultor_hr *= 0
    system.A.flux[:,:] = system.A.flux[0,:]
    system.A.flux_hr[:,:] = system.A.flux_hr[0,:]
    
    # Plot
    #system.plot_lightcurve(interactive = True)
    
    fig, ax = system.observe(stack = 1, filter = 'f770w', alpha_err = 0.5) 
    pl.show()
    
if __name__ == '__main__':
    compute(1)