#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
spitzer_example.py |github|
---------------------------

Sample observations of TRAPPIST-1 PPOs with Spitzer.

  .. plot::
     :align: center
     
     from scripts import spitzer_example
     spitzer_example._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/spitzer_example.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet.constants import *
from planetplanet import Star, Planet, System, LimbDarkenedMap, \
                         RadiativeEquilibriumMap
from planetplanet import jwst
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''
    
    '''
    
    plot()
    pl.show()

def Stacked_bc_all_filters(N = 100, albedo = 0.0, airless = True):
     '''
     `N` stacked exposures of `b` occulting `c` observed in 
     both Warm Spitzer filters

     '''
     
     if airless:
        radiancemap = RadiativeEquilibriumMap()     
     else:
        radiancemap = LimbDarkenedMap()
     
     # Instantiate the star
     mstar = 0.0802
     rstar = 0.121
     teff = (0.000524 
             * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
     star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k')

     # Instantiate `b`
     RpRs = np.sqrt(0.7266 / 100)
     r = RpRs * rstar * RSUN / REARTH
     b = Planet('b', m = 0.85, per = 1.51087081, inc = 89.65, r = r, t0 = 0,
                Omega = 0, w = 0, ecc = 0, color = 'firebrick', tnight = 40., 
                albedo = albedo, radiancemap = radiancemap, phasecurve = False)

     # Instantiate `c`
     RpRs = np.sqrt(0.687 / 100)
     r = RpRs * rstar * RSUN / REARTH
     c = Planet('c', m = 1.38, per = 2.4218233, inc = 89.67, r = r, t0 = 0,
                Omega = 0, w = 0, ecc = 0, color = 'coral', tnight = 40., 
                albedo = albedo, radiancemap = radiancemap, phasecurve = False)

     # Instantiate the system
     system = System(star, b, c, distance = 12, oversample = 10)

     # There's a triple occultation of `c` at this time
     time = np.arange(252.75, 253.50, 5 * MINUTE)

     # Compute and plot the light curve
     system.compute(time, lambda1 = 1, lambda2 = 15)
     
     # Print the temperature for reference
     print("Planet c effective temperature: %.1f K" % c.teff)
     
     # Hackishly remove the other occultations and secondary eclipses
     # so we can plot just the occultation of `c` by `b`
     a = np.argmax(time > 253.013 - 0.025)
     b = np.argmax(time > 253.013 + 0.025)
     system.c.flux[:a] = system.c.flux[0]
     system.b.flux[:a] = system.b.flux[0]
     system.c.flux[b:] = system.c.flux[0]
     system.b.flux[b:] = system.b.flux[0]
     system.c.occultor[:a] = 0
     system.c.occultor[b:] = 0
     system.b.occultor[:a] = 0
     system.b.occultor[b:] = 0
     a = np.argmax(system.time_hr > 253.013 - 0.025)
     b = np.argmax(system.time_hr > 253.013 + 0.025)
     system.c.flux_hr[:a] = system.c.flux_hr[0]
     system.b.flux_hr[:a] = system.b.flux_hr[0]
     system.c.flux_hr[b:] = system.c.flux_hr[0]
     system.b.flux_hr[b:] = system.b.flux_hr[0]
     system.c.occultor_hr[:a] = 0
     system.c.occultor_hr[b:] = 0
     system.b.occultor_hr[:a] = 0
     system.b.occultor_hr[b:] = 0

     # Let's re-center the time array for a prettier x axis
     system.A.time -= 253.013
     system.A.time_hr -= 253.013

     # Get all filters
     filters = jwst.get_spitzer_filter_wheel(warm = True)

     print("Planet c Teff : %.2f K" %system.c.teff)

     # Loop over all MIRI filters
     SNRs = np.zeros(len(filters))
     wls = np.zeros(len(filters))
     for i in range(len(filters)):

             # Observe it (ten exposures)
             np.random.seed(123)
             fig, ax = system.observe(stack = N, filter = filters[i], 
                                      instrument = 'spitzer')

             SNRs[i] = system.filter.lightcurve.event_SNRs[0]
             wls[i] = system.filter.eff_wl
             print("SNR: %.3f" % SNRs[i])

             pl.clf()
             pl.close()

     fig, ax = pl.subplots(figsize=(6,4))
     fig.subplots_adjust(left = 0.2, bottom = 0.2)
     ax.plot(wls, SNRs, "-o", color="k")
     ax.set_xlabel(r"Wavelength [$\mu$m]", fontweight = 'bold', fontsize = 14)
     ax.set_ylabel("SNR", fontweight = 'bold', fontsize = 14)

     wl_filt, dwl_filt, tputs, names = jwst.readin_miri_filters()
     ax.tick_params(labelsize=12)
     return fig, ax

def plot():
    '''
    
    '''
    
    fig1, ax1 = Stacked_bc_all_filters(N = 100000, airless = False, 
                                       albedo = 0.0)
    fig1.suptitle("100,000 stacked occultations of a blackbody at eq. temp.", 
                  fontweight = 'bold', fontsize = 10)
    fig2, ax2 = Stacked_bc_all_filters(N = 10000, airless = True, albedo = 0.0)
    fig2.suptitle("10,000 stacked occultations of airless body at eq. temp.", 
                  fontweight = 'bold', fontsize = 10)
    
    # HACK: By specifying an albedo of -38, we get an effective temperature
    # of 849 K. It's just radiation balance, so this mimics the effect of 
    # greenhouse heating...
    fig3, ax3 = Stacked_bc_all_filters(N = 10, airless = False, albedo = -38)
    fig3.suptitle("10 stacked occultations of 850 K blackbody", 
                  fontweight = 'bold', fontsize = 10)

if __name__ == '__main__':
    
    plot()
    pl.show()
