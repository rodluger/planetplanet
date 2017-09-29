#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
wasp47.py |github|
------------------

This module hosts WASP-47-specific routines.

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/wasp47.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from ..constants import *
from .ppo import Star, Planet, System
from .eyeball import LimbDarkenedMap, RadiativeEquilibriumMap
import numpy as np
import matplotlib.pyplot as pl
import os
from tqdm import tqdm

__all__ = ['Wasp47']

def Wasp47(sample = True, distance = 200, seed = None, **kwargs):
    '''
    Returns an instance of :py:obj:`planetplanet.photo.System` for the full 
    WASP-47 system. Star and planet parameters are drawn from their 
    respective prior distributions.
    
    :param bool sample: Draw a random sample from the full prior? \
           If :py:obj:`False`,returns the mean values for all parameters. \
           Default :py:obj:`True`
    :param float distance: Distance to the system in parsecs. \
           Default :py:obj:`200`
    :param int seed: Random number generator seed. Default :py:obj:`None`
    :param kwargs: Any other :py:obj:`kwargs` to be passed to \
           :py:func:`planetplanet.Star`, \
           :py:func:`planetplanet.Planet`, and :py:func:`planetplanet.System`.
    
    .. plot::
         :align: center
         
         from planetplanet.photo.wasp47 import Wasp47
         from planetplanet.constants import MINUTE
         import matplotlib.pyplot as pl
         import numpy as np
         system = Wasp47()
         system.compute(np.arange(0, 10, 1 * MINUTE))
         system.plot_lightcurve()
         pl.show()
    
    '''
    
    # Randomizer seed
    if seed is not None:
        np.random.seed(seed)
    
    # Account for the uncertainty?
    if not sample:
        N = lambda mu, sigma: mu
    else: 
        N = lambda mu, sigma: mu + sigma * np.random.randn()

    # Instantiate the star; Weiss et al. (2017)
    mstar = N(1.00, 0.05)
    rstar = N(1.12, 0.02)
    teff = N(5576, 67) # Becker et al. (2015)
    star = Star('A', m = mstar, r = rstar, teff = teff, color = 'k', **kwargs)
    
    # Parameters from Weiss et al. (2017)
    planets = [None for i in range(3)]
    names = ['e', 'b', 'd']
    
    periods = [(0.78961, 0.00001),
               (4.15912, 0.00001), 
               (9.0304, 0.0003)]
    
    # Transit times, KJD = BJD – 2454833.0
    transits = [(2146.7641, 0.0007), 
                (2149.9785, 0.0001), 
                (2155.308, 0.001)]
    
    masses = [(9.1, 1.0), 
              (358, 12), 
              (13.6, 2.0)]
                                                            
    rprs = [(0.01439, 0.00016), 
            (0.10193, 0.00018), 
            (0.02931, 0.00015)]
    
    eccentricities = [(0.03, 0.02), 
                      (0.0028, 0.0028), 
                      (0.007, 0.007)]
    
    ws = [(81.0, 146.0), 
          (51.0, 82.0), 
          (76.0, 106.0)]
    
    inclinations = [(87.0, 3.1), 
                    (89.03, 0.27), 
                    (89.36, 0.67)]
    
    # These we're just going to fix for now. We have no prior 
    # constraints on them. Let's assume the most optimistic albedos.
    albedos = [(0., 0), (0., 0), (0., 0)]
    tnights = [(100., 0), (100., 0), (100., 0)]
    Omegas = [(0., 0.5), 
              (0., 0.5), 
              (0., 0.5)]
    
    # Colors for plotting
    colors = ['firebrick', 'gold', 'cornflowerblue']
    
    # Radiance maps
    maps = [RadiativeEquilibriumMap(), LimbDarkenedMap(), LimbDarkenedMap()]
    
    # Instantiate the planets
    for i in range(3):
    
        # Period and time of transit
        per = N(*periods[i])
        t0 = N(*transits[i])
    
        # Positive mass
        m = 0
        while m <= 0:
            m = N(*masses[i])
    
        # Inclination in range [0, 90]
        inc = N(*inclinations[i])
        if inc > 90:
            inc = 180 - inc
        
        # Longitude of ascending node in degrees
        Omega = N(*Omegas[i])
        
        # Longitude of pericenter in degrees
        w = N(*ws[i])
        
        # Eccentricity
        ecc = 1
        while (ecc < 0) or (ecc >= 1):
            ecc = N(*eccentricities[i])

        # Radius from Rp / Rstar
        r = N(*rprs[i]) * rstar * RSUN / REARTH
    
        # Albedo, night side temperature, effective temperature
        albedo = N(*albedos[i])
        tnight = N(*tnights[i])
    
        # Instantiate!
        planets[i] = Planet(names[i], m = m, per = per, inc = inc, r = r, 
                            t0 = t0, Omega = Omega, w = w, ecc = ecc, 
                            color = colors[i], tnight = tnight, 
                            albedo = albedo, radiancemap = maps[i], **kwargs)

    # Return the system
    system = System(star, distance = distance, *planets, **kwargs)
    return system