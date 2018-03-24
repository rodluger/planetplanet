#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
trappist1.py |github|
---------------------

This module hosts TRAPPIST-1-specific routines.

  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/trappist1.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
    unicode_literals
from ..constants import *
from .ppo import Star, Planet, System
from . import theta
import numpy as np
import matplotlib.pyplot as pl
import os
from tqdm import tqdm
PATH = os.path.dirname(os.path.abspath(__file__))


__all__ = ['Trappist1', 'TRAPPIST_T0']


#: Time at which the mean longitudes are computed (BJD - 2,450,000)
#: This is ~ August 23 2015
TRAPPIST_T0 = 7257.93115525


def N(mu, sigma):
    """Draw a random variable from a normal distribution."""
    return mu + sigma * np.random.randn()


def Trappist1(distance=12, seed=None, flat=False, **kwargs):
    '''
    Returns an instance of :py:obj:`planetplanet.photo.System` for the full
    TRAPPIST-1 system. Star and planet parameters are drawn from their
    respective prior distributions, which are based on the observed values
    from Gillon et al. (2017), Luger et al. (2017),
    Burgasser & Mamajek (2017), and Grimm et al. (2018). Longitudes
    of ascending node are
    drawn from the :math:`\\theta` distribution derived in our paper.

    :param float distance: Distance to the system in parsecs. \
           Default :py:obj:`12`
    :param int seed: Random number generator seed. Default :py:obj:`None`
    :param kwargs: Any other :py:obj:`kwargs` to be passed to \
           :py:func:`planetplanet.Star`, \
           :py:func:`planetplanet.Planet`, and :py:func:`planetplanet.System`.

    .. plot::
         :align: center

         from planetplanet.photo.trappist1 import Trappist1
         from planetplanet.constants import MINUTE
         import matplotlib.pyplot as pl
         import numpy as np
         system = Trappist1()
         system.compute(np.arange(0, 10, 1 * MINUTE))
         system.plot_lightcurve()
         pl.show()

    '''

    # Randomizer seed
    if seed is not None:
        np.random.seed(seed)

    # Instantiate the star
    # Mass used in Grimm et al. (2018) for consistency
    mstar = 0.09
    # Radius from Burgasser & Mamajek (2017)
    rstar = N(0.121, 0.003)
    teff = (N(0.000524, 0.000034)
            * LSUN / (4 * np.pi * (rstar * RSUN) ** 2 * SBOLTZ)) ** 0.25
    star = Star('A', m=mstar, r=rstar, teff=teff, color='k', **kwargs)

    # Instantiate the planets
    planets = [None for i in range(7)]
    names = ['b', 'c', 'd', 'e', 'f', 'g', 'h']
    colors = ['firebrick', 'coral', 'gold', 'mediumseagreen', 'turquoise',
              'cornflowerblue', 'midnightblue']

    # These we're just going to fix for now. We have no prior
    # constraints on them. Let's assume the most optimistic albedos.
    albedos = [(0., 0), (0., 0), (0., 0), (0., 0),
               (0., 0), (0., 0), (0., 0)]
    tnights = [(40., 0), (40., 0), (40., 0), (40., 0),
               (40., 0), (40., 0), (40., 0)]

    # Periods from Gillon et al. (2017) and Luger et al. (2017)
    # TODO: Inconsistent!
    periods = [(1.51087081, 0.60e-6),
               (2.4218233, 0.17e-5),
               (4.049610, 0.63e-4),
               (6.099615, 0.11e-4),
               (9.206690, 0.15e-4),
               (12.35294, 0.12e-3),
               (18.767, 0.004)]

    # Depths from Gillon et al. (2017) and Luger et al. (2017)
    depths = [(0.7266, 0.0088),
              (0.687, 0.010),
              (0.367, 0.017),
              (0.519, 0.026),
              (0.673, 0.023),
              (0.782, 0.027),
              (0.346, 0.032)]

    # Compute the polar angle scatter from our MC posterior
    sig_theta = theta.sample()

    # Draw from the joint inclination distribution (Gillon et al. 2017)
    i = np.random.randint(159999)
    with open(os.path.join(PATH, "inclination.dat")) as file:
        for j, line in enumerate(file):
            if j == i:
                inclinations = np.array([float(v) for v in line.split()])
                break

    # Draw from the other distributions (Grimm et al. 2018)
    i = np.random.randint(99998)
    with open(os.path.join(PATH, "grimm.dat")) as file:
        for j, line in enumerate(file):
            if j == i:
                sample = np.array([float(v) for v in line.split()])
                sample = sample[:-1].reshape(7, 5)

    # Instantiate the planets
    for i in range(7):

        # Longitude of ascending node in degrees (from our MC posterior)
        if (i == 0):
            Omega = 0
        else:
            Omega = N(0, sig_theta)

        # Draws from posteriors
        inc = inclinations[i]

        # Perfectly flat system?
        if flat:
            inc = 90
            Omega = 0

        m = sample[i][0] * MSUN / MEARTH
        a = sample[i][1]

        # DEBUG
        """
        if i == 1:
            a *= (1 + 3.75e-5)
        elif i == 2:
            a *= (1 + 8.e-5)
        elif i == 3:
            a *= (1 + 9.2e-5)
        elif i == 4:
            a *= (1 + 1.2e-4)
        elif i == 5:
            a *= (1 + 1.5e-4)
        elif i == 6:
            a *= (1 + 1.92e-4)
        """

        ecc = sample[i][2]
        w = sample[i][3] * 180 / np.pi
        lambda0 = sample[i][4] * 180 / np.pi + Omega
        per = 2 * np.pi * np.sqrt((a * AUM) ** 3 /
                                  (G * (mstar * MSUN + m * MEARTH))) / DAYSEC

        # Radius from Rp / Rstar
        mu = np.sqrt(depths[i][0] / 100)
        sig = 0.5 * depths[i][1] / 100 / mu
        RpRs = N(mu, sig)
        r = RpRs * rstar * RSUN / REARTH

        # Albedo, night side temperature, effective temperature
        albedo = N(*albedos[i])
        tnight = N(*tnights[i])

        # Instantiate!
        planets[i] = Planet(names[i], m=m, per=per, inc=inc, r=r,
                            lambda0=lambda0, Omega=Omega, w=w, ecc=ecc,
                            color=colors[i], tnight=tnight,
                            albedo=albedo, time0=TRAPPIST_T0,
                            **kwargs)

    # Return the system
    system = System(star, distance=distance, *planets, **kwargs)
    return system
