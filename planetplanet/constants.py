#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
constants.py |github|
---------------------

Constants used throughout the code. All values are in MKS units unless otherwise specified.

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/constants.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

#: Gravitational constant
G = 6.67428e-11
#: Planck's constant
HPLANCK = 6.62607004e-34
#: Speed of light
CLIGHT = 2.998e8
#: Boltzmann constant
KBOLTZ = 1.38064852e-23
#: Stefan-Boltzmann constant
SBOLTZ = 5.670367e-8

#: Solar mass
MSUN = 1.988416e30
#: Solar luminosity
LSUN = 3.846e26
#: Solar radius
RSUN = 6.957e8
#: Astronomical Unit
AUM = 1.49598e11
#: Parsec
PARSEC = 3.086e16
#: Earth mass
MEARTH = 5.9722e24
#: Earth radius
REARTH = 6.3781e6
#: Earth insolation (solar constant)
SEARTH = 1.361e3

#: Seconds in 1 day
DAYSEC = 86400.
#: Days in 1 second
SECOND = 1. / DAYSEC
#: Days in 1 minute
MINUTE = 1. / 1440.

#: Noon UT on October 8, 2016 in BJD − 2,450,000
OCTOBER_08_2016 = 7670.

#: 1 AU in Earth radii
AUREARTH = AUM / REARTH
#: 1 solar mass in Earth masses
MSUNMEARTH = MSUN / MEARTH
#: 1 solar radius in Earth radii
RSUNREARTH = RSUN / REARTH
#: Gravitational constant in Earth units
GEARTH = G * DAYSEC ** 2 * MEARTH / REARTH ** 3

# C definitions
MDFAST = 0
NEWTON = 1
QGSL = 3
REB_INTEGRATOR_WHFAST = 1
REB_INTEGRATOR_IAS15 = 0
MAP_NONE = -1
MAP_RADIAL_DEFAULT = 0
MAP_RADIAL_CUSTOM = 1
MAP_ELLIPTICAL_DEFAULT = 2
MAP_ELLIPTICAL_CUSTOM = 3