#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
:py:mod:`wrapper.py` - C wrapper
--------------------------------

'''

from __future__ import division, print_function, absolute_import, unicode_literals
import ctypes
import numpy as np
import os
from numpy.ctypeslib import ndpointer, as_ctypes
AUREARTH = 23454.9271
MDFAST = 0
NEWTON = 1

class PLANET(ctypes.Structure):
  '''
  The class containing all the input planet parameters
  
  '''
  
  _fields_ = [("per", ctypes.c_double),
              ("inc", ctypes.c_double),
              ("ecc", ctypes.c_double),
              ("w", ctypes.c_double),
              ("a", ctypes.c_double),
              ("t0", ctypes.c_double),
              ("r", ctypes.c_double),
              ("noon", ctypes.c_double),
              ("midnight", ctypes.c_double),
              ("nlat", ctypes.c_int),
              ("x", ctypes.c_double),
              ("y", ctypes.c_double),
              ("z", ctypes.c_double)]
        
  def __init__(self, **kwargs):
    self.per = kwargs.pop('per', 1.51087081)
    inc = kwargs.pop('inc', 89.65)
    self.inc = inc * np.pi / 180.
    self.ecc = kwargs.pop('ecc', 0.)
    self.w = kwargs.pop('w', 0.)
    self.a = kwargs.pop('a', 0.01111 * AUREARTH)
    self.t0 = kwargs.pop('t0', 7322.51736)
    self.r = kwargs.pop('r', 1.086)
    self.noon = kwargs.pop('noon', 1.)
    self.midnight = kwargs.pop('midnight', 0.1)
    self.nlat = kwargs.pop('nlat', 11)
    
class SETTINGS(ctypes.Structure):
  '''
  The class that contains the model settings
  
  '''
  _fields_ = [("keptol", ctypes.c_double),
              ("maxkepiter", ctypes.c_int),
              ("kepsolver", ctypes.c_int)]
  
  def __init__(self, **kwargs):
    self.keptol = kwargs.pop('keptol', 1.e-15)
    self.maxkepiter = kwargs.pop('maxkepiter', 100)
    self.kepsolver = kwargs.pop('kepsolver', NEWTON)

# Load the C library
try:
  lib = ctypes.CDLL(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ppo.so'))
except:
  raise Exception("Can't find `ppo.so`; please run `make` to compile it.")

# Declare the C functions; user should access these through the Orbit() class below
Flux = lib.Flux
Flux.restype = ctypes.c_double
Flux.argtypes = [ctypes.c_double, ctypes.POINTER(PLANET), ctypes.POINTER(PLANET), ctypes.POINTER(SETTINGS)]

settings = SETTINGS()
t1b = PLANET(per = 1.51087081, inc = 89.65, a = 0.01111 * AUREARTH, r = 1.086, t0 = 7322.51736, midnight = 0, n = 31)
t1c = PLANET(per = 2.4218233, inc = 89.67, a = 0.01521 * AUREARTH, r = 1.056, t0 = 7282.80728, midnight = 0, n = 31)

time = np.linspace(0, 10, 100000)
flux = [Flux(t, t1b, t1c, settings) for t in time]

import matplotlib.pyplot as pl
pl.plot(time, flux)
pl.show()