#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
theta.py
--------

'''

import numpy as np
import os
import matplotlib.pyplot as pl
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from scipy.optimize import brentq

__all__ = ['sample']

# Load
PATH = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(PATH, "sigma_ecc_rho_theta_tremaine.txt")
sig_theta, rho_star, sig_ecc, P = np.loadtxt(file, unpack = True, delimiter = ',')

# Get grid points
x = np.array(sorted(list(set(sig_theta)))) * 180 / np.pi
y = np.array(sorted(list(set(rho_star))))
z = np.array(sorted(list(set(sig_ecc))))
P = P.reshape(len(x), len(y), len(z))

# Get the marginalized sigma_theta posterior
PDF = np.sum(P[:,:,0], axis = (1))
PDF /= np.max(PDF)

# Compute the cumulative distribution function
CDF = cumtrapz(PDF, x)
CDF /= np.max(CDF)

# Interpolate the CDF
try:
  f = interp1d(x[:-1], CDF, fill_value = 'extrapolate', bounds_error = False, kind = 'cubic')
except ValueError:
  f = interp1d(x[:-1], CDF, fill_value = 'extrapolate', bounds_error = False, kind = 'linear')
  
# Enforce bounds
def CDF(x):
  res = f(np.atleast_1d(x))
  res[res > 1] = 1
  res[res < 0] = 0
  if hasattr(x, '__len__'):
    return res
  else:
    return res[0]

# Draw a sample
def sample():
  y = np.random.random()
  f = lambda x: CDF(x) - y
  while np.sign(f(0)) == np.sign(f(1)):
    y = np.random.random()
    f = lambda x: CDF(x) - y
  return brentq(f, 0, 1)