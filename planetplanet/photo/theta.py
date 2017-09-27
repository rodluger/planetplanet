#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
theta.py |github|
-----------------

This module enables sampling from the :math:`\\theta` distribution for the
TRAPPIST-1 planets derived using the Monte Carlo method described in the 
paper. :math:`\\theta` is the polar angle of the angular momentum vector; we 
assume all planets in the TRAPPIST-1 system have angular momentum vectors with 
polar angle drawn from this distribution.

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/theta.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

import numpy as np
import os
import matplotlib.pyplot as pl
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d
from scipy.optimize import brentq

__all__ = ['sample', 'CDF']

# Load
PATH = os.path.dirname(os.path.abspath(__file__))
file = os.path.join(PATH, "posteriors.dat")
sig_theta, rho_star, sig_ecc, P = np.loadtxt(file, unpack = True, 
                                             delimiter = ',')

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
    f = interp1d(x[:-1], CDF, fill_value = 'extrapolate', 
                 bounds_error = False, kind = 'cubic')
except ValueError:
    f = interp1d(x[:-1], CDF, fill_value = 'extrapolate', 
                 bounds_error = False, kind = 'linear')
  
# Enforce bounds
def CDF(theta):
    '''
    Returns the value of the cumulative distribution function (CDF) of the 
    polar angle of the angular momentum vector, :math:`\\theta`.
    
    .. plot::
         :align: center
         
         from planetplanet.photo import theta
         import matplotlib.pyplot as pl
         x = np.linspace(0, 1, 1000)
         y = [theta.CDF(xi) for xi in x]
         pl.plot(x, y)
         pl.xlabel(r'$\\theta$ [deg]', fontweight = 'bold')
         pl.ylabel('Cumulative Probability', fontweight = 'bold')
         pl.show()
    
    '''
    
    res = f(np.atleast_1d(theta))
    res[res > 1] = 1
    res[res < 0] = 0
    if hasattr(theta, '__len__'):
        return res
    else:
        return res[0]

# Draw a sample
def sample():
    '''
    Draw a sample from the distribution of polar angle of the angular 
    momentum vector, :math:`\\theta`, computed using the Monte Carlo 
    technique discussed in the paper. 
    
    .. plot::
         :align: center
         
         from planetplanet.photo import theta
         import matplotlib.pyplot as pl
         x = [theta.sample() for i in range(10000)]
         pl.hist(x, bins = 50)
         pl.xlabel(r'$\\theta$ [deg]', fontweight = 'bold')
         pl.ylabel('Probability', fontweight = 'bold')
         pl.show()
         
    '''
    
    y = np.random.random()
    f = lambda x: CDF(x) - y
    while np.sign(f(0)) == np.sign(f(1)):
        y = np.random.random()
        f = lambda x: CDF(x) - y
    return brentq(f, 0, 1)