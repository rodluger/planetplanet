#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
aspect.py |github|
------------------

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/aspect.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
from scipy.optimize import brentq, minimize_scalar

def _test():
    '''
    
    '''
    
    pass

def Circle(x, r, xC, yC):
    '''
    
    '''
    
    A = r ** 2 - (x - xC) ** 2
    if hasattr(A, '__len__'):
        A[A < 0] = 0
        A = np.sqrt(A)
    else:
        if A < 0:
            A = 0
        else:
            A = np.sqrt(A)
    return yC + np.sign(r) * A

def GetRoots(a, b, xE, yE, xC, yC, r, tol = 1e-5):
    '''
    
    '''
    
    # Define some stuff
    r2 = r * r;
    a2 = a * a;
    b2 = b * b;
    a2b2 = a2 / b2;
    x0 = xE - xC;
    y0 = yE - yC;
    y2 = y0 * y0;
    x2 = x0 * x0;
    
    # Get the coefficients
    A = a2b2 - 1.;
    B = -2. * x0 * a2b2;
    C = r2 - y2 - a2 + a2b2 * x2;
    D = 4. * y2 * a2b2;
    c4 = A * A;
    c3 = 2. * A * B;
    c2 = 2. * A * C + B * B + D;
    c1 = 2. * B * C - 2. * D * x0;
    c0 = C * C - (b2 - x2) * D;
        
    # Get the real roots
    roots = [r.real + xC for r in np.roots([c4, c3, c2, c1, c0]) 
             if np.abs(r.imag) < tol]
    return roots

def Plot(theta = np.pi / 2, nz = 9, aspect = 1., ro = 0.75, xo = -0.75, yo = 0.5, tol = 1e-5):
    '''
    
    '''
    
    # Set up the figure
    fig, ax = pl.subplots(1, figsize = (6, 6))
    ax.set_aspect(1)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    
    # Plot the planet
    x = np.linspace(-1,1,1000)
    ax.plot(x, Circle(x, 1, 0, 0), 'k-', lw = 1)
    ax.plot(x, Circle(x, -1, 0, 0), 'k-', lw = 1)
    
    # Compute the ellipses (in reverse order)
    ellipses = []
    ceroots = []
    zarr = np.linspace(0, np.pi, nz + 2)[1:][::-1]
    for i, z in enumerate(zarr):

        # The ellipse
        a = np.abs(np.sin(z))
        b = max(0.001, a * np.abs(np.sin(theta))) * aspect
        xE = -np.cos(z) * np.cos(theta)
        yE = 0
        
        # Get the bounds
        A = 1 - a ** 2 / b ** 2
        B = 2 * a ** 2 / b ** 2 * xE
        C = a ** 2 - 1 - a ** 2 / b ** 2 * xE
        foo = B ** 2 - 4 * A * C
        print(a, b, foo)
        xmin = (-B - np.sqrt(foo)) / (2 * A)
        xmax = (-B + np.sqrt(foo)) / (2 * A)
        
        
        # Plot
        x = np.linspace(xmin, xmax, 10000)
        A = b ** 2 - (x - xE) ** 2
        pl.plot(x, yE + a / b * np.sqrt(A), 'k-', lw = 0.5)
        pl.plot(x, yE - a / b * np.sqrt(A), 'k-', lw = 0.5)
        
    
Plot()
pl.show()