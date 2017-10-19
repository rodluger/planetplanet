#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
layercake.py |github|
---------------------

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/layercake.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
from scipy.optimize import brentq, minimize_scalar

def _test():
    '''
    
    '''
    
    pass

def Integral(x, a, b, xE, yE, upper = True):
    '''
    
    '''
    
    sgn = 1 if upper else -1
    z = np.sqrt((b + x - xE) * (b - x + xE))
    foo = z * (x - xE) + b ** 2 * np.arctan((x - xE) / z)
    return yE * x + sgn * (a / (2 * b)) * foo

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

def Ellipse(x, a, b, xE, yE, xlimb = None, day = True):
    '''
    
    '''
    
    x = np.array(x)
    if xlimb is not None:
        if day:
            x[x < xE - xlimb] = np.nan
        else:
            x[x > xE - xlimb] = np.nan
    A = b ** 2 - (x - xE) ** 2
    if hasattr(A, '__len__'):
        A[np.isnan(A)] = 0
        A[A < 0] = 0
        A = np.sqrt(A)
    else:
        if np.isnan(A) or A < 0:
            A = 0
        else:
            A = np.sqrt(A)
    y = yE + a / b * A
    if hasattr(y, '__len__'):
        y[np.isnan(x)] = np.nan
    else:
        if np.isnan(x):
            y = np.nan
    return y

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

def Compute(theta = np.pi / 4, nz = 7, ro = 0.75, xo = -0.75, yo = 0.5, tol = 1e-5):
    '''
    
    '''
    
    # Get the C-C roots
    ccroots = []
    A0 = 0
    d = np.sqrt(xo ** 2 + yo ** 2)
    if (d < (1 + ro)):
        A = (-d + 1 - ro) * (-d - 1 + ro) * (-d + 1 + ro) * (d + 1 + ro)
        if (A >= 0):
            y = np.sqrt(A) / (2 * d)
            x = np.sqrt(1 - y ** 2)
            frac = yo / xo
            cost = 1 / np.sqrt(frac ** 2 + 1)
            sint = frac * cost
            if (xo < 0):
                cost *= -1
                sint *= -1
            x1 = x * cost + y * sint
            x2 = x * cost - y * sint
            if x1 > x2:
                x1, x2 = x2, x1
            if np.abs(Circle(x1, 1, 0, 0) - Circle(x1, ro, xo, yo)) < tol:
                root1 = (x1, Circle(x1, 1, 0, 0), np.arccos(-x1 * np.cos(theta)))
            else:
                root1 = (x1, Circle(x1, -1, 0, 0), np.arccos(-x1 * np.cos(theta)))
            if np.abs(Circle(x2, 1, 0, 0) - Circle(x2, ro, xo, yo)) < tol:
                root2 = (x2, Circle(x2, 1, 0, 0), np.arccos(-x2 * np.cos(theta)))
            else:
                root1 = (x2, Circle(x2, -1, 0, 0), np.arccos(-x2 * np.cos(theta)))
            ccroots = [root1, root2]
        
        # Total area of intersection
        if d + ro <= 1:
            A0 = np.pi * ro ** 2
        elif d + 1 <= ro:
            A0 = np.pi
        else:
            A0 = ro ** 2 * np.arccos((d ** 2 + ro ** 2 - 1) / (2 * d * ro)) \
               + np.arccos((d ** 2 + 1 - ro ** 2) / (2 * d)) \
               - 0.5 * np.sqrt(A)
            
    # Compute the ellipses (in reverse order)
    ellipses = []
    ceroots = []
    zarr = np.linspace(0, np.pi, nz + 2)[1:][::-1]
    for i, z in enumerate(zarr):

        # The ellipse
        a = np.abs(np.sin(z))
        b = max(0.001, a * np.abs(np.sin(theta)))
        xE = -np.cos(z) * np.cos(theta)
        yE = 0
        xlimb = np.cos(z) * np.sin(theta) * np.tan(theta)
        if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
            xmin = xE - b
        else:
            xmin = xE - xlimb
        if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
            xmax = xE + b
        else:
            xmax = xE - xlimb
        
        # Append to list
        ellipses.append((z, a, b, xE, yE, xlimb))
        
        # Get the C-E roots
        if xmin < xmax:
            roots = GetRoots(a, b, xE, yE, xo, yo, ro, tol = tol)
            for x in roots:
                yC1 = Circle(x, ro, xo, yo)
                yC2 = Circle(x, -ro, xo, yo)
                yE1 = Ellipse(x, a, b, xE, yE, xlimb, theta > 0)
                yE2 = Ellipse(x, -a, b, xE, yE, xlimb, theta > 0)
                if (np.abs(yC1 - yE1) < tol) or (np.abs(yC1 - yE2) < tol):
                    y = yC1
                    root = (x, y, z, a, b, xE, yE, xmin, xmax)
                    ceroots.append(root)
                elif (np.abs(yC2 - yE1) < tol) or (np.abs(yC2 - yE2) < tol):
                    y = yC2
                    root = (x, y, z, a, b, xE, yE, xmin, xmax)
                    ceroots.append(root)   
    
    # Compute the area
    print(A0)
           
    # Set up the figure
    fig, ax = pl.subplots(1, figsize = (6, 6))
    ax.set_aspect(1)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    
    # Plot the planet
    x = np.linspace(-1,1,1000)
    ax.plot(x, Circle(x, 1, 0, 0), 'k-', lw = 1)
    ax.plot(x, Circle(x, -1, 0, 0), 'k-', lw = 1)
    
    # Plot the occultor
    x = np.linspace(xo - ro, xo + ro, 1000)
    ax.plot(x, Circle(x, ro, xo, yo), 'r-', lw = 1)
    ax.plot(x, Circle(x, -ro, xo, yo), 'r-', lw = 1)     
    
    # Plot the ellipses
    for z, a, b, xE, yE, xlimb in ellipses:
        x = np.linspace(xE - b, xE + b, 1000)
        ax.plot(x, Ellipse(x, a, b, xE, yE, xlimb, theta > 0), color = 'k', ls = '-', lw = 0.5)
        ax.plot(x, Ellipse(x, -a, b, xE, yE, xlimb, theta > 0), color = 'k', ls = '-', lw = 0.5)
    
    # Plot the planet-occultor (C-C) intersections
    for x, y, z in ccroots:
        pl.plot(x, y, 'bo')
    
    # Plot the occultor-ellipse (C-E) intersections
    for x, y, *_ in ceroots:
        pl.plot(x, y, 'ro')
    
Compute()
pl.show()