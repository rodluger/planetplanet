#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
intersection.py |github|
------------------------

Interactive circle/ellipse intersection point finder. Solves a quartic polynomial
to determine the points of intersection. This method is slow because it also finds
imaginary roots, which we don't care about. Eventually :py:class:`planetplanet` will
adopt a faster solver for this.

  .. plot::
     :align: center
     
     from scripts import intersection
     import matplotlib.pyplot as pl
     intersection.Interactor()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/intersection.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
tol = 1e-5

def GetRoots(a, b, xE, yE, xC, yC, r):
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
  roots = [r.real + xC for r in np.roots([c4, c3, c2, c1, c0]) if np.abs(r.imag) < tol]
  return roots

class Interactor(object):
  '''
  
  '''
  
  def __init__(self):
    '''
  
    '''
    
    # Set up the figure
    self.fig, self.ax = pl.subplots(1, figsize = (6, 8))
    self.fig.subplots_adjust(bottom = 0.3)
    self.ax.set_aspect(1)
    self.ax.set_xlim(-2, 2)
    self.ax.set_ylim(-2, 2)
    
    # Plot the circle
    x = np.linspace(-1,1,1000)
    self.ax.plot(x, np.sqrt(1 - x**2), 'k--', lw = 1)
    self.ax.plot(x, -np.sqrt(1 - x**2), 'k--', lw = 1)
    
    # Plot the ellipse
    xE = 0
    yE = 0
    a = 1.5
    b = 0.3
    x = np.linspace(xE - b, xE + b, 1000)
    self.el1, = self.ax.plot(x, yE - a / b * np.sqrt(b ** 2 - (x - xE) ** 2), 'k-', lw = 1)
    self.el2, = self.ax.plot(x, yE + a / b * np.sqrt(b ** 2 - (x - xE) ** 2), 'k-', lw = 1)
    
    # Plot the intersections
    self.pts = [0, 0, 0, 0, 0, 0, 0, 0]
    for n in range(8):
      self.pts[n], = self.ax.plot(-99, -99, 'ro')
    n = 0
    roots = GetRoots(a, b, xE, yE, 0, 0, 1)
    for x in roots:
      yC1 = np.sqrt(1 - x ** 2)
      yC2 = -yC1
      A = b ** 2 - (x - xE) ** 2
      if A < 0:
        A = 0
      else:
        A = a / b * np.sqrt(A)
      yE1 = yE - A
      yE2 = yE + A
      if np.abs(yC1 - yE1) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC1)
        n += 1
      if np.abs(yC1 - yE2) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC1)
        n += 1
      if np.abs(yC2 - yE1) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC2)
        n += 1
      if np.abs(yC2 - yE2) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC2)
        n += 1
      
    # Sliders
    self.axa = pl.axes([0.125, 0.05, 0.775, 0.03])
    self.sla = Slider(self.axa, r'$a$', 0.01, 2, valinit = a)
    self.axb = pl.axes([0.125, 0.1, 0.775, 0.03])
    self.slb = Slider(self.axb, r'$b$', 0.01, 2, valinit = b)
    self.axx = pl.axes([0.125, 0.15, 0.775, 0.03])
    self.slx = Slider(self.axx, r'$x$', -1, 1, valinit = xE)
    self.axy = pl.axes([0.125, 0.2, 0.775, 0.03])
    self.sly = Slider(self.axy, r'$y$', -1, 1, valinit = yE)
    
    # Init
    self.sla.on_changed(self.update)
    self.slb.on_changed(self.update)
    self.slx.on_changed(self.update)
    self.sly.on_changed(self.update)
    self.update(None)
    
    pl.show()
  
  def update(self, val):
    '''
  
    '''
  
    # Get the values
    a = self.sla.val
    b = self.slb.val
    xE = self.slx.val
    yE = self.sly.val
  
    # Update the ellipse
    x = np.linspace(xE - b, xE + b, 1000)
    A = b ** 2 - (x - xE) ** 2
    A[A < 0] = 0
    A = np.sqrt(A)
    self.el1.set_xdata(x)
    self.el1.set_ydata(yE - a / b * A)
    self.el2.set_xdata(x)
    self.el2.set_ydata(yE + a / b * A)
    
    # Get the roots
    for n in range(8):
      self.pts[n].set_xdata(-99)
    n = 0
    roots = GetRoots(a, b, xE, yE, 0, 0, 1)
    for x in roots:
      yC1 = np.sqrt(1 - x ** 2)
      yC2 = -yC1
      A = b ** 2 - (x - xE) ** 2
      if A < 0:
        A = 0
      else:
        A = a / b * np.sqrt(A)
      yE1 = yE - A
      yE2 = yE + A
      if np.abs(yC1 - yE1) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC1)
        n += 1
      if np.abs(yC1 - yE2) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC1)
        n += 1
      if np.abs(yC2 - yE1) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC2)
        n += 1
      if np.abs(yC2 - yE2) < tol:
        self.pts[n].set_xdata(x)
        self.pts[n].set_ydata(yC2)
        n += 1

if __name__ == '__main__':       
  Interactor()
  pl.show()