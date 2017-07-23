#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py
----------

Interactive eyeball planet visualizer.

'''

import numpy as np
import matplotlib.pyplot as pl
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.widgets import Slider

def ZenithColor(z):
  '''
  Color of a given zenith angle slice
  
  '''
  
  return pl.get_cmap('inferno')(0.3 + 0.3 * (np.cos(z) + 1))

def PlotEyeball(theta = np.pi / 3, nz = 11, dphi = 0, dlam = 0, fig = None):
  '''

  '''
  
  # The coordinates of the substellar point
  x_ss = -np.cos(theta + dphi) * np.cos(dlam)
  y_ss = np.sin(dlam)
  
  # The rotation angle that makes the planet symmetric about the x-axis
  alpha = np.arctan2(y_ss, -x_ss)
  
  # Compute the new effective theta 
  if theta + dphi < -np.pi:
    theta = 2 * np.pi - np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  elif theta + dphi < 0:
    theta = -np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  elif theta + dphi > np.pi:
    theta = -np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  else:
    theta = np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  
  # Set up the floating axis
  if fig is None:
    fig = pl.figure(figsize = (6,6))
  tr = Affine2D().rotate_deg(-alpha * 180 / np.pi)
  x = 1. / (np.abs(np.cos(alpha)) + np.abs(np.sin(alpha)))
  grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=(-x, x, -x, x), tick_formatter1 = None, tick_formatter2 = None)
  ax_orig = floating_axes.FloatingSubplot(fig, 111, grid_helper = grid_helper)
  ax_orig.axis["bottom"].set_visible(False)
  ax_orig.axis["top"].set_visible(False)
  ax_orig.axis["left"].set_visible(False)
  ax_orig.axis["right"].set_visible(False)
  fig.add_subplot(ax_orig)
  ax = ax_orig.get_aux_axes(tr)
  
  # Plot the occulted body
  r = 1
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x, y, color = 'k', zorder = 98, lw = 1)
  ax.plot(x, -y, color = 'k', zorder = 98, lw = 1)

  # Plot the zenith angle ellipses
  zarr = np.linspace(0, np.pi, nz + 2)
  for i, z in enumerate(zarr[1:]):

    # The ellipse
    a = r * np.abs(np.sin(z))
    b = max(0.001, a * np.abs(np.sin(theta)))
    xE = -r * np.cos(z) * np.cos(theta)
    xlimb = r * np.cos(z) * np.sin(theta) * np.tan(theta)
    if ((theta > 0) and (b < xlimb)) or ((theta <= 0) and (b > xlimb)):
      xmin = xE - b
    else:
      xmin = xE - xlimb
    if ((theta > 0) and (b > -xlimb)) or ((theta <= 0) and (b < -xlimb)):
      xmax = xE + b
    else:
      xmax = xE - xlimb
        
    # Plot it
    x = np.linspace(xE - b, xE + b, 1000)
    if theta > 0:
      x[x < xE - xlimb] = np.nan
    elif theta > -np.pi / 2:
      x[x > xE - xlimb] = np.nan
    A = b ** 2 - (x - xE) ** 2
    A[A < 0] = 0
    y = (a / b) * np.sqrt(A)
    if np.abs(np.cos(z)) < 1e-5:
      style = dict(color = 'k', ls = '--', lw = 0.5, zorder = 9999)
    else:
      style = dict(color = 'k', ls = '-', lw = 0.5, zorder = 9999)
    ax.plot(x, y, **style)
    ax.plot(x, -y, **style)
    
    # Fill the ellipses
    if theta < 0:
      ax.fill_between(x, -y, y, color = ZenithColor(z), zorder = int(100 * z))
    else:
      ax.fill_between(x, -y, y, color = ZenithColor(z), zorder = int(-100 * z))
  
    # Fill the ellipses that are cut off by the limb
    if theta < 0:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x_, -y_, y_, color = ZenithColor(zarr[i]), zorder = int(-100 * z))
    else:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x_, -y_, y_, color = ZenithColor(z), zorder = int(-100 * z))
    
  return fig, ax_orig, ax

class Interact(object):
  '''
  
  '''
  
  def __init__(self):
    '''
  
    '''
  
    self.fig = pl.figure(figsize = (6,6))
    self.fig.subplots_adjust(bottom = 0.25)
    self.axtheta = pl.axes([0.3, 0.05, 0.44, 0.03])
    self.theta = Slider(self.axtheta, r'$\theta$', -180., 180., valinit = 90.)
    self.axpsi = pl.axes([0.3, 0.1, 0.44, 0.03])
    self.psi = Slider(self.axpsi, r'$\psi$', -90, 90., valinit = 0.)
    self.axlam = pl.axes([0.3, 0.15, 0.44, 0.03])
    self.lam = Slider(self.axlam, r'$\lambda$', -90, 90., valinit = 0.)
    self.theta.on_changed(self.update)
    self.psi.on_changed(self.update)
    self.lam.on_changed(self.update)
    self.update(90.)
    pl.show()
  
  def update(self, val):
    '''
  
    '''
  
    # Remove all axes except sliders
    for ax in self.fig.get_axes():
      if ax not in [self.axtheta, self.axpsi, self.axlam]:
        ax.remove()
  
    # Plot the planet
    PlotEyeball(fig = self.fig, theta = self.theta.val * np.pi / 180, dphi = self.psi.val * np.pi / 180, dlam = self.lam.val * np.pi / 180)

if __name__ == '__main__':
  Interact()