#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py
----------

Interactive eyeball planet visualizer.

'''

import numpy as np
np.seterr(invalid = 'ignore')
import matplotlib.pyplot as pl
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.widgets import Slider

def ZenithColor(z):
  '''
  Color of a given zenith angle slice
  
  '''
  
  return pl.get_cmap('inferno')(0.3 + 0.3 * (np.cos(z) + 1))

def DrawEyeball(x0 = 0, y0 = 0, r = 1, theta = np.pi / 3, nz = 11, dphi = 0, dlam = 0, 
                occultors = [], fig = None):
  '''

  '''
  
  # The coordinates of the substellar point
  x_ss = -np.cos(theta + dphi) * np.cos(dlam)
  y_ss = np.sin(dlam)
  
  # The rotation angle that makes the planet symmetric about the x-axis
  alpha = np.arctan2(y_ss, -x_ss)
  
  # Compute the new effective theta 
  if theta + dphi < -np.pi:
    theta_eff = 2 * np.pi - np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  elif theta + dphi < 0:
    theta_eff = -np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  elif theta + dphi > np.pi:
    theta_eff = -np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  else:
    theta_eff = np.arccos(-x_ss * np.cos(alpha) + y_ss * np.sin(alpha))
  
  # Set up the floating axis
  if fig is None:
    fig = pl.figure(figsize = (6,6))
  tr = Affine2D().rotate_deg(-alpha * 180 / np.pi)
  x = 1. / (np.abs(np.cos(alpha)) + np.abs(np.sin(alpha)))
  grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=(-x, x, -x, x))
  ax_orig = floating_axes.FloatingSubplot(fig, 111, grid_helper = grid_helper)
  ax_orig.axis["bottom"].set_visible(False)
  ax_orig.axis["top"].set_visible(False)
  ax_orig.axis["left"].set_visible(False)
  ax_orig.axis["right"].set_visible(False)
  fig.add_subplot(ax_orig)
  ax = ax_orig.get_aux_axes(tr)
  
  # Plot the occultors. Note that we need to transform
  # their position vectors since we're in a rotated frame.
  for occultor in occultors:
    xo = occultor['x']
    yo = occultor['y']
    ro = occultor['r']
    zo = occultor.get('zorder', 1)
    co = occultor.get('color', 'gray')
    xo, yo = xo * np.cos(alpha) - yo * np.sin(alpha), \
             yo * np.cos(alpha) + xo * np.sin(alpha)
    x = np.linspace(-ro, ro, 1000)
    y = np.sqrt(ro ** 2 - x ** 2)
    ax.plot(xo + x, yo + y, color = 'k', zorder = zo, lw = 1)
    ax.plot(xo + x, yo - y, color = 'k', zorder = zo, lw = 1)
    ax.fill_between(xo + x, yo - y, yo + y, color = co, zorder = zo - 0.01, lw = 1)
    
  # Plot the occulted body
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x0 + x, y0 + y, color = 'k', zorder = 0, lw = 1)
  ax.plot(x0 + x, y0 - y, color = 'k', zorder = 0, lw = 1)

  # Plot the zenith angle ellipses
  zarr = np.linspace(0, np.pi, nz + 2)
  for i, z in enumerate(zarr[1:]):

    # The ellipse
    a = r * np.abs(np.sin(z))
    b = max(0.001, a * np.abs(np.sin(theta_eff)))
    xE = -r * np.cos(z) * np.cos(theta_eff)
    xlimb = r * np.cos(z) * np.sin(theta_eff) * np.tan(theta_eff)
    if ((theta_eff > 0) and (b < xlimb)) or ((theta_eff <= 0) and (b > xlimb)):
      xmin = xE - b
    else:
      xmin = xE - xlimb
    if ((theta_eff > 0) and (b > -xlimb)) or ((theta_eff <= 0) and (b < -xlimb)):
      xmax = xE + b
    else:
      xmax = xE - xlimb
        
    # Plot it
    x = np.linspace(xE - b, xE + b, 1000)
    if theta_eff > 0:
      x[x < xE - xlimb] = np.nan
    elif theta_eff > -np.pi / 2:
      x[x > xE - xlimb] = np.nan
    A = b ** 2 - (x - xE) ** 2
    A[A < 0] = 0
    y = (a / b) * np.sqrt(A)
    if np.abs(np.cos(z)) < 1e-5:
      style = dict(color = 'k', ls = '--', lw = 0.5, zorder = 0)
    else:
      style = dict(color = 'k', ls = '-', lw = 0.5, zorder = 0)
    ax.plot(x0 + x, y0 + y, **style)
    ax.plot(x0 + x, y0 - y, **style)
    
    # Fill the ellipses
    if theta_eff < 0:
      ax.fill_between(x0 + x, y0 - y, y0 + y, color = ZenithColor(zarr[i+1]), zorder = 0.5 * (z / np.pi - 1))
    else:
      ax.fill_between(x0 + x, y0 - y, y0 + y, color = ZenithColor(zarr[i]), zorder = 0.5 * (-z / np.pi - 1))
  
    # Fill the ellipses that are cut off by the limb
    if theta_eff < 0:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x0 + x_, y0 - y_, y0 + y_, color = ZenithColor(zarr[i]), zorder = 0.5 * (-z / np.pi - 1))
    else:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x0 + x_, y0 - y_, y0 + y_, color = ZenithColor(zarr[i]), zorder = 0.5 * (-z / np.pi - 1))
  
  return fig, ax
  
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
    DrawEyeball(fig = self.fig, theta = self.theta.val * np.pi / 180, dphi = self.psi.val * np.pi / 180, dlam = self.lam.val * np.pi / 180)

if __name__ == '__main__':
  
  Interact()
  