#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py
----------

'''

import numpy as np
np.seterr(invalid = 'ignore')
import matplotlib.pyplot as pl
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.widgets import Slider

def ZenithColor(z, cmap = 'inferno'):
  '''
  Color of a given zenith angle slice
  
  '''
  
  return pl.get_cmap(cmap)(0.2 + 0.4 * (np.cos(z) + 1))

def Draw(x0 = 0, y0 = 0, r = 1, theta = np.pi / 3, nz = 11, alpha = 0, occultors = [], cmap = 'inferno', fig = None, pos = None):
  '''

  '''

  # The rotation transformation, Equation (E6) in the paper
  xy = lambda x, y: (x * np.cos(alpha) + y * np.sin(alpha), y * np.cos(alpha) - x * np.sin(alpha))
  
  # Set up the floating axis
  if fig is None:
    fig = pl.figure(figsize = (6,6))
  tr = Affine2D().rotate_deg(alpha * 180 / np.pi)
  x = 1. / (np.abs(np.cos(alpha)) + np.abs(np.sin(alpha)))
  scale = max([r] + [occultor['r'] for occultor in occultors])
  x *= scale
  grid_helper = floating_axes.GridHelperCurveLinear(tr, extremes=(-x, x, -x, x))
  ax_orig = floating_axes.FloatingSubplot(fig, 111, grid_helper = grid_helper)
  if pos is not None:
    ax_orig.set_position(pos)
  ax_orig.axis["bottom"].set_visible(False)
  ax_orig.axis["top"].set_visible(False)
  ax_orig.axis["left"].set_visible(False)
  ax_orig.axis["right"].set_visible(False)
  ax_orig.patch.set_alpha(0)
  fig.add_subplot(ax_orig)
  ax = ax_orig.get_aux_axes(tr)
  
  # Plot the occultors. Note that we need to transform
  # their position vectors since we're in a rotated frame.
  occ = [None for i in occultors]
  for i, occultor in enumerate(occultors):
    xo = occultor['x']
    yo = occultor['y']
    ro = occultor['r']
    zo = occultor.get('zorder', 1)
    co = occultor.get('color', 'lightgray')
    ao = occultor.get('alpha', 1)
    xo, yo = xy(xo, yo)
    occ[i] = pl.Circle((xo, yo), ro, color = co, ec = 'k', alpha = ao, zorder = zo, clip_on = False)
    ax.add_artist(occ[i])
    
  # Plot the occulted body
  x = np.linspace(-r, r, 1000)
  y = np.sqrt(r ** 2 - x ** 2)
  ax.plot(x0 + x, y0 + y, color = 'k', zorder = 0, lw = 1, clip_on = False)
  ax.plot(x0 + x, y0 - y, color = 'k', zorder = 0, lw = 1, clip_on = False)

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
      style = dict(color = 'k', ls = '--', lw = 0.5, zorder = 0, clip_on = False)
    else:
      style = dict(color = 'k', ls = '-', lw = 0.5, zorder = 0, clip_on = False)
    if scale < 3 * r:
      ax.plot(x0 + x, y0 + y, **style)
      ax.plot(x0 + x, y0 - y, **style)
    
    # Fill the ellipses
    if theta < 0:
      ax.fill_between(x0 + x, y0 - y, y0 + y, color = ZenithColor(zarr[i+1], cmap = cmap), zorder = 0.5 * (z / np.pi - 1), clip_on = False)
    else:
      ax.fill_between(x0 + x, y0 - y, y0 + y, color = ZenithColor(zarr[i], cmap = cmap), zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
  
    # Fill the ellipses that are cut off by the limb
    if theta < 0:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x0 + x_, y0 - y_, y0 + y_, color = ZenithColor(zarr[i], cmap = cmap), zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
    else:
      x_ = np.linspace(-r, xE - xlimb, 1000)
      y_ = np.sqrt(r ** 2 - x_ ** 2)
      ax.fill_between(x0 + x_, y0 - y_, y0 + y_, color = ZenithColor(zarr[i], cmap = cmap), zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
  
  return fig, ax, occ, xy
  
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
  
    # Get the angles
    theta = self.theta.val * np.pi / 180
    dpsi = self.psi.val * np.pi / 180
    dlambda = self.lam.val * np.pi / 180
  
    # The coordinates of the substellar point
    x_ss = -np.cos(theta + dpsi) * np.cos(dlambda)
    y_ss = np.sin(dlambda)

    # The rotation angle that makes the planet symmetric about the x-axis
    alpha = -np.arctan2(y_ss, -x_ss)

    # Compute the new effective theta
    if theta + dpsi < -np.pi:
      theta = 2 * np.pi - np.arccos(-x_ss * np.cos(alpha) - y_ss * np.sin(alpha))
    elif theta + dpsi < 0:
      theta = -np.arccos(-x_ss * np.cos(alpha) - y_ss * np.sin(alpha))
    elif theta + dpsi > np.pi:
      theta = -np.arccos(-x_ss * np.cos(alpha) - y_ss * np.sin(alpha))
    else:
      theta = np.arccos(-x_ss * np.cos(alpha) - y_ss * np.sin(alpha))
    
    # Plot the planet
    Draw(fig = self.fig, theta = theta, alpha = alpha)