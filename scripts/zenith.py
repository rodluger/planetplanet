#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
zenith.py |github|
------------------

Plots an interactive heatmap of the zenith angle for a given point on a sphere
oriented at an angle `theta` away from the observer. This is used
to calculate the radiance map of occulted planets.

  .. plot::
     :align: center
     
     from scripts import zenith
     zenith._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/zenith.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import numpy as np
import matplotlib.pyplot as pl
from matplotlib.widgets import Slider
cmap = pl.get_cmap('inferno_r')

def _test():
    '''
    
    '''

    plot()
    pl.show()

def ZenithAngle(x, y, r, theta):
    '''
    Compute the zenith angle.
    
    '''
    
    # Normalize
    x = x / r
    y = y / r
    x2 = x * x
    y2 = y * y
    
    # This is a solution to a quadratic equation in z = sin(za) **    2 
    z = 0.5 * ((1 - 2 * x2 - y2) * np.cos(2 * theta) + 2 * x \
        * np.sqrt(1 - x2 - y2) * np.sin(2 * theta) + y2 + 1)
    
    # Where are we relative to the terminator?
    xterm = np.sin(theta) * np.sqrt(np.abs(1 - y2))
    
    # Solve for the zenith angle
    if np.abs(theta) <= np.pi / 2:
        if (x <= xterm):
            return np.arcsin(np.sqrt(z)) 
        else:
            return np.pi - np.arcsin(np.sqrt(z))
    else:    
        if (x >= -xterm):
            return np.arcsin(np.sqrt(z)) 
        else:
            return np.pi - np.arcsin(np.sqrt(z))

def plot():
    '''
    
    '''
                
    fig, ax = pl.subplots(1)
    ax.axis('off')
    fig.subplots_adjust(bottom = 0.2)

    z = np.zeros((100, 100)) * np.nan
    img = pl.imshow(z, cmap = cmap, vmin = 0, vmax = 180., 
                    extent = (-1, 1, -1, 1))
    x = np.linspace(-0.99,0.99,1000)
    ax.plot(x, np.sqrt(0.99 ** 2 - x ** 2), 'k-', lw = 2)
    ax.plot(x, -np.sqrt(0.99 ** 2 - x ** 2), 'k-', lw = 2)
    ax.set_xlim(-1.1,1.1)
    ax.set_ylim(-1.1,1.1)

    axslider = pl.axes([0.3, 0.05, 0.44, 0.03])
    slider = Slider(axslider, r'$\theta$', -180., 180., valinit = 45.)

    def update(val):
        theta = slider.val
        for i, x in enumerate(np.linspace(-1,1,100)):
            for j, y in enumerate(np.linspace(-1,1,100)):
                if (x ** 2 + y ** 2 <= 1):
                    z[j,i] = ZenithAngle(x, y, 1, theta * np.pi / 180) \
                             * 180 / np.pi
        img.set_data(z)
        fig.canvas.draw_idle()
    
    slider.on_changed(update)
    update(45.)

if __name__ == '__main__':

    plot()
    pl.show()