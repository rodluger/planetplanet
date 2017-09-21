#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
orbits_wasp47.py |github|
-------------------------

Plots the orbital paths of the WASP-47 planets.

  .. plot::
     :align: center
     
     from scripts import orbits
     orbits._test()

  .. role:: raw-html(raw)
     :format: html
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/orbits_wasp47.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import numpy as np
import matplotlib.pyplot as pl

def _test():
    '''
    
    '''
    
    plot()
    pl.show()

def plot():
    '''
    
    '''
    
    # Setup
    AURJUP = 2092.51204
    rstar = 11.188 / AURJUP
    semis = np.array([0.052, 0.088, 0.0173])
    radii = np.array([1.17, 0.331, 0.167])
    incs = np.array([89.02, 89.22, 86.2]) * np.pi / 180
    colors = ['firebrick', 'gold', 'mediumseagreen']
    labels = ['b', 'd', 'e']
    fig = pl.figure(figsize = (12,5))

    # Plot the star
    x = np.linspace(-rstar, rstar, 1000)
    pl.fill_between(x, -np.sqrt(rstar ** 2 - x ** 2) * AURJUP, 
                    np.zeros_like(x), color = 'sandybrown')
    pl.plot(x, -np.sqrt(rstar ** 2 - x ** 2) * AURJUP, color = 'k', lw = 0.5)
    pl.fill_between(x, np.sqrt(rstar ** 2 - x ** 2) * AURJUP, 
                    np.zeros_like(x), color = 'sandybrown', zorder = 99)
    pl.plot(x, np.sqrt(rstar ** 2 - x ** 2) * AURJUP, 
            color = 'k', lw = 0.5, zorder = 99)

    # Plot the planet orbits
    for i in range(len(labels)):

        # Contours
        x = np.linspace(-semis[i], semis[i], 1000)
        a = semis[i]
        b = semis[i] * np.cos(incs[i])
        y = (b / a) * np.sqrt(a ** 2 - x ** 2) * AURJUP
        pl.plot(x, y, lw = 1, color = colors[i])
        pl.plot(x, -y, lw = 1, color = colors[i])
    
        # Dummy line for legend
        pl.plot(x, y + 999, lw = 2, color = colors[i], label = labels[i])
    
        # Fill
        pl.fill_between(x, y - radii[i], y + radii[i], color = colors[i], 
                        alpha = 0.5)
        pl.fill_between(x, -y - radii[i], -y + radii[i], color = colors[i], 
                        alpha = 0.5)

    pl.xlabel('x [AU]', fontsize = 16, fontweight = 'bold')
    pl.ylabel(r'y [R$_\mathrm{Jup}$]', fontsize = 16, fontweight = 'bold')
    pl.ylim(-5, 5)
    pl.legend(ncol = 2, loc = 'lower left', frameon = False)
    
    return fig, pl.gca()
    
if __name__ == '__main__':
    fig, _ = plot()
    fig.savefig('orbits_wasp47.pdf', bbox_inches = 'tight')