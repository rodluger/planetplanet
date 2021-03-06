#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
lightcurve.py |github|
----------------------

Computes a light curve of the TRAPPIST-1 system over ten days, with
orbital parameters drawn at random from their prior distributions.
All transits, secondary eclipses, planet-planet occultations, and mutual
transits are shown. Click on an individual event to see it in detail.
Then click on "Play" to animate the event.

  .. plot::
     :align: center
     
     from scripts import lightcurve
     lightcurve._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/lightcurve.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from planetplanet import Trappist1
from planetplanet.constants import *
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np

def _test():
    '''

    '''

    plot()
    pl.show()

def plot(interactive = True):
    '''
    
    '''
    
    # Instantiate the Trappist-1 system
    system = Trappist1(sample = True, phasecurve = True, 
                       nbody = True)

    # Get the occultation light curves over 10 random days
    tstart = OCTOBER_08_2016 + np.random.random() * 365
    time = np.arange(tstart, tstart + 10., MINUTE)
    system.compute(time)
    
    # Normalize the time array
    for body in system.bodies:
        body.time -= body.time[0]
        body.time_hr -= body.time_hr[0]
    
    if interactive:
    
        system.plot_lightcurve(interactive = True)
        
    else:
    
        fig, ax = system.plot_lightcurve(interactive = False)

        # Appearance
        ax.set_xlim(0, 10.)
        ax.set_ylim(0.9815, 1.005)
        for body in system.bodies:
            ax.plot([-1.001e3, -1.000e3], [1, 1], color = body.color, 
                    label = body.name, lw = 4)
        ax.annotate("Occultations by", xy = (0.175, 0.92), ha = 'left', 
                    va = 'bottom', fontweight = 'bold', fontsize = 10, 
                    xycoords = "axes fraction")
        ax.legend(loc = (0.325, 0.9), ncol = 8, frameon = False, 
                  handlelength = 1)
        ax.set_xlabel("Time [days]", fontweight = 'bold', fontsize = 16)
        ax.set_ylabel("Normalized Flux", fontweight = 'bold', fontsize = 16)
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontsize(12)
        ax.get_xaxis().set_major_locator(MaxNLocator(10))
    
        return fig, ax
  
if __name__ == '__main__':
  
    # Generate the PDF
    fig, ax = plot(interactive = False)
    fig.savefig("lightcurve.pdf", bbox_inches = 'tight')
    pl.close()

    # Interactive mode
    plot(interactive = True) 
    pl.show()