#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
triple_transit.py |github|
--------------------------

Computes and plots a hypothetical triple mutual transit event, where three large 
planets transit the star and occult each other simultaneously.

  .. plot::
     :align: center
     
     from scripts import triple_transit
     triple_transit._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/triple_transit.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                unicode_literals
from planetplanet import Planet, Star, System, DrawEyeball, LimbDarkenedMap
import matplotlib
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
import numpy as np
np.random.seed(1234)

def _test():
    '''
    
    '''
    
    plot()

def u1(lam):
    '''
    A really silly linear limb darkening law with a linear
    wavelength dependence.
    
    '''
    
    lam = np.atleast_1d(lam)
    result = 0.5 * (1 - (lam - 5) / 10) + 0.5
    result[lam < 5] = 1.
    result[lam > 15] = 0.5
    return result

def plot():
    '''
    
    '''
    
    # Instantiate the star
    star = Star('Star', m = 0.1, r = 0.1, nz = 31, color = 'k', 
                limbdark = [u1])

    # Planet b
    b = Planet('b', m = 1, per = 2, inc = 90.4, r = 2., t0 = 0, 
               nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

    # Planet c
    c = Planet('c', m = 1, per = 8, inc = 90., r = 2., t0 = 0.0005, 
               nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

    # Planet d
    d = Planet('d', m = 1, per = 32, inc = 89.94, r = 2., t0 = 0.002, 
               nz = 1, Omega = 0, w = 0., ecc = 0., phasecurve = False)

    # System
    system = System(star, b, c, d)

    # Get the occultation light curves
    time = np.linspace(-0.06, 0.06, 1000)
    system.compute(time)

    # Set up the figure
    fig = pl.figure(figsize = (7, 7))
    fig.subplots_adjust(left = 0.175)

    # Plot three different wavelengths (first, mid, and last)
    axlc = pl.subplot2grid((60, 5), (15, 0), colspan = 5, rowspan = 40)
    axlc.plot(star.time * 1440, star.flux[:, 0] / star.flux[0, 0], 'b-', 
              label = r"$" + 
                '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[0])) 
                + r"\ \mu\mathrm{m}$")
    axlc.plot(star.time * 1440, 
              star.flux[:, star.flux.shape[-1] // 2] 
              / star.flux[0, star.flux.shape[-1] // 2], 'g-', 
              label = r"$" + '{:.4s}'.format('{:0.2f}'.format(
                1e6 * star.wavelength[star.flux.shape[-1] // 2])) 
                + r"\ \mu\mathrm{m}$")
    axlc.plot(star.time * 1440, star.flux[:, -1] / star.flux[0, -1], 'r-', 
              label = r"$" + 
                '{:.4s}'.format('{:0.2f}'.format(1e6 * star.wavelength[-1])) 
                + r"\ \mu\mathrm{m}$")
    axlc.set_xlabel('Time [minutes]', fontweight = 'bold', fontsize = 14)
    axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 14)
    axlc.get_yaxis().set_major_locator(MaxNLocator(4))
    axlc.get_xaxis().set_major_locator(MaxNLocator(8))
    axlc.margins(0, None)
    for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
        tick.set_fontsize(12)
    axlc.legend(loc = 'lower right', fontsize = 12)

    # Plot the images
    t = [300, 400, 500, 600, 700]

    x0 = 0.535
    dx = 0.15
    px = [x0 - 2 * dx, x0 - dx, x0, x0 + dx, x0 + 2 * dx]
    for n in range(5):
    
        # Convert them into a list of dicts
        occ_dict = []
        for i, occultor in enumerate([b, c, d]):
            occ_dict.append(dict(x = occultor.x_hr[t[n]] / star._r, 
                                 y = occultor.y_hr[t[n]] / star._r, 
                                 r = occultor._r / star._r, 
                                 zorder = i + 1, alpha = 1))
    
        # Draw the eyeball planet and the occultors
        DrawEyeball(px[n], 0.85, 0.0425, LimbDarkenedMap(), theta = np.pi / 2, 
                    nz = 31, gamma = 0, occultors = occ_dict, cmap = 'inferno', 
                    fig = fig, draw_ellipses = False, teff = star.teff, 
                    limbdark = star.limbdark)

    # Arrows
    axlc.annotate("", xy = (star.time[t[0]] * 1440, 1.006), xycoords = "data", 
                  xytext = (-80, 63), textcoords = "offset points", 
                  clip_on = False, 
                  arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

    axlc.annotate("", xy = (star.time[t[1]] * 1440, 1.006), xycoords = "data", 
                  xytext = (-40, 63), textcoords = "offset points", 
                  clip_on = False, 
                  arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

    axlc.annotate("", xy = (star.time[t[2]] * 1440, 1.006), xycoords = "data", 
                  xytext = (0, 63), textcoords = "offset points", 
                  clip_on = False, 
                  arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

    axlc.annotate("", xy = (star.time[t[3]] * 1440, 1.006), xycoords = "data", 
                  xytext = (40, 63), textcoords = "offset points", 
                  clip_on = False, 
                  arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

    axlc.annotate("", xy = (star.time[t[4]] * 1440, 1.006), xycoords = "data", 
                  xytext = (80, 63), textcoords = "offset points", 
                  clip_on = False, 
                  arrowprops = dict(arrowstyle = '-', alpha = 0.5, lw = 1))

    # Save it if we're running this as a script
    if __name__ == '__main__':
        fig.savefig('triple.pdf', bbox_inches = 'tight')

    # Animate!
    fig2, axlc, axxz, axim = system.plot_occultation('A', -0.05, nz = 51, 
                                                     draw_ellipses = False, 
                                                     draw_terminator = False)
    axlc.set_xlabel('Time [days]', fontsize = 10, fontweight = 'bold')
    
    # Show both plots
    pl.show()

if __name__ == '__main__':

    plot()