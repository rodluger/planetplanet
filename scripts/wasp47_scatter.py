#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
wasp47_scatter.py |github|
--------------------------

Computes all occultations that occur in the TRAPPIST-1 system over a 
3 year time period for a random draw from the prior. Plots each
occultation as a circle in a top-view of the system; the circle size, 
transparency, and color indicate the duration, SNR, and 
occulting body, respectively.

  .. plot::
     :align: center
     
     from scripts import scatter
     scatter._test()

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/scripts/scatter.py"><i class="fa fa-github" aria-hidden="true"></i></a>`


'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import planetplanet
from planetplanet import Wasp47, jwst
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''

    '''

    Run(save = False)

def scatter_plot(system, tstart, tend, dt = 1 * MINUTE, sz = 0.2,
                 cadence = 10, filter = 'f770w'):
    '''
    Compute all occultations between `tstart` and `tend` and plot an 
    occultation scatter plot like the one in the paper.

    :param float tstart: The integration start time (BJD − 2,450,000)
    :param float tend: The integration end time (BJD − 2,450,000)
    :param float dt: The time resolution in days. Occultations shorter \
           than this will not be registered.
    :param float sz: The size scaling for the occultation circles. \
           Default `0.2`
    
    '''
    
    # Reset
    system._reset()
    time = np.arange(tstart, tend, dt)
        
    # Compute the wavelength grid for all MIRI filters
    lambda1 = 6
    lambda2 = 9
    R = 100
    wav = [lambda1]
    while(wav[-1] < lambda2):
        wav.append(wav[-1] + wav[-1] / R)
    wavelength = np.array(wav)
    
    # Compute all limb darkening coefficients
    for body in system.bodies:
        body.u = [None for ld in body.limbdark]
        for n, ld in enumerate(body.limbdark):
            if callable(ld):
                body.u[n] = ld(wavelength)
            elif not hasattr(ld, '__len__'):
                body.u[n] = ld * np.ones_like(wavelength)
            else:
                raise Exception("Limb darkening coefficients must be "
                                + "provided as a list of scalars or "
                                + "as a list of functions.")
        body.u = np.array(body.u)
        
        # HACK: Disable phase curves. The code would take *way*
        # too long to run, and they don't affect these statistics.
        body.phasecurve = False

    # Convert from microns to meters
    wavelength *= 1e-6
                
    # No oversampling
    time_hr = np.array(time)
    
    # Continuum flux
    system._continuum = np.zeros(len(time_hr) * len(wavelength))
    
    # Allocate memory
    system._malloc(len(time_hr), len(wavelength))
    
    # Call the light curve routine
    err = system._Flux(len(time_hr), np.ctypeslib.as_ctypes(time_hr), 
                     len(wavelength), 
                     np.ctypeslib.as_ctypes(wavelength), 
                     np.ctypeslib.as_ctypes(system._continuum), 
                     len(system.bodies), system._ptr_bodies, 
                     system.settings)
    assert err <= 0, "Error in C routine `Flux` (%d)." % err

    # Loop over all bodies and plot each occultation event as a circle
    nppo = 0
    figp, axp = pl.subplots(1, figsize = (8,8))
    figp.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)
    axp.axis('off')
    for bi, body in enumerate(system.bodies[1:]):

        # Identify the different events
        inds = np.where(body.occultor > 0)[0]
        difs = np.where(np.diff(inds) > 1)[0]

        # Plot the orbit outline
        f = np.linspace(0, 2 * np.pi, 1000)
        r = body.a * (1 - body.ecc ** 2) / (1 + body.ecc * np.cos(f))
        x = r * np.cos(body._w + f) - r * np.sin(body._w + f) \
              * np.cos(body._inc) * np.sin(body._Omega)
        z = r * np.sin(body._w + f) * np.sin(body._inc)
        axp.plot(x, z, 'k-', lw = 1, alpha = 0.05)
        n = np.argmin(1e10 * (x < 0) + np.abs(z))
        axp.annotate(body.name, xy = (x[n], z[n]), color = 'k', 
                     alpha = 0.2, fontweight = 'bold',
                     fontsize = 8, zorder = -99, ha = 'center', 
                     va = 'center')

        # Loop over individual ones
        plot_secondary = True
        for i in inds[difs]:

            # Loop over all possible occultors
            for occ in range(len(system.bodies)):

                # Is body `occ` occulting?
                if (body.occultor[i] & 2 ** occ):

                    # Note that `i` is the last index of the occultation
                    duration = np.argmax(body.occultor[:i][::-1] 
                                         & 2 ** occ == 0)
                    
                    # Indices of occultation
                    idx = range(i - duration, i + 1)
                    if len(idx) < 2:
                        continue
                    if (occ == 0) and not plot_secondary:
                        continue
                    
                    # Get SNR
                    w = jwst.get_miri_filter_wheel()
                    f = w[np.argmax([f.name.lower() == filter.lower() for f in w])]
                    f.compute_lightcurve(time[idx], body.flux[idx,:], 
                                         system.continuum[idx,:], 
                                         system.wavelength, 
                                         stack = 1,
                                         atel = 25., 
                                         thermal = True,
                                         quiet = True)                        
                    Nplan = f.lightcurve.Nsys
                    Nback = f.lightcurve.Nback
                    Nstar = f.lightcurve.Ncont
                    Nplan = Nplan[0] - Nplan
                    snr = np.sqrt(np.sum((Nplan) ** 2 / (Nstar + Nback)))
                    noise = np.nanmedian(np.sqrt(Nstar + Nback) / np.median(Nstar))

                    # Alpha proportional to the SNR
                    alpha = max(0.1, min(0.9, snr / 5))

                    # Size = duration in minutes * sz
                    ms = 20 * min(100, (duration * dt * 1440 * sz))

                    # If the occultor is the star, plot it only once
                    if (occ == 0):
                        point = axp.scatter(body.x[i], body.z[i], marker = 'o', 
                                 color = system.colors[occ], 
                                 alpha = alpha, s = ms, edgecolor = 'none',
                                 zorder = 100, picker = True)
                        
                        point.body = body
                        point.idx = idx
                        point.snr = snr
                        point.occultor = system.bodies[occ]
                        point.noise = noise
                        plot_secondary = False
                    else:
                        point = axp.scatter(body.x[i], body.z[i], marker = 'o', 
                                 color = system.colors[occ], 
                                 alpha = alpha, s = ms, edgecolor = 'none',
                                 zorder = bi, picker = True)
                        
                        point.body = body
                        point.idx = idx
                        point.snr = snr
                        point.occultor = system.bodies[occ]
                        point.noise = noise
                        nppo += 1

            # Check for mutual transits
            if system.bodies[0].occultor[i]:

                # Get all bodies currently occulting the star
                occultors = []
                for occ in range(1, len(system.bodies)):
                    if (system.bodies[0].occultor[i] & 2 ** occ):
                        occultors.append(occ)

                # Check if any of these occult each other
                for occ1 in occultors:
                    for occ2 in occultors:
                        if system.bodies[occ1].occultor[i] & 2 ** occ2:
                            axp.scatter(system.bodies[occ1].x[i], 
                                     system.bodies[occ1].z[i], marker = 'x',
                                     color = system.colors[occ2], 
                                     alpha = 1, zorder = 100, s = 20,
                                     picker = True)

    # Legend 1: Occultor names/colors
    axl1 = pl.axes([0.025, 0.775, 0.2, 0.2])
    axl1.axis('off')
    axl1.set_xlim(-0.5, 1.5)
    axl1.set_ylim(-len(system.bodies) // 2 - 1, 1.5)
    axl1.annotate('Occultations by', xy = (0.5, 1), ha = 'center', 
                  va = 'center', fontweight = 'bold')
    for j, body in enumerate(system.bodies):
        if j < len(system.bodies) // 2:
            x, y = (0, -j)
        else:
            x, y = (0.825, len(system.bodies) // 2 - j)
        axl1.plot(x, y, 'o', color = system.colors[j], ms = 6, alpha = 1, 
                  markeredgecolor = 'none')
        axl1.annotate(body.name, xy = (x + 0.1, y), xycoords = 'data',
                      ha = 'left', va = 'center', color = system.colors[j])

    # Legend 2: Size/duration
    axl2 = pl.axes([0.775, 0.775, 0.2, 0.2])
    axl2.axis('off')
    axl2.set_xlim(-1, 1)
    axl2.set_ylim(-3, 1.5)
    axl2.annotate('Duration', xy = (0., 1), ha = 'center', va = 'center', 
                  fontweight = 'bold')
    for j, duration in enumerate([10, 30, 60]):
        ms = 20 * min(100, (duration * sz))
        axl2.scatter(-0.65, -0.75 * j + 0.2, marker = 'o', color = 'k', s = ms, 
                  alpha = 0.65, edgecolor = 'none')
        axl2.annotate('%d minutes' % duration, 
                      xy = (-0.3, -0.75 * j + 0.2), xycoords = 'data',
                      ha = 'left', va = 'center', color = 'k')

    # Legend 3: Transparency/SNR
    axl3 = pl.axes([0.025, 0.025, 0.2, 0.2])
    axl3.axis('off')
    axl3.set_xlim(-0.5, 1.5)
    axl3.set_ylim(-2, 1.5)
    axl3.annotate('SNR', xy = (0.5, 0.65), ha = 'center', 
                  va = 'center', fontweight = 'bold')
    for j, snr in enumerate([0.0, 1.0, 2.0]):
        alpha = max(0.1, min(0.9, snr / 5))
        axl3.plot(-0.15, -0.75 * j, 'o', color = 'k', ms = 8, 
                  alpha = alpha, markeredgecolor = 'none')
        if j == 0:
            l = '<0.1'
        else:
            l = ' %.1f' % snr
        ann = axl3.annotate(l, xy = (-0.05, -0.75 * j), 
                            xycoords = 'data', ha = 'left', va = 'center', 
                            color = 'k')
        ann.set_fontname('Andale Mono')
    for j, snr in enumerate([3.0, 4.0, 5.0]):
        alpha = max(0.1, min(0.9, snr / 5))
        axl3.plot(0.675, -0.75 * j, 'o', color = 'k', ms = 8, 
                  alpha = alpha, markeredgecolor = 'none')
        if j == 2:
            l = '>5.0'
        else:
            l = ' %.1f' % snr
        ann = axl3.annotate(l, xy = (0.775, -0.75 * j), 
                            xycoords = 'data', ha = 'left', va = 'center', 
                            color = 'k')
        ann.set_fontname('Andale Mono')
    
    # Observer direction
    axp.annotate("To observer", xy = (0.5, -0.1), 
                 xycoords = "axes fraction", xytext = (0, 30),
                 ha = 'center', va = 'center', annotation_clip = False, 
                 color = 'cornflowerblue', textcoords = "offset points", 
                 arrowprops = dict(arrowstyle = "-|>", 
                                   color = 'cornflowerblue'))
    
    # Log
    if not system.settings.quiet:
        print("There were %d PPOs between t = %.2f and t = %.2f." % 
              (nppo, tstart, tend))
    
    # Interactivity
    def onpick(event):
    
        # Get the event info
        point = event.artist
        body = point.body
        idx = point.idx
        snr = point.snr
        occultor = point.occultor
        noise = point.noise
        
        # Plot the occultation
        fig, ax = pl.subplots(1, figsize = (8, 6))
        fig.subplots_adjust()
        sz = idx[-1] - idx[0]
        a = max(0, idx[0] - 2 * sz)
        b = min(len(system.time), idx[-1] + 2 * sz)
        b -= (b - a) % cadence
        fnorm = body.flux[a:b,len(system.wavelength) // 2] / system.continuum[a:b,len(system.wavelength) // 2]
        fnorm += (1 - fnorm[0])
        time = system.time[a:b]
        ax.plot(time, fnorm, color = 'b', lw = 1)
        
        # Downbin to cadence
        time = np.nanmean(time.reshape(-1, cadence), axis = 1)
        fnorm = np.nanmean(fnorm.reshape(-1, cadence), axis = 1)
        noise /= np.sqrt(cadence)
        
        # Plot the simulated observation
        obs = fnorm + np.random.randn(len(fnorm)) * noise
        ax.errorbar(time, obs, yerr = noise, fmt = "o", c = "k", ms = 2,
                    alpha = 0.4, zorder = 10, lw = 1)
        
        ax.set_xlabel('Time [days]', fontweight = 'bold')
        ax.set_ylabel('Normalized Flux', fontweight = 'bold')
        ax.set_title('%s occulted by %s, SNR = %.2f' % (body.name, occultor.name, snr))
        
        pl.show()
        

    figp.canvas.mpl_connect('pick_event', onpick)
    
    return figp

def Run(save = True):
    '''
    Computes and plots the scatter plot for the WASP-47
    system seen with JWST MIRI at 15 microns over the
    course of 3 years.
    
    :param bool eyeball: Assume the radiative equilibrium surface map? \
           Default :py:obj:`True`. If :py:obj:`False`, uses the limb-darkened
           surface map.
    
    '''
    

    # Instantiate the Wasp-47 system
    # Plot all occultations from October 2016 to October 2019
    system = Wasp47(sample = True, nbody = True, seed = 45)

    fig = scatter_plot(system, OCTOBER_08_2016, OCTOBER_08_2016 + 20) # debug 365 * 3)
    
    if save:
        fig.savefig("scatter.pdf", bbox_inches = 'tight')
    
    pl.show()

if __name__ == '__main__':
    
    Run()