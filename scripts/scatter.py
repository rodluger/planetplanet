#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
scatter.py |github|
-------------------

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
from planetplanet import Trappist1, jwst
from planetplanet.constants import *
import matplotlib.pyplot as pl
import numpy as np

def _test():
    '''

    '''

    Run(save = False)

def scatter_plot(system, tstart, tend, dt = 0.001, sz = 0.2):
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

    # Compute the wavelength grid. We are hard-coding the
    # 15 micron JWST MIRI band here.
    lambda1 = 12.5
    lambda2 = 17.5
    R = 50
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

        # Simulate an observation w/ JWST at 15 microns
        # Same syntax as in `observe()`
        w = jwst.get_miri_filter_wheel()
        filter = w[np.argmax([f.name.lower() == 'f1500w' for f in w])]
        filter.compute_lightcurve(time, body.flux,
                                  system.continuum,
                                  system.wavelength,
                                  stack = 1,
                                  atel = 25.,
                                  thermal = True,
                                  quiet = True)

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

        # Total body photons
        total_body_photons = np.nanmedian(filter.lightcurve.Nsys)

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

                    # Planet, background, and star photons
                    Nplan = filter.lightcurve.Nsys[idx]
                    Nback = filter.lightcurve.Nback[idx]
                    Nstar = filter.lightcurve.Ncont[idx]

                    # Compute the number of photons *missing*
                    # NOTE: There was a BUG in the previous version,
                    # where we did
                    # >>> Nplan = np.nanmedian(Nplan) - Nplan
                    # which gets the wrong baseline for the planet's
                    # continuum. This led to low SNR in the previous
                    # versions of the plots!
                    Nplan = total_body_photons - Nplan

                    # Compute signal of and noise on the event
                    # in parts per million
                    norm = 1.e6 / np.sum(Nstar + Nback)
                    signal = norm * np.sum(np.fabs(Nplan))
                    noise = norm * np.sqrt(np.sum(Nstar + Nback))

                    # Compute the actual SNR on event. Note that this
                    # is NOT the sum of the signals divided by the sum
                    # of the noises! We need to add the SNR of each
                    # *datapoint* individually in quadrature.
                    snr = np.sqrt(np.sum((Nplan) ** 2 / (Nstar + Nback)))

                    # Alpha proportional to the SNR
                    alpha = max(0.05, min(0.95, snr / 3))

                    # Size = duration in minutes * sz
                    ms = min(100, (duration * dt * 1440 * sz))

                    # If the occultor is the star, plot it only once
                    if (occ == 0):
                        if plot_secondary:
                            axp.plot(body.x[i], body.z[i], 'o',
                                     color = system.colors[occ],
                                     alpha = alpha, ms = ms,
                                     markeredgecolor = 'none',
                                     zorder = 100)
                            plot_secondary = False
                    else:
                        axp.plot(body.x[i], body.z[i], 'o',
                                 color = system.colors[occ],
                                 alpha = alpha, ms = ms,
                                 markeredgecolor = 'none',
                                 zorder = bi)
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
                            axp.plot(system.bodies[occ1].x[i],
                                     system.bodies[occ1].z[i], 'x',
                                     color = system.colors[occ2],
                                     alpha = 1, zorder = 100, ms = 20)

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
        ms = min(100, (duration * sz))
        axl2.plot(-0.65, -0.75 * j + 0.2, 'o', color = 'k', ms = ms,
                  alpha = 0.65, markeredgecolor = 'none')
        axl2.annotate('%d minutes' % duration,
                      xy = (-0.3, -0.75 * j + 0.2), xycoords = 'data',
                      ha = 'left', va = 'center', color = 'k')

    # Legend 3: Transparency/SNR
    axl3 = pl.axes([0.025, 0.025, 0.2, 0.2])
    axl3.axis('off')
    axl3.set_xlim(-0.5, 1.5)
    axl3.set_ylim(-2, 1.5)
    axl3.annotate('S/N', xy = (0.5, 0.65), ha = 'center',
                  va = 'center', fontweight = 'bold')
    for j, snr in enumerate([0, 0.5, 1.0, 1.5]):
        alpha = max(0.05, min(0.95, snr / 3))
        axl3.plot(-0.15, -0.65 * j, 'o', color = 'k', ms = 8,
                  alpha = alpha, markeredgecolor = 'none', clip_on = False)
        if j == 0:
            l = '<0.1'
        else:
            l = ' %.1f' % snr
        ann = axl3.annotate(l, xy = (-0.05, -0.65 * j),
                            xycoords = 'data', ha = 'left', va = 'center',
                            color = 'k', clip_on = False)
        ann.set_fontname('Andale Mono')
    for j, snr in enumerate([2.0, 2.5, 3.0]):
        alpha = max(0.05, min(0.95, snr / 3))
        axl3.plot(0.675, -0.65 * j, 'o', color = 'k', ms = 8,
                  alpha = alpha, markeredgecolor = 'none')
        if j == 2:
            l = '>3.0'
        else:
            l = ' %.1f' % snr
        ann = axl3.annotate(l, xy = (0.775, -0.65 * j),
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

    return figp

def Run(eyeball = True, save = True):
    '''
    Computes and plots the scatter plot for the TRAPPIST-1
    system seen with JWST MIRI at 15 microns over the
    course of 3 years.

    :param bool eyeball: Assume the radiative equilibrium surface map? \
           Default :py:obj:`True`. If :py:obj:`False`, uses the limb-darkened
           surface map.

    '''

    # Choose the correct radiance map
    if eyeball:
        radiancemap = planetplanet.RadiativeEquilibriumMap()
    else:
        radiancemap = planetplanet.LimbDarkenedMap()

    # Instantiate the Trappist-1 system
    # Plot all occultations from October 2016 to October 2019
    system = Trappist1(sample = True, nbody = True, radiancemap = radiancemap,
                       seed = 44)

    fig = scatter_plot(system, OCTOBER_08_2016, OCTOBER_08_2016 + 365 * 3)

    if save:
        fig.savefig("scatter.pdf", bbox_inches = 'tight')

    pl.show()

if __name__ == '__main__':

    Run(eyeball = True)
