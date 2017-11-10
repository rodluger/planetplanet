#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
ppo.py |github|
---------------

The main Python interface to the C photodynamics code. Allows users to compute
and plot planet-planet occultation light curves, as well as transits, secondary
eclipses, phase curves, mutual transits, planet-moon events, and more.


  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/ppo.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
from ..constants import *
from ..detect import jwst
from .eyeball import DrawEyeball, GetAngles
from .structs import *
import ctypes
import numpy as np
np.seterr(invalid = 'ignore')
import os, shutil
import sysconfig
from numpy.ctypeslib import ndpointer, as_ctypes
import matplotlib
import matplotlib.pyplot as pl
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import Button
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap
from scipy.integrate import quad
from tqdm import tqdm
from six import string_types
import numba
import rebound

__all__ = ['System']

# Find system suffix and import the shared library
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
    suffix = ".so"
dn = os.path.dirname
libppo = ctypes.cdll.LoadLibrary(os.path.join(dn(dn(dn(
                                 os.path.abspath(__file__)))),
                                 "libppo" + suffix))

def _colorline(ax, x, y, color = (0, 0, 0), **kwargs):
    '''
    Plots the curve `y(x)` with linearly increasing alpha.
    Adapted from `http://nbviewer.jupyter.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb`_.

    '''

    # A bit hacky... But there doesn't seem to be
    # an easy way to get the hex code for a named color...
    if isinstance(color, string_types):
        if color.startswith("#"):
            hex = color[1:]
        else:
            if len(color) == 1:
                if color == 'k':
                    color = 'black'
                elif color == 'r':
                    color = 'red'
                elif color == 'b':
                    color = 'blue'
                elif color == 'g':
                    color = 'green'
                elif color == 'y':
                    color = 'yellow'
                elif color == 'w':
                    color = 'white'
                else:
                    # ?!
                    color = 'black'
            hex = matplotlib.colors.cnames[color.lower()][1:]
        r, g, b = tuple(int(hex[i:i+2], 16) / 255. for i in (0, 2, 4))
    else:
        r, g, b = color

    colors = [(r, g, b, i) for i in np.linspace(0, 1, 3)]
    cmap = LinearSegmentedColormap.from_list('alphacmap', colors, N = 1000)
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, array = np.linspace(0.0, 1.0, len(x)),
                        cmap = cmap, norm = pl.Normalize(0.0, 1.0), **kwargs)
    ax.add_collection(lc)

    return lc

class _Animation(object):
    '''
    An animation class for occultation movies.

    '''

    def __init__(self, t, fig, axim, tracker, occ, pt_xz, pt_zy, body, bodies,
                 occultors, interval = 50, gifname = None, xy = None,
                 quiet = False):
        '''

        '''

        self.t = t
        self.fig = fig
        self.axim = axim
        self.tracker = tracker
        self.occ = occ
        self.pt_xz = pt_xz
        self.pt_zy = pt_zy
        self.body = body
        self.bodies = bodies
        self.occultors = occultors
        self.xy = xy
        self.pause = True
        self.animation = animation.FuncAnimation(self.fig, self.animate,
                                                 frames = 100,
                                                 interval = interval,
                                                 repeat = True)

        # Save?
        if gifname is not None:
            self.pause = False
            if not gifname.endswith('.gif'):
                gifname += '.gif'
            if not quiet:
                print("Saving %s..." % gifname)
            self.animation.save(gifname, writer = 'imagemagick',
                                fps = 20, dpi = 150)
            self.pause = True

        # Toggle button
        self.axbutton = pl.axes([0.185, 0.12, 0.1, 0.03])
        self.button = Button(self.axbutton, 'Play', color = 'lightgray')
        self.button.on_clicked(self.toggle)

    def toggle(self, event):
        '''
        Pause/play the animation. Unfortunately I haven't been able to
        figure out how to freeze the animation so that it resumes at the
        same frame it was paused on...

        '''

        if self.pause:
            self.button.label.set_text('Pause')
        else:
            self.button.label.set_text('Play')
        self.pause ^= True

    def animate(self, j):
        '''
        Play frame `j` of the animation.

        '''

        if not self.pause:

            # Normalize the time index
            j = int(j * len(self.t) / 100.)

            # Time tracker
            self.tracker.set_xdata(self.bodies[0].time_hr[self.t[j]])

            # Occultor images
            x0 = self.body.x_hr[self.t[j]]
            y0 = self.body.y_hr[self.t[j]]
            for k, occultor in enumerate(self.occultors):
                if self.xy is None:
                    xo, yo = occultor.x_hr[self.t[j]] - x0, \
                             occultor.y_hr[self.t[j]] - y0
                else:
                    xo, yo = self.xy(occultor.x_hr[self.t[j]] - x0,
                                     occultor.y_hr[self.t[j]] - y0)
                self.occ[k].center = (xo / self.body._r, yo / self.body._r)

                # BODY orbits
                for k, b in enumerate(self.bodies):
                    self.pt_xz[k].set_xdata(b.x_hr[self.t[j]])
                    self.pt_xz[k].set_ydata(b.z_hr[self.t[j]])
                    self.pt_zy[k].set_xdata(b.z_hr[self.t[j]])
                    self.pt_zy[k].set_ydata(b.y_hr[self.t[j]])

class System(object):
    '''

    The planetary system class. This is the main interface to the
    photodynamical core. Instantiate with all bodies in the system and the
    desired settings, passed as :py:obj:`kwargs`.

    :param bodies: Any number of :py:func:`Planet`, :py:func:`Moon`, or \
           :py:func:`Star` instances comprising all the bodies in the system. \
           The first body is assumed to be the primary.
    :param bool nbody: Uses the :py:obj:`REBOUND` N-body code to compute \
           orbits. Default :py:obj:`True`
    :param float keptol: Kepler solver tolerance. Default `1.e-15`
    :param int maxkepiter: Maximum number of Kepler solver iterations. \
           Default `100`
    :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
    :param float timestep: Timestep in days for the N-body solver. \
           Default `0.01`
    :param bool adaptive: Adaptive grid for limb-darkened bodies? \
           Default :py:obj:`True`
    :param bool quiet: Suppress output? Default :py:obj:`False`
    :param float mintheta: Absolute value of the minimum phase angle in \
           degrees. Below this angle, elliptical boundaries of constant \
           surface brightness on the planet surface are treated as vertical \
           lines. Default `1.`
    :param int maxvertices: Maximum number of vertices allowed in the area \
           computation. Default `999`
    :param int maxfunctions: Maximum number of functions allowed in the area \
           computation. Default `999`
    :param int oversample: Oversampling factor for each exposure. Default `1`
    :param float distance: Distance to the system in parsecs. Default `10.`
    :param bool circleopt: Solve the simpler quadratic problem for \
           circle-ellipse intersections when the axes of the ellipse are \
           equal to within :math:`10^{-10}`? Default :py:obj:`True`
    :param bool batmanopt: Use the :py:mod:`batman` algorithm to compute \
           light curves of radially symmetric bodies? This can significantly \
           speed up the code. Default :py:obj:`True`
    :param str integrator: The N-body integrator \
           (:py:obj:`whfast` | :py:obj:`ias15`) to use. \
           Default :py:obj:`ias15`

    '''

    def __init__(self, *bodies, **kwargs):
        '''
        Initialize the class.

        '''

        # Initialize
        self.bodies = bodies
        self.settings = SETTINGS(**kwargs)
        self._reset()

    def _reset(self):
        '''
        Applies all settings, makes bodies accessible as properties, resets
        flags, and computes some preliminary orbital information for the
        system.

        '''

        # Move params set by the user over to the settings class
        for param in self.settings.params:
            if hasattr(self, param):
                setattr(self.settings, param, getattr(self, param))
                delattr(self, param)

        # Make planets accessible as properties
        self.primary = self.bodies[0]
        for body in self.bodies:
            setattr(self, body.name, body)
            # Force an update of `tperi0`
            body.ecc = body._ecc
        self._names = np.array([p.name for p in self.bodies])
        self.colors = [b.color for b in self.bodies]

        # Evaluate the body host names
        for i, body in enumerate(self.bodies):
            if body.host is None:
                # Host is CM of all interior bodies
                body.host = CM(*self.bodies[:i])
                body._host = -1
            elif type(body.host) is int:
                # Host is an index
                body._host = body.host
                body.host = self.bodies[body.host]
            elif isinstance(body.host, string_types):
                # Host is another body in the system
                body.host = getattr(self, body.host)
                body._host = np.argmax(self._names == body.host.name)
            elif body.host in self.bodies:
                # Host is another body in the system
                body._host = np.argmax(self._names == body.host.name)
            else:
                try:
                  assert body.host.body_type == 'cm', "Error!"
                except:
                  # Something went wrong
                  raise ValueError("Invalid `host` setting for body `%s`."
                                     % body.name)

        # Compute the semi-major axis for each body (in Earth radii)
        for body in self.bodies:
            body._computed = False
            body.a = ((body.per * DAYSEC) ** 2 * G
                      * (body.host._m + body._m) * MEARTH
                      / (4 * np.pi ** 2)) ** (1. / 3.) / REARTH

        # Reset continuum
        self._continuum = np.empty((0,), dtype = 'float64')

        # Reset animations
        self._animations = []

        # Reset flag
        self._computed = False

    def _malloc(self, nt, nw, noccultors = 0, noccultations = 0):
        '''
        Allocate memory for all the C arrays.

        '''

        # Check that the first body is a star
        assert self.bodies[0].body_type == 'star', \
               'The first body must be a :py:class:`Star`.'

        # Force cartesian coords for the star.
        self.bodies[0].cartesian = True

        # Check that if there are N stars, the first N bodies are stars
        nstars = np.sum([int(b.body_type == 'star') for b in self.bodies])
        assert np.all([b.body_type == 'star' for b in self.bodies[:nstars]]), \
               'Stars must be passed as the first arguments to `System`,' \
               'before any planets or moons.'
        self.settings.nstars = nstars
        self.nstars = nstars

        # Check that if there are multiple stars, Keplerian solver is off
        if nstars > 1:
            assert self.settings.nbody, "N-body integrator must be selected" \
                   "for systems with multiple stars."

        # If there's more than one star, disable the RadiativeEquilibriumMap
        if 'star' in [b.body_type for b in self.bodies[1:]]:
            for b in self.bodies[1:]:
                assert b.radiancemap.maptype not in [MAP_ELLIPTICAL_DEFAULT,
                    MAP_ELLIPTICAL_CUSTOM], "Elliptically-symmetric radiance "\
                    "maps not implemented for multiple-star systems."

        # Initialize the C interface
        self._Orbits = libppo.Orbits
        self._Orbits.restype = ctypes.c_int
        self._Orbits.argtypes = [ctypes.c_int,
                                 ctypes.ARRAY(ctypes.c_double, nt),
                                 ctypes.c_int,
                                 ctypes.POINTER(ctypes.POINTER(BODY)),
                                 SETTINGS]

        self._Flux = libppo.Flux
        self._Flux.restype = ctypes.c_int
        self._Flux.argtypes = [ctypes.c_int,
                               ctypes.ARRAY(ctypes.c_double, nt),
                               ctypes.c_int,
                               ctypes.ARRAY(ctypes.c_double, nw),
                               ctypes.ARRAY(ctypes.c_double, nw * nt),
                               ctypes.c_int,
                               ctypes.POINTER(ctypes.POINTER(BODY)),
                               SETTINGS]

        # Are we just checking for occultations?
        if (noccultors > 0) and (noccultations > 0):
            self._NextOccultation = libppo.NextOccultation
            self._NextOccultation.restype = ctypes.c_int
            self._NextOccultation.argtypes = \
                                [ctypes.c_int,
                                ctypes.ARRAY(ctypes.c_double, nt),
                                ctypes.c_int,
                                ctypes.POINTER(ctypes.POINTER(BODY)),
                                SETTINGS,
                                ctypes.c_int,
                                ctypes.c_int,
                                ctypes.ARRAY(ctypes.c_int, noccultors),
                                ctypes.c_int,
                                ctypes.ARRAY(ctypes.c_double, noccultations),
                                ctypes.ARRAY(ctypes.c_int, noccultations),
                                ctypes.ARRAY(ctypes.c_double, noccultations)]

        # Allocate memory for all the arrays
        for body in self.bodies:
            body.time = np.zeros(nt)
            body._time = np.ctypeslib.as_ctypes(body.time)
            body.wavelength = np.zeros(nw)
            body._wavelength = np.ctypeslib.as_ctypes(body.wavelength)
            body.x = np.zeros(nt)
            body.x[0] = body.x0
            body._x = np.ctypeslib.as_ctypes(body.x)
            body.y = np.zeros(nt)
            body.y[0] = body.y0
            body._y = np.ctypeslib.as_ctypes(body.y)
            body.z = np.zeros(nt)
            body.z[0] = body.z0
            body._z = np.ctypeslib.as_ctypes(body.z)
            body.vx = np.zeros(nt)
            body.vx[0] = body.vx0
            body._vx = np.ctypeslib.as_ctypes(body.vx)
            body.vy = np.zeros(nt)
            body.vy[0] = body.vy0
            body._vy = np.ctypeslib.as_ctypes(body.vy)
            body.vz = np.zeros(nt)
            body.vz[0] = body.vz0
            body._vz = np.ctypeslib.as_ctypes(body.vz)
            body.occultor = np.zeros(nt, dtype = 'int32')
            body._occultor = np.ctypeslib.as_ctypes(body.occultor)
            body.total_flux = np.zeros(nw)
            body._total_flux = np.ctypeslib.as_ctypes(body.total_flux)

            # HACK: The fact that flux is 2d is a nightmare for ctypes. We will
            # treat it as a 1d array within C and keep track of the row/column
            # indices by hand...
            body.flux = np.zeros((nt, nw))
            body._flux1d = body.flux.reshape(-1)
            body._flux = np.ctypeslib.as_ctypes(body._flux1d)

            # HACK: Same for the limb darkening coefficients
            body._u1d = body.u.reshape(-1)
            body._u = np.ctypeslib.as_ctypes(body._u1d)

            # Dimensions
            body.nu = len(body.u)
            body.nt = nt
            body.nw = nw

        # A pointer to a pointer to `BODY`. This is an array of `n`
        # `BODY` instances, passed by reference. The contents can all be
        # accessed through `bodies`
        # NOTE: Before I subclassed BODY, this used to be
        # >>> self._ptr_bodies = (ctypes.POINTER(BODY) * \
        # >>> len(self.bodies))(*[ctypes.pointer(p) for p in self.bodies])
        # I now cast the `Planet`, `Star`, and `Moon` instances as `BODY`
        # pointers, as per https://stackoverflow.com/a/37827528
        self._ptr_bodies = (ctypes.POINTER(BODY) * len(self.bodies)) \
                           (*[ctypes.cast(ctypes.byref(p),
                           ctypes.POINTER(BODY)) for p in self.bodies])

    def compute(self, time, lambda1 = 5, lambda2 = 15, R = 100):
        '''
        Compute the full system light curve over a given `time` array between
        wavelengths `lambda1` and `lambda2` at resolution `R`. This method
        runs the photodynamical core, populates all bodies with their
        individual light curves, and flags all occultation events.

        :param array_like time: The times at which to evaluate the light \
               curve (BJD − 2,450,000)
        :param float lambda1: The start point of the wavelength grid in \
               microns. Default `5`
        :param float lambda2: The end point of the wavelength grid in \
               microns. Default `15`
        :param float R: The spectrum resolution, \
               :math:`R = \\frac{\lambda}{\Delta\lambda}`. Default `100`

        '''

        # Reset
        self._reset()

        # Compute the wavelength grid
        wav = [lambda1]
        while(wav[-1] < lambda2):
            wav.append(wav[-1] + wav[-1] / R)
        wavelength = np.array(wav)

        # Compute all limb darkening coefficients
        for body in self.bodies:
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

        # Convert from microns to meters
        wavelength *= 1e-6

        # Oversample the time array to generate light curves
        # observed w/ finite exposure time.
        # Ensure `oversample` is odd.
        if self.settings.oversample % 2 == 0:
            oversample = self.settings.oversample + 1
        else:
            oversample = self.settings.oversample
        if oversample == 1:
            time_hr = np.array(time)
        else:
            tmid = (time[:-1] + time[1:]) / 2.
            tmid = np.concatenate(([tmid[0] - (tmid[1] - tmid[0])], tmid,
                                   [tmid[-1] + (tmid[-1] - tmid[-2])]))
            time_hr = [np.linspace(tmid[i], tmid[i + 1], oversample)
                       for i in range(len(tmid) - 1)]
            time_hr = np.concatenate(time_hr)

        # Continuum flux
        self._continuum = np.zeros(len(time_hr) * len(wavelength))

        # Allocate memory
        self._malloc(len(time_hr), len(wavelength))

        # Call the light curve routine
        err = self._Flux(len(time_hr), np.ctypeslib.as_ctypes(time_hr),
                         len(wavelength), np.ctypeslib.as_ctypes(wavelength),
                         np.ctypeslib.as_ctypes(self._continuum),
                         len(self.bodies), self._ptr_bodies, self.settings)
        assert err <= 0, "Error in C routine `Flux` (%d)." % err
        self._computed = True

        # Downbin the continuum flux to the original time array
        self._continuum_hr = np.array(self._continuum)
        if self.settings.oversample > 1:
            self._continuum = np.mean(self._continuum.reshape(-1, oversample,
                                      len(wavelength)), axis = 1)

        # Downbin other arrays to original time array
        for body in self.bodies:

            # Store the oversampled light curve
            body.time_hr = np.array(body.time)
            body.flux_hr = np.array(body.flux)
            body.occultor_hr = np.array(body.occultor)
            body.x_hr = np.array(body.x)
            body.y_hr = np.array(body.y)
            body.z_hr = np.array(body.z)
            body.vx_hr = np.array(body.vx)
            body.vy_hr = np.array(body.vy)
            body.vz_hr = np.array(body.vz)

            # Store the binned light curve
            if self.settings.oversample > 1:

                # Revert to original time array
                body.time = time

                # Average the flux over the exposure
                body.flux = np.mean(body.flux.reshape(-1, oversample,
                                                      len(wavelength)),
                                                      axis = 1)

                # Take the XYZ position at the bin center
                body.x = body.x[oversample // 2::oversample]
                body.y = body.y[oversample // 2::oversample]
                body.z = body.z[oversample // 2::oversample]
                body.vx = body.vx[oversample // 2::oversample]
                body.vy = body.vy[oversample // 2::oversample]
                body.vz = body.vz[oversample // 2::oversample]

                # Get all bodies that occult at some point over the exposure
                # The operation `bitwise_or.reduce` is *magical*
                body.occultor = np.bitwise_or.reduce(body.occultor.reshape(-1,
                                                     oversample), axis = 1)

        # Loop over all bodies to check for occultations
        for body in self.bodies:

            # Set the flag
            body._computed = True

            # Loop over all possible occultors
            for occ in range(len(self.bodies)):

                # Indices of occultation
                inds = np.where(body.occultor_hr & 2 ** occ)[0]

                # Store each one individually
                if len(inds):
                    inds = np.split(inds, 1 + np.where(np.diff(inds) > 1)[0])

                    # Loop over each occultation event
                    for i in inds:

                        # Add a little padding
                        di = (i[-1] - i[0]) // 2
                        ia = max(0, i[0] - di)
                        ib = min(len(body.time_hr), i[-1] + di + 1)
                        i = np.array(list(range(ia, i[0])) + list(i)
                                     + list(range(i[-1] + 1, ib)))

                        # Store it
                        body._inds.append(i)

    def compute_orbits(self, time):
        '''
        Run the dynamical code to compute the positions of all bodies over a
        given `time` array. This method does not compute light curves, but it
        does check for occultations.

        :param array_like time: The times at which to store the body \
               positions (BJD − 2,450,000)

        '''

        # Reset
        self._reset()
        for body in self.bodies:
            body.u = np.array([], dtype = float)

        self._malloc(len(time), 1)

        # Call the light curve routine
        err = self._Orbits(len(time), np.ctypeslib.as_ctypes(time),
                           len(self.bodies), self._ptr_bodies, self.settings)
        assert err <= 0, "Error in C routine `Orbits` (%d)." % err

    def observe(self, save = None, filter = 'f1500w', stack = 1,
                instrument = 'jwst', alpha_err = 0.7, figsize = (12,4),
                time_unit = 'BJD − 2,450,000'):
        '''
        Run telescope observability calculations for a system that has had its
        lightcurve computed. Calculates a synthetic noised lightcurve in the
        user specified filter.

        :param bool save: Save a text file and a plot. Default :py:obj:`None`
        :param filter: Filter name or :py:func:`Filter` object. \
               Default `'f1500w'`
        :type filter: str or :py:func:`Filter`
        :param int stack: Number of exposures to stack. Default `1`
        :param str instrument: Telescope instrument to use. Default `'jwst'`

        '''

        # Have we computed the light curves?
        assert self._computed, "Please run `compute()` first."

        # Get the filter "wheel"
        if instrument.lower() == 'jwst':

            wheel = jwst.get_miri_filter_wheel()
            atel = 25.0
            thermal = True

        elif instrument.lower() == 'ost':

            # Using JWST filters for OST for now
            wheel = jwst.get_miri_filter_wheel()
            # create new filter for 50 microns
            filt50 = jwst.create_tophat_filter(47.5, 52.5, dlam=0.1,
                                               Tput=0.3, name="OST50")
            # Will have mirror between 8-15m, let's say: 13.5m
            # This gives a telescope area of 144 m^2
            atel = np.pi * (13.5 / 2.) ** 2
            # No thermal noise
            thermal = False

        elif instrument.lower() == 'spitzer':

            wheel = jwst.get_spitzer_filter_wheel()
            atel = np.pi * (0.85 / 2.) ** 2
            thermal = False

        else:

            raise ValueError("Invalid instrument.")

        # Get the filter
        if type(filter).__name__ == "Filter":

            self.filter = filter

        elif filter.lower() in [f.name.lower() for f in wheel]:

            self.filter = wheel[np.argmax([f.name.lower() == filter.lower()
                                for f in wheel])]

        else:

            raise ValueError("Invalid filter.")

        # Check that filter wavelength width is fully contained within
        # the photodynamical wavelength grid
        assert (self.filter.eff_wl - self.filter.eff_dwl/2.) > \
               np.min(self.wavelength), "Wavelength grid does not extend short"\
               " enough for filter."
        assert (self.filter.eff_wl + self.filter.eff_dwl/2.) < \
               np.max(self.wavelength), "Wavelength grid does not extend long"\
               " enough for filter."

        # Compute lightcurve in this filter
        self.filter.compute_lightcurve(self.time, self.flux, self.continuum,
                                       self.wavelength, stack = stack,
                                       time_hr = self.time_hr,
                                       flux_hr = self.flux_hr,
                                       atel = atel, thermal = thermal,
                                       quiet = self.settings.quiet)

        # Setup plot
        fig, ax = pl.subplots(figsize = figsize)
        ax.set_title(r"%s" % self.filter.name)
        ax.set_ylabel("Relative Flux", fontweight = 'bold', fontsize = 10)
        ax.set_xlabel("Time [%s]" % time_unit, fontweight = 'bold',
                      fontsize = 10)

        # Plot lightcurve
        self.filter.lightcurve.plot(ax0 = ax, alpha_err = alpha_err)

        # Determine average precision attained in lightcurve
        ppm = np.mean(self.filter.lightcurve.sig
                      / self.filter.lightcurve.obs) * 1e6
        if not self.settings.quiet:
            print("Average lightcurve precision: %.3f ppm" %ppm)

        """
        Determine SNR on each event:
        """

        # Mask in and out of occultation times
        inmask = np.logical_or.reduce([planet.occultor > 0
                                       for planet in self.bodies])
        outmask = ~inmask
        arr = np.arange(len(inmask))

        # Break lightcurve mask into event chunks
        n = 0
        prev = False
        events = []
        for i in range(len(inmask)):

            if inmask[i]:

                if not prev:
                    n += 1
                    events.append([i])

                else:

                    events[n-1].append(i)

            prev = inmask[i]

        # Determine SNR on event detections
        SNRs = []
        signal = []
        noise = []
        for i in range(n):

            # In-event indices
            mask = events[i]

            # System photons over event
            Nsys = self.filter.lightcurve.Nsys[mask]

            # Background photons over event
            Nback = self.filter.lightcurve.Nback[mask]

            # Continuum photons over event
            Ncont = self.filter.lightcurve.Ncont[mask]

            # Compute signal of and noise on the event
            S = np.sum(np.fabs(Ncont - Nsys)) / np.sum(Nsys + Nback) * 1e6
            N = np.sqrt(np.sum(Nsys + Nback)) / np.sum(Nsys + Nback) * 1e6
            signal.append(S)
            noise.append(N)

            # Compute the actual SNR on event. Note that this is NOT
            # the sum of the signals divided by the sum of the noises!
            # We need to add the SNR of each *datapoint* individually
            # in quadrature.
            SNR = np.sqrt(np.sum((Ncont - Nsys) ** 2 / (Nsys + Nback)))
            SNRs.append(SNR)

        # Annotate event SNRs on each event
        for i in range(n):
                k = np.argmin(self.filter.lightcurve.obs[np.array(events[i])])
                j = events[i][k]
                ax.annotate(r"$%.1f \sigma$" % SNRs[i],
                            xy = (self.filter.lightcurve.time[
                                       int(np.mean(events[i]))],
                                  self.filter.lightcurve.obs[j]
                                  - 1 * self.filter.lightcurve.sig[j]),
                            ha = 'center',
                            va = 'center', color = "black",
                            fontweight = 'bold',
                            fontsize = 10, xytext = (0, -15),
                            textcoords = 'offset points')

        # Save SNR and ampl as attributes in Lightcurve
        self.filter.lightcurve.event_SNRs = SNRs
        self.filter.lightcurve.event_signal = signal
        self.filter.lightcurve.event_noise = noise

        # Save to disk?
        if save is not None:

            # Compose data array to save
            data = np.array([self.filter.lightcurve.time,
                             self.filter.lightcurve.obs,
                             self.filter.lightcurve.sig]).T

            # Did the user provide a file name?
            if save is True:

                # No
                np.savetxt("jwst_lc_%s.txt" % self.filter.name,
                           data, fmt = str("%.6e"),
                           header = "time [days]            flux" +
                                    "                 error",
                           comments = "")
                fig.savefig("jwst_lc_%s.png" % self.filter.name,
                            bbox_inches="tight")

            else:

                # Yes
                np.savetxt("%s.txt" % save,
                           data, fmt = str("%.6e"),
                           header = "time [days]            flux"+
                                    "                 error",
                           comments = "")
                fig.savefig("%s.png" % save, bbox_inches="tight")

        # Adjust for on-screen viewing
        fig.subplots_adjust(bottom = 0.2)

        return fig, ax

    def next_occultation(self, occulted, occultors = None,
                         tstart = OCTOBER_08_2016,
                         tend = OCTOBER_08_2016 + 100,
                         dt = 0.001, noccultations = 1):
        '''
        Computes the time of the next occultation of body `occulted`.

        :param occulted: The occulted body, passed as a string corresponding \
               to the body name, the body instance, or the index of the body.
        :param occultors: The occultor(s), passed as a list of strings, body \
               instances, or indices. Default is to consider occultations by \
               all bodies in the system.
        :param float tstart: Time at which to start searching for \
               occultations (BJD − 2,450,000). Default `8000.` \
               (12:00:00 UT October 8, 2016)
        :param float tend: Time at which to end the search if fewer than \
               `noccultations` occultations have been found \
               (BJD − 2,450,000). Default `8100.`
        :param float dt: The search timestep in days. Occultations shorter \
               than this value will be missed. Default `0.001` \
               (about 2 minutes)
        :param int noccultations: The number of occultations to search for. \
               Once this many occultations have been found, halts the N-body \
               integration and returns. Default `1`
        :returns: Arrays corresponding to the times of the occultations \
               (BJD − 2,450,000), the occulting bodies, and the durations \
               (in minutes)

        '''

        assert tend > tstart, "The end time must be larger than the start time!"

        # Convert `occulted` to an index
        if isinstance(occulted, string_types):
            occulted = np.argmax([b.name == occulted for b in self.bodies])
        elif occulted in self.bodies:
            occulted = np.argmax(np.array(self.bodies) == occulted)
        else:
            try:
                if occulted < len(self.bodies):
                    pass
                else:
                    raise TypeError("")
            except TypeError:
                raise ValueError("Invalid value for `occulted`.")

        # If not set, default to occultations by all bodies
        if occultors is None:
            occultors = np.arange(len(self.bodies))

        # Convert `occultors` to an array of indices
        occultors = np.atleast_1d(occultors)
        for i, occultor in enumerate(occultors):
            if isinstance(occultor, string_types) or \
              (type(occultor) is np.str_):
                occultors[i] = np.argmax([b.name == occultor
                                          for b in self.bodies])
            elif occultor in self.bodies:
                occultors[i] = np.argmax(np.array(self.bodies) == occultor)
            else:
                try:
                    if occultor < len(self.bodies):
                        pass
                    else:
                        raise TypeError("")
                except TypeError:
                    raise ValueError("Invalid value for `occultors`.")
        occultors = np.array(occultors, dtype = 'int32')
        noccultors = len(occultors)

        # Get the time array
        time = np.arange(tstart, tend, dt)

        # Construct empty occultation arrays
        occultation_times = np.zeros(noccultations, dtype = 'float64') * np.nan
        occultation_inds = np.zeros(noccultations, dtype = 'int32')
        occultation_durs = np.zeros(noccultations, dtype = 'float64')

        # Reset and allocate memory
        self._reset()
        for body in self.bodies:
            body.u = np.array([], dtype = float)
        self._malloc(len(time), 1, noccultors = noccultors,
                     noccultations = noccultations)

        # Final check: currently only supported for the N-body integrator.
        # TODO: Make this work for the Keplerian integrator.
        assert self.settings.nbody, "Currently, `next_occultation` works " + \
                "only with the N-Body solver. Please set `nbody = True`."

        # Call the orbit routine
        err = self._NextOccultation(len(time), np.ctypeslib.as_ctypes(time),
                                    len(self.bodies),
                                    self._ptr_bodies, self.settings, occulted,
                                    noccultors,
                                    np.ctypeslib.as_ctypes(occultors),
                                    noccultations,
                                    np.ctypeslib.as_ctypes(occultation_times),
                                    np.ctypeslib.as_ctypes(occultation_inds),
                                    np.ctypeslib.as_ctypes(occultation_durs))

        times = []
        occultors = []
        durations = []
        for t, i, d in zip(occultation_times, occultation_inds,
                           occultation_durs):
            if not np.isnan(t):
                if not self.settings.quiet:
                    print(self.bodies[occulted].name + " occulted by " +
                          self.bodies[i].name + " at t = %8.3f days" % t +
                          " lasting %5.2f minutes." % (d / MINUTE))
                times.append(t)
                occultors.append(self.bodies[i])
                durations.append(d / MINUTE)

        return times, occultors, durations

    def plot_occultation(self, body, time, wavelength = 15., interval = 50,
                         gifname = None, spectral = False,
                         time_unit = 'BJD − 2,450,000',
                         trail_dur = 1., full_lightcurve = False,
                         **kwargs):
        '''
        Plots and animates an occultation event.

        :param body: The occulted body
        :type body: :py:class:`BODY`
        :param float time: The approximate time of the occultation event \
               (BJD − 2,450,000)
        :param float wavelength: The wavelength in microns at which to plot \
               the light curve. Must be within the wavelength grid. \
               Default `15`
        :param int interval: The interval between frames in the animation in \
               ms. Default `50`
        :param str gifname: If specified, saves the occultation animation as \
               a :py:obj:`gif` in the current directory. Default :py:obj:`None`
        :param bool spectral: Plot the light curve at different wavelengths? \
               If :py:obj:`True`, plots the first, middle, and last \
               wavelength in the wavelength grid. If :py:obj:`False`, plots \
               the specified `wavelength`. Default :py:obj:`True`
        :param float trail_dur: Duration of the orbital trails in days when \
               plotting the orbits. The plotting code actually calls \
               :py:mod:`rebound` to compute these. Default `1`
        :param bool full_lightcurve: Plot the full system light curve? If \
               :py:obj:`True`, the light curve will include the contribution \
               of any other bodies being occulted at the given time. Default \
               :py:obj:`False`, in which case only the flux from :py:obj:`body`\
               is plotted.
        :param kwargs: Any additional keyword arguments to be passed to \
               :py:func:`plot_image`

        :returns: :py:obj:`(fig, ax1, ax2, ax3)`, a figure instance and \
                  its three axes

        .. plot::
             :align: center

             from scripts import occultation
             import matplotlib.pyplot as pl
             occultation.plot()
             pl.show()

        '''

        # Have we computed the light curves?
        assert self._computed, "Please run `compute()` first."

        # Get the wavelength index
        assert (wavelength >= self.bodies[0].wavelength[0]) and \
               (wavelength >= self.bodies[0].wavelength[-1]), \
               "Wavelength value outside of computed grid."
        w = np.argmax(1e-6 * wavelength <= self.bodies[0].wavelength)

        # Check file name
        if gifname is not None:
            if gifname.endswith(".gif"):
                gifname = gifname[:-4]

        # Get the occulted body
        if isinstance(body, string_types):
            p = np.argmax(self._names == body)
            body = self.bodies[p]
        else:
            p = np.argmax(self.bodies == body)

        # Get the indices of the occultation
        tind = np.argmin(np.abs(body.time_hr - time))
        if len(body._inds) == 0:
            if not self.settings.quiet:
                print("No occultation occurs at the given time.")
            return None
        iind = np.argmax([tind in inds for inds in body._inds])
        if (iind == 0) and not (tind in body._inds[0]):
            if not self.settings.quiet:
                print("No occultation occurs at the given time.")
            return None
        t = body._inds[iind]

        if not self.settings.quiet:
            print("Plotting the occultation...")

        # Set up the figure
        fig = pl.figure(figsize = (7, 8))
        fig.subplots_adjust(left = 0.175)
        axlc = pl.subplot2grid((5, 3), (3, 0), colspan = 3, rowspan = 2)
        ti = t[0]
        tf = t[-1]

        # Are we including the flux from other bodies?
        if full_lightcurve:
            flux = self.flux_hr
            total_flux = np.sum([b.total_flux for b in self.bodies], axis = 0)
        else:
            flux = body.flux_hr
            total_flux = body.total_flux

        # Get the duration in minutes. We find the occultor(s)
        # at the exact time requested by the user, and compute how
        # many cadences this occultation lasts. Note that this kind
        # of assumes a constant cadence and doesn't work that well
        # for mutual transits and other weird events.
        occultor = body.occultor_hr[tind]
        ncad = np.count_nonzero(body.occultor_hr[ti:tf] == occultor)
        duration = ncad * np.nanmedian(np.diff(body.time_hr[ti:tf])) * 1440.

        # Plot three different wavelengths (first, mid, and last)
        if spectral:
            fluxb = flux[t, 0] / total_flux[0]
            fluxg = flux[t, flux.shape[-1] // 2] \
                    / total_flux[flux.shape[-1] // 2]
            fluxr = flux[t, -1] / total_flux[-1]

            # Plot
            axlc.plot(body.time_hr[t], fluxb, 'b-',
                      label = r"$" + '{:.4s}'.format('{:0.2f}'.format(
                              1e6 * body.wavelength[0])) + r"\ \mu\mathrm{m}$")
            axlc.plot(body.time_hr[t], fluxg, 'g-',
                      label = r"$" + '{:.4s}'.format('{:0.2f}'.format(
                              1e6 * body.wavelength[flux.shape[-1] // 2]))
                              + r"\ \mu\mathrm{m}$")
            axlc.plot(body.time_hr[t], fluxr, 'r-',
                      label = r"$" + '{:.4s}'.format('{:0.2f}'.format(
                              1e6 * body.wavelength[-1]))
                              + r"\ \mu\mathrm{m}$")
            axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold',
                            fontsize = 10)
            axlc.legend(loc = 'lower right', fontsize = 8)
        else:
            fluxw = flux[t, w] / total_flux[w]
            axlc.plot(body.time_hr[t], fluxw, 'k-')
            label = r"$" + '{:.4s}'.format('{:0.2f}'.format(
                    1e6 * body.wavelength[w])) + r"\ \mu\mathrm{m}$"
            axlc.set_ylabel(r'Normalized Flux (%s)' % label,
                            fontweight = 'bold', fontsize = 10)

        axlc.set_xlabel('Time [%s]' % time_unit, fontweight = 'bold',
                        fontsize = 10)
        axlc.get_yaxis().set_major_locator(MaxNLocator(4))
        axlc.get_xaxis().set_major_locator(MaxNLocator(4))
        tracker = axlc.axvline(body.time_hr[ti], color = 'k',
                               alpha = 0.5, lw = 1, ls = '--')
        for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
            tick.set_fontsize(8)
        axlc.ticklabel_format(useOffset = False)
        if body.time_hr[ti] > 1e4:
            for label in axlc.get_xmajorticklabels():
                label.set_rotation(30)

        # Sort occultors by z-order (occultor closest to observer last)
        bits = np.bitwise_or.reduce(body.occultor_hr[t])
        occultors = []
        for b in range(len(self.bodies)):
            if (bits & 2 ** b):
                occultors.append(b)
        zorders = [-self.bodies[o].z_hr[ti] for o in occultors]
        occultors = [o for (z,o) in sorted(zip(zorders, occultors))]

        # Our orbital plot axes
        axxz = pl.subplot2grid((5, 4), (0, 0), colspan = 2, rowspan = 2)
        axzy = pl.subplot2grid((5, 4), (0, 2), colspan = 2, rowspan = 2)

        # Plot the locations of the bodies
        pt_xz = [None for b in self.bodies]
        pt_zy = [None for b in self.bodies]
        for bi, b in enumerate(self.bodies):
            pt_xz[bi], = axxz.plot(b.x_hr[ti], b.z_hr[ti], 'o', color = 'lightgrey',
                                   alpha = 1, markeredgecolor = b.color,
                                   zorder = 99, clip_on = False,
                                   ms = 5, markeredgewidth = 2)
            pt_zy[bi], = axzy.plot(b.z_hr[ti], b.y_hr[ti], 'o', color = 'lightgrey',
                                   alpha = 1, markeredgecolor = b.color,
                                   zorder = 99, clip_on = False,
                                   ms = 5, markeredgewidth = 2)

        # Plot the orbits. We will run rebound
        # backwards in time to get the orbital
        # trails. This avoids having to deal with
        # the case when the occultation is at the
        # beginning of the timeseries!
        if not self.settings.quiet:
            print("Computing orbits for plotting...")
        per = np.max([b.per for b in self.bodies])
        tlong = np.arange(0, trail_dur, self.settings.timestep)
        sim = rebound.Simulation()
        sim.G = GEARTH
        for b in self.bodies:
            sim.add(m = b._m, x = b.x_hr[tf], y = b.y_hr[tf],
                    z = b.z_hr[tf], vx = -b.vx_hr[tf],
                    vy = -b.vy_hr[tf], vz = -b.vz_hr[tf])
        x = np.zeros((len(self.bodies), len(tlong)))
        y = np.zeros((len(self.bodies), len(tlong)))
        z = np.zeros((len(self.bodies), len(tlong)))
        for i, time in enumerate(tlong):
            sim.integrate(time, exact_finish_time = 1)
            for q in range(len(self.bodies)):
                x[q,i] = sim.particles[q].x
                y[q,i] = sim.particles[q].y
                z[q,i] = sim.particles[q].z
        maxx = maxy = maxz = 0.
        for j, b in enumerate(self.bodies):

            # Top view
            _colorline(axxz, x[j][::-1], z[j][::-1], lw = 1,
                       clip_on = False, color = b.color)

            # Side view
            _colorline(axzy, z[j][::-1], y[j][::-1], lw = 1,
                       clip_on = False, color = b.color)

            # Axis limits
            maxx = max(maxx, np.max(x[j][::-1]), -np.min(x[j][::-1]))
            maxy = max(maxy, np.max(y[j][::-1]), -np.min(y[j][::-1]))
            maxz = max(maxz, np.max(z[j][::-1]), -np.min(z[j][::-1]))

        # Symmetrical limits
        axxz.set_ylim(-maxz, maxz)
        axxz.set_xlim(-maxx, maxx)
        axzy.set_ylim(-maxy, maxy)
        axzy.set_xlim(-maxz, maxz)
        axxz.set_aspect('equal')
        axzy.set_aspect('equal')
        axxz.axis('off')
        axzy.axis('off')

        # Legend
        axlg = pl.axes([0.2, 0.4, 0.8, 0.4])
        axlg.axis('off')
        axlg.annotate(r"$z$", xy=(0, 0.5), xytext=(0, 20),
                      fontsize = 8,
                      xycoords = 'axes fraction',
                      textcoords = 'offset points',
                      ha = 'center', va = 'center',
                      arrowprops=dict(arrowstyle="<|-", lw = 0.5, color = 'k'))
        axlg.annotate(r"$x$", xy=(-0.005, 0.5075), xytext=(18, 0),
                      fontsize = 8,
                      xycoords = 'axes fraction',
                      textcoords = 'offset points',
                      ha = 'center', va = 'center',
                      arrowprops=dict(arrowstyle="<|-", lw = 0.5, color = 'k'))
        axlg.annotate(r"$y$", xy=(0.5, 0.5), xytext=(0, 20),
                      fontsize = 8,
                      xycoords = 'axes fraction',
                      textcoords = 'offset points',
                      ha = 'center', va = 'center',
                      arrowprops=dict(arrowstyle="<|-", lw = 0.5, color = 'k'))
        axlg.annotate(r"$z$", xy=(0.495, 0.5075), xytext=(18, 0),
                      fontsize = 8,
                      xycoords = 'axes fraction',
                      textcoords = 'offset points',
                      ha = 'center', va = 'center',
                      arrowprops=dict(arrowstyle="<|-", lw = 0.5, color = 'k'))
        axlg.annotate(r"$y$", xy=(0, 0.05), xytext=(0, 20),
                      fontsize = 8,
                      xycoords = 'axes fraction',
                      textcoords = 'offset points',
                      ha = 'center', va = 'center',
                      arrowprops=dict(arrowstyle="<|-", lw = 0.5, color = 'k'))
        axlg.annotate(r"$x$", xy=(-0.005, 0.0575), xytext=(18, 0),
                      fontsize = 8,
                      xycoords = 'axes fraction',
                      textcoords = 'offset points',
                      ha = 'center', va = 'center',
                      arrowprops=dict(arrowstyle="<|-", lw = 0.5, color = 'k'))

        # Plot the image
        axim, occ, xy = self.plot_image(ti, body, occultors, fig = fig,
                                        wavelength = wavelength, **kwargs)
        bodies = [self.bodies[o] for o in occultors] + [body]
        xmin = min(np.concatenate([o.x_hr[t] - body.x_hr[t] - 1.1 * o._r
                                   for o in bodies]))
        xmax = max(np.concatenate([o.x_hr[t] - body.x_hr[t] + 1.1 * o._r
                                   for o in bodies]))
        ymin = min(np.concatenate([o.y_hr[t] - body.y_hr[t] - 1.1 * o._r
                                   for o in bodies]))
        ymax = max(np.concatenate([o.y_hr[t] - body.y_hr[t] + 1.1 * o._r
                                   for o in bodies]))
        axim.set_xlim(xmin, xmax)
        axim.set_ylim(ymin, ymax)
        axim.axis('off')
        axim.set_aspect('equal')

        # The title
        if len(occultors) == 1:
            axxz.annotate("%s occulted by %s" % (body.name,
                          self.bodies[occultors[0]].name), xy = (0.55, 0.95),
                          xycoords = "figure fraction", ha = 'center',
                          va = 'center',
                          fontweight = 'bold', fontsize = 12, clip_on = False)
        else:
            axxz.annotate("%s occulted by %s" % (body.name,
                          ", ".join([occultor.name for occultor in
                          [self.bodies[o] for o in occultors]])),
                          xy = (0.55, 0.95),
                          xycoords = "figure fraction", ha = 'center',
                          va = 'center',
                          fontweight = 'bold', fontsize = 12, clip_on = False)
        axxz.annotate("Duration: %.2f minutes" % duration,
                      xy = (0.55, 0.92), ha = 'center', va = 'center',
                      xycoords = 'figure fraction', fontsize = 10,
                      style = 'italic', clip_on = False)

        # Animate!
        if gifname is not None:
            tmp = '%s.%03d.gif' % (gifname, len(self._animations) + 1)
        else:
            tmp = None
        self._animations.append(_Animation(range(ti, tf), fig, axim, tracker,
                                occ, pt_xz, pt_zy, body, self.bodies,
                                [self.bodies[o] for o in occultors],
                                interval = interval, gifname = tmp,
                                quiet = self.settings.quiet, xy = xy))

        return fig, axlc, axxz, axim

    def plot_image(self, t, occulted, occultor = None, wavelength = 15.,
                   fig = None, figx = 0.535, figy = 0.5, figr = 0.05,
                   **kwargs):
        '''
        Plots an image of the `occulted` body and its occultor(s) at a
        given index of the `time_hr` array `t`.

        .. plot::
             :align: center

             from planetplanet.photo.eyeball import DrawEyeball
             import matplotlib.pyplot as pl
             DrawEyeball(0.535, 0.5, 0.25, theta = 1., gamma = 1.,
                         occultors = [{'x': 0.2, 'y': 0.4, 'r': 0.5}],
                         cmap = 'inferno')
             pl.show()

        :param int t: The index of the occultation in the high resolution \
               time array `time_hr`
        :param occulted: The occulted body instance
        :type occultor: :py:class:`BODY` or :py:obj:`list`
        :param occultor: The occultor(s). Default :py:obj:`None`
        :type occultor: :py:class:`BODY` or :py:obj:`list`
        :param float wavelength: The wavelength in microns at which to plot \
               the light curve. Must be within the wavelength grid. \
               Default `15`
        :param fig: The figure on which to plot the image
        :type fig: :py:class:`matplotlib.Figure`
        :param float figx: The x coordinate of the image in figure units. \
               Default `0.535`
        :param float figy: The y coordinate of the image in figure units. \
               Default `0.5`
        :param float figr: The radius of the image in figure units.
               Default `0.05`

        '''

        # Set up the plot
        if fig is None:
            fig = pl.gcf()

        # Get the occulted body
        rp = occulted._r
        x0 = occulted.x_hr[t]
        y0 = occulted.y_hr[t]
        z0 = occulted.z_hr[t]
        vx0 = occulted.vx_hr[t]
        vy0 = occulted.vy_hr[t]
        vz0 = occulted.vz_hr[t]

        # Get the angles
        theta, gamma = GetAngles(x0, y0, z0, vx0, vy0, vz0,
                                 Lambda = occulted._Lambda,
                                 Phi = occulted._Phi)

        # Get the occultor
        if occultor is None:
            occultor = []
            for b in range(len(self.bodies)):
                if (occulted.occultor[t] & 2 ** b):
                    occultor.append(b)
        elif not hasattr(occultor, '__len__'):
            occultor = [occultor]

        # If they are BODY instances, turn them into indices
        for i, occ in enumerate(occultor):
            if occ in self.bodies:
                occultor[i] = np.argmax([b == occ for b in self.bodies])

        # Convert them into a list of dicts
        occ_dict = []
        for i, occultor in enumerate(occultor):
            occultor = self.bodies[occultor]
            occ_dict.append(dict(x = (occultor.x_hr[t] - x0) / rp,
                                 y = (occultor.y_hr[t] - y0) / rp,
                                 r = occultor._r / rp, zorder = i + 1,
                                 alpha = 1,
                                 color = occultor.color))

        # Draw the eyeball planet and the occultor(s)
        fig, ax, occ, xy = DrawEyeball(figx, figy, figr, occulted.radiancemap,
                                       theta = theta, gamma = gamma,
                                       occultors = occ_dict,
                                       cmap = occulted.cmap,
                                       fig = fig, wavelength = wavelength,
                                       teff = occulted.teff,
                                       limbdark = occulted.limbdark,
                                       color = occulted.color,
                                       **kwargs)

        return ax, occ, xy

    def plot_lightcurve(self, wavelength = 15., interactive = True):
        '''

        Plot the light curve of the system after running :py:func:`compute`.

        .. plot::
             :align: center

             from scripts import lightcurve
             import matplotlib.pyplot as pl
             lightcurve.plot()
             pl.show()

        :param float wavelength: The wavelength in microns at which to plot \
               the light curve. Must be within the wavelength grid. \
               Default `15`
        :param bool interactive: Interactive (clickable) plot? Default \
               :py:obj:`True`. If :py:obj:`False`, returns a :py:obj:`fig` \
               and an :py:obj:`axis` instance.
        :returns: :py:obj:`(fig, ax)` if \
                  :py:obj:`interactive` = :py:obj:`False`

        '''

        # Have we computed the light curves?
        assert self._computed, "Please run `compute()` first."

        if not self.settings.quiet:
            print("Plotting the light curve...")

        # Plot
        fig, ax = pl.subplots(1, figsize = (12, 4))
        assert (wavelength >= self.bodies[0].wavelength[0]) and \
               (wavelength >= self.bodies[0].wavelength[-1]), \
               "Wavelength value outside of computed grid."
        w = np.argmax(1e-6 * wavelength <= self.bodies[0].wavelength)
        flux = self.flux[:,w] / np.nanmedian(self.flux[:,w])
        ax.plot(self.time, flux, 'k-', lw = 1)

        # Plot the hi resolution light curve
        flux_hr = self.flux_hr[:,w] / np.nanmedian(self.flux_hr[:,w])
        ax.plot(self.time_hr, flux_hr, 'k-', lw = 1, alpha = 0.3, picker = 10)
        fig.canvas.mpl_connect('pick_event', self._onpick)

        # Appearance
        ax.set_xlabel('Time [BJD − 2,450,000]', fontweight = 'bold',
                      fontsize = 10)
        ax.set_ylabel(r'Normalized Flux @ %.1f$\mathbf{\mu}$m' % wavelength,
                      fontweight = 'bold', fontsize = 10)
        ax.get_yaxis().set_major_locator(MaxNLocator(4))
        ax.get_xaxis().set_major_locator(MaxNLocator(4))
        for tick in ax.get_xticklabels() + ax.get_yticklabels():
            tick.set_fontsize(8)
        ax.ticklabel_format(useOffset = False)

        # Limits
        ymax = np.nanmax(flux)
        ymin = np.nanmin(flux)
        yrng = ymax - ymin
        ax.set_ylim(ymin - 0.2 * yrng, ymax + 0.2 * yrng)
        ax.margins(0, None)

        # Label all of the events. Store the list of annotations
        # if the user wants to tweak their positions
        self._lc_annotations = []
        nann = np.zeros(100)
        events = []

        for body in self.bodies:

            for i, t in enumerate(body._inds):

                tstart = t[0] + np.argmax(body.occultor_hr[t] > 0)
                tend = t[0] + len(body.time_hr[t]) - \
                              np.argmax(body.occultor_hr[t][::-1] > 0)
                tmid = (tstart + tend) // 2
                occultors = []

                for b in range(len(self.bodies)):
                    for ti in t:
                        if (body.occultor_hr[ti] & 2 ** b):
                            occultors.append(b)

                occultors = list(set(occultors))
                time = body.time_hr[tmid]
                tmin = body.time_hr[t][np.argmin(flux_hr[t])]
                fmin = np.min(flux_hr[t])

                for n, occultor in enumerate([self.bodies[o]
                                              for o in occultors]):

                    # Keep track of events so we don't double count them
                    event_id = [body.name, occultor.name, tmid]

                    if event_id not in events:

                        # This keeps track of how many annotations we're
                        # trying to squeeze into the same time bin. Add a
                        # vertical shift to avoid overlaps
                        idx = int(100 * (tmin - body.time_hr[0]) /
                                 (body.time_hr[-1] - body.time_hr[0]))

                        if idx < len(nann):
                            dy = -15 * (1 + nann[idx])
                            nann[idx] += 1
                            self._lc_annotations.append(ax.annotate("%s"
                                                        % body.name,
                                                        xy = (tmin, fmin),
                                                        ha = 'center',
                                                        va = 'center',
                                                        color = occultor.color,
                                                        fontweight = 'bold',
                                                        fontsize = 10,
                                                        xytext = (0, dy),
                                                        textcoords = \
                                                         'offset points',
                                                        clip_on = True))
                            events.append(event_id)

        # Interactive?
        if interactive:
            pl.show()
        else:
            return fig, ax

    def plot_orbits(self, bodies = 'all', trail_dur = 1.):
        '''
        Plot the orbits of a given set of bodies over `nper` periods.

        .. plot::
             :align: center

             import matplotlib.pyplot as pl
             import numpy as np
             import planetplanet as pp
             T1 = pp.Trappist1()
             T1.compute_orbits(np.linspace(0,10,10000))
             fig, ax = T1.plot_orbits()
             pl.show()

        :param planets: The planet(s) to plot. Can be a list of strings \
               corresponding to their names or a list of \
               :py:obj:`Planet <planetplanet.photo.structs.Planet>` instances.\
               Default is to plot all planets in the system
        :param float trail_dur: The duration of the orbital trails in days. \
               Default `1`
        :returns: A figure and an axis instance

        '''

        # Have we computed the orbits:
        try:
            self.time
        except AttributeError:
            raise Exception("Please run `compute()` or `compute_orbits()` " +
                            "first.")

        # Get planets
        if (bodies == 'all') or (bodies is None):
            bodies = self.bodies
        elif not hasattr(bodies, '__len__'):
            bodies = [bodies]

        # Set up
        fig, ax = pl.subplots(2,2, figsize = (8,8))
        fig.subplots_adjust(wspace = 0.2, hspace = 0.2)
        arrs = np.array([])

        for i, body in enumerate(bodies):

            # Ensure we have a planet instance
            if isinstance(body, string_types):
                body = self.bodies[np.argmax([b.name == body
                                   for b in self.bodies])]

            # Trim the arrays
            N = np.argmax(body.time > body.time[-1] - trail_dur)
            x = body.x[N:]
            y = body.y[N:]
            z = body.z[N:]

            # x-z
            _colorline(ax[0,0], x, z, lw = 1, color = body.color)
            ax[0,0].plot(x[-1], z[-1], '.', color = body.color)

            # x-y
            _colorline(ax[1,0], x, y, lw = 1, color = body.color)
            ax[1,0].plot(x[-1], y[-1], '.', color = body.color)

            # z-y
            _colorline(ax[1,1], z, y, lw = 1, color = body.color)
            ax[1,1].plot(z[-1], y[-1], '.', color = body.color)

            # Append to running array
            arrs = np.hstack((arrs, x, y, z))

        # Appearance
        ax[0,0].set_ylabel(r'z [R$_\oplus$]', fontsize = 14,
                           fontweight = 'bold')
        ax[1,0].set_xlabel(r'x [R$_\oplus$]', fontsize = 14,
                           fontweight = 'bold')
        ax[1,0].set_ylabel(r'y [R$_\oplus$]', fontsize = 14,
                           fontweight = 'bold')
        ax[1,1].set_xlabel(r'z [R$_\oplus$]', fontsize = 14,
                           fontweight = 'bold')
        ax[0,1].axis('off')

        ax[0,0].get_shared_x_axes().join(ax[0,0], ax[1,0])
        ax[1,0].get_shared_y_axes().join(ax[1,0], ax[1,1])

        lim = 1.2 * np.max((np.abs(np.min(arrs)), np.abs(np.max(arrs))))
        for axis in [ax[0,0], ax[1,0], ax[1,1]]:
            for tick in axis.get_xticklabels() + axis.get_yticklabels():
                tick.set_fontsize(8)
            axis.set_xlim(-lim, lim)
            axis.set_ylim(-lim, lim)

        return fig, ax

    def _onpick(self, event):
        '''
        Picker event for interactive light curves

        '''

        index = event.ind[len(event.ind) // 2]
        for body in self.bodies:
            for occultation in body._inds:
                if index in occultation:
                    self.plot_occultation(body.name, body.time_hr[index])
        pl.show()

    @property
    def flux(self):
        '''
        The total flux of the system computed on a grid of time and wavelength.

        '''

        return np.sum([b.flux for b in self.bodies], axis = 0)

    @property
    def continuum(self):
        '''
        The total *continuum* flux of the system computed on a
        grid of time and wavelength. This is the same as :py:attr:`flux`,
        but without any occultations/transits.

        '''

        return self._continuum.reshape(self.bodies[0].flux.shape)

    @property
    def time(self):
        '''
        Time array in days, BJD − 2,450,000

        '''

        return self.bodies[0].time

    @property
    def flux_hr(self):
        '''
        The total flux of the system computed on a grid of high
        resolution time and wavelength.

        '''

        return np.sum([b.flux_hr for b in self.bodies], axis = 0)

    @property
    def time_hr(self):
        '''
        High-resolution time array in days, BJD − 2,450,000

        '''

        return self.bodies[0].time_hr

    @property
    def wavelength(self):
        '''
        Wavelength array in microns.

        '''

        return self.bodies[0].wavelength * 1.e6
