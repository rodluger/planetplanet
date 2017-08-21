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

from __future__ import division, print_function, absolute_import, unicode_literals
from ..constants import *
from ..detect import jwst
from .eyeball import DrawEyeball
from .maps import NullMap
from .structs import *
import ctypes
import numpy as np
np.seterr(invalid = 'ignore')
import os, shutil
import sysconfig
from numpy.ctypeslib import ndpointer, as_ctypes
import matplotlib.pyplot as pl
import matplotlib.animation as animation
from matplotlib.ticker import MaxNLocator
from matplotlib.widgets import Button
from scipy.integrate import quad
from tqdm import tqdm
import numba

__all__ = ['System']

# Find system suffix and import the shared library
suffix = sysconfig.get_config_var('EXT_SUFFIX')
if suffix is None:
  suffix = ".so"
dn = os.path.dirname
libppo = ctypes.cdll.LoadLibrary(os.path.join(dn(dn(dn(os.path.abspath(__file__)))), "libppo" + suffix))

class _Animation(object):
  '''
  An animation class for occultation movies.

  '''

  def __init__(self, t, fig, axim, tracker, occ, ptb, body, bodies, occultors,
               interval = 50, gifname = None, xy = None, quiet = False):
    '''

    '''

    self.t = t
    self.fig = fig
    self.axim = axim
    self.tracker = tracker
    self.occ = occ
    self.ptb = ptb
    self.body = body
    self.bodies = bodies
    self.occultors = occultors
    self.xy = xy
    self.pause = True
    self.animation = animation.FuncAnimation(self.fig, self.animate, frames = 100,
                                             interval = interval, repeat = True)

    # Save?
    if gifname is not None:
      self.pause = False
      if not gifname.endswith('.gif'):
        gifname += '.gif'
      if not quiet:
        print("Saving %s..." % gifname)
      self.animation.save(gifname, writer = 'imagemagick', fps = 20, dpi = 150)
      self.pause = True

    # Toggle button
    self.axbutton = pl.axes([0.185, 0.12, 0.1, 0.03])
    self.button = Button(self.axbutton, 'Play', color = 'lightgray')
    self.button.on_clicked(self.toggle)

  def toggle(self, event):
    '''
    Pause/play the animation. Unfortunately I haven't been able to figure out
    how to freeze the animation so that it resumes at the same frame it was paused
    on...

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
        r = occultor._r

        if self.xy is None:
          xo, yo = occultor.x_hr[self.t[j]] - x0, occultor.y_hr[self.t[j]] - y0
        else:
          xo, yo = self.xy(occultor.x_hr[self.t[j]] - x0, occultor.y_hr[self.t[j]] - y0)

        self.occ[k].center = (xo / self.body._r, yo / self.body._r)

      # BODY orbits
      for k, b in enumerate(self.bodies):
        self.ptb[k].set_xdata(b.x_hr[self.t[j]])
        self.ptb[k].set_ydata(b.z_hr[self.t[j]])

class System(object):
  '''

  The planetary system class. This is the main interface to the photodynamical core.
  Instantiate with all bodies in the system and the desired settings, passed as :py:obj:`kwargs`.

  :param bodies: Any number of :py:func:`Planet`, :py:func:`Moon`, or :py:func:`Star` instances \
         comprising all the bodies in the system. The first body is assumed to be the primary.
  :param bool nbody: Uses the :py:obj:`REBOUND` N-body code to compute orbits. Default :py:obj:`False`
  :param float keptol: Kepler solver tolerance. Default `1.e-15`
  :param int maxkepiter: Maximum number of Kepler solver iterations. Default `100`
  :param str kepsolver: Kepler solver (`newton` | `mdfast`). Default `newton`
  :param float timestep: Timestep in days for the N-body solver. Default `0.01`
  :param bool adaptive: Adaptive grid for limb-darkened bodies? Default :py:obj:`True`
  :param bool quiet: Suppress output? Default :py:obj:`False`
  :param float mintheta: Absolute value of the minimum phase angle in degrees. Below this \
         angle, elliptical boundaries of constant surface brightness on the planet surface are \
         treated as vertical lines. Default `1.`
  :param int maxvertices: Maximum number of vertices allowed in the area computation. Default `999`
  :param int maxfunctions: Maximum number of functions allowed in the area computation. Default `999`
  :param int oversample: Oversampling factor for each exposure. Default `1`
  :param float distance: Distance to the system in parsecs. Default `10.`
  :param bool circleopt: Solve the simpler quadratic problem for circle-ellipse intersections when \
         the axes of the ellipse are equal to within :math:`10^{-10}`? Default :py:obj:`True`
  :param bool batmanopt: Use the :py:mod:`batman` algorithm to compute light curves of radially \
         symmetric bodies? This can significantly speed up the code. Default :py:obj:`True`
  :param str integrator: The N-body integrator (:py:obj:`whfast` | :py:obj:`ias15`) to use. Default :py:obj:`whfast`

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
    Applies all settings, makes bodies accessible as properties, resets flags,
    and computes some preliminary orbital information for the system.

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

    # Compute the semi-major axis for each planet/moon (in Earth radii)
    for body in self.bodies:
      body._computed = False
      body._host = 0
      if body.body_type == 'planet':
        body.a = ((body.per * DAYSEC) ** 2 * G * (self.primary._m + body._m) * MEARTH / (4 * np.pi ** 2)) ** (1. / 3.) / REARTH
      elif body.body_type == 'moon':
        if type(body.host) is str:
          body.host = getattr(self, body.host)
        body._host = np.argmax(self._names == body.host.name)
        body.a = ((body.per * DAYSEC) ** 2 * G * (body.host._m + body._m) * MEARTH / (4 * np.pi ** 2)) ** (1. / 3.) / REARTH

    # Reset animations
    self._animations = []

    # Reset flag
    self._computed = False

  def _malloc(self, nt, nw):
    '''
    Allocate memory for all the C arrays.

    '''

    # Initialize the C interface
    self._Orbits = libppo.Orbits
    self._Orbits.restype = ctypes.c_int
    self._Orbits.argtypes = [ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nt),
                             ctypes.c_int, ctypes.POINTER(ctypes.POINTER(BODY)),
                             SETTINGS]

    self._Flux = libppo.Flux
    self._Flux.restype = ctypes.c_int
    self._Flux.argtypes = [ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nt),
                           ctypes.c_int, ctypes.ARRAY(ctypes.c_double, nw),
                           ctypes.c_int, ctypes.POINTER(ctypes.POINTER(BODY)),
                           SETTINGS]

    # Allocate memory for all the arrays
    for body in self.bodies:
      body.time = np.zeros(nt)
      body._time = np.ctypeslib.as_ctypes(body.time)
      body.wavelength = np.zeros(nw)
      body._wavelength = np.ctypeslib.as_ctypes(body.wavelength)
      body.x = np.zeros(nt)
      body._x = np.ctypeslib.as_ctypes(body.x)
      body.y = np.zeros(nt)
      body._y = np.ctypeslib.as_ctypes(body.y)
      body.z = np.zeros(nt)
      body._z = np.ctypeslib.as_ctypes(body.z)
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
      
      # Radiance map
      if body.radiancemap is None:
        body._custommap = 0
        body._radiancemap = NullMap().ctypes
      else:
        assert type(body.radiancemap) is numba.ccallback.CFunc, "Parameter `radiancemap` must be a NUMBA C callback function."
        assert body.radiancemap.ctypes.argtypes == (ctypes.c_double, ctypes.c_double), "The `radiancemap` callback must accept two C doubles."
        assert body.radiancemap.ctypes.restype == ctypes.c_double, "The `radiancemap` callback must return a C double."
        body._custommap = 1
        body._radiancemap = body.radiancemap.ctypes
      
      # Dimensions
      body.nu = len(body.u)
      body.nt = nt
      body.nw = nw

    # A pointer to a pointer to `BODY`. This is an array of `n` `BODY` instances,
    # passed by reference. The contents can all be accessed through `bodies`
    self._ptr_bodies = (ctypes.POINTER(BODY) * len(self.bodies))(*[ctypes.pointer(p) for p in self.bodies])

  def compute(self, time, lambda1 = 5, lambda2 = 15, R = 100):
    '''
    Compute the full system light curve over a given `time` array between
    wavelengths `lambda1` and `lambda2` at resolution `R`. This method runs the
    photodynamical core, populates all bodies with their individual light curves,
    and flags all occultation events.

    :param array_like time: The times at which to evaluate the light curve in days
    :param float lambda1: The start point of the wavelength grid in microns. Default `5`
    :param float lambda2: The end point of the wavelength grid in microns. Default `15`
    :param float R: The spectrum resolution, :math:`R = \\frac{\lambda}{\Delta\lambda}`. Default `100`

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
          raise Exception("Limb darkening coefficients must be provided as a list of scalars or as a list of functions.")
      body.u = np.array(body.u)

    # Convert from microns to meters
    wavelength *= 1e-6

    # Oversample the time array to generate light curves observed w/ finite exposure time.
    # Ensure `oversample` is odd.
    if self.settings.oversample % 2 == 0:
      oversample = self.settings.oversample + 1
    else:
      oversample = self.settings.oversample
    if oversample == 1:
      time_hr = np.array(time)
    else:
      tmid = (time[:-1] + time[1:]) / 2.
      tmid = np.concatenate(([tmid[0] - (tmid[1] - tmid[0])], tmid, [tmid[-1] + (tmid[-1] - tmid[-2])]))
      time_hr = [np.linspace(tmid[i], tmid[i + 1], oversample) for i in range(len(tmid) - 1)]
      time_hr = np.concatenate(time_hr)

    # Allocate memory
    self._malloc(len(time_hr), len(wavelength))

    # Call the light curve routine
    err = self._Flux(len(time_hr), np.ctypeslib.as_ctypes(time_hr), len(wavelength),
                     np.ctypeslib.as_ctypes(wavelength), len(self.bodies),
                     self._ptr_bodies, self.settings)
    assert err <= 0, "Error in C routine `Flux` (%d)." % err
    self._computed = True

    # Downbin to original time array
    for body in self.bodies:

      # Store the oversampled light curve
      body.time_hr = np.array(body.time)
      body.flux_hr = np.array(body.flux)
      body.occultor_hr = np.array(body.occultor)
      body.x_hr = np.array(body.x)
      body.y_hr = np.array(body.y)
      body.z_hr = np.array(body.z)

      # Store the binned light curve
      if self.settings.oversample > 1:

        # Revert to original time array
        body.time = time

        # Average the flux over the exposure
        body.flux = np.mean(body.flux.reshape(-1, oversample, len(wavelength)), axis = 1)

        # Take the XYZ position at the bin center
        body.x = body.x[oversample // 2::oversample]
        body.y = body.y[oversample // 2::oversample]
        body.z = body.z[oversample // 2::oversample]

        # Get all bodies that occult at some point over the exposure
        # The operation `bitwise_or.reduce` is *magical*
        body.occultor = np.bitwise_or.reduce(body.occultor.reshape(-1, oversample), axis = 1)

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
            i = np.array(list(range(ia, i[0])) + list(i) + list(range(i[-1] + 1, ib)))

            # Store it
            body._inds.append(i)

  def compute_orbits(self, time):
    '''
    Run the dynamical code to compute the positions of all bodies over a given `time`
    array. This method does not compute light curves, but it does check for occultations.

    :param array_like time: The times at which to store the body positions in days

    '''

    # Reset
    self._reset()
    for body in self.bodies:
      body.u = np.array([], dtype = float)

    self._malloc(len(time), 1)

    # Call the light curve routine
    err = self._Orbits(len(time), np.ctypeslib.as_ctypes(time), len(self.bodies),
                       self._ptr_bodies, self.settings)
    assert err <= 0, "Error in C routine `Orbits` (%d)." % err

  def observe(self, save = None, filter = 'f1500w', stack = 1, instrument = 'jwst'):
    '''
    Run telescope observability calculations for a system that has had its lightcurve
    computed. Calculates a synthetic noised lightcurve in the user specified filter. 

    :param bool save: Save a text file and a plot. Default :py:obj:`None`
    :param filter: Filter name or :py:func:`Filter` object. Default `'f1500w'`
    :type filter: str or :py:func:`Filter`
    :param int stack: Number of exposures to stack. Default `1`
    :param str instrument: Telescope instrument to use. Default `'jwst'`
    
    '''

    # Have we computed the light curves?
    assert self._computed, "Please run `compute()` first."

    # Get the filter "wheel"
    if instrument.lower() == 'jwst':
      wheel = jwst.get_miri_filter_wheel()
      atel  = 25.0
      thermal = True
    elif instrument.lower() == 'ost':
      # Using JWST filters for OST for now
      wheel = jwst.get_miri_filter_wheel()
      # create new filter for 50 microns
      filt50 = jwst.create_tophat_filter(47.5, 52.5, dlam=0.1, Tput=0.3, name="OST50")
      # Will have mirror between 8-15m, let's say: 12m
      atel = 144.0
      # No thermal noise
      thermal = False
    elif instrument.lower() == 'spitzer':
      wheel = jwst.get_spitzer_filter_wheel()
      atel = np.pi * (0.85/2.)**2
      thermal = False
    else:
      raise ValueError("Invalid instrument.")

    # Get the filter
    if type(filter).__name__ == "Filter":
      self.filter = filter
    elif filter.lower() in [f.name.lower() for f in wheel]:
      self.filter = wheel[np.argmax([f.name.lower() == filter.lower() for f in wheel])]
    else:
      raise ValueError("Invalid filter.")

    # Compute lightcurve in this filter
    self.filter.compute_lightcurve(self.flux, self.time, self.wavelength, stack = stack,
                                   time_hr = self.time_hr, flux_hr = self.flux_hr,
                                   atel = atel, thermal = thermal)

    # Setup plot
    fig, ax = pl.subplots(figsize = (12,4))
    ax.set_title(r"%s" % self.filter.name)
    ax.set_ylabel("Relative Flux", fontweight = 'bold', fontsize = 10)
    ax.set_xlabel("Time [days]", fontweight = 'bold', fontsize = 10)

    # Plot lightcurve
    self.filter.lightcurve.plot(ax0 = ax)

    # Determine average precision attained in lightcurve
    ppm = np.mean(self.filter.lightcurve.sig/self.filter.lightcurve.obs) * 1e6
    print("Average lightcurve precision: %.3f ppm" %ppm)

    """
    Determine SNR on each event:
    """
    # Use polynomial fit to continuum
    method = "poly"

    # Mask in and out of occultation times
    inmask = np.logical_or.reduce([planet.occultor > 0 for planet in self.bodies])
    outmask = ~inmask
    arr = np.arange(len(inmask))

    # Determine continuum flux if there were no event
    if method == "interp":
        # Interpolation
        x = self.filter.lightcurve.time[inmask]
        xp = self.filter.lightcurve.time[outmask]
        fp = self.filter.lightcurve.obs[outmask]
        vals = np.interp(x, xp, fp)
    elif method == "poly":
        # Poly
        deg = 3
        # Compute for plot
        x = self.filter.lightcurve.time[inmask]
        xp = self.filter.lightcurve.time[outmask]
        yp = self.filter.lightcurve.obs[outmask]
        zp = np.polyfit(xp,yp,deg)
        p1 = np.poly1d(zp)
        vals = p1(x)
        # Compute for photons
        xp = self.filter.lightcurve.time[outmask]
        yp = self.filter.lightcurve.Nphot[outmask]
        zp = np.polyfit(xp,yp,deg)
        p2 = np.poly1d(zp)

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
        mask = events[i]
        # System photons over event
        Nsys = self.filter.lightcurve.Nphot[mask]
        # Background photons over event
        Nback = self.filter.lightcurve.Nback[mask]
        if method == "interp":
            # Interpolated continuum photons over event
            Ncont = np.interp(self.filter.lightcurve.time[mask],
                                 self.filter.lightcurve.time[outmask],
                                 self.filter.lightcurve.Nphot[outmask])
        elif method == "poly":
            # Polynomial fit to continuum over event
            Ncont = p2(self.filter.lightcurve.time[mask])
        # Compute signal of and noise on the event
        signal.append(np.sum(np.fabs(Ncont - Nsys))/np.sum(Nsys) * 1e6)
        noise.append(np.sqrt(np.sum(Nsys + Nback))/np.sum(Nsys) * 1e6)
        #imax = np.argmax(Ncont - Nsys)
        #amp.append(np.sum(Ncont - Nsys)/np.sum(Ncont)*1e6)
        # Compute SNR on event
        SNR = np.sum(np.fabs(Ncont - Nsys)) / np.sqrt(np.sum(Nsys + Nback))
        SNRs.append(SNR)

    # Annotate event SNRs on each event
    for i in range(n):
        k = np.argmin(self.filter.lightcurve.obs[np.array(events[i])])
        j = events[i][k]
        ax.annotate(r"$%.1f \sigma$" %SNRs[i], xy = (self.filter.lightcurve.time[int(np.mean(events[i]))],
                   self.filter.lightcurve.obs[j]-1*self.filter.lightcurve.sig[j]), ha = 'center',
                   va = 'center', color = "black", fontweight = 'bold',
                   fontsize = 10, xytext = (0, -15), textcoords = 'offset points')

    # Save SNR and ampl as attributes in Lightcurve
    self.filter.lightcurve.event_SNRs = SNRs
    self.filter.lightcurve.event_signal = signal
    self.filter.lightcurve.event_noise = noise

    # Save to disk?
    if save is not None:

      # Compose data array to save
      data = np.array([self.filter.lightcurve.time, self.filter.lightcurve.obs, self.filter.lightcurve.sig]).T

      # Did the user provide a file name?
      if save is True:

        # No
        np.savetxt("jwst_lc_%s.txt" % self.filter.name,
                   data, fmt = str("%.6e"),
                   header = "time [days]      flux         error",
                   comments = "")
        fig.savefig("jwst_lc_%s.png" % self.filter.name, bbox_inches="tight")

      else:

        # Yes
        np.savetxt("%s.txt" % save,
                   data, fmt = str("%.6e"),
                   header = "time [days]      flux         error",
                   comments = "")
        fig.savefig("%s.png" % save, bbox_inches="tight")

    # Adjust for on-screen viewing
    fig.subplots_adjust(bottom = 0.2)

    return fig, ax

  def scatter_plot(self, tstart, tend, dt = 0.001):
    '''

    Compute all occultations between `tstart` and `tend` and plot an occultation
    scatter plot like the one in the paper.

    .. plot::
       :align: center

       from planetplanet.photo import Trappist1
       import matplotlib.pyplot as pl
       system = Trappist1(sample = True, nbody = True)
       system.scatter_plot(0, 365 * 3)
       pl.show()

    :param float tstart: The integration start time in days
    :param float tend: The integration end time in days
    :param float dt: The time resolution in days. Occultations shorter than this \
           will not be registered.

    '''

    # Reset
    self._reset()
    for body in self.bodies:
      body.u = np.array([], dtype = float)
    time = np.arange(tstart, tend, dt)
    self._malloc(len(time), 1)

    # Call the orbit routine
    err = self._Orbits(len(time), np.ctypeslib.as_ctypes(time), len(self.bodies),
                       self._ptr_bodies, self.settings)
    assert err <= 0, "Error in C routine `Orbits` (%d)." % err

    # Loop over all bodies and plot each occultation event as a circle
    nppo = 0
    figp, axp = pl.subplots(1, figsize = (8,8))
    figp.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.1, top = 0.9)
    axp.axis('off')
    for body in self.bodies[1:]:

      # Identify the different events
      inds = np.where(body.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]

      # Plot the orbit outline
      f = np.linspace(0, 2 * np.pi, 1000)
      r = body.a * (1 - body.ecc ** 2) / (1 + body.ecc * np.cos(f))
      x = r * np.cos(body._w + f) - r * np.sin(body._w + f) * np.cos(body._inc) * np.sin(body._Omega)
      z = r * np.sin(body._w + f) * np.sin(body._inc)
      axp.plot(x, z, 'k-', lw = 1, alpha = 0.05)
      n = np.argmin(1e10 * (x < 0) + np.abs(z))
      axp.annotate(body.name, xy = (x[n], z[n]), color = 'k', alpha = 0.2, fontweight = 'bold',
                   fontsize = 8, zorder = -99, ha = 'center', va = 'center')

      # Loop over individual ones
      plot_secondary = True
      for i in inds[difs]:

        # Loop over all possible occultors
        for occ in range(len(self.bodies)):

          # Is body `occ` occulting?
          if (body.occultor[i] & 2 ** occ):

            # Note that `i` is the last index of the occultation
            duration = np.argmax(body.occultor[:i][::-1] & 2 ** occ == 0)

            # Compute the minimum impact parameter
            idx = range(i - duration, i + 1)
            impact = np.min(np.sqrt((self.bodies[occ].x[idx] - body.x[idx]) ** 2 +
                                    (self.bodies[occ].y[idx] - body.y[idx]) ** 2)) / (self.bodies[occ]._r + body._r)

            # Transparency proportional to the impact parameter
            alpha = 0.8 * (1 - impact) + 0.01

            # Size = duration in minutes / 3
            ms = duration * dt * 1440 / 3

            # If the occultor is the star, plot it only once
            if (occ == 0):
              if plot_secondary:
                axp.plot(body.x[i], body.z[i], 'o', color = self.colors[occ], alpha = alpha, ms = ms, markeredgecolor = 'none')
                plot_secondary = False
            else:
              axp.plot(body.x[i], body.z[i], 'o', color = self.colors[occ], alpha = alpha, ms = ms, markeredgecolor = 'none')
              nppo += 1

        # Check for mutual transits
        if self.bodies[0].occultor[i]:

          # Get all bodies currently occulting the star
          occultors = []
          for occ in range(1, len(self.bodies)):
            if (self.bodies[0].occultor[i] & 2 ** occ):
              occultors.append(occ)

          # Check if any of these occult each other
          for occ1 in occultors:
            for occ2 in occultors:
              if self.bodies[occ1].occultor[i] & 2 ** occ2:
                axp.plot(self.bodies[occ1].x[i], self.bodies[occ1].z[i], 'x',
                         color = self.colors[occ2], alpha = 1, zorder = 100,
                         ms = 20)

    # Legend 1: Occultor names/colors
    axl1 = pl.axes([0.025, 0.775, 0.2, 0.2])
    axl1.axis('off')
    axl1.set_xlim(-0.5, 1.5)
    axl1.set_ylim(-len(self.bodies) // 2 - 1, 1.5)
    axl1.annotate('Occultations by', xy = (0.5, 1), ha = 'center', va = 'center', fontweight = 'bold')
    for j, body in enumerate(self.bodies):
      if j < len(self.bodies) // 2:
        x, y = (0, -j)
      else:
        x, y = (0.825, len(self.bodies) // 2 - j)
      axl1.plot(x, y, 'o', color = self.colors[j], ms = 6, alpha = 1, markeredgecolor = 'none')
      axl1.annotate(body.name, xy = (x + 0.1, y), xycoords = 'data',
                    ha = 'left', va = 'center', color = self.colors[j])

    # Legend 2: Size/duration
    axl2 = pl.axes([0.775, 0.775, 0.2, 0.2])
    axl2.axis('off')
    axl2.set_xlim(-1, 1)
    axl2.set_ylim(-3, 1.5)
    axl2.annotate('Duration', xy = (0., 1), ha = 'center', va = 'center', fontweight = 'bold')
    for j, duration in enumerate([10, 30, 60]):
      ms = duration / 3.
      axl2.plot(-0.65, -0.75 * j + 0.2, 'o', color = 'k', ms = ms, alpha = 0.65, markeredgecolor = 'none')
      axl2.annotate('%d minutes' % duration, xy = (-0.3, -0.75 * j + 0.2), xycoords = 'data',
                    ha = 'left', va = 'center', color = 'k')

    # Legend 3: Transparency/impact parameter
    axl3 = pl.axes([0.025, 0.025, 0.2, 0.2])
    axl3.axis('off')
    axl3.set_xlim(-0.5, 1.5)
    axl3.set_ylim(-2, 1.5)
    axl3.annotate('Impact parameter', xy = (0.5, 0.65), ha = 'center', va = 'center', fontweight = 'bold')
    for j, impact in enumerate([0, 0.2, 0.4]):
      alpha = 0.8 * (1 - impact) + 0.01
      axl3.plot(-0.15, -0.75 * j, 'o', color = 'k', ms = 8, alpha = alpha, markeredgecolor = 'none')
      axl3.annotate('%.1f' % impact, xy = (-0.05, -0.75 * j), xycoords = 'data',
                    ha = 'left', va = 'center', color = 'k')
    for j, impact in enumerate([0.6, 0.8, 1.]):
      alpha = 1 - impact
      axl3.plot(0.675, -0.75 * j, 'o', color = 'k', ms = 8, alpha = alpha, markeredgecolor = 'none')
      axl3.annotate('%.1f' % impact, xy = (0.775, -0.75 * j), xycoords = 'data',
                    ha = 'left', va = 'center', color = 'k')

    # Observer direction
    axp.annotate("To observer", xy = (0.5, -0.1), xycoords = "axes fraction", xytext = (0, 30),
                 ha = 'center', va = 'center', annotation_clip = False, color = 'cornflowerblue',
                 textcoords = "offset points", arrowprops=dict(arrowstyle = "-|>", color = 'cornflowerblue'))

    # Log
    if not self.settings.quiet:
      print("There were %d PPOs between t = %.2f and t = %.2f." % (nppo, tstart, tend))

    return figp

  def histogram(self, tstart, tend, dt = 0.0001):
    '''

    Computes statistical properties of planet-planet occultations that occur over a given interval.
    Computes the frequency of occultations as a function of orbital phase, duration,
    and impact parameter, as well as the fully marginalized occultation frequency for
    each planet in the system. Occultations by the star are not included, nor are occultations
    occuring behind the star, which are not visible to the observer.

    .. plot::
       :align: center

       from scripts import hist
       import matplotlib.pyplot as pl
       figs = hist.Plot()
       for fig in figs[:-1]:
         pl.close(fig)
       pl.show()

    :param float tstart: The integration start time in days
    :param float tend: The integration end time in days
    :param float dt: The time resolution in days. Occultations shorter than this \
           will not be registered.

    :returns: A list of :py:obj:`(phase angle, impact parameter, duration)` tuples for each planet in \
              the system. The phase angle is measured in degrees and the duration is measured in days.

    .. warning:: This routine computes the **orbital phase angle**, which is measured from **transit**. \
              This is different from the mean longitude by :math:`\pi/2`

    '''

    # Reset
    self._reset()
    for body in self.bodies:
      body.u = np.array([], dtype = float)
    time = np.arange(tstart, tend, dt)
    self._malloc(len(time), 1)

    # Call the orbit routine
    err = self._Orbits(len(time), np.ctypeslib.as_ctypes(time), len(self.bodies),
                       self._ptr_bodies, self.settings)
    assert err <= 0, "Error in C routine `Orbits` (%d)." % err

    # A histogram of the distribution of phases, impact parameters, and durations
    hist = [[] for body in self.bodies[1:]]
    for k, body in enumerate(self.bodies[1:]):

      # Identify the different planet-planet events
      inds = np.where(body.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]

      # Loop over individual ones
      for i in inds[difs]:

        # Loop over possible occultors
        for occ in range(1, len(self.bodies)):

          # Is body `occ` occulting (but not behind the star)?
          if (body.occultor[i] & 2 ** occ) and (body.occultor[i] & 1 == 0):

            # Note that `i` is the last index of the occultation
            duration = np.argmax(body.occultor[:i][::-1] & 2 ** occ == 0)
            if duration > 0:

              # Orbital phase, **measured from transit**
              # At transit, phase = 0; at secondary, phase = 180.
              phase = np.arctan2(body.x[i], -body.z[i]) * 180 / np.pi

              # Compute the minimum impact parameter
              idx = range(i - duration, i + 1)
              impact = np.min(np.sqrt((self.bodies[occ].x[idx] - body.x[idx]) ** 2 +
                                      (self.bodies[occ].y[idx] - body.y[idx]) ** 2)) / (self.bodies[occ]._r + body._r)

              # Convert duration to log
              duration = np.log10(duration * dt * 1440)

              # Running list
              hist[k].append((phase, impact, duration))

      # Make into array
      hist[k] = np.array(hist[k])

    return hist

  def next_occultation(self, tstart, occulted, min_duration = 10, max_impact = 0.5, occultor = None, maxruns = 100, dt = 0.001):
    '''
    Computes the time of the next occultation of body `occulted`.

    :param float tstart: The time at which to start the occultation search (days)
    :param occulted: The occulted body
    :type occulted: :py:class:`BODY`
    :param occultor: The occultor(s). If :py:obj:`None`, occultations by any body are considered. Default :py:obj:`None`
    :type occultor: :py:class:`BODY` or :py:obj:`list`
    :param float min_duration: Minimum occultation duration in **minutes**. Default `10`
    :param float max_impact: The maximum impact parameter. Default `0.5`
    :param int maxruns: Maximum number of 100-day runs to search for the occultation. Default `100`
    :param float dt: The integration timestep

    :returns: The time of the occultation. If no occultation was found, returns :py:class:`np.nan`

    '''

    # Quiet?
    quiet = self.settings.quiet
    self.settings.quiet = True
    if quiet:
      iterator = lambda x: range(x)
    else:
      iterator = lambda x: tqdm(range(x))

    # Convert occultors to indices if necessary
    if occultor is None:
      occultor = list(range(1, len(self.bodies)))
    elif not hasattr(occultor, '__len__'):
      occultor = [occultor]
    for i, occ in enumerate(occultor):
      if occ in self.bodies:
        occultor[i] = np.argmax([b == occ for b in self.bodies])

    # Loop until we find one
    for n in iterator(maxruns):

      # Compute the orbits, 100 days at a time
      time = np.arange(tstart + n * 100., tstart + (n + 1) * 100., dt)
      self.compute_orbits(time)

      # Identify the different planet-planet events
      inds = np.where(occulted.occultor > 0)[0]
      difs = np.where(np.diff(inds) > 1)[0]

      # Loop over individual ones
      for i in inds[difs]:

        # Loop over possible planet occultors
        for occ in occultor:

          # Is body `occ` occulting (but not behind the star)?
          if (occulted.occultor[i] & 2 ** occ) and ((occ == 0) or (occulted.occultor[i] & 1 == 0)):

            # Note that `i` is the last index of the occultation
            duration = np.argmax(occulted.occultor[:i][::-1] & 2 ** occ == 0)
            if duration * dt * 1440 >= min_duration:

              # Compute the minimum impact parameter
              idx = range(i - duration, i + 1)
              b = np.sqrt((self.bodies[occ].x[idx] - occulted.x[idx]) ** 2 + (self.bodies[occ].y[idx] - occulted.y[idx]) ** 2) / (self.bodies[occ]._r + occulted._r)
              ind = np.argmin(b)
              if b[ind] <= max_impact:
                if not quiet:
                  print("The next occultation of %s by %s is at t = %.2f days." % (occulted.name, self.bodies[occ].name, time[idx[ind]]))
                self.settings.quiet = quiet
                return time[idx[ind]]

    # Nothing found...
    if not quiet:
      print("No occultation of %s by %s found." % (occulted.name, self.bodies[occ].name))
    self.settings.quiet = quiet
    return np.nan

  def plot_occultation(self, body, time, wavelength = 15., interval = 50, gifname = None, spectral = False, **kwargs):
    '''
    Plots and animates an occultation event.
    
    :param body: The occulted body
    :type body: :py:class:`BODY`
    :param float time: The time of the occultation event in days
    :param float wavelength: The wavelength in microns at which to plot the light curve. \
           Must be within the wavelength grid. Default `15`
    :param int interval: The interval between frames in the animation in ms. Default `50`
    :param str gifname: If specified, saves the occultation animation as a :py:obj:`gif` in the current directory. Default :py:obj:`None`
    :param bool spectral: Plot the light curve at different wavelengths? If :py:obj:`True`, plots the first, middle, \
           and last wavelength in the wavelength grid. If :py:obj:`False`, plots the specified `wavelength`. Default :py:obj:`True`
    :param kwargs: Any additional keyword arguments to be passed to :py:func:`plot_image`

    :returns: :py:obj:`(fig, ax1, ax2, ax3)`, a figure instance and its three axes

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
    assert (wavelength >= self.bodies[0].wavelength[0]) and (wavelength >= self.bodies[0].wavelength[-1]), "Wavelength value outside of computed grid."
    w = np.argmax(1e-6 * wavelength <= self.bodies[0].wavelength)
    
    # Check file name
    if gifname is not None:
      if gifname.endswith(".gif"):
        gifname = gifname[:-4]

    # Get the occulted body
    if type(body) is str:
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
    
    # Plot three different wavelengths (first, mid, and last)
    if spectral:
      fluxb = body.flux_hr[t, 0] / body.total_flux[0]
      fluxg = body.flux_hr[t, body.flux_hr.shape[-1] // 2] / body.total_flux[body.flux.shape[-1] // 2]
      fluxr = body.flux_hr[t, -1] / body.total_flux[-1]

      # Add a baseline?
      if not (body.phasecurve or body.body_type == 'star'):
        fluxg += 1
        fluxr += 1
        fluxb += 1

      # Plot
      axlc.plot(body.time_hr[t], fluxb, 'b-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[0])) + r"\ \mu\mathrm{m}$")
      axlc.plot(body.time_hr[t], fluxg, 'g-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[body.flux.shape[-1] // 2])) + r"\ \mu\mathrm{m}$")
      axlc.plot(body.time_hr[t], fluxr, 'r-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[-1])) + r"\ \mu\mathrm{m}$")
    else:
      
      fluxw = body.flux_hr[t, w] / body.total_flux[w]
      if not (body.phasecurve or body.body_type == 'star'):
        fluxw += 1
      axlc.plot(body.time_hr[t], fluxw, 'k-', label = r"$" + '{:.4s}'.format('{:0.2f}'.format(1e6 * body.wavelength[w])) + r"\ \mu\mathrm{m}$")

    axlc.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    axlc.set_ylabel(r'Normalized Flux', fontweight = 'bold', fontsize = 10)
    axlc.get_yaxis().set_major_locator(MaxNLocator(4))
    axlc.get_xaxis().set_major_locator(MaxNLocator(4))
    tracker = axlc.axvline(body.time_hr[ti], color = 'k', alpha = 0.5, lw = 1, ls = '--')
    for tick in axlc.get_xticklabels() + axlc.get_yticklabels():
      tick.set_fontsize(8)
    axlc.legend(loc = 'lower right', fontsize = 8)
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

    # Plot the orbits of all bodies except moons
    # NOTE: These are plotted from a top-down view for each planet
    # individually. This places non-coplanar planets in the same
    # orbital plane so that their orbits can be better visualized.
    # In the future, we should show a top-down and edge-on view
    # in fixed reference frames.
    axxz = pl.subplot2grid((5, 3), (0, 0), colspan = 3, rowspan = 2)
    f = np.linspace(0, 2 * np.pi, 1000)
    for j, b in enumerate(self.bodies):
      if b.body_type == 'moon':
        continue
      if j == p:
        style = dict(color = 'r', alpha = 1, ls = '-', lw = 1)
      elif j in occultors:
        style = dict(color = 'k', alpha = 1, ls = '-', lw = 1)
      else:
        style = dict(color = 'k', alpha = 0.1, ls = '--', lw = 1)
      r = b.a * (1 - b.ecc ** 2) / (1 + b.ecc * np.cos(f))
      x = r * np.cos(b._w + f) - r * np.sin(b._w + f) * np.cos(b._inc) * np.sin(b._Omega)
      z = r * np.sin(b._w + f) * np.sin(b._inc)
      axxz.plot(x, z, **style)

    # Plot the locations of the bodies
    ptb = [None for b in self.bodies]
    for bi, b in enumerate(self.bodies):
      if b == body:
        ptb[bi], = axxz.plot(b.x_hr[ti], b.z_hr[ti], 'o', color = 'r', alpha = 1, markeredgecolor = 'k', zorder = 99)
      elif bi in occultors:
        ptb[bi], = axxz.plot(b.x_hr[ti], b.z_hr[ti], 'o', color = 'lightgrey', alpha = 1, markeredgecolor = 'k', zorder = 99)
      else:
        ptb[bi], = axxz.plot(b.x_hr[ti], b.z_hr[ti], 'o', color = '#dddddd', alpha = 1, markeredgecolor = '#999999', zorder = 99)

    # Appearance
    axxz.set_ylim(-max(np.abs(axxz.get_ylim())), max(np.abs(axxz.get_ylim())))
    axxz.set_xlim(-max(np.abs(axxz.get_xlim())), max(np.abs(axxz.get_xlim())))
    axxz.set_aspect('equal')
    axxz.axis('off')

    # Plot the image
    axim, occ, xy = self.plot_image(ti, body, occultors, fig = fig, wavelength = wavelength, **kwargs)
    bodies = [self.bodies[o] for o in occultors] + [body]
    xmin = min(np.concatenate([o.x_hr[t] - body.x_hr[t] - 1.1 * o._r for o in bodies]))
    xmax = max(np.concatenate([o.x_hr[t] - body.x_hr[t] + 1.1 * o._r for o in bodies]))
    ymin = min(np.concatenate([o.y_hr[t] - body.y_hr[t] - 1.1 * o._r for o in bodies]))
    ymax = max(np.concatenate([o.y_hr[t] - body.y_hr[t] + 1.1 * o._r for o in bodies]))
    axim.set_xlim(xmin, xmax)
    axim.set_ylim(ymin, ymax)
    axim.axis('off')
    axim.set_aspect('equal')

    # The title
    if len(occultors) == 1:
      axxz.annotate("%s occulted by %s" % (body.name, self.bodies[occultors[0]].name), xy = (0.5, 1.25),
                       xycoords = "axes fraction", ha = 'center', va = 'center',
                       fontweight = 'bold', fontsize = 12)
    else:
      axxz.annotate("%s occulted by %s" % (body.name,
                      ", ".join([occultor.name for occultor in [self.bodies[o] for o in occultors]])),
                       xy = (0.5, 1.25),
                       xycoords = "axes fraction", ha = 'center', va = 'center',
                       fontweight = 'bold', fontsize = 12)
    axxz.annotate("Duration: %.2f minutes" % ((body.time_hr[tf] - body.time_hr[ti]) * 1440.),
                     xy = (0.5, 1.1), ha = 'center', va = 'center', xycoords = 'axes fraction',
                     fontsize = 10, style = 'italic')

    # Animate!
    if gifname is not None:
      tmp = '%s.%03d.gif' % (gifname, len(self._animations) + 1)
    else:
      tmp = None
    self._animations.append(_Animation(range(ti, tf), fig, axim, tracker, occ, ptb, body,
                            self.bodies, [self.bodies[o] for o in occultors],
                            interval = interval, gifname = tmp, quiet = self.settings.quiet,
                            xy = xy))

    return fig, axlc, axxz, axim

  def plot_image(self, t, occulted, occultor = None, wavelength = 15., fig = None, figx = 0.535, figy = 0.5, figr = 0.05, **kwargs):
    '''
    Plots an image of the `occulted` body and its occultor(s) at a given index of the `time_hr` array `t`.

    .. plot::
       :align: center

       from planetplanet.photo.eyeball import DrawEyeball
       import matplotlib.pyplot as pl
       DrawEyeball(0.535, 0.5, 0.25, theta = 1., gamma = 1.,
                   occultors = [{'x': 0.2, 'y': 0.4, 'r': 0.5}], cmap = 'inferno')
       pl.show()

    :param int t: The index of the occultation in the high resolution time array `time_hr`
    :param occulted: The occulted body instance
    :type occultor: :py:class:`BODY` or :py:obj:`list`
    :param occultor: The occultor(s). Default :py:obj:`None`
    :type occultor: :py:class:`BODY` or :py:obj:`list`
    :param float wavelength: The wavelength in microns at which to plot the light curve. \
           Must be within the wavelength grid. Default `15`
    :param fig: The figure on which to plot the image
    :type fig: :py:class:`matplotlib.Figure`
    :param float figx: The x coordinate of the image in figure units. Default `0.535`
    :param float figy: The y coordinate of the image in figure units. Default `0.5`
    :param float figr: The radius of the image in figure units. Default `0.05`

    '''
    
    # Set up the plot
    if fig is None:
      fig = pl.gcf()

    # Get the occulted body
    rp = occulted._r
    x0 = occulted.x_hr[t]
    y0 = occulted.y_hr[t]
    z0 = occulted.z_hr[t]

    # Get the angles
    if (occulted.nu == 0) and not (occulted.body_type in ['planet', 'moon'] and occulted.airless == False):

      x = x0 * np.cos(occulted._Omega) + y0 * np.sin(occulted._Omega)
      y = y0 * np.cos(occulted._Omega) - x0 * np.sin(occulted._Omega)
      z = z0
      r = np.sqrt(x ** 2 + y ** 2 + z ** 2)

      # Coordinates of the hotspot in a frame where the planet is
      # at x, y, z = (0, 0, r), at full phase
      xprime = occulted._r * np.cos(occulted._Phi) * np.sin(occulted._Lambda)
      yprime = occulted._r * np.sin(occulted._Phi)
      zprime = r - occulted._r * np.cos(occulted._Phi) * np.cos(occulted._Lambda)

      # Transform to the rotated sky plane
      rxz = np.sqrt(x ** 2 + z ** 2)
      xstar = ((z * r) * xprime - (x * y) * yprime + (x * rxz) * zprime) / (r * rxz)
      ystar = (rxz * yprime + y * zprime) / r
      zstar = (-(x * r) * xprime - (y * z) * yprime + (z * rxz) * zprime) / (r * rxz)

      # Transform back to the true sky plane
      xstar, ystar = xstar * np.cos(occulted._Omega) - ystar * np.sin(occulted._Omega), \
                     ystar * np.cos(occulted._Omega) + xstar * np.sin(occulted._Omega)
      x = x0
      y = y0

      # Distance from planet center to hotspot
      d = np.sqrt((xstar - x) ** 2 + (ystar - y) ** 2)

      # Get the rotation and phase angles
      gamma = np.arctan2(ystar - y, xstar - x) + np.pi
      if (zstar - z) <= 0:
        theta = np.arccos(d / occulted._r)
      else:
        theta = -np.arccos(d / occulted._r)

    else:

      theta = np.pi / 2
      gamma = 0

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
      occ_dict.append(dict(x = (occultor.x_hr[t] - x0) / rp, y = (occultor.y_hr[t] - y0) / rp, r = occultor._r / rp, zorder = i + 1, alpha = 1))

    # Draw the eyeball planet and the occultor(s)
    fig, ax, occ, xy = DrawEyeball(figx, figy, figr, theta = theta, gamma = gamma,
                                   occultors = occ_dict, cmap = 'inferno', fig = fig, 
                                   radiancemap = occulted.radiancemap, 
                                   wavelength = wavelength, **kwargs)

    return ax, occ, xy

  def plot_lightcurve(self, wavelength = 15.):
    '''

    Plot the light curve of the system after running :py:func:`compute`.

    .. plot::
       :align: center

       from scripts import lightcurve
       import matplotlib.pyplot as pl
       lightcurve.plot()
       pl.show()

    :param float wavelength: The wavelength in microns at which to plot the light curve. \
           Must be within the wavelength grid. Default `15`

    :returns: :py:obj:`(fig, ax)`

    '''

    # Have we computed the light curves?
    assert self._computed, "Please run `compute()` first."

    if not self.settings.quiet:
      print("Plotting the light curve...")

    # Plot
    fig, ax = pl.subplots(1, figsize = (12, 4))
    assert (wavelength >= self.bodies[0].wavelength[0]) and (wavelength >= self.bodies[0].wavelength[-1]), "Wavelength value outside of computed grid."
    w = np.argmax(1e-6 * wavelength <= self.bodies[0].wavelength)
    flux = self.flux[:,w] / np.nanmedian(self.flux[:,w])
    ax.plot(self.time, flux, 'k-', lw = 1)

    # Plot the hi resolution light curve
    flux_hr = self.flux_hr[:,w] / np.nanmedian(self.flux_hr[:,w])
    ax.plot(self.time_hr, flux_hr, 'k-', lw = 1, alpha = 0.3, picker = 10)
    fig.canvas.mpl_connect('pick_event', self._onpick)

    # Appearance
    ax.set_xlabel('Time [days]', fontweight = 'bold', fontsize = 10)
    ax.set_ylabel(r'Normalized Flux @ %.1f$\mathbf{\mu}$m' % wavelength, fontweight = 'bold', fontsize = 10)
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
        tend = t[0] + len(body.time_hr[t]) - np.argmax(body.occultor_hr[t][::-1] > 0)
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
        for n, occultor in enumerate([self.bodies[o] for o in occultors]):
          # Keep track of events so we don't double count them
          event_id = [body.name, occultor.name, tmid]
          if event_id not in events:
            # This keeps track of how many annotations we're trying to squeeze
            # into the same time bin. Add a vertical shift to avoid overlaps
            idx = int(100 * (tmin - body.time_hr[0]) / (body.time_hr[-1] - body.time_hr[0]))
            dy = -15 * (1 + nann[idx])
            nann[idx] += 1
            self._lc_annotations.append(ax.annotate("%s" % body.name, xy = (tmin, fmin), ha = 'center',
                                        va = 'center', color = occultor.color, fontweight = 'bold',
                                        fontsize = 10, xytext = (0, dy), textcoords = 'offset points', clip_on = True))
            events.append(event_id)

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
  def time(self):
    '''
    Time array in days.

    '''

    return self.bodies[0].time

  @property
  def flux_hr(self):
    '''
    The total flux of the system computed on a grid of high resolution time and wavelength.

    '''

    return np.sum([b.flux_hr for b in self.bodies], axis = 0)

  @property
  def time_hr(self):
    '''
    High-resolution time array in days.

    '''

    return self.bodies[0].time_hr

  @property
  def wavelength(self):
    '''
    Wavelength array in microns.

    '''

    return self.bodies[0].wavelength * 1.e6