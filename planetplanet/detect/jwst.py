#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
jwst.py |github|
----------------

A Python interface for simulating time-series filter photometry with the James
Webb Space Telescope (JWST).

  .. role:: raw-html(raw)
     :format: html

  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/detect/jwst.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from __future__ import division, print_function, absolute_import, \
                       unicode_literals
import numpy as np
import matplotlib as mpl
from matplotlib import rc
import matplotlib.pyplot as plt
import os
import astropy.units as u

__all__ = [
    "Filter",
    "estimate_eclipse_snr",
    "create_tophat_filter",
    "get_spitzer_filter_wheel",
    "get_miri_filter_wheel"
]

HERE = os.path.abspath(os.path.split(__file__)[0])
MIRI_FILTERS = "filter_files/miri_filter.csv"
MIRI_PATH = os.path.join(HERE,MIRI_FILTERS)

def readin_miri_filters(f=MIRI_PATH):
    '''
    Function to read-in MIRI photometic filters

    :param str f: File location/name of JWST/MIRI photometric filter file. \
                  Default is relative to file.

    :returns: array_like wl_filt: Wavelength grid for filter throughputs \
              [:math:`\mu \mathrm{m}`]
    :returns: array_like dwl_filt: Wavelength bin widths for filter \
              throughputs [:math:`\mu \mathrm{m}`]
    :returns: array_like tputs: Filter throughputs
    :returns: list tputs: List of MIRI filter names
    '''
    import pandas as pd

    # Read in MIRI filter response functions
    df_filt = pd.read_csv(f)
    vals = np.array(df_filt.values[1:,:], dtype=float)

    # Wavelength grid
    wl_filt = vals[:,0]

    # Delta wavelengths
    dwl_filt = wl_filt[1:] - wl_filt[:-1]
    dwl_filt = np.hstack([dwl_filt, dwl_filt[-1]])

    # Filter throughputs
    tputs = vals[:,1:]

    # Filter names
    names = np.array(df_filt.columns[1:], dtype=str)

    return wl_filt, dwl_filt, tputs, names

def get_spitzer_filter_wheel(warm = True):
    """
    Read-in and instantiate list of Spitzer IRAC filters as
    :py:func:`Filter` objects

    :param bool warm: If :py:obj:`True` uses only Warm Spitzer filters

    :returns: A list of :py:func:`Filter` objects
    """

    # Use only warm Spitzer filters?
    if warm:
        N = 2
    else:
        N = 4

    # Read-in Spitzer IRAC filters
    data1 = np.genfromtxt(os.path.join(HERE, "filter_files",
                          "080924ch1trans_full.txt"), skip_header=3)
    data2 = np.genfromtxt(os.path.join(HERE, "filter_files",
                          "080924ch2trans_full.txt"), skip_header=3)
    data3 = np.genfromtxt(os.path.join(HERE, "filter_files",
                          "080924ch3trans_full.txt"), skip_header=3)
    data4 = np.genfromtxt(os.path.join(HERE, "filter_files",
                          "080924ch4trans_full.txt"), skip_header=3)

    spitzer_names = ["IRAC3.6", "IRAC4.5", "IRAC5.8", "IRAC8.0"]
    spitzer_data = [data1, data2, data3, data4]

    # Construct filter wheel
    wheel = [Filter(name = spitzer_names[i], wl=spitzer_data[i][:,0],
                    throughput=spitzer_data[i][:,1])
             for i in range(N)]

    return wheel

################################################################################

class Filter(object):
    '''
    A photometruc filter class. This is an interface for filter
    photometry calculations.

    :param str name: Name of the filter
    :param array_like wl: Wavelength grid for filter [:math:`\mu \mathrm{m}`]
    :param array_like throughput: Filter throughput
    :param array_like dwl: Wavelength width grid for filter \
           [:math:`\mu \mathrm{m}`]
    :param float eff_wl: Effective filter wavelenth [:math:`\mu \mathrm{m}`]
    :param float eff_dwl: Effective filter width [:math:`\mu \mathrm{m}`]
    '''
    def __init__(self, name=None, wl=None, throughput=None, dwl=None,
                 eff_wl=None, eff_dwl=None):
        '''
        Instantiate the class.
        '''
        self.name = name
        self.wl = wl
        self.throughput = throughput
        self._dwl = dwl
        self._eff_wl = eff_wl
        self._eff_dwl = eff_dwl

        if self._dwl is None:
            self._dwl = self._calc_dwl()
        if self._eff_wl is None:
            self._eff_wl = self._calc_eff_wl()
        if self._eff_dwl is None:
            self._eff_dwl = self._calc_eff_dwl()

    @property
    def dwl(self):
        return self._dwl

    @dwl.setter
    def dwl(self, value):
        # if dwl is reset, then recalculate eff_wl and eff_dwl
        self._dwl = value
        self._eff_wl = self._calc_eff_wl()
        self._eff_dwl = self._calc_eff_dwl()

    @property
    def eff_wl(self):
        return self._eff_wl

    @eff_wl.setter
    def eff_wl(self, value):
        self._eff_wl = value

    @property
    def eff_dwl(self):
        return self._eff_dwl

    @eff_dwl.setter
    def eff_dwl(self, value):
        self._eff_dwl = value

    def _calc_dwl(self):
        dwl = self.wl[1:] - self.wl[:-1]
        dwl = np.hstack([dwl, dwl[-1]])
        return dwl

    def _calc_eff_wl(self):
        return sum(self.wl*self.throughput*self.dwl)/sum(self.throughput*self.dwl)

    def _calc_eff_dwl(self):
        return sum(self.dwl*self.throughput)/np.max(self.throughput)

    def plot(self, ax=None):
        '''
        Plots the filter throughput curve.

        :param ax: An axis instance
        :type ax: :py:obj:`axis`

        '''
        import matplotlib.pyplot as plt
        if ax is None:
            fig, axi = plt.subplots(figsize=(10,6))
            axi.set_xlabel("Wavelength")
            axi.set_ylabel("Throughput")
        else:
            axi = ax

        axi.plot(self.wl, self.throughput)
        #axi.errorbar(self.eff_wl, np.max(self.throughput)/2,
        # xerr=self.eff_dwl/2, fmt="o", c="k", ms=5)

        plt.show()

    def photon_rate(self, lam, flux, atel=25.0, dlam=None):
        '''
        Compute the photon count rate registered by the detector.

        :param array_like lam: High-res wavelengths [:math:`\mu \mathrm{m}`]
        :param array_like flux: Spectral flux density \
               [:math:`\mathrm{W/m}^2 / \mu \mathrm{m}`]
        :param float atel: Telescope collecting area [:math:`\mathrm{m}^2`]. \
               Default is `25`
        :param array_like dlam: Delta-wavelength grid. Default is \
               :py:obj:`None` and is calculated

        :returns array_like cphot: Photon count rate [:math:`\mathrm{s}^{-1}`]

        '''

        # (Speed of light) * (Planck's constant)
        hc = 1.986446e-25  # h*c (kg*m**3/s**2)

        if dlam is None:
            dlam = lam[1:] - lam[:-1]
            dlam = np.hstack([dlam, dlam[-1]])

        # interpolate filter throughout to HR grid
        T = np.interp(lam, self.wl, self.throughput)

        cphot = np.sum(flux*dlam*(lam*1e-6)/hc*T*atel, axis=flux.ndim-1)

        return cphot

    def convolve(self, lam, flux):
        '''
        Convolve flux with normalized filter throughput.

        :param array_like lam: High-res wavelength grid [:math:`\mu \mathrm{m}`]
        :param array_like flux:  High-res flux grid \
               [:math:`\mathrm{W/m}^2 / \mu \mathrm{m}`]

        :returns array_like F: Flux convolved with normalized throughput

        '''

        # interpolate filter throughout to HR grid
        T = np.interp(lam, self.wl, self.throughput)

        # Convolve with normalized throughput
        F = np.sum(flux * T) / np.sum(T)

        return F

    def compute_lightcurve(self, time, flux, continuum, lam, stack = 1,
                           atel = 25., thermal = True, time_hr = None,
                           flux_hr = None, quiet = False):
        """
        Computes an observed lightcurve in the :py:func:`Filter`.

        :param array_like time: Time grid [days]
        :param array_like flux: Observed flux grid (`time` by `lam`) \
               [:math:`\mathrm{W/m}^2 / \mu \mathrm{m}`]
        :param array_like continuum: Observed *continuum* flux grid (no \
               occultations, shape `time` by `lam`) \
               [:math:`\mathrm{W/m}^2 / \mu \mathrm{m}`]
        :param array_like lam: Wavelength [:math:`\mu \mathrm{m}`]
        :param int stack: Number of exposures to stack. Default `1`
        :param float atel: Telescope collecting area [:math:`\mathrm{m}^2`]. \
               Default `25`
        :param bool thermal: Compute thermal noise. Default :py:obj:`True`
        :param array_like time_hr: High-res time grid [Days]. \
               Default :py:obj:`None`
        :param array_like flux_hr: High-res flux grid. Default :py:obj:`None`

        """

        if not quiet:
            print("Computing observed light curve in %s filter..." % self.name)

        Ntime = len(time)
        Nlam = len(lam)
        tmin = np.min(time)
        tmax = np.max(time)

        # Low-res exposure time array
        dtlr = time[1:] - time[:-1]
        dtlr = np.hstack([dtlr, dtlr[-1]])

        # Hi-res exposure time array
        if time_hr is not None:
            zeros = np.where(np.diff(time_hr) == 0)[0]
            time_hr = np.delete(time_hr, zeros)
            flux_hr = np.delete(flux_hr, zeros, axis = 0)
            dthr = time_hr[1:] - time_hr[:-1]
            dthr = np.hstack([dthr, dthr[-1]])

        # Calculate jwst background flux
        if thermal:
            Fback = jwst_background(lam)
        else:
            Fback = np.zeros_like(lam)

        # Exposure time [s]
        tint = dtlr * 3600. * 24

        # Calculate SYSTEM photons
        Nsys = stack * tint * self.photon_rate(lam, flux[:,:], atel = atel)

        # Calculate CONTINUUM photons
        Ncont = stack * tint * self.photon_rate(lam, continuum[:,:],
                                                atel = atel)

        # Hi-res light curve
        if time_hr is not None:
            tint_hr = dthr * 3600. * 24
            Nsys_hr = stack * tint_hr * self.photon_rate(lam, flux_hr[:,:],
                                                         atel = atel)
            norm_hr = np.median(Nsys_hr)
        else:
            Nsys_hr = None
            norm_hr = None

        # Calculate BACKGROUND photons
        Nback = stack * tint * self.photon_rate(lam, Fback, atel = atel)

        # Signal-to-noise
        SNR = Nsys / np.sqrt(Nsys + Nback)

        # Generate synthetic data points
        norm = np.median(Nsys)
        sig = np.sqrt(Nsys + Nback) / norm
        obs = random_draw(Nsys / norm, sig)

        # Create lightcurve object to hold observed quantities
        self.lightcurve = Lightcurve(time = time,
                                     Nsys = Nsys,
                                     Ncont = Ncont,
                                     Nback = Nback,
                                     SNR = SNR,
                                     obs = obs,
                                     sig = sig,
                                     norm = norm,
                                     tint = tint,
                                     stack = stack,
                                     time_hr = time_hr,
                                     Nsys_hr = Nsys_hr,
                                     norm_hr = norm_hr)

        return

###############################################################################

class Lightcurve(object):
    """
    A lightcurve class to contain all outputs from a snythetic
    lightcurve observation.

    :param array_like time: Observed time grid [days]
    :param array_like Nsys: Number of photons from `System`
    :param array_like Ncont: Number of continuum photons from `System`
    :param array_like Nback: Number of photons from background
    :param array_like SNR: Signal-to-noise on `System`
    :param array_like obs: Observed photon signal (normalized)
    :param array_like sig: 1-sigma errors on signal (normalized)
    :param float norm: Normalization constant (median of lightcurve)
    :param array_like tint: Integration time [mins]

    """
    def __init__(self, time = None, Nsys = None, Ncont = None, Nback = None,
                 SNR = None, obs = None, sig = None, norm = None, tint = None,
                 stack = None, time_hr = None, Nsys_hr = None,
                 norm_hr = None):
        '''

        '''

        self.time = time
        self.Nsys = Nsys
        self.Ncont = Ncont
        self.Nback = Nback
        self.SNR = SNR
        self.obs = obs
        self.sig = sig
        self.norm = norm
        self.tint = tint
        self.stack = stack
        self.time_hr = time_hr
        self.Nsys_hr = Nsys_hr
        self.norm_hr = norm_hr

    def plot(self, ax0=None, title="", alpha_err = 0.7):
        '''
        Plots the synthetic lightcurve.

        :param ax0: User provided plot :py:obj:`axis`. Default :py:obj:`None`
        :type ax0: :py:obj:`axis`
        :param str title: Plot title. Defult ""

        '''

        # Create new fig if axis is not user provided
        if ax0 is None:
            fig, ax = plt.subplots(figsize=(16,6))
            ax.set_title(r"%s" %title)
            ax.set_ylabel("Relative Flux")
            ax.set_xlabel("Time [days]")
        else:
            ax = ax0

        # Plot
        ax.plot(self.time, self.Nsys / self.norm, label='Binned', zorder = 11,
                 alpha=0.75, lw = 1.5, color = 'b')

        if (self.time_hr is not None) and (self.Nsys_hr is not None):
          ax.plot(self.time_hr, self.Nsys_hr /  self.norm_hr, zorder = 11,
                  alpha=0.75, lw = 1, label = 'Unbinned', color = 'g')

        ax.errorbar(self.time, self.obs, yerr=self.sig, fmt="o", c="k", ms=2,
                    alpha=alpha_err, zorder=10, lw = 1)
        if self.stack > 1:
            ax.text(0.02, 0.95, r"$\Delta t = %.1f$ mins ($\times$ %d)"
                    % (self.tint[0]/60., self.stack),
                    ha="left", va="top", transform = ax.transAxes,
                    fontsize=12)
        else:
            ax.text(0.02, 0.95, r"$\Delta t = %.1f$ mins"
                    % (self.tint[0]/60.),
                    ha="left", va="top", transform = ax.transAxes,
                    fontsize=12)

        ax.legend(loc = 'upper right')
        ax.ticklabel_format(useOffset = False)
        ax.margins(0, 0.15)

        if ax0 is None:
            fig.subplots_adjust(bottom=0.2)

        return ax0

################################################################################

def planck(temp, wav):
    '''
    Planck blackbody function.

    :param float or array_like temp: Temperature [K]
    :param float or array_like wav: Wavelength [:math:`\mu \mathrm{m}`]

    :returns array_like B_lambda: Planck function \
             [:math:`\mathrm{W/m}^2 / \mu \mathrm{m}`]

    '''

    h = 6.62607e-34       # Planck constant (J * s)
    c = 2.998e8           # Speed of light (m / s)
    k = 1.3807e-23        # Boltzmann constant (J / K)
    wav = wav * 1e-6

    # Returns B_lambda [W/m^2/um/sr]
    return 1e-6 * (2. * h * c**2) / (wav**5) / \
           (np.exp(h * c / (wav * k * temp)) - 1.0)

def plot_miri_filters(ax, wl, filters, names, ylim=(0.0,1.0), leg=True):
    """
    Plot JWST/MIRI photometric filter bands.

    :param :py:obj:`Axis` ax: Matplotlib Axis
    :param array_like wl: Wavelengths
    :param array_like filters: Filter throughput curves
    :param list names: Names of filters
    :param tuple ylim: y-axis limits. Default `(0.0, 1.0)`
    :param bool leg: Add a legend to the plot. Default :py:obj:`True`

    """

    ax1 = ax.twinx()
    for i in range(len(names)):
        # Plot filter throughputs
        ax1.fill_between(wl, filters[:,i], label=names[i], alpha=0.3)
        ax1.set_ylabel("Filter Throughput", rotation=270, labelpad=25)
    if leg:
        leg=ax1.legend(loc=2, fontsize=16, ncol=1)
        leg.get_frame().set_alpha(0.0)
    ax1.set_ylim(ylim)

def jwst_background(wl):
    """
    Calculates the JWST thermal background from (Glasse et al. 2015)

    :param array_like wl: Wavelength grid [:math:`\mu \mathrm{m}`]

    :returns array_like Fback: Spectral flux density \
             [:math:`\mathrm{W/m}^2 / \mu \mathrm{m}`]
    """

    # Solid angle
    omega_background = np.pi*(.42/206265.*wl/10.)**2.

    # Number of background grey body components (see table 1 of Glasse et al.):
    nback = 7
    emiss = [4.2e-14,4.3e-8,3.35e-7,9.7e-5,1.72e-3,1.48e-2,1.31e-4]
    tb = [5500.,270.,133.8,71.0,62.0,51.7,86.7]

    # Set up a vector for computing the background flux emissivity:
    Fback = np.zeros_like(wl)

    # Sum over the background components (insert 1 steradian for omega):
    for i in range(nback):
        Fback += emiss[i] * planck(tb[i], wl) * omega_background

    return Fback

def estimate_eclipse_snr(tint = 36.4*60., nout = 4.0, lammin = 1.0,
                         lammax = 30.0, Nlam = 10000, Tstar = 2560.,
                         Tplan = 400., Rs = 0.12, Rp = 1.086, d = 12.2,
                         atel = 25.0, verbose=True, plot=True, thermal = True,
                         filters = 'MIRI'):
    """
    Estimate the signal-to-noise on the detection of secondary eclipses in
    JWST/MIRI photometric filters.

    :param float tint: Exposure time [s]. Default `36.4*60.`
    :param float nout: Out-of-transit time observed [transit durations]. \
           Default `4`
    :param float lammin: Wavelength minimum [:math:`\mu \mathrm{m}`]. \
           Default `1.0`
    :param float lammax: Wavelength maximum [:math:`\mu \mathrm{m}`]. \
           Default `30.0`
    :param int Nlam: Number of wavelengths. Default `10000`
    :param float Tstar: Stellar effective temperature [K]. Default `2560`
    :param float Tplan: Planet equilibrium temperature [K]. Default `400`
    :param float Rs: Stellar radius [solar radii]. Default `0.12`
    :param float Rp: Planet radius [Earth radii]. Default `1.086`
    :param float d: System distance [pc]. Default `12.2`
    :param float atel: Telescope collecting area [:math:`\mathrm{m}^2`]. \
           Default `25`
    :param verbose: Print things. Default `True`
    :param plot: Make a plot. Detault `True`
    :param thermal: Include thermal noise. Default `True`
    :param filters: User provided filter names or list of filters. Default \
            `MIRI`
    :type filters: str or list of :py:func:`Filter` objects

    """

    # Integration time is eclipse duration for planet b
    # (assume same duration as transit):
    #tint = 36.4 * 60. # sec
    # Number of out-of-eclipse transit durations observed (with planet + star):
    #nout = 4.0

    # Generate high-res wavelength grid
    lam = np.linspace(lammin, lammax, Nlam)
    dlam = lam[1:] - lam[:-1]
    dlam = np.hstack([dlam, dlam[-1]])

    # Calculate BB intensities for the star and planet [W/m^2/um/sr]
    Bstar = planck(Tstar, lam)
    Bplan = planck(Tplan, lam)

    # solid angle in steradians
    omega_star = np.pi*(Rs*u.Rsun.in_units(u.km)/(d*u.pc.in_units(u.km)))**2.
    omega_planet = np.pi*(Rp*u.Rearth.in_units(u.km)/
                         (d*u.pc.in_units(u.km)))**2.

    # fluxes at earth [W/m^2/um]
    Fstar = Bstar * omega_star
    Fplan = Bplan * omega_planet
    Fback = jwst_background(lam)

    # Decide on filters to use
    if type(filters) is list:

        # Use provided list of filters
        wheel = filters
        nfilt = len(filters)

    elif filters == "MIRI":

        # Read-in MIRI filters
        wl_filt, dwl_filt, tputs, fnames = readin_miri_filters()
        nfilt = len(fnames)

        # Construct Filter "wheel"
        wheel = [
            Filter(
                name=fnames[i], wl=wl_filt, throughput=tputs[:,i]
            )
            for i in range(nfilt)
            ]

    if plot:
        fig, ax = plt.subplots(figsize=(11,7))
        ax.set_xlabel(r"Wavelength [$\mu$m]")
        ax.set_ylabel("SNR")
        if filters == "MIRI": plot_miri_filters(ax, wl_filt, tputs, fnames)

    # Loop over filters
    for i in range(nfilt):

        # Count STELLAR photons
        Nphot_star = tint * wheel[i].photon_rate(lam, Fstar, atel = atel)

        # Count PLANET photons
        Nphot_planet = tint * wheel[i].photon_rate(lam, Fplan, atel = atel)

        # Count BACKGROUND photons
        if thermal:
            Nphot_bg = tint * wheel[i].photon_rate(lam, Fback, atel = atel)
        else:
            Nphot_bg = np.zeros_like(Nphot_planet)

        # Calculate SNR on planet photons
        SNR = Nphot_planet / np.sqrt((1+1./nout)*Nphot_star
                             + 1./nout*Nphot_planet+(1+1./nout)*Nphot_bg)

        # Optionally print
        if verbose:
            print("Filter: ",wheel[i].name," Wavelength: ",wheel[i].eff_wl)
            print("Photons from star: %.3e" %(Nphot_star))
            print("Photons from planet: %.3e" %(Nphot_planet))
            print("Photons from background: %.3e" %(Nphot_bg))
            print("S/N: %.3f" %(SNR))
            print("------------")

        # Optionally plot
        if plot:
            ax.plot(wheel[i].eff_wl, SNR, "o", c="k")

    if plot:
        return fig, ax

def fake_time_func(t, factor=0.01):
    """
    Returns a function of time for testing.

    :param array_like t: Time grid
    :param float factor: Relative variability factor

    :returns array_like f: A function of time
    """
    f = np.sin(t/2)**2 * 0.5*np.sin(t/10.)
    return 1.0 + factor*(f/np.max(f))


def create_fake_data(Nlam=4000, Ntime=4000, lammin=1.0, lammax=30.0, tmin=0.0,
                     tmax=4.0, Tstar = 2560., Tplan=400., Rs=0.12, Rp=1.086,
                     d=12.2, tfact=0.01):
    """
    Creates a fake flux dataset as a function of wavelength and
    time for testing.

    :param int Nlam: Number of wavelengths. Default `4000`
    :param int Ntime: Number of times. Default `4000`
    :param float lammin: Wavelength min [:math:`\mu \mathrm{m}`]. Defult `1.0`
    :param float lammax: Wavelength max [:math:`\mu \mathrm{m}`]. Defult `30.0`
    :param float tmin: Time min [days]. Default `0.0`
    :param float tmax: Time max [days]. Default `4.0`
    :param float Tstar: Stellar effective temperature [K]. Defult `2560`
    :param float Tplan: Planet equilibrium temperature [K]. Defult `400`
    :param float Rs: Stellar radius [solar radii]. Default `0.12`
    :param float Rp: Planet radius [Earth radii]. Default `1.086`
    :param float d: System distance [pc]. Default `12.2`
    :param float tfact: Relative variability factor. Default `0.01`

    :returns tuple ttup: (`time`, `delta_time`)
    :returns tuple ltup: (`lam`, `delta_lam`)
    :returns tuple data: Flux as a function of time and wavelength
    """

    """
    1. WAVELENGTH (using Planck function)
    """

    # Generate high-res wavelength grid
    lam = np.linspace(lammin, lammax, Nlam)
    dlam = lam[1:] - lam[:-1]
    dlam = np.hstack([dlam, dlam[-1]])

    # Calculate BB intensities for the star and planet
    Bstar = planck(Tstar, lam)
    Bplan = planck(Tplan, lam)

    # solid angle in steradians
    omega_star = np.pi*(Rs*u.Rsun.in_units(u.km)/(d*u.pc.in_units(u.km)))**2.
    omega_planet = np.pi*(Rp*u.Rearth.in_units(u.km)/
                         (d*u.pc.in_units(u.km)))**2.
    omega_background = np.pi*(.42/206265.*lam/10.)**2.

    # fluxes at earth [W/m^2/um]
    Fstar = Bstar * omega_star
    Fplan = Bplan * omega_planet

    """
    2. TIME (using fake function)
    """

    # Generate high-res time grid [hours]
    t = np.linspace(tmin, tmax, Ntime)
    dt = t[1:] - t[:-1]
    dt = np.hstack([dt, dt[-1]])

    # Generate fake time-dependent function
    ftime = fake_time_func(t*24, factor=tfact)

    # Create 2d dataset
    data = np.outer(ftime, Fstar)

    return (t, dt), (lam, dlam), data

def gen_lr_grid(lomin, lomax, Nlow):
    '''
    Generate a low-resolution grid and bin widths

    :param float lomin: Min value
    :param float lomax: Max value
    :param int Nlow: Number of final points

    :returns array_like lr: New evenly-spaced grid
    :returns array_like dlr: New bin widths
    '''
    # Use linspace to create an evenly spaced grid
    lr = np.linspace(lomin,lomax,Nlow+1)
    # Bin widths
    dlr = lr[1:] - lr[:-1]
    # Find xmean of neighboring bins
    lr = (lr[:-1] + lr[1:])/2.
    return lr, dlr

def downbin_series(yhr, xhr, xlr, dxlr=None):
    """
    Re-bin series to lower resolution using `sp.binned_statistic`

    :param array_like yhr: Series to be degraded
    :param array_like xhr: High-res x grid
    :param array_like xlr: Low-res x grid
    :param array_like dxlr: Low-res x width grid. Default :py:obj:`None`

    :returns: array_like ylr: Low-res series
    """
    from scipy.stats import binned_statistic

    if dxlr is None:
        # Calculate low-res
        #dxlr = xlr
        ValueError("Please supply dlam in downbin_spec()")

    # Reverse ordering if wl vector is decreasing with index
    if len(xlr) > 1:
        if xhr[0] > xhr[1]:
            xhr = np.array(xhr[::-1])
            spec = np.array(yhr[::-1])
        if xlr[0] > xlr[1]:
            xlr = np.array(xlr[::-1])
            dxlr = np.array(dxlr[::-1])

    # Calculate bin edges
    LRedges = np.hstack([xlr - 0.5*dxlr, xlr[-1]+0.5*dxlr[-1]])

    # Call scipy.stats.binned_statistic()
    ylr = binned_statistic(xhr, yhr, statistic="mean", bins=LRedges)[0]

    return ylr

def random_draw(val, sig):
    """
    Draw fake data points from model `val` with errors `sig`

    :param float or array_like val: Mean of Gaussian to sample from
    :param float or array_like sig: Standard dev. of Gaussian to sample from

    :returns float or array_like: Randomly sampled fake data points
    """
    if type(val) is np.ndarray:
        return val + np.random.randn(len(val))*sig
    elif (type(val) is float) or (type(val) is int):
        return val + np.random.randn(1)[0]*sig

def get_miri_filter_wheel():
    """
    Create a `list` of MIRI :py:func:`Filter` objects

    :returns list wheel: `list` of MIRI :py:func:`Filter` objects
    """

    # Read-in MIRI filters
    wl_filt, dwl_filt, tputs, fnames = readin_miri_filters()
    nfilt = len(fnames)

    # Construct Filter "wheel"
    wheel = [
                Filter(name=fnames[i], wl=wl_filt, throughput=tputs[:,i])
                for i in range(nfilt)
            ]

    return wheel

def create_tophat_filter(lammin, lammax, dlam=0.1, Tput=0.3, name="custom"):
    """
    Create a tophat :py:func:`Filter`

    :param float lammin: Wavelength minimum [:math:`\mu \mathrm{m}`]
    :param float lammax: Wavelength maximum [:math:`\mu \mathrm{m}`]
    :param float dlam: Wavelength resolution [:math:`\mu \mathrm{m}`]. \
           Default `0.1`
    :param float Tput: Filter throughput. Default `0.3`
    :param str name: Name of filter. Default "custom"

    :returns filt: New custom tophat :py:func:`Filter` object
    :type filt: :py:func:`Filter`

    """

    # Number of wavelength points
    N = int(round((lammax - lammin) / dlam))

    # Create wavelength grid
    lam = np.arange(lammin-2*dlam, lammax+2*dlam, dlam)

    # Create throughput curve
    Tputs = np.ones_like(lam) * Tput

    # Set the first and last 2 elements to zero
    Tputs[0] = 0.0
    Tputs[1] = 0.0
    Tputs[-1] = 0.0
    Tputs[-2] = 0.0

    # Construct filter
    filt = Filter(name=name, wl=lam, throughput=Tputs)

    return filt
