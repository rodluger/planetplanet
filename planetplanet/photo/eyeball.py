#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
eyeball.py |github|
-------------------

Code for plotting and visualizing "eyeball" planets.

  .. role:: raw-html(raw)
     :format: html
     
  .. |github| replace:: :raw-html:`<a href = "https://github.com/rodluger/planetplanet/blob/master/planetplanet/photo/eyeball.py"><i class="fa fa-github" aria-hidden="true"></i></a>`

'''

from .maps import RadiativeEquilibriumMap, LimbDarkenedMap
from ..constants import *
import numpy as np
np.seterr(invalid = 'ignore')
import matplotlib.pyplot as pl
from matplotlib.transforms import Affine2D
import mpl_toolkits.axisartist.floating_axes as floating_axes
from matplotlib.widgets import Slider

__all__ = ['DrawEyeball', 'DrawOrbit', 'GetAngles']

def rodrigues(v, k, theta):
    '''
    The Rodrigues rotation formula in 3D, given a vector `v`, a unit vector `k`
    normal to the plane of rotation, and the angle of rotation `theta` 
    in radians.
    
    '''
        
    return v * np.cos(theta) + np.cross(k, v) * np.sin(theta) \
           + k * np.dot(k, v) * (1 - np.cos(theta))

def ucross(u, v):
    '''
    The cross product of vectors `u` and `v`, normalized to unity.
    
    '''
    
    res = np.cross(u, v)
    res /= np.sqrt(np.sum(res ** 2))
    return res

def GetAngles(x, y, z, vx, vy, vz, Lambda = 0., Phi = 0.):
    '''
    Computes the eyeball angles :math:`\\theta` and :math:`\gamma` given the 
    Cartesian orbital elements and the hotspot offset angles :math:`\Lambda` 
    and :math:`\Phi`.
    
    :param float x: The **x** component of the planet's position vector
    :param float y: The **y** component of the planet's position vector
    :param float z: The **z** component of the planet's position vector
    :param float vx: The **x** component of the planet's velocity vector
    :param float vy: The **y** component of the planet's velocity vector
    :param float vz: The **z** component of the planet's velocity vector
    :param float Lambda: The longitudinal hotspot offset in radians. \
           Default `0`
    :param float Phi: The latitudinal hotspot offset in radians. Default `0`
    :returns: :math:`\\theta` and :math:`\gamma`, the eyeball phase angles, \
              in radians
    
    '''
    
    # The position vector of the center of the planet and its magnitude
    r = np.array([x, y, z])
    R = np.sqrt(np.sum(r ** 2))
    
    # The velocity vector of the planet and its magnitude
    v = np.array([vx, vy, vz])
    V = np.sqrt(np.sum(v ** 2))

    # The position vector of the substellar point,
    # relative to the planet center, normalized to 1
    rstar = -r / R
    
    # Vector normal to the longitudinal plane
    vlon = ucross(r / R, v / V)
        
    # Apply the longitudinal offset
    rstar = rodrigues(rstar, vlon, Lambda)

    # Vector normal to the latitudinal plane
    vlat = ucross(vlon, rstar)
    
    # Apply the latitudinal offset
    rstar = rodrigues(rstar, vlat, Phi)

    # The Cartesian coordinates of the hotspot,
    # relative to the planet center, normalized to 1
    xstar, ystar, zstar = rstar

    # Projected distance from planet center to hotspot, normalized to 1
    d = min(1, np.sqrt(xstar ** 2 + ystar ** 2))
    
    # Get the rotation and phase angles
    gamma = np.arctan2(ystar, xstar) + np.pi
    if zstar <= 0:
        theta = np.arccos(d)
    else:
        theta = -np.arccos(d)

    return theta, gamma
  
def LimbDarkenedFlux(lam, z, teff = 2500, limbdark = [1.]):
    '''
    
    '''
    
    # Evaluate the limb darkening function if necessary 
    u = [None for ld in limbdark]
    for n, ld in enumerate(limbdark):
        if callable(ld):
            u[n] = ld(lam)
        elif not hasattr(ld, '__len__'):
            u[n] = ld
        else:
            raise Exception("Limb darkening coefficients must be provided " 
                          + "as a list of scalars or as a list of functions.")
    limbdark = u
    
    # Convert to m
    lam *= 1e-6
        
    # Compute the normalization term, Equation (E5)
    norm = 0
    for i, u in enumerate(limbdark):
        norm += u / ((i + 2) * (i + 3))
    norm = 1 - 2 * norm
    a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
    b = HPLANCK * CLIGHT / (lam * KBOLTZ * teff)
    B0 = a / (np.exp(b) - 1.) / norm
    
    # Initialize
    flux = B0
    cosz = np.cos(z)

    # Loop over the coefficient order
    for i, u in enumerate(limbdark):

        # The Taylor expansion is in (1 - mu)
        x = (1 - cosz) ** (i + 1)
    
        # Compute the wavelength-dependent intensity
        flux -= u * B0 * x
    
    return flux

def RadiativeEquilibriumFlux(lam, z, albedo = 0.3, tnight = 40., 
                             irrad = SEARTH):
    '''
    
    '''

    # Convert to m
    lam *= 1e-6

    # Compute the temperature
    if (z < np.pi / 2):
        temp = ((irrad * np.cos(z) * (1 - albedo)) / SBOLTZ) ** 0.25
        if (temp < tnight):
            temp = tnight
    else:
        temp = tnight

    # Compute the radiance
    a = 2 * HPLANCK * CLIGHT * CLIGHT / (lam * lam * lam * lam * lam)
    b = HPLANCK * CLIGHT / (lam * KBOLTZ * temp)
    return a / (np.exp(b) - 1.)

def DrawEyeball(x0 = 0.5, y0 = 0.5, r = 0.5, 
                radiancemap = RadiativeEquilibriumMap(), 
                theta = np.pi / 3, nz = 31, gamma = 0, occultors = [], 
                cmap = 'inferno', fig = None, draw_terminator = False, 
                draw_ellipses = False, rasterize = False,
                cpad = 0.2, limbdark = [1.], teff = 2500., wavelength = 15.,
                color = None):
    '''
    Creates a floating axis and draws an "eyeball" planet at given 
    phase and rotation angles.
    
    .. plot::
         :align: center
         
         from planetplanet import DrawEyeball
         from planetplanet.photo.maps import RadiativeEquilibriumMap
         import matplotlib.pyplot as pl
         fig = pl.figure(figsize = (3,3))
         DrawEyeball(radiancemap = RadiativeEquilibriumMap(), fig = fig)
         pl.show()
    
    :param float x0: The `x` position of the center of the planet in \
           figure coordinates
    :param float y0: The `y` position of the center of the planet in \
           figure coordinates
    :param float r: The radius of the planet in figure coordinates
    :param float theta: The phase angle of the eyeball in radians. \
           Default :math:`\pi/3`
    :param float gamma: The rotation angle of the eyeball in radians. \
           Default `0`
    :param int nz: The number of zenith angle wedges. Default `11`
    :param array_like occultors: A list of :py:obj:`dict` instances with \
           information on each of the occultors to draw. \
    Each dictionary must have keywords `x`, `y`, and `r`, corresponding to \
           the xy position and radius of the \
    occultor, respectively. These are defined relative to the center of the \
           planet and scaled so that the \
    radius of the planet is unity. Optional keywords are `zorder`, `color`, \
           and `alpha`. Default :py:obj:`[]`
    :param str cmap: The name of the :py:class:`matplotlib` colormap. \
           Default `inferno`
    :param fig: The figure object in which to create the axis. \
           Default :py:obj:`None`
    :type fig: :py:class:`matplotlib.Figure`
    :param bool draw_terminator: Draw the terminator ellipse outline? \
           Default :py:obj:`False`
    :param bool draw_ellipses: Draw the zenith angle ellipse outlines? \
           Default :py:obj:`False`
    :param bool rasterize: Rasterize the image? Default :py:obj:`False`
    :param str color: Occulted body outline color. Default :py:obj:`None`
    
    :return: **fig**, **ax**, **occ**, **xy**. These are the figure and \
             floating axis instances, a lits of :py:obj:`Circle` \
             instances corresponding to each of the occultors, and a function \
             **xy** that performs the rotation \
             transformation into the axis-symmetric eyeball frame
    
    '''
    
    # Check the symmetry
    if (radiancemap.maptype == MAP_RADIAL_DEFAULT) or \
       (radiancemap.maptype == MAP_RADIAL_CUSTOM):
        theta = np.pi / 2
        gamma = 0
        
    # The rotation transformation, Equation (E6) in the paper
    xy = lambda x, y: (x * np.cos(gamma) + y * np.sin(gamma), 
                       y * np.cos(gamma) - x * np.sin(gamma))
    
    # Set up the floating axis
    if fig is None:
        fig = pl.figure(figsize = (6,6))
    tr = Affine2D().rotate_deg(gamma * 180 / np.pi)
    x = 1. / (np.abs(np.cos(gamma)) + np.abs(np.sin(gamma)))
    scale = max([1] + [occultor['r'] for occultor in occultors])
    x *= scale
    grid_helper = floating_axes.GridHelperCurveLinear(tr, 
                                                      extremes=(-x, x, -x, x))
    ax_orig = floating_axes.FloatingSubplot(fig, 111, 
                                            grid_helper = grid_helper)
    ax_orig.set_position([x0 - r, y0 - r, 2 * r, 2 * r])
    ax_orig.axis["bottom"].set_visible(False)
    ax_orig.axis["top"].set_visible(False)
    ax_orig.axis["left"].set_visible(False)
    ax_orig.axis["right"].set_visible(False)
    ax_orig.patch.set_alpha(0)
    fig.add_subplot(ax_orig)
    ax = ax_orig.get_aux_axes(tr)
    if rasterize:
        ax_orig.set_rasterization_zorder(9999)
        ax.set_rasterization_zorder(9999)
    
    # Plot the occultors. Note that we need to transform
    # their position vectors since we're in a rotated frame.
    occ = [None for i in occultors]
    for i, occultor in enumerate(occultors):
        xo = occultor['x']
        yo = occultor['y']
        ro = occultor['r']
        zo = occultor.get('zorder', 1)
        ec = occultor.get('color', 'k')
        ao = occultor.get('alpha', 1)
        xo, yo = xy(xo, yo)
        occ[i] = pl.Circle((xo, yo), ro, color = 'lightgrey',
                           ec = ec, lw = 2,
                           alpha = ao, zorder = zo, clip_on = False)
        ax.add_artist(occ[i])
        
    # Plot the occulted body
    x = np.linspace(-1, 1, 1000)
    y = np.sqrt(1 - x ** 2)
    if color is not None:
        ax.plot(x, y, color = color, zorder = 0, lw = 2, clip_on = False)
        ax.plot(x, -y, color = color, zorder = 0, lw = 2, clip_on = False)
        
    # Get the radiance map. If it's one of the default maps,
    # we need to call their special Python implementations defined
    # above. Otherwise we use the `ctypes` method of the map.
    # We will normalize it so that we can use it as a colormap.
    if radiancemap.maptype == MAP_RADIAL_DEFAULT:
        # Function wrapper to use correct limb darkening parameters
        func = lambda lam, z: LimbDarkenedFlux(lam, z, teff = teff, 
                                               limbdark = limbdark)
    elif radiancemap.maptype == MAP_ELLIPTICAL_DEFAULT:
        # TODO: Technically we should pass the albedo and night side 
        # temperature to the function here But the albedo doesn't matter 
        # for drawing because of the normalization; the night side
        # temperature could matter, but only if it's comparable to the 
        # dayside temperature, otherwise the nightside is dark no matter what.
        func = RadiativeEquilibriumFlux
    else:
        func = radiancemap.ctypes
    zarr = np.linspace(0, np.pi, 100)
    rarr = [func(wavelength, za) for za in zarr]
    rmax = np.max(rarr)
    rmin = np.min(rarr)
    rrng = rmax - rmin
    rmax += cpad * rrng
    rmin -= cpad * rrng
    if rrng > 0:
        color = lambda z: pl.get_cmap(cmap)((func(wavelength, z) - rmin) 
                                            / (rmax - rmin))
    else:
        color = lambda z: pl.get_cmap(cmap)(1 - cpad)

    # Plot the zenith angle ellipses
    zarr = np.linspace(0, np.pi, nz + 2)
    for i, z in enumerate(zarr[1:]):

        # The ellipse
        a = np.abs(np.sin(z))
        b = max(0.001, a * np.abs(np.sin(theta)))
        xE = -np.cos(z) * np.cos(theta)
        xlimb = np.cos(z) * np.sin(theta) * np.tan(theta)
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
            # This is the terminator
            if draw_terminator and scale < 3:
                style = dict(color = 'k', ls = '--', lw = 0.5, 
                             zorder = 0, clip_on = False)
                ax.plot(x, y, **style)
                ax.plot(x, -y, **style)
        else:
            # These are the ellipse boundaries
            if draw_ellipses and scale < 3:
                style = dict(color = 'k', ls = '-', lw = 0.5, 
                             zorder = 0, clip_on = False)
                ax.plot(x, y, **style)
                ax.plot(x, -y, **style)
            
        # Fill the ellipses
        if theta < 0:
            ax.fill_between(x, -y, y, color = color(zarr[i+1]), 
                            zorder = 0.5 * (z / np.pi - 1), clip_on = False)
        else:
            ax.fill_between(x, -y, y, color = color(zarr[i]), 
                            zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
    
        # Fill the ellipses that are cut off by the limb
        if theta < 0:
            x_ = np.linspace(-1, xE - xlimb, 1000)
            y_ = np.sqrt(1 - x_ ** 2)
            ax.fill_between(x_, -y_, y_, color = color(zarr[i]), 
                            zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
        else:
            x_ = np.linspace(-1, xE - xlimb, 1000)
            y_ = np.sqrt(1 - x_ ** 2)
            ax.fill_between(x_, -y_, y_, color = color(zarr[i]), 
                            zorder = 0.5 * (-z / np.pi - 1), clip_on = False)
    
    return fig, ax, occ, xy

def DrawOrbit(radiancemap = RadiativeEquilibriumMap(), inc = 70., Omega = 0., 
              ecc = 0., w = 0., Phi = 0., Lambda = 0., nphases = 20, size = 1, 
              draw_orbit = True, draw_orbital_vectors = True, 
              plot_phasecurve = False, label_phases = False, 
              figsize = (8, 8), **kwargs):
    '''
    Draw an "eyeball" planet's orbit on the sky, illustrating the changing 
    phases over the course of the orbit.
    
    .. plot::
         :align: center
     
         from planetplanet import DrawOrbit
         import matplotlib.pyplot as pl
         DrawOrbit(Omega = 45., size = 2, figsize = (6, 6))
         pl.show()
    
    :param float inc: The orbital inclination in degrees. Default `70.`
    :param float Omega: The longitude of ascending node in degrees. \
           Default `0.`
    :param float ecc: The orbital eccentricity. Default `0.`
    :param float w: The longitude of pericenter in degrees. Default `0.`
    :param float Phi: The latitudinal hot spot offset in degrees. Default `0.`
    :param float Lambda: The longitudinal hot spot offset in degrees. \
           Default `0.`
    :param int nphases: The number of planets to draw at different phases. \
           Default `20`
    :param float size: The size of the planets in arbitrary units. Default `1.`
    :param bool draw_orbit: Draw the orbit outline? Default :py:obj:`True`
    :param bool draw_orbital_vectors: Draw the orbital radial vectors? \
           Default :py:obj:`True`
    :param bool plot_phasecurve: Compute and plot the phase curve for one \
           orbit? Default :py:obj:`False`
    :param bool label_phases: Label each of the phases? Default :py:obj:`False`
    :param tuple figsize: The size of the figure in inches. \
           Default :py:obj:`(8, 8)`
    :param dict kwargs: Any other :py:obj:`kwargs` to be passed to \
           :py:func:`DrawEyeball`
    
    :returns: `fig`, `axes`, and optionally `figphase`, `axphase`; these are \
              all the figure and axes objects generated by the function
    
    '''
        
    # Convert the hotspot angles to radians
    Lambda *= np.pi / 180
    Phi *= np.pi / 180
    
    # Get the orbital elements over a full orbit of the planet
    # We are assuming a period of 10 days, but it doesn't matter for the plot!
    # We make the star tiny so that secondary eclipse is negligible
    from . import Star, Planet, System
    star = Star('A', r = 0.01)
    b = Planet('b', per = 10., inc = inc, Omega = Omega, t0 = 0, ecc = ecc, 
               w = w, Phi = Phi, Lambda = Lambda, airless = True, 
               phasecurve = True)
    system = System(star, b, mintheta = 0.001)
    time = np.linspace(-5, 5, 1000)
    if plot_phasecurve:
        system.compute(time)
    else:
        system.compute_orbits(time)

    # Plot stuff
    fig, ax = pl.subplots(1, figsize = figsize)
    ax.margins(0.1, 0.1)
    
    # Phase curve
    if plot_phasecurve:
        figphase, axphase = pl.subplots(1, figsize = (8, 2))
        figphase.subplots_adjust(bottom = 0.3)
        axphase.plot(np.linspace(0, 1, len(b.time)), b.flux[:,0] 
                     / (np.nanmax(b.flux[:,0])), 'k-')
        axphase.set_xlabel('Orbital phase', fontweight = 'bold', fontsize = 12)
        axphase.set_ylabel('Relative flux', fontweight = 'bold', fontsize = 12)
            
    # Orbit outline
    if draw_orbit:
        ax.plot(b.x, b.y, 'k-', alpha = 0.5)
    
    # Adjust the figure dimensions so the aspect ratio is unity
    left = 0.125
    right = 0.9
    xmin, xmax = ax.get_xlim()
    if (xmax - xmin) < 2:
        xmax = 1
        xmin = -1
    ymin, ymax = ax.get_ylim()
    if (ymax - ymin) < 2:
        ymax = 1
        ymin = -1
    width = right - left
    height = width * (ymin - ymax) / (xmin - xmax)
    bottom = 0.5 - height / 2
    top = 0.5 + height / 2
    fig.subplots_adjust(left = left, right = right, bottom = bottom, top = top)
    ax.axis('off')

    # Get the indices of the images we'll plot, sorted by zorder
    inds = np.array(list(range(0, 1000, 1000 // nphases)), dtype = int)
    inds = inds[np.argsort([-b.z[i] for i in inds])]

    # Plot images at different phases
    axes = [ax]
    for i in inds:
        
        # Get the eyeball angles
        theta, gamma = GetAngles(b.x[i], b.y[i], b.z[i], b.vx[i], b.vy[i], 
                                 b.vz[i], Lambda = Lambda, Phi = Phi)

        # Plot the radial vector
        if draw_orbital_vectors:
            ax.plot([0, b.x[i]], [0, b.y[i]], 'k-', alpha = 0.5, lw = 1)
    
        # Get the figure coordinates of the point
        disp_coords = ax.transData.transform((b.x[i], b.y[i]))
        xf, yf = fig.transFigure.inverted().transform(disp_coords)
    
        # Draw the planet
        _, tmp, _, _ = DrawEyeball(xf, yf, 0.015 * size, radiancemap, 
                       theta = theta, gamma = gamma, fig = fig, **kwargs)
        axes.append(tmp)
        
        # Indicate the orbital phase
        if label_phases:
            dx = b.x[i] #/ r
            dy = b.y[i] #/ r
            dr = np.sqrt(dx ** 2 + dy ** 2)
            dx *= 16 * size / dr
            dy *= 16 * size / dr
            tmp.annotate("%.2f" % (i / 1000.), xy = (0, 0), xytext = (dx, dy), 
                         xycoords = 'data', textcoords = 'offset points', 
                         fontsize = 8, ha = 'center', va = 'center',
                         zorder = 10000)

    if plot_phasecurve:
        return fig, axes, figphase, axphase
    else:
        return fig, axes
 
class Interact(object):
    '''
    Generates an interactive "eyeball" planet viewer, where the user can 
    change the orbital phase and the latitudinal/longitudinal hot spot 
    offset angles.
    
    .. plot::
         :align: center
     
         from planetplanet.photo.eyeball import Interact
         Interact()
    
    :param dict kwargs: Any :py:obj:`kwargs` to be passed to \
           :py:func:`DrawEyeball`
    
    '''
    
    def __init__(self, **kwargs):
        '''
    
        '''
    
        self.kwargs = kwargs
        self.fig = pl.figure(figsize = (6,6))
        self.fig.subplots_adjust(bottom = 0.25)
        self.axtheta = pl.axes([0.3, 0.05, 0.44, 0.03])
        self.theta = Slider(self.axtheta, r'$\theta$', -180., 180., 
                            valinit = 90.)
        self.axphi = pl.axes([0.3, 0.1, 0.44, 0.03])
        self.phi = Slider(self.axphi, r'$\Phi$', -90, 90., valinit = 0.)
        self.axlam = pl.axes([0.3, 0.15, 0.44, 0.03])
        self.lam = Slider(self.axlam, r'$\Lambda$', -90, 90., valinit = 0.)
        self.theta.on_changed(self._update)
        self.phi.on_changed(self._update)
        self.lam.on_changed(self._update)
        self._update(90.)
        pl.show()
    
    def _update(self, val):
        '''
    
        '''
    
        # Remove all axes except sliders
        for ax in self.fig.get_axes():
            if ax not in [self.axtheta, self.axphi, self.axlam]:
                ax.remove()
    
        # Get the angles
        theta = self.theta.val * np.pi / 180
        Lambda = self.lam.val * np.pi / 180
        Phi = self.phi.val * np.pi / 180
    
        # The coordinates of the substellar point
        x_ss = -np.cos(theta + Lambda) * np.cos(Phi)
        y_ss = np.sin(Phi)

        # The rotation angle that makes the planet symmetric about the x-axis
        gamma = -np.arctan2(y_ss, -x_ss)

        # Compute the new effective theta
        if theta + Lambda < -np.pi:
            theta = 2 * np.pi - np.arccos(-x_ss * np.cos(gamma) - y_ss 
                    * np.sin(gamma))
        elif theta + Lambda < 0:
            theta = -np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
        elif theta + Lambda > np.pi:
            theta = -np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
        else:
            theta = np.arccos(-x_ss * np.cos(gamma) - y_ss * np.sin(gamma))
        
        # Plot the planet
        DrawEyeball(0.525, 0.6, 0.3, RadiativeEquilibriumMap(), 
                    fig = self.fig, theta = theta, gamma = gamma, 
                    **self.kwargs)