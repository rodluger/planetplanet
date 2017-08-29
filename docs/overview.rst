Overview
========

The :py:obj:`planetplanet` software package was designed to compute light curves of
planet-planet occultations (PPOs) in exoplanet systems. During a PPO, a planet
occults (transits) the disk of another planet in the same system, blocking its thermal
(and reflected) light, which can be measured photometrically by a distant observer.
We developed this package with the TRAPPIST-1 planetary system in mind, but :py:obj:`planetplanet`
is generally applicable to any exoplanet system.

Because of the general way in which :py:obj:`planetplanet` computes occultation light curves,
it natively supports the computation of transit light curves, secondary eclipse light curves,
and planetary phase curves, as well as occultations of planets by moons and mutual transits of planets
across the face of their host star. Check out our `paper <PPOs.pdf>`_ for more information.

This documentation explains how to `install <install.html>`_ the code and use it to compute
light curves for TRAPPIST-1 and other exoplanet systems. Check out the `tutorial <quickstart.html>`_
for a quick guide on how to use the code and the `scripts page <scripts.html>`_ for additional
examples and code to reproduce several of the figures in the paper. Finally, the `API section <api.html>`_
contains detailed descriptions of the Python and C :py:obj:`planetplanet` backends.