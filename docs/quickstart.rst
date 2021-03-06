Tutorial
========

In this section I'll go over some of the :py:obj:`planetplanet` basics, like instantiating
a planetary system, computing orbits and light curves, and doing some simple plotting. Users are
encouraged to also check out the `scripts <scripts.html>`_ page for a collection of examples,
including several of the scripts used to plot the figures in the paper. Those scripts cover most
of the functionality of :py:obj:`planetplanet`.

Instantiating a planetary system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After `installing <install.html>`_, you can easily import :py:obj:`planetplanet` in Python:

.. code-block:: python

   >>> import planetplanet as pp

The simplest thing you can do is to instantiate a one-planet system and compute some transit and
secondary eclipse light curves:

.. code-block:: python

   >>> star = pp.Star('star', m = 0.1, r = 0.1)
   >>> planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)

Here we've instantiated a star with mass 0.1 :math:`\mathrm{M_\odot}` and radius 0.1 :math:`\mathrm{R_\odot}`
-- corresponding to an M7 (ish) dwarf. We've also given it the name "star", which we will use to access
it later. We instantiated the planet with mass 1 :math:`\mathrm{M_\oplus}`, radius 1 :math:`\mathrm{R_\oplus}` --
i.e., an Earth-size planet. We gave it an inclination of 89.5 degrees and specified that the time of transit :py:obj:`t0`
is at :math:`t = 4`. Note that by default, :py:class:`Star <planetplanet.photo.structs.Star>` parameters
are specified in Solar units, while :py:class:`Planet <planetplanet.photo.structs.Planet>` parameters are specified
in Earth units. Check out the :py:class:`Star <planetplanet.photo.structs.Star>` and :py:class:`Planet <planetplanet.photo.structs.Planet>`
classes for a description of all the available keywords and their default settings.

Next, we instantiate a :py:class:`System <planetplanet.photo.ppo.System>` object, which is the main gateway to the C
photodynamical routines:

.. code-block:: python

   >>> system = pp.System(star, planet, distance = 10)

We're telling :py:obj:`planetplanet` that this planetary system is composed of a star and a planet, and that it's at a distance
of 10 parsecs. Note that the :py:class:`Star <planetplanet.photo.structs.Star>` instance must always come first; currently, only
one :py:class:`Star <planetplanet.photo.structs.Star>` is allowed per system, as :py:obj:`planetplanet` can't yet handle binary
systems (but stay tuned!). Check out the :py:class:`System <planetplanet.photo.ppo.System>` class for a list of all available parameters.

Computing orbits and light curves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Now let's compute the light curve over the span of ten days, at a cadence of 1.44 minutes:

.. code-block:: python

   >>> import numpy as np
   >>> time = np.arange(0, 10, 0.001)
   >>> system.compute(time)
   Computing orbits with the Kepler solver...
   [==================================================] 100% 1ms
   Computing occultation light curves...
   Done!

Several things just happened. First, :py:obj:`planetplanet` computed the orbital solution for the system over the given
time array using a Keplerian solver and stored the planet's Cartesian coordinates in the :py:obj:`x`, :py:obj:`y`, and
:py:obj:`z` attributes:

.. code-block:: python

   >>> planet.x
   array([-383.81118951, -384.2744517 , -384.73602829, ...,
            -2.78460764,   -1.85641188,   -0.92820798])
   >>> planet.y
   array([ 1.93374349,  1.92672441,  1.91969687, ...,
          -3.86741063, -3.86745305, -3.86747849])
   >>> planet.z
   array([ 221.58505598,  220.78074891,  219.97547339, ...,
          -443.16136416, -443.16622404, -443.16913998])

.. warning:: :py:obj:`planetplanet` uses a **left-handed** Cartesian coordinate system. This is somewhat unconventional, \
             but it is convenient in that the **x** axis points to the right on the sky, the **y** axis points up, and \
             the **z** axis points *into* the sky. The observer is thus always at z = :math:`-\infty`. \
             Prograde orbits proceed counter-clockwise when looking down the **y** \
             axis. In practice this choice doesn't matter, since the absolute sense of the orbit (and whether a planet is to the left \
             or to the right of the star) cannot usually be established from photometry.

We can view the orbit by running

.. code-block:: python

   >>> system.plot_orbits()

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   time = np.arange(0, 10, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time)
   system.plot_orbits()
   pl.show()

The code also computed light curves for all occultation events. These are stored in the :py:obj:`flux`
attributes of each of the bodies and in the :py:class:`System <planetplanet.photo.ppo.System>` instance:

.. code-block:: python

   >>> system.flux
   array([[  9.01190198e-15,   8.68430010e-15,   8.36835988e-15, ...]])
   >>> star.flux
   array([[  9.01190198e-15,   8.68430010e-15,   8.36835988e-15, ...]])
   >>> planet.flux
   array([[  0.,  0.,  0., ...]])

The :py:obj:`system.flux` attribute is the sum of the light curves of all bodies in the system. Note
that :py:obj:`flux` is a two-dimensional array:

.. code-block:: python

   >>> system.flux.shape
   (10000, 112)

The first axis is the time axis (we computed stuff over 10,000 cadences); the second axis is the
wavelength axis. Though we didn't specify it, :py:obj:`planetplanet` computed the light curve
for the system over a grid of wavelengths:

.. code-block:: python

   >>> pl.plot(system.wavelength, system.flux[0, :])

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   time = np.arange(0, 10, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time)
   pl.plot(system.wavelength, system.flux[0, :])
   pl.xlabel(r'Wavelength [$\mu$ m]', fontweight = 'bold')
   pl.ylabel(r'Flux [W / m$^2$]', fontweight = 'bold')
   pl.show()

What we see is the Rayleigh-Jeans tail of the stellar flux. If we poke around in the docs, we see that
the default effective temperature for a :py:obj:`Star <planetplanet.photo.structs.Star>` instance is 5577 K,
so this is the Sun's blackbody spectrum. By default, :py:obj:`planetplanet` computes light curves in the range
5 - 15 :math:`\mu\mathrm{m}` at a resolution `R = 100`. We can change these values when we call :py:obj:`compute()`:

.. code-block:: python

   >>> system.compute(time, lambda1 = 0.1, lambda2 = 15, R = 1000)
   >>> pl.plot(system.wavelength, system.flux[0, :])
   >>> pl.xscale('log')

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   time = np.arange(0, 10, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time, lambda1 = 0.1, lambda2 = 15, R = 1000)
   pl.plot(system.wavelength, system.flux[0, :])
   pl.xscale('log')
   pl.xlabel(r'Wavelength [$\mu$ m]', fontweight = 'bold')
   pl.ylabel(r'Flux [W / m$^2$]', fontweight = 'bold')
   pl.show()

Now let's look at the system light curve at a given wavelength over the entire time array:

.. code-block:: python

   >>> system.plot_lightcurve(wavelength = 15)

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   time = np.arange(0, 10, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time, lambda1 = 0.1, lambda2 = 15, R = 1000)
   system.plot_lightcurve(wavelength = 15)
   pl.show()

Three transits and three secondary eclipses are clearly visible. Clicking on one of the
events brings up an interactive window. This is what you get when you click on a transit:

.. image:: /transit.gif
   :width: 400px
   :align: center

At the top, we see the orbital configuration at the time of transit (observer toward the
bottom); in the middle, we see an animation of the event; and at the bottom, the event
light curve. Note that the image of the star reflects the limb darkening parameters :py:obj:`planetplanet` assumed,
`u1 = 1` and `u2 = 0` (see :py:obj:`Star <planetplanet.photo.structs.Star>`). By default,
limb darkening coefficients can be specified up to any order as a Taylor expansion in the
quantity :math:`(1 - \mu) = (1 - \cos\phi)`, where :math:`\phi` is the viewing angle:

.. math::

   B_\lambda(\mu) = B_\lambda^0 \left[ 1 - \sum_{i=1}^{n} u_i(\lambda) (1 - \mu)^i \right]

Different limb darkening parameters can be specified when instantiating a
:py:obj:`Star <planetplanet.photo.structs.Star>` object via the :py:obj:`limbdark`
keyword. Note that :py:obj:`planetplanet` also allows users to specify wavelength-dependent limb
darkening coefficients. See `this script <scripts/mutual_transit.html>`_ for an example.

Computing phase curves
~~~~~~~~~~~~~~~~~~~~~~

We can also plot the phase curve for the planet in the examples above, but we would have
needed to specify `phasecurve = True` when instantiating it. Alternatively, we can just
set the attribute directly:

.. code-block:: python

   >>> planet.phasecurve = True
   >>> system.compute(time, lambda1 = 0.1, lambda2 = 15, R = 1000)
   >>> system.plot_lightcurve(wavelength = 10)

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   planet.phasecurve = True
   time = np.arange(0, 10, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time, lambda1 = 0.1, lambda2 = 15, R = 1000)
   system.plot_lightcurve(wavelength = 15)
   pl.show()

By default, planets are instantiated with the :py:obj:`RadiativeEquilibriumMap <planetplanet.photo.maps.RadiativeEquilibirumMap>`
surface map with an :py:class:`albedo <planetplanet.photo.structs.Planet>` of `0.3` and a nightside temperature
:py:class:`tnight <planetplanet.photo.structs.Planet>` of `40` K. The latter two can be passed as keywords to the
:py:class:`Planet <planetplanet.photo.structs.Planet>` class or specified directly by setting the respective attributes. The
:py:obj:`RadiativeEquilibriumMap <planetplanet.photo.maps.RadiativeEquilibirumMap>` surface map computes radiances assuming
instant re-readiation, which is valid for planets with atmospheres that have negligible thermal inertia and negligible recirculation
(i.e., the airless planet limit). These are "eyeball" planets, which look like `this <scripts/eyeball_orbit.html>`_.
Alternatively, users may specify

.. code-block:: python

   >>> planet.radiancemap = pp.LimbDarkenedMap()

to treat the planet as a limb-darkened body, whose emission is always symmetric about the center of the planet disk, regardless of
the orbital phase or viewing angle. In this case, no phase curve is visible:

.. code-block:: python

   >>> planet.phasecurve = True
   >>> system.plot_lightcurve(wavelength = 10)

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   planet.phasecurve = True
   planet.radiancemap = pp.LimbDarkenedMap()
   time = np.arange(0, 10, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time, lambda1 = 0.1, lambda2 = 15, R = 1000)
   system.plot_lightcurve(wavelength = 15)
   pl.show()

Because of the generalized integration scheme in :py:obj:`planetplanet`, users can also specify custom surface maps, provided they
are radially symmetric about the hotspot (which need not point toward the star!). Check out this `script <scripts/custom_map.html>`_.

.. plot::
   :align: center

   from scripts import custom_map
   import matplotlib.pyplot as pl
   custom_map.view_planet()
   pl.show()

Computing planet-planet occultations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Finally, let's look at how we compute planet-planet occultation (PPO) light curves for the TRAPPIST-1 system. The
:py:mod:`trappist1 <planetplanet.photo.trappist1>` module contains utilities for instantiating the TRAPPIST-1
planetary system:

.. code-block:: python

   >>> system = pp.Trappist1(sample = True, seed = 543210)

Note that I specified :py:obj:`sample = True`, meaning we will draw at random from the prior on all the orbital
parameters (if :py:obj:`sample = False`, the mean values are used for all parameters). The prior is currently informed
by the observational constraints on the system from Gillon et al. (2017),
Luger et al. (2017), and Burgasser & Mamajek (2017), as well as on the Monte Carlo simulations in Luger, Lustig-Yaeger
and Agol (2017) for the planets' mutual inclinations. First, let's look at the orbits:

.. code-block:: python

   >>> time = np.arange(0, 10, 0.001)
   >>> system.compute(time)
   >>> system.plot_orbits()

.. plot::
   :align: center

   import planetplanet as pp
   import matplotlib.pyplot as pl
   import numpy as np
   time = np.arange(0, 10, 0.001)
   system = pp.Trappist1(sample = True)
   system.compute(time)
   system.plot_orbits()
   pl.show()

The plot in the lower left-hand corner is the view from Earth. We can also plot the full system light curve as before:

.. code-block:: python

   >>> system.plot_lightcurve(wavelength = 15)

.. plot::
   :align: center

   import planetplanet as pp
   import matplotlib.pyplot as pl
   import numpy as np
   time = np.arange(0, 10, 0.001)
   system = pp.Trappist1(sample = True, seed = 543210)
   system.compute(time)
   system.plot_lightcurve(wavelength = 15)
   pl.show()

All the occultations are labeled: the label text indicates the occulted body, while the label color indicates the
occultor (black, red, orange, yellow, green, aqua, light blue, dark blue, for the star and each of the seven planets,
respectively). Transits are thus labeled with "A" (for the star) and secondary eclipses are labeled with the planet
name and colored black. There are several interesting features in the light curve: a simultaneous secondary eclipse of
b and c at 1.6 days, a simultaneous transit of b and e at 5.3 days, and two prominent planet-planet occultations of c by
b at 1.3 and 1.87 days. If we click on the one at 1.87 days, we can see its light curve:

.. image:: /ppo.gif
   :width: 400px
   :align: center

Note that by default we display the flux normalized to the planet's total emission at full phase. The
:py:obj:`flux` attribute of planet c contains the actual flux in :math:`\mathrm{W/m^2}` if that's what
you need.

Hunting for occultations
~~~~~~~~~~~~~~~~~~~~~~~~

The last thing we'll go over here is how to use :py:obj:`planetplanet` to predict when PPOs occur. We
make this easy via the :py:obj:`next_occultation() <planetplanet.photo.System.next_occultation>` method,
which returns the times (and durations) of the next :py:obj:`N` occultations of a given body in the system.
This `example script <scripts/next_occultation.html>`_ shows how to do this. You can easily sort the results
to find the longest upcoming occultation, which for a given instance of the system can reveal fun events like
this prograde-retrograde occultation of TRAPPIST-1c by TRAPPIST-1b lasting *5 hours*!

.. image:: /retro.gif
   :width: 400px
   :align: center

Simulating observations
~~~~~~~~~~~~~~~~~~~~~~~

Now we will cover how to compute synthetic observations of :py:obj:`planetplanet`
light curves with the James Webb Space Telescope (JWST)  Mid-Infrared Instrument
(MIRI) Imager. After running :py:obj:`compute()` on a :py:class:`System <planetplanet.photo.ppo.System>`
instance, simply call :py:obj:`observe()`:

.. code-block:: python

   >>> fig, ax = system.observe(stack = 1, filter = 'f1500w', alpha_err = 0.5)
   Computing observed light curve in F1500W filter...
   Average lightcurve precision: 222.090 ppm

.. plot::
   :align: center

   import planetplanet as pp
   import numpy as np
   import matplotlib.pyplot as pl
   star = pp.Star('star', m = 0.1, r = 0.1)
   planet = pp.Planet('planet', m = 1., r = 1., per = 3., inc = 89.5, t0 = 4.)
   planet.phasecurve = True
   time = np.arange(0, 9.5, 0.001)
   system = pp.System(star, planet, distance = 10)
   system.compute(time, lambda1 = 0.1, lambda2 = 20, R = 1000)
   fig, ax = system.observe(stack = 1, filter = 'f1500w', alpha_err = 0.5)
   pl.show()

.. warning:: Make sure that :py:obj:`compute()` was run with a wavelength range
             that fully covers any filters you may want to use. Otherwise an
             `AssertionError` will be thrown when you run :py:obj:`observe()`.

We can explore all the available MIRI filters and plot their throughput curves:

.. code-block:: python

   >>> miri_filters = pp.jwst.get_miri_filter_wheel()
   >>> for filter in miri_filters: print("%s : %.1f um" %(filter.name, filter.eff_wl))
   F560W : 5.6 um
   F770W : 7.7 um
   F1000W : 10.0 um
   F1130W : 11.3 um
   F1280W : 12.8 um
   F1500W : 15.1 um
   F1800W : 18.0 um
   F2100W : 20.8 um
   F2550W : 25.4 um
   >>> miri_filters[4].plot()

.. plot::
   :align: center

   import planetplanet as pp
   miri_filters = pp.jwst.get_miri_filter_wheel()
   miri_filters[4].plot()

:py:class:`Filter <planetplanet.detect.jwst.Filter>` objects can also be passed
to :py:obj:`observe()`:

.. code-block:: python

   >>> fig, ax = system.observe(filter = miri_filters[4])

After :py:obj:`observe()` has run, the :py:class:`system <planetplanet.photo.ppo.System>`
will have the :py:class:`Filter <planetplanet.detect.jwst.Filter>` as an attribute
with an instantiated :py:class:`Lightcurve <planetplanet.detect.jwst.Lightcurve>`
containing various useful quantities. For instance, the number of photons detected
from the system `Nsys` and the number of background photons `Nback` are given by:

.. code-block:: python

   >>> system.filter.lightcurve.Nsys
   array([ 23620244.09945266,  23620230.07335656,  23620216.03934677, ...,
        23613770.63008701,  23613761.01311478,  23613751.41684124])
   >>> system.filter.lightcurve.Nback
   array([ 3892306.56855529,  3892306.56855529,  3892306.56855529, ...,
        3892306.56855313,  3892306.56856005,  3892306.56856005])

which are a function of `time`:

.. code-block:: python

   >>> system.filter.lightcurve.time
   array([  0.00000000e+00,   1.00000000e-03,   2.00000000e-03, ...,
         9.49700000e+00,   9.49800000e+00,   9.49900000e+00])


More examples
~~~~~~~~~~~~~

That's all for the tutorial (for now), though we'll keep adding features to :py:obj:`planetplanet` and
posting updated info here. Make sure to check out the `scripts <scripts.html>`_ page for more examples
and some advanced usage information. And there's always the `API <api.html>`_ if you're feeling adventurous
or would like to adapt :py:mod:`planetplanet` for your needs!
