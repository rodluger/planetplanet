<div align="center">
<img src="https://rodluger.github.io/planetplanet/_images/title.gif" width="400px">
</img>
<br/><br/>
<p><a href="https://travis-ci.org/rodluger/planetplanet"><img src="https://travis-ci.org/rodluger/planetplanet.svg?branch=master"/></a>
<a href="http://dx.doi.org/10.5281/zenodo.997391"><img src="https://img.shields.io/badge/doi-zenodo-568AB8.svg?style=flat"/></a>
<a href="https://raw.githubusercontent.com/rodluger/planetplanet/master/LICENSE?token=AI5FKxGMJTv55h2EE_AuXW2gofnIaRDeks5Zm0unwA%3D%3D"><img src="https://img.shields.io/badge/license-GPL-a2a2a2.svg?style=flat"/></a>
<a href="https://rodluger.github.io/planetplanet/PPOs.pdf"><img src="https://img.shields.io/badge/read-the_paper-fd7709.svg?style=flat"/></a>
<a href="https://rodluger.github.io/planetplanet/index.html"><img src="https://img.shields.io/badge/read-the_docs-AF5891.svg?style=flat"/></a>
</p>
</div>

# Overview
`planetplanet` is a general photodynamical code for modeling exoplanet transits, secondary eclipses, phase curves, and exomoons, as well as eclipsing binaries, circumbinary planets, and more. The code was originally developed to model planet-planet occultation (PPO) light curves for the TRAPPIST-1 system. During a PPO, a planet
occults (transits) the disk of another planet in the same planetary system, blocking its thermal
(and reflected) light, which can be measured photometrically by a distant observer.

`planetplanet` is coded in C and wrapped in a user-friendly Python interface. Once installed, generating light curves is as easy as

```python
import planetplanet as pp
import numpy as np
import matplotlib.pyplot as pl

star = pp.Star('A', m = 0.1, r = 0.1, limbdark = [0.4, 0.26])
planet = pp.Planet('b', m = 1, r = 1, t0 = 0., per = 3.)
system = pp.System(star, planet)
time = np.arange(-1., 1., 0.0001)
system.compute(time)
system.plot_occultation('A', time = 0)
pl.show()
```

<div align="center">
<img src="https://rodluger.github.io/misc/transit.gif" alt="Exoplanet transit light curve" width="500px">
</div>

Please check out the [documentation](https://rodluger.github.io/planetplanet/index.html) or read the [paper](https://rodluger.github.io/planetplanet/PPOs.pdf) for more information.

# Installation
The `planetplanet` code is now `pip`-installable:

```
pip install planetplanet
```

Alternatively, to install from source:

```
git clone git@github.com:rodluger/planetplanet.git
cd planetplanet
git submodule init && git submodule update
python setup.py develop
```

Note that you may need to install the [GNU Scientific Library](https://www.gnu.org/software/gsl/). On a Mac, it's as simple as

```
brew install gsl
```

# Just for fun
Here's a an example of a planet-planet occultation [**[code]**](https://github.com/rodluger/planetplanet/blob/master/scripts/occultation.py):

<div align="center">
<img src="https://rodluger.github.io/misc/ppo.gif" alt="Planet-planet occultation" width="500px">
</div>

And here's a wacky example of a transit of a circumbinary exomoon [**[code]**](https://github.com/rodluger/planetplanet/blob/master/scripts/circumbinary_exomoon.py):

<div align="center">
<img src="https://rodluger.github.io/misc/cbexomoon.gif" alt="Circumbinary exomoon" width="500px">
</div>
