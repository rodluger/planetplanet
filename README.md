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

A general photodynamical code for modeling exoplanet transits, secondary eclipses, phase curves, and exomoons, as well as eclipsing binaries, circumbinary planets, and more. The code was originally developed to model planet-planet occultation (PPO) light curves for the TRAPPIST-1 system. During a PPO, a planet
occults (transits) the disk of another planet in the same planetary system, blocking its thermal
(and reflected) light, which can be measured photometrically by a distant observer.

Here's a wacky example of what the code can do: a transits of a circumbinary exomoon!

<div align="center">
<img src="https://rodluger.github.io/misc/cbexomoon.gif" alt="Circumbinary exomoon" width="500px">
</div>

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

Please check out the [documentation](https://rodluger.github.io/planetplanet/ndex.html) or read the [paper](https://rodluger.github.io/planetplanet/PPOs.pdf) for more information. If you find something is amiss, please submit an [issue](https://github.com/rodluger/planetplanet/issues)!
