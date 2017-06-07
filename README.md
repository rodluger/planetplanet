# planetplanet

Generates planet-planet occultation light curves in Python. To install:

```
git clone git@github.com:rodluger/planetplanet.git
git submodule init && git submodule update
make -C planetplanet/photo
```

<p align='center'><img src="img/eyeball.png" width="800"/></p>

## TODO

- Normalize stellar limb darkening coefficients and compute planet flux self-consistently
- Add limb-darkened planet option
- Add nightside heat redistribution option
- Speed up N-body code with adaptive grid
- Compute orbital phase angle in the general case (not edge-on limit)