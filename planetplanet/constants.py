# Constants
G = 6.67428e-11
HPLANCK = 6.62607004e-34
CLIGHT = 2.998e8
KBOLTZ = 1.38064852e-23
SBOLTZ = 5.670367e-8

# Units
MSUN = 1.988416e30
LSUN = 3.846e26
RSUN = 6.957e8
AUM = 1.49598e11
PARSEC = 3.086e16
MEARTH = 5.9722e24
REARTH = 6.3781e6
SEARTH = 1.361e3

# Time
DAYSEC = 86400.
SECOND = 1. / DAYSEC
MINUTE = 1. / 1440.

# Scalings
AUREARTH = AUM / REARTH
MSUNMEARTH = MSUN / MEARTH
RSUNREARTH = RSUN / REARTH
GEARTH = G * DAYSEC ** 2 * MEARTH / REARTH ** 3

# C definitions
MDFAST = 0
NEWTON = 1
REB_INTEGRATOR_WHFAST = 1
REB_INTEGRATOR_IAS15 = 0