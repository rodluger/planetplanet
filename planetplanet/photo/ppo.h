/**
@file ppo.h
@brief Main header for the photodynamical routines.
*/
#include <stdio.h>
#include <math.h>

// Models
#define MDFAST                  0                                                     /**< Use the Murray & Dermott fast Kepler solver */
#define NEWTON                  1                                                     /**< Use the standard Newton Kepler solver */
#define QGSL                    3                                                     /**< Use the GSL complex polynomial solver */

// Surface maps
#define MAP_NONE               -1
#define MAP_RADIAL_DEFAULT      0
#define MAP_RADIAL_CUSTOM       1
#define MAP_ELLIPTICAL_DEFAULT  2
#define MAP_ELLIPTICAL_CUSTOM   3

// Errors
#define ERR_NONE                0                                                     /**< No error occurred */
#define ERR_NOT_IMPLEMENTED     1                                                     /**< Function/option not yet implemented */
#define ERR_KEPLER              2                                                     /**< Error in the Kepler solver; probably didn't converge */
#define ERR_INPUT               3                                                     /**< Bad input value */
#define ERR_OOB                 -1                                                    /**< Warning: out of bounds */

// Constants
#define PI                      acos(-1.)                                             /**< Plain old pi */
#define BIGG                    6.67428e-11                                           /**< Gravitational constant (mks) */
#define DAYSEC                  86400.                                                /**< Number of seconds in one day */
#define KBOLTZ                  1.38064852e-23                                        /**< Boltzmann constant in W s / K */
#define SBOLTZ                  5.670367e-8                                           /**< Stefan-Boltzmann constant in W / m^2 / K^4 */
#define SEARTH                  1.361e3                                               /**< Solar constant in W / m^2 */
#define HPLANCK                 6.62607004e-34                                        /**< Planck constant in m^2 kg / s */
#define CLIGHT                  2.998e8                                               /**< Speed of light in m / s */
#define REARTH                  6.3781e6                                              /**< Radius of the Earth in m */
#define MICRON                  1e-6                                                  /**< Meters in 1 micron */
#define PARSEC                  3.086e16                                              /**< Meters in 1 parsec */
#define MEARTH                  5.9722e24                                             /**< Mass of Earth in kg */
#define GEARTH                  (BIGG * DAYSEC * DAYSEC * MEARTH / (REARTH * REARTH * REARTH)) /**< Graviational constant in Earth units */

// Settings
#define MAXIM                   1.e-2                                                 /**< Maximum magnitude of the imaginary component for root to be treated as real */
#define SMALL                   1.e-10                                                /**< Tolerance in the geometry routines */
#define TINY                    1.e-15                                                /**< Tolerance in the geometry routines */
#define MINCRESCENT            (-70. * PI / 180)                                      /**< Numerical issues pop up when the dayside crescent is too thin, so for phase angles smaller than this we use more zenith slices */
#define CRESCENTNZ              31                                                    /**< Number of zenith slices in the tiny crescent limit */

/**
A radiance map function of the wavelength and the zenith angle. 

*/
typedef double (*RADIANCEMAP)(double, double);

/**
Struct containing all the information pertaining to a body (star, planet, or moon) in the system.

*/
typedef struct {
  double m;                                                                           /**< Body mass in Earth masses */
  double per;                                                                         /**< Orbital period in days */
  double inc;                                                                         /**< Orbital inclination in radians */
  double ecc;                                                                         /**< Orbital eccentricity */
  double w;                                                                           /**< Longitude of pericenter in radians */
  double Omega;                                                                       /**< Longitude of ascending node in radians */
  double a;                                                                           /**< Orbital semi-major axis in Earth radii */
  double tperi0;                                                                      /**< Time of pericenter passage in days */
  double r;                                                                           /**< Body radius in Earth radii */
  double albedo;                                                                      /**< Body albedo */
  double teff;                                                                        /**< Effective temperature in K */
  double tnight;                                                                      /**< Night side temperature in K */
  int phasecurve;                                                                     /**< Compute the phase curve for this body? */
  double Lambda;                                                                      /**< Longitudinal hotspot offset in radians */
  double Phi;                                                                         /**< Latitudinal hotspot offset in radians */
  int host;                                                                           /**< The index of this body's host (host star if planet; host planet if moon) */
  int nu;                                                                             /**< Number of limb darkening coefficients (per wavelength) */
  int nz;                                                                             /**< Number of zenith angle slices */
  int nt;                                                                             /**< Size of time array */
  int nw;                                                                             /**< Size of wavelength array */
  double *u;                                                                          /**< Limb darkening coefficient grid */
  double *time;                                                                       /**< Time array */
  double *wavelength;                                                                 /**< Wavelength array */
  double *x;                                                                          /**< The Cartesian x position on the sky (right positive) */
  double *y;                                                                          /**< The Cartesian y position on the sky (up positive) */
  double *z;                                                                          /**< The Cartesian z position on the sky (into sky positive) */
  int *occultor;                                                                      /**< The array of occultor bit flags */
  double *flux;                                                                       /**< The grid of observed flux from this body in time/wavelength */
  double *total_flux;                                                                 /**< The total unocculted flux of this body at full phase */
  int maptype;                                                                        /**< Which radiance map to use? */
  RADIANCEMAP radiancemap;                                                            /**< A function that returns a radiance map for the body's surface */
} BODY;

/**
Struct containing all the settings for computing orbits and light curves for the system.

*/
typedef struct {
  int nbody;                                                                          /**< Use N-Body solver? */
  int integrator;                                                                     /**< Which N-Body solver to use */
  double keptol;                                                                      /**< Kepler solver tolerance */
  int maxkepiter;                                                                     /**< Maximum iterations in Kepler solver */
  int kepsolver;                                                                      /**< Which Kepler solver to use */
  double timestep;                                                                    /**< N-Body timestep in days */
  int adaptive;                                                                       /**< Adaptive zenith angle grid for limb-darkened bodies? */
  int circleopt;                                                                      /**< Treat zenith angle slices as circles in the limb-darkened limit? No reason to set this to 0 */
  int batmanopt;                                                                      /**< Use BATMAN algorithm to speed up limb-darkened occultation light curves? */
  int quarticsolver;                                                                  /**< Which quartic solver to use? */
  int quiet;                                                                          /**< Suppress output? */
  double mintheta;                                                                    /**< Minimum absolute value of the phase angle below which it is assumed to be constant to prevent numerical errors */
  int maxvertices;                                                                    /**< Maximum number of vertices (for memory allocation) */
  int maxfunctions;                                                                   /**< Maximum number of functions (for memory allocation) */
  double distance;                                                                    /**< Distance to the system in parsecs */
} SETTINGS;

/**
A generic ellipse class.

*/
typedef struct {
  int circle;                                                                         /**< Is this a circle? */
  double r;                                                                           /**< Radius (if it's a circle) */
  double a;                                                                           /**< Semi-major axis (if it's an ellipse) */
  double b;                                                                           /**< Semi-minor axis (if it's an ellipse) */
  double x0;                                                                          /**< x coordinate of the center of the ellipse */
  double y0;                                                                          /**< y coordinate of the center of the ellipse */
  double xmin;                                                                        /**< Leftmost point of the ellipse */
  double xmax;                                                                        /**< Rightmost point of the ellipse */
} ELLIPSE;

/**
A generic ellipse function.

*/
typedef double (*CURVE)(double, ELLIPSE*);

/**
A generic ellipse integral.

*/
typedef double (*INTEGRAL)(double, double, ELLIPSE*, int*);

/**
A container for the properties, functions, and integrals of an ellipse.

*/
typedef struct {
  double y;                                                                           /**< The value of the function at a given point */
  ELLIPSE *ellipse;                                                                   /**< The ellipse class corresponding to this function */
  CURVE curve;                                                                        /**< The curve function */
  INTEGRAL integral;                                                                  /**< The integral function */
} FUNCTION;

// Global functions
double Blackbody(double lambda, double T);
int NBody(int np, BODY **body, SETTINGS settings);
int Kepler(int np, BODY **body, SETTINGS settings);
void OccultedFlux(double r, int no, double x0[no], double y0[no], double ro[no], double theta, double tnight, double teff, double distance, double mintheta, int maxvertices, int maxfunctions, int adaptive, int circleopt, int batmanopt, int quarticsolver, int nu, int nz, int nw, double u[nu * nw], double lambda[nw], double flux[nw], int maptype, RADIANCEMAP radiancemap, int quiet, int *iErr);
void UnoccultedFlux(double r, double theta, double tnight, double teff, double distance, double mintheta, int maxvertices, int maxfunctions, int adaptive, int circleopt, int batmanopt, int quarticsolver, int nu, int nz, int nw, double u[nu * nw], double lambda[nw], double flux[nw], int maptype, RADIANCEMAP radiancemap, int quiet, int *iErr);
int Orbits(int nt, double time[nt], int np, BODY **body, SETTINGS settings);
int Flux(int nt, double time[nt], int nw, double wavelength[nw], int np, BODY **body, SETTINGS settings);