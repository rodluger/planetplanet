#include <stdio.h>
#include <math.h>

// Models
#define MDFAST                  9
#define NEWTON                  10

// Errors
#define ERR_NONE                0                                                     // We're good!
#define ERR_NOT_IMPLEMENTED     1                                                     // Function/option not yet implemented
#define ERR_KEPLER              3                                                     // Error in the Kepler solver; probably didn't converge
#define ERR_BAD_ECC             5                                                     // Bad value for eccentricity
#define ERR_RADIUS              9                                                     // Bad input radius
#define ERR_INC                 10                                                    // Bad inclination
#define ERR_PER                 13                                                    // Bad period
#define ERR_RHOS_ARS            14                                                    // Must specify either rhos or aRs!
#define ERR_RHOS                15                                                    // Bad rhos
#define ERR_ECC_W               16                                                    // Bad eccentricity/omega

// Arrays
#define ARR_M                   2
#define ARR_E                   3
#define ARR_F                   4
#define ARR_R                   5
#define ARR_X                   6
#define ARR_Y                   7
#define ARR_Z                   8
#define ARR_B                   9

// Constants
#define PI                      acos(-1.)
#define G                       6.672e-8
#define DAYSEC                  86400.

// Structs
typedef struct {
  double inc;
  double esw;
  double ecw;
  double per;
  double ecc;
  double w;
  double a;
  double t0;
} PLANET;

typedef struct {
  int npts;
  double *M;
  double *E;
  double *f;
  double *r;
  double *x;
  double *y;
  double *z;
  double *b; 
} ARRAYS;

typedef struct {
  double keptol;
  int maxkepiter;
  int computed;
  int kepsolver;
} SETTINGS;

// Functions
double sgn(double x);
double TrueAnomaly(double E, double ecc);
double EccentricAnomalyFast(double dMeanA, double dEcc, double tol, int maxiter);
double EccentricAnomaly(double dMeanA, double dEcc, double tol, int maxiter);
int Compute(int npts, double *time, PLANET *planet, SETTINGS *settings, ARRAYS *arr);
void dbl_free(double *ptr);