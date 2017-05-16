#include <stdio.h>
#include <math.h>

// Models
#define MDFAST                  0
#define NEWTON                  1

// Errors
#define ERR_NONE                0                                                     // We're good!
#define ERR_NOT_IMPLEMENTED     1                                                     // Function/option not yet implemented
#define ERR_KEPLER              2                                                     // Error in the Kepler solver; probably didn't converge
#define ERR_INPUT               3                                                     // Bad input value

// Constants
#define PI                      acos(-1.)
#define G                       6.672e-8
#define DAYSEC                  86400.

// Structs
typedef struct {
  double per;
  double inc;
  double ecc;
  double w;
  double a;
  double t0;
  double r;
  double noon;
  double midnight;
  int nlat;
  double x;
  double y;
  double z;
} PLANET;

typedef struct {
  double keptol;
  int maxkepiter;
  int kepsolver;
} SETTINGS;

// Functions
double sgn(double x);
double TrueAnomaly(double E, double ecc);
double EccentricAnomalyFast(double dMeanA, double dEcc, double tol, int maxiter);
double EccentricAnomaly(double dMeanA, double dEcc, double tol, int maxiter);
int OrbitXYZ(double time, PLANET *planet, SETTINGS *settings);