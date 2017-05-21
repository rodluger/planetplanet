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
#define KBOLTZ                  1.38064852e-23                                        // W s / K
#define SBOLTZ                  5.670367e-8                                           // W / m^2 / K^4
#define SEARTH                  1.361e3                                               // W / m^2
#define HPLANCK                 6.62607004e-34                                        // m^2 kg / s
#define CLIGHT                  2.998e8                                               // m / 2
#define REARTH                  6371000.                                              // m

// Settings
#define MAXVERTICES             200
#define MAXFUNCTIONS            200
#define DTOL1                   1.e-10
#define DTOL2                   1.e-15

// Structs
typedef struct {
  double per;
  double inc;
  double ecc;
  double w;
  double a;
  double t0;
  double r;
  double albedo;
  double irrad;
  int phasecurve;
  int nl;
  int t;
  int nt;
  int nw;
  double *time;
  double *wavelength;
  double *x;
  double *y;
  double *z;
  int *occultor;
  double **flux;
} PLANET;

typedef struct {
  double keptol;
  int maxkepiter;
  int kepsolver;
  double polyeps1;
  double polyeps2;
  int maxpolyiter;
} SETTINGS;

typedef struct {
  int circle;
  double r;
  double a;
  double b;
  double x0;
  double y0;
  double xmin;
  double xmax;
} ELLIPSE;

typedef double (*CURVE)(double, ELLIPSE*);

typedef double (*INTEGRAL)(double, double, ELLIPSE*);

typedef struct {
  double y;
  ELLIPSE *ellipse;
  CURVE curve;
  INTEGRAL integral;
} FUNCTION;

// Functions
int OrbitXYZ(int np, PLANET **planet, SETTINGS settings);
void OccultedFlux(double r, double x0, double y0, double ro, double theta, double albedo, double irrad, double polyeps1, double polyeps2, int maxpolyiter, int nlat, int nlam, double lambda[nlam], double flux[nlam]);
void UnoccultedFlux(double r, double theta, double albedo, double irrad, double polyeps1, double polyeps2, int maxpolyiter, int nlat, int nlam, double lambda[nlam], double flux[nlam]);
int Flux(int nt, double time[nt], int nw, double wavelength[nw], int np, PLANET **planet, SETTINGS settings);