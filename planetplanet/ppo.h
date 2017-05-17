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
int OrbitXYZ(double time, PLANET *planet, SETTINGS *settings);
double SurfaceBrightness(double lat, double noon, double midnight, int n);
double OccultedFlux(double r, double x0, double y0, double ro, double theta, double noon, double midnight, int n);
double UnoccultedFlux(double r, double theta, double noon, double midnight, int n);
double Flux(double time, PLANET *planet1, PLANET *planet2, SETTINGS *settings);