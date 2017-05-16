#include <stdio.h>
#include <math.h>

// Constants
#define PI                      acos(-1.)
#define MAXVERTICES             100
#define MAXFUNCTIONS            100
#define DTOL1                   1.e-10
#define DTOL2                   1.e-15

// Structs
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

double SurfaceBrightness(double lat, double noon, double midnight, int n);
double DeltaFlux(double r, double x0, double y0, double ro, double theta, double noon, double midnight, int n);
double Flux(double r, double theta, double noon, double midnight, int n);