#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include "complex.h"

void GetRoots(double a, double b, double x0, double y0, double r, double roots[2]){
  /*
  
  */
  
  int i, j;
  double A, B, C, D;
  fcomplex c[5];
  fcomplex croots[5];
  double r2 = r * r;
  double a2 = a * a;
  double b2 = b * b;
  double a2b2 = a2 / b2;
  double y2 = y0 * y0;
  double x2 = x0 * x0;
  
  // Get the coefficients
  A = a2b2 - 1.;
  B = -2. * x0 * a2b2;
  C = r2 - y2 - a2 + a2b2 * x2;
  D = 4. * y2 * a2b2;
  c[4].r = A * A;
  c[3].r = 2. * A * B;
  c[2].r = 2. * A * C + B * B + D;
  c[1].r = 2. * B * C - 2. * D * x0;
  c[0].r = C * C - (b2 - x2) * D;
  
  // Zero out the rest of the elements
  for (i = 0; i < 5; i++) c[i].i = 0.;
  for (i = 0; i < 5; i++) {
    croots[i].r = NAN;
    croots[i].i = NAN;
  }
  
  // Solve the quartic w/ polishing
  zroots(c, 4, croots, 1);
  
  // Get the real roots (up to 2)
  roots[0] = NAN;
  roots[1] = NAN;
  j = 0;
  for (i = 1; i < 5; i++) {
    if (!(isnan(croots[i].r)) && (fabs(croots[i].i) < DTOL1)) {
      roots[j] = croots[i].r;
      j += 1;
    }
    if (j == 2) break;
  }

}

double Latitude(double x, double y, double x0, double y0, double r, double theta) {
  /*
  
  */
  
  int i;
  double alpha, beta, gamma, c1, c2, c3, a, b, c;
  double xterm, xlimb, xE, yE, dx, xmin, xmax;
  double z[2];
  
  // Compute the latitude. We will solve
  // a quadratic equation for z = sin(lat) **  2  
  alpha = (x - x0) / (r * sin(theta));
  beta = 1 / tan(theta);
  gamma = (y - y0) / r;
  c1 = 4 * alpha * alpha * beta * beta;
  c2 = 1 + beta * beta;
  c3 = alpha * alpha + gamma * gamma + beta * beta;
  b = (c1 - 2 * c2 * c3) / (c2 * c2);
  c = (c3 * c3 - c1) / (c2 * c2);
  z[0] = (-b + sqrt(b * b - 4 * c)) / 2;
  z[1] = (-b - sqrt(b * b - 4 * c)) / 2;
  
  // Where is the terminator for this value of `y`?
  if (theta <= 0)
    xterm = x0 - fabs(sin(theta)) * sqrt(r * r - (y - y0) * (y - y0));
  else
    xterm = x0 + fabs(sin(theta)) * sqrt(r * r - (y - y0) * (y - y0));

  // Are we on the dayside?
  if (x <= xterm) {
  
    // We have two possible solutions. But only one is on the
    // observer's side of the occulted planet. TODO: Speed this up?
    for (i = 0; i < 2; i++) {
      if ((z[i] >= 0) && (z[i] <= 1)) {
        a = r * sqrt(z[i]);
        b = a * fabs(sin(theta));
        xE = x0 - r * sqrt(1 - z[i]) * cos(theta);
        yE = y0;
        dx = (b / a) * sqrt(a * a - (y - yE) * (y - yE));
        xlimb = r * sqrt(1 - z[i]) * sin(theta) * tan(theta);
        if (((theta > 0) && (b < xlimb)) || ((theta <= 0) && (b > xlimb)))
          xmin = xE - b;
        else
          xmin = xE - xlimb;
        if (((theta > 0) && (b > -xlimb)) || ((theta <= 0) && (b < -xlimb)))
          xmax = xE + b;
        else
          xmax = xE - xlimb;
        if ((x >= xmin - DTOL1) && (x <= xmax + DTOL1)) {
          if (((fabs(x - (xE + dx)) < DTOL1) || (fabs(x - (xE - dx)) < DTOL1)))
            return asin(sqrt(z[i]));
        }
      }
    }
    
  // Or the nightside?
  } else {
    
    // We have two possible solutions. But only one is on the
    // observer's side of the occulted. TODO: Speed this up?
    for (i = 0; i < 2; i++) {
      if ((z[i] >= 0) && (z[i] <= 1)) {
        a = r * sqrt(z[i]);
        b = a * fabs(sin(theta));
        xE = x0 + r * sqrt(1 - z[i]) * cos(theta);
        yE = y0;
        dx = (b / a) * sqrt(a * a - (y - yE) * (y - yE));
        xlimb = -r * sqrt(1 - z[i]) * sin(theta) * tan(theta);
        if (((theta > 0) && (b < xlimb)) || ((theta <= 0) && (b > xlimb)))
          xmin = xE - b;
        else
          xmin = xE - xlimb;
        if (((theta > 0) && (b > -xlimb)) || ((theta <= 0) && (b < -xlimb)))
          xmax = xE + b;
        else
          xmax = xE - xlimb;
        if ((x >= xmin - DTOL1) && (x <= xmax + DTOL1)) {
          if (((fabs(x - (xE + dx)) < DTOL1) || (fabs(x - (xE - dx)) < DTOL1)))
            return PI - asin(sqrt(z[i]));
        }
      }
    }
  
  }
  
  return NAN;
}

double SurfaceBrightness(double lat, double noon, double midnight, int n) {
  /*
  
  */
  
  int i;
  double l;
  
  for (i = 0; i < n + 2; i++) {
    if (lat < PI / 2.) {
      l = i * PI / (n + 1);
      if (l > lat) {
        l -= PI / (n + 1);
        break;
      }
    }
    else {
      l = i * PI / (n + 1);
      if (l > lat) {
        break;
      }
    }
  }
  
  return 0.5 * (noon - midnight) * (cos(l) + 1) + midnight;
  
}

double curve(double x, FUNCTION function) {
  return function.curve(x, function.ellipse);
}

double integral(double x0, double x1, FUNCTION function) {
  return function.integral(x0, x1, function.ellipse);
}

double fupper(double x, ELLIPSE *ellipse) {
  double A;
  if (ellipse->circle) {
    A = ellipse->r * ellipse->r - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < DTOL2) A = 0;
    if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 + sqrt(A);
  } else {
    A = ellipse->b * ellipse->b - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < DTOL2) A = 0;
    if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 + (ellipse->a / ellipse->b) * sqrt(A);
  }
}

double flower(double x, ELLIPSE *ellipse) {
  double A;
  if (ellipse->circle) {
    A = ellipse->r * ellipse->r - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < DTOL2) A = 0;
    if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 - sqrt(A);
  } else {
    A = ellipse->b * ellipse->b - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < DTOL2) A = 0;
    if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 - (ellipse->a / ellipse->b) * sqrt(A);
  }
}

double iupper(double xL, double xR, ELLIPSE *ellipse) {
  double yL,yR,zL,zR,FL,FR;
  if ((xL > ellipse->xmax) || (xL < ellipse->xmin)) return NAN;
  if ((xR > ellipse->xmax) || (xR < ellipse->xmin)) return NAN;
  yL = xL - ellipse->x0;
  yR = xR - ellipse->x0;
  if (ellipse->circle) {
    zL = sqrt((ellipse->r - yL) * (ellipse->r + yL));
    zR = sqrt((ellipse->r - yR) * (ellipse->r + yR));
    FL = 0.5 * (yL * zL + ellipse->r * ellipse->r * atan(yL / zL)) + ellipse->y0 * xL;
    FR = 0.5 * (yR * zR + ellipse->r * ellipse->r * atan(yR / zR)) + ellipse->y0 * xR;
  } else {
    zL = sqrt((ellipse->b - yL) * (ellipse->b + yL));
    zR = sqrt((ellipse->b - yR) * (ellipse->b + yR));
    FL = (ellipse->a / (2 * ellipse->b)) * (yL * zL + ellipse->b * ellipse->b * atan(yL / zL)) + ellipse->y0 * xL;
    FR = (ellipse->a / (2 * ellipse->b)) * (yR * zR + ellipse->b * ellipse->b * atan(yR / zR)) + ellipse->y0 * xR;
  }
  return FR - FL;
}

double ilower(double xL, double xR, ELLIPSE *ellipse) {
  double yL,yR,zL,zR,FL,FR;
  if ((xL > ellipse->xmax) || (xL < ellipse->xmin)) return NAN;
  if ((xR > ellipse->xmax) || (xR < ellipse->xmin)) return NAN;
  yL = xL - ellipse->x0;
  yR = xR - ellipse->x0;
  if (ellipse->circle) {
    zL = sqrt((ellipse->r - yL) * (ellipse->r + yL));
    zR = sqrt((ellipse->r - yR) * (ellipse->r + yR));
    FL = -0.5 * (yL * zL + ellipse->r * ellipse->r * atan(yL / zL)) + ellipse->y0 * xL;
    FR = -0.5 * (yR * zR + ellipse->r * ellipse->r * atan(yR / zR)) + ellipse->y0 * xR;
  } else {
    zL = sqrt((ellipse->b - yL) * (ellipse->b + yL));
    zR = sqrt((ellipse->b - yR) * (ellipse->b + yR));
    FL = -(ellipse->a / (2 * ellipse->b)) * (yL * zL + ellipse->b * ellipse->b * atan(yL / zL)) + ellipse->y0 * xL;
    FR = -(ellipse->a / (2 * ellipse->b)) * (yR * zR + ellipse->b * ellipse->b * atan(yR / zR)) + ellipse->y0 * xR;
  }
  return FR - FL;
}

int dblcomp( const void* a, const void* b) {
  /*
  
  */
  
  double dbl_a = * ( (double*) a );
  double dbl_b = * ( (double*) b );
  if ( dbl_a == dbl_b ) return 0;
  else if ( dbl_a < dbl_b ) return -1;
  else return 1;
}

int funcomp( const void* a, const void* b) {
  /*
  
  */
  
  FUNCTION *fun_a = (FUNCTION *)a;
  FUNCTION *fun_b = (FUNCTION *)b;
  if ( fun_a->y == fun_b->y ) return 0;
  else if ( fun_a->y < fun_b->y ) return -1;
  else return 1;
}

void AddLatitudeSlice(double latitude, double r, double x0, double y0, double ro, double theta, double *vertices, int *v, FUNCTION *functions, int *f){
  /*
  
  */    
  
  double xlimb, x, y;
  double roots[2];
  double r2 = ro * ro + DTOL2;
  ELLIPSE *ellipse;
  ellipse = malloc(sizeof(ELLIPSE)); 
  
  // Construct the ellipse
  ellipse->circle = 0;
  ellipse->a = r * fabs(sin(latitude));
  ellipse->b = ellipse->a * fabs(sin(theta));
  ellipse->x0 = x0 - r * cos(latitude) * cos(theta);
  ellipse->y0 = y0;
  
  // The x position of the intersection with the occulted planet
  // limb, relative to the occulted planet center
  xlimb = r * cos(latitude) * sin(theta) * tan(theta);

  // Ellipse x minimum
  if (((theta > 0) && (ellipse->b < xlimb)) || ((theta <= 0) && (ellipse->b > xlimb))) {
    ellipse->xmin = ellipse->x0 - ellipse->b;
    if (ellipse->xmin * ellipse->xmin + ellipse->y0 * ellipse->y0 < r2)
      vertices[(*v)++] = ellipse->xmin;    
  } else {
    ellipse->xmin = ellipse->x0 - xlimb;
  }
  
  // Ellipse x maximum
  if (((theta > 0) && (ellipse->b > -xlimb)) || ((theta <= 0) && (ellipse->b < -xlimb))) {
    ellipse->xmax = ellipse->x0 + ellipse->b;
    if (ellipse->xmax * ellipse->xmax + ellipse->y0 * ellipse->y0 < r2)
      vertices[(*v)++] = ellipse->xmax;    
  } else {
    ellipse->xmax = ellipse->x0 - xlimb;
  }
  
  // Are the limb vertices inside the occultor? 
  if (ellipse->b >= fabs(xlimb)) {
    x = ellipse->x0 - xlimb;
    y = fupper(x, ellipse);
    if (x * x + y * y < r2)
      vertices[(*v)++] = x;
    y = flower(x, ellipse);
    if (x * x + y * y < r2)
      vertices[(*v)++] = x;
  }
  
  // Ellipse-occultor vertices
  GetRoots(ellipse->a, ellipse->b, ellipse->x0, ellipse->y0, ro, roots);

  if (!(isnan(roots[0])) && (((theta > 0) && (roots[0] > ellipse->x0 - xlimb)) || ((theta <= 0) && (roots[0] < ellipse->x0 - xlimb)))) {
    vertices[(*v)++] = roots[0];
  }
  if (!(isnan(roots[1])) && (((theta > 0) && (roots[1] > ellipse->x0 - xlimb)) || ((theta <= 0) && (roots[1] < ellipse->x0 - xlimb)))) {
    vertices[(*v)++] = roots[1];
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = ellipse;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = ellipse;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  
}

void AddOccultor(double r, double x0, double y0, double ro, double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  double r2 = r * r + DTOL2;
  ELLIPSE *occultor;
  occultor = malloc(sizeof(ELLIPSE)); 
  
  // Initialize the occultor (always at the origin)
  occultor->r = ro;
  occultor->x0 = 0.;
  occultor->y0 = 0.;
  occultor->circle = 1;
    
  // Occultor x minimum
  occultor->xmin = -occultor->r;
  if ((occultor->xmin - x0) * (occultor->xmin - x0) + (y0 * y0) < r2)
    vertices[(*v)++] = occultor->xmin;    
  
  // Occultor x maximum
  occultor->xmax = occultor->r;
  if ((occultor->xmax - x0) * (occultor->xmax - x0) + (y0 * y0) < r2)
    vertices[(*v)++] = occultor->xmax;    

  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = occultor;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = occultor;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  
}

void AddOcculted(double r, double x0, double y0, double ro, double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  double r2 = ro * ro + DTOL2;
  double d, A, x, y, cost, sint;
  ELLIPSE *occulted;
  occulted = malloc(sizeof(ELLIPSE));
  
  // Initialize the occulted planet
  occulted->r = r;
  occulted->x0 = x0;
  occulted->y0 = y0;
  occulted->circle = 1;
  
  // Occulted planet x minimum
  occulted->xmin = x0 - r;
  if (occulted->xmin * occulted->xmin + occulted->y0 * occulted->y0 < r2)
    vertices[(*v)++] = occulted->xmin;    
  
  // Occulted planet x maximum
  occulted->xmax = x0 + r;
  if (occulted->xmax * occulted->xmax + occulted->y0 * occulted->y0 < r2)
    vertices[(*v)++] = occulted->xmax;    
  
  // Vertices of intersection with the occultor
  // Adapted from http://mathworld.wolfram.com/Circle-CircleIntersection.html
  d = sqrt(x0 * x0 + y0 * y0);
  if (d < (ro + r)) {
    A = (-d + r - ro) * (-d - r + ro) * (-d + r + ro) * (d + r + ro);
    if (A >= 0) {
      y = sqrt(A) / (2 * d);
      x = (d * d - r * r + ro * ro) / (2 * d);
      cost = -x0 / d;
      sint = -y0 / d;
      vertices[(*v)++] = -x * cost + y * sint;
      vertices[(*v)++] = -x * cost - y * sint;
    }
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = occulted;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = occulted;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  
}

double OccultedFlux(double r, double x0, double y0, double ro, double theta, double noon, double midnight, int n) {
  /*
  
  */
  
  int i, j;
  int b = 0;
  int v = 0;
  int f = 0;
  double s;
  double lat;
  double xL, xR, x, y, area;
  double flux = 0.;
  double r2 = r * r + DTOL2;
  double ro2 = ro * ro + DTOL2;
  double vertices[MAXVERTICES];
  FUNCTION functions[MAXFUNCTIONS];
  FUNCTION boundaries[MAXFUNCTIONS];
  
  // Avoid the singular point
  if (fabs(theta) < 1e-2)
    theta = 1e-2;
  
  // Generate all the shapes and get their vertices and curves
  AddOccultor(r, x0, y0, ro, vertices, &v, functions, &f);   
  AddOcculted(r, x0, y0, ro, vertices, &v, functions, &f); 
  for (i = 0; i < n; i++) {
    lat = (i + 1) * PI / (n + 1);
    AddLatitudeSlice(lat, r, x0, y0, ro, theta, vertices, &v, functions, &f);
  }
    
  // Sort the vertices
  qsort(vertices, v, sizeof(double), dblcomp);
    
  // Loop over all vertices
  for (i = 0; i < v - 1; i++) {
    
    // Get the integration limits, with a
    // small perturbation for stability
    xL = vertices[i] + DTOL2;
    xR = vertices[i+1] - DTOL2;

    // Check if they are identical
    if (xR <= xL + DTOL2)
      continue;
        
    // Bisect the limits. Find the boundary functions that are
    // finite valued at this point and sort their integrals 
    // in order of increasing function value.
    x = 0.5 * (xL + xR);
    
    // Loop over all functions
    b = 0;
    for (j = 0; j < f; j++) {
      y = curve(x, functions[j]);
      if (!(isnan(y)) && (x * x + y * y < ro2) && ((x - x0) * (x - x0) + (y - y0) * (y - y0) < r2)) {
        boundaries[b] = functions[j];
        boundaries[b++].y = y;
      }
    }

    // Sort boundaries based on y
    qsort(boundaries, b, sizeof(FUNCTION), funcomp);

    // Loop over all regions bounded by xL and xR
    for (j = 0; j < b - 1; j++) {
      // The area of each region is just the difference of successive integrals
      area = integral(xL, xR, boundaries[j + 1]) - integral(xL, xR, boundaries[j]);

      // Get the latitude of the midpoint and the corresponding surface brightness
      y = 0.5 * (boundaries[j + 1].y + boundaries[j].y);
      lat = Latitude(x, y, x0, y0, r, theta);
      s = SurfaceBrightness(lat, noon, midnight, n);
      flux += s * area;
      
    } 
    
  }
  
  // Free all of the ellipse instances
  for (j = 0; j < f; j+=2) {
    free(functions[j].ellipse);
  }  

  return flux;
   
}

double UnoccultedFlux(double r, double theta, double noon, double midnight, int n) {
  /*
  
  */
  
  return OccultedFlux(r, 0, 0, 2 * r, theta, noon, midnight, n);
}