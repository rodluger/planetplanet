#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include "complex.h"

double EyeballTemperature(double lat, double irrad, double albedo) {
  /*
  
  */
  
  if (lat < PI / 2)
    return pow((irrad * cos(lat) * (1 - albedo)) / SBOLTZ, 0.25);
  else
    return 0;
}

double Blackbody(double lambda, double T) {
  /*
  
  */
  
  // Planck's law
  double a = 2 * HPLANCK * CLIGHT * CLIGHT / (lambda * lambda * lambda * lambda * lambda);
  double b = HPLANCK * CLIGHT / (lambda * KBOLTZ * T);
  return a / (exp(b) - 1);
}

void GetRoots(double a, double b, double x0, double y0, double r, double polyeps1, double polyeps2, int maxpolyiter, double roots[2]) {
  /*
  
  */
  
  int i, j;
  double A, B, C, D;
  dcomplex c[5];
  dcomplex croots[5];
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
  zroots(c, 4, croots, 1, polyeps1, polyeps2, maxpolyiter);
  
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
  double xE, yE;
  double z[2];
  double l[2];
  double d;
  double la, lb, da, db;
  double tmp, tmp1, tmp2;
  double sintheta, costheta, cotantheta;
  double solution;

  // Are we dealing with circles?
  if ((fabs(fabs(theta) - PI / 2)) < 1.e-5) {
    
    // Trivial!
    d = sqrt((x - x0) * (x - x0) + (y - y0) * (y - y0));
    solution = asin(d / r);
    if (theta < 0)
      solution = PI - solution;
    return solution;
    
  }
  
  // A little trig
  sintheta = sin(theta);
  costheta = cos(theta);
  cotantheta = costheta / sintheta;
  
  // Compute the latitude. We will solve
  // a quadratic equation for z = sin(lat) **  2  
  alpha = (x - x0) / (r * sintheta);
  beta = cotantheta;
  gamma = (y - y0) / r;
  c1 = 4 * alpha * alpha * beta * beta;
  c2 = 1 + beta * beta;
  c3 = alpha * alpha + gamma * gamma + beta * beta;
  b = (c1 - 2 * c2 * c3) / (c2 * c2);
  c = (c3 * c3 - c1) / (c2 * c2);
  tmp = b * b - 4 * c;
  if (tmp < 0) tmp = 0;
  z[0] = (-b + sqrt(tmp)) / 2;
  z[1] = (-b - sqrt(tmp)) / 2;
  
  // Find the two solutions that thread the point
  for (i = 0; i < 2; i++) {
    
    // A: Compute the distance to the desired (x,y) point
    la = asin(sqrt(z[i]));
    a = r * fabs(sin(la));
    b = a * fabs(sintheta);
    xE = x0 - r * cos(la) * costheta;
    yE = y0;
    tmp = (a / b) * sqrt(fabs(b * b - (x - xE) * (x - xE)));
    tmp1 = fabs(y - (yE - tmp));
    tmp2 = fabs(y - (yE + tmp));
    if (tmp1 < tmp2) da = tmp1;
    else da = tmp2;
    
    // B: Compute the distance to the desired (x,y) point
    lb = PI - la;
    a = r * fabs(sin(lb));
    b = a * fabs(sintheta);
    xE = x0 - r * cos(lb) * costheta;
    yE = y0;
    tmp = (a / b) * sqrt(fabs(b * b - (x - xE) * (x - xE)));
    tmp1 = fabs(y - (yE - tmp));
    tmp2 = fabs(y - (yE + tmp));
    if (tmp1 < tmp2) db = tmp1;
    else db = tmp2;
    
    // Get the closest solution
    if (da < db)
      l[i] = la;
    else
      l[i] = lb;
    
  }
    
  // Only one of these solutions is on the observer's side
  if (theta < 0) {
    if (l[0] > l[1])
      solution = l[0];
    else
      solution = l[1];
  } else {
    if (l[0] < l[1])
      solution = l[0];
    else
      solution = l[1];
  }
  
  return solution;
  
}

void SurfaceIntensity(double albedo, double irrad, int nu, double u[nu], double lmin, double lmax, int nlat, int nlam, double lambda[nlam], double B[nlat + 1][nlam]) {
  /*
  Returns the blackbody intensity at the center of each latitude slice 
  evaluated at a given array of wavelengths.
  
  */
  
  int i, j;
  double latitude, coslat;
  double T = 0.;
  
  // Loop over each slice
  for (i = 0; i < nlat + 1; i++) {
    
    // Get the latitude at the *center* of this slice
    latitude = (i + 0.5) / (nlat + 1.) * (lmax - lmin) + lmin;
    
    // Get its blackbody temperature
    if (nu == 0) {
      
      // No limb darkening
      T = EyeballTemperature(latitude, irrad, albedo);
    
    } else {
    
      // Polynomial limb darkening
      T = 0;
      coslat = cos(latitude);
      for (j = 0; j < nu; j++) {
        T += u[j] * pow(1 - coslat, j);
      }
      
    }

    // Loop over each wavelength bin
    for (j = 0; j < nlam; j++) {
        
      // Get the blackbody intensity (W / m^2 / m / sr)
      B[i][j] = Blackbody(lambda[j], T);
      
    }
  } 
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

void AddLatitudeSlice(double latitude, double r, double x0, double y0, double ro, double theta, double polyeps1, double polyeps2, int maxpolyiter, double *vertices, int *v, FUNCTION *functions, int *f){
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
  // limb, relative to the ellipse center
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
  GetRoots(ellipse->a, ellipse->b, ellipse->x0, ellipse->y0, ro, polyeps1, polyeps2, maxpolyiter, roots);

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

void OccultedFlux(double r, double x0, double y0, double ro, double theta, double albedo, double irrad, double polyeps1, double polyeps2, int maxpolyiter, int nu, int nlat, int nlam, double u[nu], double lambda[nlam], double flux[nlam]) {
  /*
  
  */
  
  int i, j, k, m;
  int b = 0;
  int v = 0;
  int f = 0;
  double d, lmin, lmax, lat;
  double xL, xR, x, y, area;
  double r2 = r * r + DTOL2;
  double ro2 = ro * ro + DTOL2;
  double vertices[MAXVERTICES];
  FUNCTION functions[MAXFUNCTIONS];
  FUNCTION boundaries[MAXFUNCTIONS];
  double B[nlat + 1][nlam];
  
  // Avoid the singular point
  if (fabs(theta) < 1e-2)
    theta = 1e-2;
  
  // If we're doing limb darkening, theta must be pi/2 (full phase)
  // If the occultor is completely within the occulted body (usually
  // the case during transit), only compute latitude circles that
  // intersect the occultor
  if (nu > 0) {
    theta = PI / 2.;
    d = sqrt(x0 * x0 + y0 * y0);
    if (d < r - ro) {
      lmin = asin((d - ro) / r) - 1e-5;
      lmax = asin((d + ro) / r) + 1e-5;
    } else {
      lmin = 0 - 1e-5;
      lmax = PI + 1e-5;
    }
  } else {
    lmin = 0 - 1e-5;
    lmax = PI + 1e-5;
  }
  
  // Zero out the flux
  for (m = 0; m < nlam; m++)
    flux[m] = 0.;

  // Pre-compute the surface intensity in each latitude slice
  SurfaceIntensity(albedo, irrad, nu, u, lmin, lmax, nlat, nlam, lambda, B);
    
  // Generate all the shapes and get their vertices and curves
  AddOccultor(r, x0, y0, ro, vertices, &v, functions, &f);   
  AddOcculted(r, x0, y0, ro, vertices, &v, functions, &f); 
  for (i = 0; i < nlat; i++) {
    lat = (i + 1.) / (nlat + 1.) * (lmax - lmin) + lmin;
    AddLatitudeSlice(lat, r, x0, y0, ro, theta, polyeps1, polyeps2, maxpolyiter, vertices, &v, functions, &f);
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
      
      // Get the latitude of the midpoint and the index of the latitude
      // grid *above* this latitude
      y = 0.5 * (boundaries[j + 1].y + boundaries[j].y);
      lat = Latitude(x, y, x0, y0, r, theta);      
      k = (int) ((lat - lmin) / ((lmax - lmin) / (nlat + 1.)));
      
      // Sanity check
      if ((k < 0) || (k > nlat)) {
        printf("Invalid index for latitude grid! Aborting.\n");
        abort();
      }
      
      // Multiply by the area of the latitude slice to get the flux at each wavelength
      for (m = 0; m < nlam; m++)
        flux[m] += B[k][m] * area;

    } 
    
  }
  
  // Scale flux by REARTH**2
  // Convert to W m^2 / m^2 / um / sr
  for (m = 0; m < nlam; m++)
    flux[m] *= REARTH * REARTH * 1e-6;  
  
  // Free all of the ellipse instances
  for (j = 0; j < f; j+=2) {
    free(functions[j].ellipse);
  }  
   
}

void UnoccultedFlux(double r, double theta, double albedo, double irrad, double polyeps1, double polyeps2, int maxpolyiter, int nu, int nlat, int nlam, double u[nu], double lambda[nlam], double flux[nlam]) {
  /*
  
  */
  
  OccultedFlux(r, 0, 0, 2 * r, theta, albedo, irrad, polyeps1, polyeps2, maxpolyiter, nu, nlat, nlam, u, lambda, flux);
  
}