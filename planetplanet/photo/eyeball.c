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

void GetRoots(double a, double b, double xE, double yE, double xC, double yC, double r, double polyeps1, double polyeps2, int maxpolyiter, double roots[2]) {
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
  double x0 = xE - xC;
  double y0 = yE - yC;
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
      roots[j] = croots[i].r + xC;
      j += 1;
    }
    if (j == 2) break;
  }

}

double Latitude(double x, double y, double r, double theta) {
  /*
  
  */
    
  int i;
  double alpha, beta, gamma, c1, c2, c3, a, b, c;
  double xE;
  double z[2];
  double l[2];
  double d;
  double la, lb, da, db;
  double tmp, tmp1, tmp2;
  double sintheta, costheta, cotantheta;
  double solution;

  // Are we dealing with circles?
  if ((fabs(fabs(theta) - PI / 2)) < DTOL1) {
    
    // Trivial!
    d = sqrt(x * x + y * y);
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
  alpha = x / (r * sintheta);
  beta = cotantheta;
  gamma = y / r;
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
    xE = -r * cos(la) * costheta;
    tmp = (a / b) * sqrt(fabs(b * b - (x - xE) * (x - xE)));
    tmp1 = fabs(y + tmp);
    tmp2 = fabs(y - tmp);
    if (tmp1 < tmp2) da = tmp1;
    else da = tmp2;
    
    // B: Compute the distance to the desired (x,y) point
    lb = PI - la;
    a = r * fabs(sin(lb));
    b = a * fabs(sintheta);
    xE = -r * cos(lb) * costheta;
    tmp = (a / b) * sqrt(fabs(b * b - (x - xE) * (x - xE)));
    tmp1 = fabs(y + tmp);
    tmp2 = fabs(y - tmp);
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

void AddLatitudeSlice(double latitude, double r, int no, double x0[no], double y0[no], double ro[no], double theta, double polyeps1, double polyeps2, int maxpolyiter, double *vertices, int *v, FUNCTION *functions, int *f){
  /*
  
  */    
  
  int i;
  double xlimb, x, y;
  double roots[2];
  double ro2[no];
  for (i = 0; i < no; i++) ro2[i] = ro[i] * ro[i] + DTOL1;
  ELLIPSE *ellipse;
  ellipse = malloc(sizeof(ELLIPSE)); 
  
  // Construct the ellipse
  ellipse->circle = 0;
  ellipse->a = r * fabs(sin(latitude));
  ellipse->b = ellipse->a * fabs(sin(theta));
  ellipse->x0 = -r * cos(latitude) * cos(theta);
  ellipse->y0 = 0;
  
  // The x position of the intersection with the occulted planet
  // limb, relative to the ellipse center
  xlimb = r * cos(latitude) * sin(theta) * tan(theta);

  // Ellipse x minimum
  if (((theta > 0) && (ellipse->b < xlimb)) || ((theta <= 0) && (ellipse->b > xlimb))) {
    ellipse->xmin = ellipse->x0 - ellipse->b;
    // Is this point inside __at least__ one occultor?
    for (i = 0; i < no; i++) {
      if ((ellipse->xmin - x0[i]) * (ellipse->xmin - x0[i]) + y0[i] * y0[i] < ro2[i]) {
        vertices[(*v)++] = ellipse->xmin;
        break;    
      }
    }
  } else {
    ellipse->xmin = ellipse->x0 - xlimb;
  }
  
  // Ellipse x maximum
  if (((theta > 0) && (ellipse->b > -xlimb)) || ((theta <= 0) && (ellipse->b < -xlimb))) {
    ellipse->xmax = ellipse->x0 + ellipse->b;
    // Is this point inside __at least__ one occultor?
    for (i = 0; i < no; i++) {
      if ((ellipse->xmax - x0[i]) * (ellipse->xmax - x0[i]) + y0[i] * y0[i] < ro2[i]) {
        vertices[(*v)++] = ellipse->xmax;
        break;    
      }
    }
  } else {
    ellipse->xmax = ellipse->x0 - xlimb;
  }
  
  // Are the limb vertices inside __at least__ one occultor? 
  if (ellipse->b >= fabs(xlimb)) {
    x = ellipse->x0 - xlimb;
    y = fupper(x, ellipse);
    for (i = 0; i < no; i++) {
      if ((x - x0[i]) * (x - x0[i]) + (y - y0[i]) * (y - y0[i]) < ro2[i]) {
        vertices[(*v)++] = x;
        break;    
      }
    }
    y = flower(x, ellipse);
    for (i = 0; i < no; i++) {
      if ((x - x0[i]) * (x - x0[i]) + (y - y0[i]) * (y - y0[i]) < ro2[i]) {
        vertices[(*v)++] = x;
        break;    
      }
    }
  }
  
  // Ellipse-occultor vertices
  for (i = 0; i < no; i++) {
    
    // Solve the quartic
    GetRoots(ellipse->a, ellipse->b, ellipse->x0, ellipse->y0, x0[i], y0[i], ro[i], polyeps1, polyeps2, maxpolyiter, roots);

    if (!(isnan(roots[0])) && (((theta > 0) && (roots[0] > ellipse->x0 - xlimb)) || ((theta <= 0) && (roots[0] < ellipse->x0 - xlimb)))) {
      vertices[(*v)++] = roots[0];
    }
    if (!(isnan(roots[1])) && (((theta > 0) && (roots[1] > ellipse->x0 - xlimb)) || ((theta <= 0) && (roots[1] < ellipse->x0 - xlimb)))) {
      vertices[(*v)++] = roots[1];
    }
    
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = ellipse;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = ellipse;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  
}

void AddOccultors(double r, int no, double x0[no], double y0[no], double ro[no], double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  int i;
  double r2 = r * r + DTOL2;
  ELLIPSE *occultor[no];
   
  for (i = 0; i < no; i++) {
    
    // Allocate memory
    occultor[i] = malloc(sizeof(ELLIPSE));
    
    // Initialize the occultor
    occultor[i]->r = ro[i];
    occultor[i]->x0 = x0[i];
    occultor[i]->y0 = y0[i];
    occultor[i]->circle = 1;
    
    // Occultor x minimum
    occultor[i]->xmin = x0[i] - occultor[i]->r;
    if (occultor[i]->xmin * occultor[i]->xmin + occultor[i]->y0 * occultor[i]->y0 < r2)
      vertices[(*v)++] = occultor[i]->xmin;    
  
    // Occultor x maximum
    occultor[i]->xmax = x0[i] + occultor[i]->r;
    if (occultor[i]->xmax * occultor[i]->xmax + occultor[i]->y0 * occultor[i]->y0 < r2)
      vertices[(*v)++] = occultor[i]->xmax;    

    // Finally, the boundary curves and their integrals
    functions[*f].ellipse = occultor[i];
    functions[*f].curve = fupper;
    functions[(*f)++].integral = iupper;
    functions[*f].ellipse = occultor[i];
    functions[*f].curve = flower;
    functions[(*f)++].integral = ilower;
  
  }
  
}

void AddOcculted(double r, int no, double x0[no], double y0[no], double ro[no], double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  int i;
  double ro2[no];
  for (i = 0; i < no; i++) ro2[i] = ro[i] * ro[i] + DTOL1;
  double d, A, x, y, cost, sint;
  ELLIPSE *occulted;
  occulted = malloc(sizeof(ELLIPSE));
  
  // Initialize the occulted planet at the origin
  occulted->r = r;
  occulted->x0 = 0.;
  occulted->y0 = 0.;
  occulted->circle = 1;
  
  // Occulted planet x minimum
  occulted->xmin = -r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((occulted->xmin - x0[i]) * (occulted->xmin - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = occulted->xmin;
      break;    
    }
  }  
  
  // Occulted planet x maximum
  occulted->xmax = r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((occulted->xmax - x0[i]) * (occulted->xmax - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = occulted->xmax;
      break;    
    }
  }  

  // Vertices of intersection with the occultors
  // Adapted from http://mathworld.wolfram.com/Circle-CircleIntersection.html
  for (i = 0; i < no; i++) {
    d = sqrt(x0[i] * x0[i] + y0[i] * y0[i]);
    if (d < (ro[i] + r)) {
      A = (-d + r - ro[i]) * (-d - r + ro[i]) * (-d + r + ro[i]) * (d + r + ro[i]);
      if (A >= 0) {
        y = sqrt(A) / (2 * d);
        x = (d * d - r * r + ro[i] * ro[i]) / (2 * d);
        cost = -x0[i] / d;
        sint = -y0[i] / d;
        vertices[(*v)++] = -x * cost + y * sint + x0[i];
        vertices[(*v)++] = -x * cost - y * sint + x0[i];
      }
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

void OccultedFlux(double r, int no, double x0[no], double y0[no], double ro[no], double theta, double albedo, double irrad, double polyeps1, double polyeps2, int maxpolyiter, int adaptive, int nu, int nlat, int nlam, double u[nu], double lambda[nlam], double flux[nlam]) {
  /*
  
  */
  
  int i, j, k, m;
  int b = 0;
  int v = 0;
  int f = 0;
  double d, lmin, lmax, lat;
  double xL, xR, x, y, area;
  double r2 = r * r + DTOL1;
  double tmp, dmin, dmax;
  double ro2[no];
  for (i = 0; i < no; i++) ro2[i] = ro[i] * ro[i] + DTOL1;
  double vertices[MAXVERTICES];
  FUNCTION functions[MAXFUNCTIONS];
  FUNCTION boundaries[MAXFUNCTIONS];
  double B[nlat + 1][nlam];
  
  // Avoid the singular point
  if (fabs(theta) < MINTHETA)
    theta = MINTHETA;
  
  // Zero out the flux
  for (m = 0; m < nlam; m++) 
    flux[m] = 0.;
  
  // If we're doing limb darkening, theta must be pi/2 (full phase)
  if (nu > 0)
    theta = PI / 2.;
  
  // Latitude grid bounds
  lmin = -DTOL1;
  lmax = PI + DTOL1;
  
  // TODO! Adaptive grid is broken!!!
  
  // Adaptive grid? Let's only compute latitude circles that intersect the occultors
  if ((adaptive) && (nu > 0)) {
  
    dmin = 1.;
    dmax = 0.;
    
    // Loop over each occultor
    for (i = 0; i < no; i++) {
    
      // Distance to occultor center
      d = sqrt(x0[i] * x0[i] + y0[i] * y0[i]);
      
      // Distance to far side
      tmp = (d + ro[i]) / r;
      if (tmp > dmax) dmax = tmp;
      
      // Distance to near side
      tmp = (d - ro[i]) / r;
      if (tmp < dmin) dmin = tmp;
    
    }
    
    // Compute min and max latitudes
    if (dmax >= 1 - DTOL1)
      lmax = PI + DTOL1;
    else
      lmax = asin(dmax) + DTOL1;
    if (dmin <= DTOL1)
      lmin = -DTOL1;
    else
      lmin = asin(dmin) - DTOL1;  

  }
  
  // Pre-compute the surface intensity in each latitude slice
  SurfaceIntensity(albedo, irrad, nu, u, lmin, lmax, nlat, nlam, lambda, B);
    
  // Generate all the shapes and get their vertices and curves
  AddOccultors(r, no, x0, y0, ro, vertices, &v, functions, &f);   
  AddOcculted(r, no, x0, y0, ro, vertices, &v, functions, &f); 
  for (i = 0; i < nlat; i++) {
    lat = (i + 1.) / (nlat + 1.) * (lmax - lmin) + lmin;
    AddLatitudeSlice(lat, r, no, x0, y0, ro, theta, polyeps1, polyeps2, maxpolyiter, vertices, &v, functions, &f);
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
      
      // Check that it's inside the planet
      if (!(isnan(y)) && (x * x + y * y < r2)) {
      
        // Check that it's inside __at least__ one occultor
        for (k = 0; k < no; k++) {
          if ((x - x0[k]) * (x - x0[k]) + (y - y0[k]) * (y - y0[k]) < ro2[k]) {
            
            // Add the function to the list of boundaries
            boundaries[b] = functions[j];
            boundaries[b++].y = y;
            break;
          
          }
        
        }
      
      }
      
    }

    // Sort boundaries based on y
    qsort(boundaries, b, sizeof(FUNCTION), funcomp);
        
    // Loop over all regions bounded by xL and xR
    for (j = 0; j < b - 1; j++) {
    
      // The area of each region is just the difference of successive integrals
      area = integral(xL, xR, boundaries[j + 1]) - integral(xL, xR, boundaries[j]);
      
      // Get the latitude of the midpoint
      y = 0.5 * (boundaries[j + 1].y + boundaries[j].y);
      if (x * x + y * y > r * r) {
        // We are off-planet!
        continue;
      }
      lat = Latitude(x, y, r, theta);
      
      // Get the index of the latitude grid *above* this latitude
      k = (int) ((lat - lmin) / ((lmax - lmin) / (nlat + 1.)));
            
      // Sanity check
      if ((k < 0) || (k > nlat)) {
        printf("Invalid index for latitude grid.\n");
        printf("k = %d; nlat = %d\n", k, nlat);
        printf("Latitude = %.3f; min, max = %.3f, %.3f\n", lat * 180/PI, lmin * 180/PI, lmax * 180/PI);
        printf("x, y, r, theta = %.3f, %.3f, %.3f, %.3f\n", x, y, r, theta);
        printf("Aborting.\n");        
        abort();
      }
      
      // Multiply by the area of the latitude slice to get the flux at each wavelength
      for (m = 0; m < nlam; m++)
        flux[m] += B[k][m] * area;
    } 
    
  }
  
  // Scale flux by REARTH**2 and convert to W m^2 / m^2 / um / sr
  for (m = 0; m < nlam; m++) 
    flux[m] *= REARTH * REARTH * MICRON;
  
  // Free all of the ellipse instances
  for (j = 0; j < f; j+=2) {
    free(functions[j].ellipse);
  }  
   
}

void UnoccultedFlux(double r, double theta, double albedo, double irrad, double polyeps1, double polyeps2, int maxpolyiter, int adaptive, int nu, int nlat, int nlam, double u[nu], double lambda[nlam], double flux[nlam]) {
  /*
  
  */
  
  double x0[1] = {0};
  double y0[1] = {0};
  double ro[1] = {2 * r};
  
  // Hack: compute the occulted flux with a single huge occultor
  OccultedFlux(r, 1, x0, y0, ro, theta, albedo, irrad, polyeps1, polyeps2, maxpolyiter, adaptive, nu, nlat, nlam, u, lambda, flux);
    
}