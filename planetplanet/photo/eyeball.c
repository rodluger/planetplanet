#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ppo.h"
#include "complex.h"

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

double EyeballTemperature(double lat, double irrad, double albedo, double tnight) {
  /*
  
  */
  
  double temp;
  
  if (lat < PI / 2) {
    temp = pow((irrad * cos(lat) * (1 - albedo)) / SBOLTZ, 0.25);
    if (temp > tnight)
      return temp;
    else
      return tnight;
  } else {
    return tnight;
  }
}

double Blackbody(double lambda, double T) {
  /*
  
  */
  
  // Planck's law
  double a = 2 * HPLANCK * CLIGHT * CLIGHT / (lambda * lambda * lambda * lambda * lambda);
  double b = HPLANCK * CLIGHT / (lambda * KBOLTZ * T);
  return a / (exp(b) - 1);
}

void GetRoots(double a, double b, double xE, double yE, double xC, double yC, double r, double polyeps1, double polyeps2, int maxpolyiter, double roots[4]) {
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
  
  // Get the real roots (up to 4)
  roots[0] = NAN;
  roots[1] = NAN;
  roots[2] = NAN;
  roots[3] = NAN;
  j = 0;
  for (i = 1; i < 5; i++) {
    if (!(isnan(croots[i].r)) && (fabs(croots[i].i) < MAXIM)) {
      roots[j] = croots[i].r + xC;
      j += 1;
    }
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
  if ((fabs(fabs(theta) - PI / 2)) < SMALL) {
    
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

void SurfaceIntensity(double albedo, double irrad, double tnight, double teff, int nlat, double latgrid[nlat], int nlam, double lambda[nlam], int nu, double u[nu * nlam], double B[nlat + 1][nlam]) {
  /*
  Returns the blackbody intensity at the center of each latitude slice 
  evaluated at a given array of wavelengths.
  
  */
  
  int i, j, k;
  double latitude, coslat;
  double T = 0.;
  double I1[nlam];
  double x;
  
  // Loop over each slice
  for (i = 0; i < nlat + 1; i++) {
    
    // Get the latitude halfway between the `latgrid` points.
    // In the regular grid, B[0] is the intensity of the region between the 
    // sub-stellar point and latgrid[0] (the location of the first latitude ellipse); 
    // B[1] is the intensity between latgrid[1] and latgrid[2], and so forth.
    // B[nlat] is the intensity between latgrid[nlat - 1] and the anti-stellar point.
    if (i == 0) {
      latitude = latgrid[0] - 0.5 * (latgrid[1] - latgrid[0]);
      if (latitude < 0) latitude = 0;
    } else if (i == nlat) {
      latitude = latgrid[nlat - 1] + 0.5 * (latgrid[nlat - 1] - latgrid[nlat - 2]);
      if (latitude > PI) latitude = PI;
    } else
      latitude = 0.5 * (latgrid[i] + latgrid[i - 1]);

    // Get its blackbody temperature
    if (teff == 0) {
      
      // No limb darkening
      T = EyeballTemperature(latitude, irrad, albedo, tnight);

      // Loop over each wavelength bin
      for (j = 0; j < nlam; j++) {
        
        // Get the blackbody intensity (W / m^2 / m / sr)
        B[i][j] = Blackbody(lambda[j], T);
      
      }
    
    } else {
    
      // Polynomial limb darkening
      
      // TODO: We are assuming that the intensity at the center
      // of the body is the blackbody intensity at its effective
      // temperature. This is not strictly true, but probably OK
      // for modest limb darkening. We should eventually find the
      // correct normalization.
      for (j = 0; j < nlam; j++) {
        I1[j] = Blackbody(lambda[j], teff);
        B[i][j] = I1[j];
      }
      
      // Pre-compute mu
      coslat = cos(latitude);
      
      // Loop over the coefficient order
      for (k = 0; k < nu; k++) {
      
        // The Taylor expansion is in (1 - mu)
        x = pow(1 - coslat, k + 1);
        
        // Compute the wavelength-dependent intensity
        for (j = 0; j < nlam; j++) {
          B[i][j] -= u[nlam * k + j] * I1[j] * x;
        } 
      } 
    }
  }   
}

double curve(double x, FUNCTION function) {
  return function.curve(x, function.ellipse);
}

double integral(double x0, double x1, FUNCTION function, int *oob) {
  return function.integral(x0, x1, function.ellipse, oob);
}

double fupper(double x, ELLIPSE *ellipse) {
  double A;
  if (ellipse->circle) {
    A = ellipse->r * ellipse->r - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 + sqrt(A);
  } else {
    A = ellipse->b * ellipse->b - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 + (ellipse->a / ellipse->b) * sqrt(A);
  }
}

double flower(double x, ELLIPSE *ellipse) {
  double A;
  if (ellipse->circle) {
    A = ellipse->r * ellipse->r - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 - sqrt(A);
  } else {
    A = ellipse->b * ellipse->b - (x - ellipse->x0) * (x - ellipse->x0);
    if (fabs(A) < TINY) A = 0;
    else if ((x > ellipse->xmax) || (x < ellipse->xmin)) return NAN;
    return ellipse->y0 - (ellipse->a / ellipse->b) * sqrt(A);
  }
}

double iupper(double xL, double xR, ELLIPSE *ellipse, int *oob) {
  /*
  
  */
  
  double yL,yR,zL,zR,FL,FR;
  
  // Catch NaNs silently
  if ((xL > ellipse->xmax) || (xL < ellipse->xmin)) *oob = 1;
  if ((xR > ellipse->xmax) || (xR < ellipse->xmin)) *oob = 1;
  
  yL = xL - ellipse->x0;
  yR = xR - ellipse->x0;
  if (ellipse->circle) {
    zL = sqrt((ellipse->r - yL) * (ellipse->r + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->r - yR) * (ellipse->r + yR));
    if (isnan(zR)) zR = 0;
    FL = 0.5 * (yL * zL + ellipse->r * ellipse->r * atan(yL / zL)) + ellipse->y0 * xL;
    FR = 0.5 * (yR * zR + ellipse->r * ellipse->r * atan(yR / zR)) + ellipse->y0 * xR;
  } else {
    zL = sqrt((ellipse->b - yL) * (ellipse->b + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->b - yR) * (ellipse->b + yR));
    if (isnan(zR)) zR = 0;
    FL = (ellipse->a / (2 * ellipse->b)) * (yL * zL + ellipse->b * ellipse->b * atan(yL / zL)) + ellipse->y0 * xL;
    FR = (ellipse->a / (2 * ellipse->b)) * (yR * zR + ellipse->b * ellipse->b * atan(yR / zR)) + ellipse->y0 * xR;
  }
  return FR - FL;
}

double ilower(double xL, double xR, ELLIPSE *ellipse, int *oob) {
  /*
  
  */
  
  double yL,yR,zL,zR,FL,FR;
  
  // Catch NaNs silently
  if ((xL > ellipse->xmax) || (xL < ellipse->xmin)) *oob = 1;
  if ((xR > ellipse->xmax) || (xR < ellipse->xmin)) *oob = 1;
    
  yL = xL - ellipse->x0;
  yR = xR - ellipse->x0;
  if (ellipse->circle) {
    zL = sqrt((ellipse->r - yL) * (ellipse->r + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->r - yR) * (ellipse->r + yR));
    if (isnan(zR)) zR = 0;
    FL = -0.5 * (yL * zL + ellipse->r * ellipse->r * atan(yL / zL)) + ellipse->y0 * xL;
    FR = -0.5 * (yR * zR + ellipse->r * ellipse->r * atan(yR / zR)) + ellipse->y0 * xR;
  } else {
    zL = sqrt((ellipse->b - yL) * (ellipse->b + yL));
    if (isnan(zL)) zL = 0;
    zR = sqrt((ellipse->b - yR) * (ellipse->b + yR));
    if (isnan(zR)) zR = 0;
    FL = -(ellipse->a / (2 * ellipse->b)) * (yL * zL + ellipse->b * ellipse->b * atan(yL / zL)) + ellipse->y0 * xL;
    FR = -(ellipse->a / (2 * ellipse->b)) * (yR * zR + ellipse->b * ellipse->b * atan(yR / zR)) + ellipse->y0 * xR;
  }
  return FR - FL;
}

void AddLatitudeSlice(double latitude, double r, int no, double x0[no], double y0[no], double ro[no], 
                      double theta, double polyeps1, double polyeps2, int maxpolyiter, int maxvertices,
                      int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f){
  /*
  
  */    
  
  int i, j;
  double xlimb, x, y;
  double roots[4];
  double ro2[no];
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += TINY * ro2[i];
  }
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
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
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
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
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
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
        break;    
      }
    }
    y = flower(x, ellipse);
    for (i = 0; i < no; i++) {
      if ((x - x0[i]) * (x - x0[i]) + (y - y0[i]) * (y - y0[i]) < ro2[i]) {
        vertices[(*v)++] = x;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
        break;    
      }
    }
  }
  
  // Ellipse-occultor vertices
  for (i = 0; i < no; i++) {
    
    // Is the occulted planet entirely within the occultor?
    if (sqrt(x0[i] * x0[i] + y0[i] * y0[i]) + r + SMALL < ro[i]) continue;
    
    // Solve the quartic
    GetRoots(ellipse->a, ellipse->b, ellipse->x0, ellipse->y0, x0[i], y0[i], ro[i], polyeps1, polyeps2, maxpolyiter, roots);
    
    for (j = 0; j < 4; j++) {
      if (!(isnan(roots[j])) && (((theta > 0) && (roots[j] > ellipse->x0 - xlimb)) || ((theta <= 0) && (roots[j] < ellipse->x0 - xlimb)))) {
        vertices[(*v)++] = roots[j];
        if (*v > maxvertices) {
            printf("ERROR: Maximum number of vertices exceeded.\n");
            abort();
          }
      }
    }
    
  }
  
  // Finally, the boundary curves and their integrals
  functions[*f].ellipse = ellipse;
  functions[*f].curve = fupper;
  functions[(*f)++].integral = iupper;
  functions[*f].ellipse = ellipse;
  functions[*f].curve = flower;
  functions[(*f)++].integral = ilower;
  if (*f > maxfunctions) {
    printf("ERROR: Maximum number of functions exceeded.\n");
    abort();
  }
  
}

void AddOccultors(double r, int no, double x0[no], double y0[no], double ro[no], int maxvertices, int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  int i;
  double r2 = r * r;
  r2 += TINY * r2;
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
    if (occultor[i]->xmin * occultor[i]->xmin + occultor[i]->y0 * occultor[i]->y0 < r2) {
      vertices[(*v)++] = occultor[i]->xmin;  
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }  
    }
    
    // Occultor x maximum
    occultor[i]->xmax = x0[i] + occultor[i]->r;
    if (occultor[i]->xmax * occultor[i]->xmax + occultor[i]->y0 * occultor[i]->y0 < r2) {
      vertices[(*v)++] = occultor[i]->xmax;    
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
    }

    // Finally, the boundary curves and their integrals
    functions[*f].ellipse = occultor[i];
    functions[*f].curve = fupper;
    functions[(*f)++].integral = iupper;
    functions[*f].ellipse = occultor[i];
    functions[*f].curve = flower;
    functions[(*f)++].integral = ilower;
    if (*f > maxfunctions) {
      printf("ERROR: Maximum number of functions exceeded.\n");
      abort();
    }
  
  }
  
}

void AddOcculted(double r, int no, double x0[no], double y0[no], double ro[no], int maxvertices, int maxfunctions, double *vertices, int *v, FUNCTION *functions, int *f) {
  /*
  
  */ 

  int i;
  double ro2[no];
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += TINY * ro2[i];
  }
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
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
      break;    
    }
  }  
  
  // Occulted planet x maximum
  occulted->xmax = r;
  // Is this point inside __at least__ one occultor?
  for (i = 0; i < no; i++) {
    if ((occulted->xmax - x0[i]) * (occulted->xmax - x0[i]) + y0[i] * y0[i] < ro2[i]) {
      vertices[(*v)++] = occulted->xmax;
      if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
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
        vertices[(*v)++] = -x * cost + y * sint;
        vertices[(*v)++] = -x * cost - y * sint;
        if (*v > maxvertices) {
          printf("ERROR: Maximum number of vertices exceeded.\n");
          abort();
        }
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
  if (*f > maxfunctions) {
    printf("ERROR: Maximum number of functions exceeded.\n");
    abort();
  }
  
}

void OccultedFlux(double r, int no, double x0[no], double y0[no], double ro[no], double theta, double albedo, 
                  double irrad, double tnight, double teff, double polyeps1, double polyeps2, int maxpolyiter, double mintheta, int maxvertices,
                  int maxfunctions, int adaptive, int nu, int nlat, int nlam, double u[nu], double lambda[nlam], 
                  double flux[nlam], int quiet, int *iErr) {
  /*
  
  */
  
  int i, j, k, m;
  int oob;
  int b = 0;
  int v = 0;
  int f = 0;
  int good;
  double lmin, lmax, lat;
  double xL, xR, x, y, area;
  double r2 = r * r;
  r2 += SMALL * r2;
  double d, dmin, dmax;
  double ro2[no];
  for (i = 0; i < no; i++) {
    ro2[i] = ro[i] * ro[i];
    ro2[i] += SMALL * ro2[i];
  }
  double vertices[maxvertices];
  FUNCTION functions[maxfunctions];
  FUNCTION boundaries[maxfunctions];
  double B[nlat * no + 1][nlam];
  double latgrid[nlat * no];
  *iErr = ERR_NONE;
  
  // Avoid the singular point
  if (fabs(theta) < mintheta)
    theta = mintheta;

  // Zero out the flux
  for (m = 0; m < nlam; m++) 
    flux[m] = 0.;
  
  // If we're doing limb darkening, theta must be pi/2 (full phase)
  if (teff > 0)
    theta = PI / 2.;
  
  // Generate all the shapes and get their vertices and curves
  AddOccultors(r, no, x0, y0, ro, maxvertices, maxfunctions, vertices, &v, functions, &f);     
  AddOcculted(r, no, x0, y0, ro, maxvertices, maxfunctions, vertices, &v, functions, &f); 
  
  // Compute the latitude grid
  if (adaptive && (teff > 0)) {
  
    // Adaptive latitude grid
    // Loop over each occultor
    for (i = 0; i < no; i++) {
    
      // Distance to occultor center
      d = sqrt(x0[i] * x0[i] + y0[i] * y0[i]);
      
      // Distance to far side
      dmax = (d + ro[i]) / r + SMALL;
      if (dmax >= 1 - SMALL) dmax = 1;
      lmax = asin(dmax);
          
      // Distance to near side
      dmin = (d - ro[i]) / r - SMALL;
      if (dmin <= SMALL) dmin = 0;
      lmin = asin(dmin);

      // Add to latitude grid
      for (j = 0; j < nlat; j++) {
        latgrid[nlat * i + j] = (j + 1.) / (nlat + 1.) * (lmax - lmin) + lmin;
      }
      
      // Adjust the boundaries
      latgrid[nlat * i] = lmin;
      latgrid[nlat * (i + 1) - 1] = lmax;
      
    }

    // Sort the grid
    qsort(latgrid, nlat * no, sizeof(double), dblcomp);
    
  } else {
  
    // Linearly spaced latitude grid
    lmin = 0 - SMALL;
    lmax = PI + SMALL;
    for (i = 0; i < nlat * no; i++) {
      latgrid[i] = (i + 1.) / (nlat * no + 1.) * (lmax - lmin) + lmin;
    }
    
  }

  // Add the ellipses  
  for (i = 0; i < nlat * no; i++) {
    AddLatitudeSlice(latgrid[i], r, no, x0, y0, ro, theta, polyeps1, polyeps2, maxpolyiter, maxvertices, maxfunctions, vertices, &v, functions, &f);
  }
  // Pre-compute the surface intensity in each latitude slice
  SurfaceIntensity(albedo, irrad, tnight, teff, nlat * no, latgrid, nlam, lambda, nu, u, B);
  
  // Sort the vertices
  qsort(vertices, v, sizeof(double), dblcomp);
      
  // Loop over all vertices
  for (i = 0; i < v - 1; i++) {

    // Get the integration limits, with a
    // small perturbation for stability
    xL = vertices[i] + TINY;
    xR = vertices[i+1] - TINY;
    
    // Check if they are identical
    if (xR <= xL + TINY)
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
    
      // Get the midpoint
      y = 0.5 * (boundaries[j + 1].y + boundaries[j].y);
      
      // Is it in at least one occultor?
      good = 0;
      for (k = 0; k < no; k++) {
        if ((x - x0[k]) * (x - x0[k]) + (y - y0[k]) * (y - y0[k]) < ro2[k]) {
          good = 1;
          break;
        }
      }
      if (!good) continue;
      
      // The area of each region is just the difference of successive integrals
      oob = 0;
      area = integral(xL, xR, boundaries[j + 1], &oob) - integral(xL, xR, boundaries[j], &oob);      
      if (oob) *iErr = ERR_OOB;
      
      // Get the latitude of the midpoint
      lat = Latitude(x, y, r, theta);
      
      // Get the index `k` of the latitude grid *above* this latitude.
      // B[k] is the intensity of this region.
      for (k = 0; k < nlat * no; k++) {
        if (latgrid[k] > lat)
          break;
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

void UnoccultedFlux(double r, double theta, double albedo, double irrad, double tnight, double teff, double polyeps1, 
                    double polyeps2, int maxpolyiter, double mintheta, int maxvertices, int maxfunctions, int adaptive, 
                    int nu, int nlat, int nlam, double u[nu], double lambda[nlam], double flux[nlam], int quiet, int *iErr) {
  /*
  
  */
  
  double x0[1] = {0};
  double y0[1] = {0};
  double ro[1] = {2 * r};
  
  // Hack: compute the occulted flux with a single huge occultor
  OccultedFlux(r, 1, x0, y0, ro, theta, albedo, irrad, tnight, teff, polyeps1, polyeps2, maxpolyiter, mintheta, maxvertices, maxfunctions, adaptive, nu, nlat, nlam, u, lambda, flux, quiet, iErr);
    
}