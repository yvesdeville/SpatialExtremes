#include "header.h"

/* These functions estimate the spatial dependence parameters by using
   a least square optimization. More precisely, it compares the
   expected extremal coefficients to the "observed" ones.
*/

void fitcovmat2d(double *cov11, double *cov12, double *cov22,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist, pen = 1.0;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  pen += mahalDistFct(distVec, *nPairs, cov11, cov12, cov22,
		     mahalDist);

  for (i=0;i<*nPairs;i++)
    *ans += pen * R_pow_di((2 * pnorm(mahalDist[i] / 2, 0.0, 1.0, 1, 0)
			    - extcoeff[i]) / weights[i], 2);

  return;
}


void fitcovmat3d(double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist, pen = 1.0;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  pen += mahalDistFct3d(distVec, *nPairs, cov11, cov12, cov13, 
			cov22, cov23, cov33, mahalDist);

  for (i=0;i<*nPairs;i++)
    *ans += pen * R_pow_di((2 * pnorm(mahalDist[i] / 2, 0.0, 1.0, 1, 0) -
			    extcoeff[i]) / weights[i], 2);

  return;
}

void fitcovariance(int *covmod, double *sill, double *range, double *smooth,
		   int *nPairs, double *dist, double *extcoeff,
		   double *weights, double *ans){
  
  int i;
  double *rho, pen = 1.0;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  switch (*covmod){
  case 1:
    pen += whittleMatern(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    pen += cauchy(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    pen += powerExp(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  }

  for (i=0;i<*nPairs;i++)
    *ans += pen * R_pow_di((1 + sqrt(1 - 0.5 * (rho[i] + 1)) -
			    extcoeff[i]) / weights[i], 2);

  return;
}
