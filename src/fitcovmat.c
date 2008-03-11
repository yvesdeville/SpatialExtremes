#include "header.h"

void fitcovmat2d(double *cov11, double *cov12, double *cov22,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  *ans = -mahalDistFct(distVec, *nPairs, cov11, cov12, cov22,
		       mahalDist);

  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++)
    *ans += R_pow_di((2 * pnorm(mahalDist[i] / 2, 0.0, 1.0, 1, 0)
		      - extcoeff[i]) / weights[i], 2);
  return;
}


void fitcovmat3d(double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33,
		 int *nPairs, double *distVec, double *extcoeff,
		 double *weights, double *ans){
  
  int i;
  double *mahalDist;

  mahalDist = (double *)R_alloc(*nPairs, sizeof(double));

  *ans = -mahalDistFct3d(distVec, *nPairs, cov11, cov12, cov13, 
			 cov22, cov23, cov33, mahalDist);

  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++)
    *ans += R_pow_di((2 * pnorm(mahalDist[i] / 2, 0.0, 1.0, 1, 0) -
		      extcoeff[i]) / weights[i], 2);
  
  return;
}

void fitcovariance(int *covmod, double *sill, double *range, double *smooth,
		   int *nPairs, double *dist, double *extcoeff,
		   double *weights, double *ans){
  
  int i;
  double *rho;

  rho = (double *)R_alloc(*nPairs, sizeof(double));

  switch (*covmod){
  case 1:
    *ans = -whittleMatern(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    *ans = -cauchy(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    *ans = -powerExp(dist, *nPairs, *sill, *range, *smooth, rho);
    break;
  }

  if (*ans != 0.0)
    return;

  for (i=0;i<*nPairs;i++)
    *ans += R_pow_di((1 + sqrt(1 - 0.5 * (rho[i] + 1)) -
		      extcoeff[i]) / weights[i], 2);

  return;
}
