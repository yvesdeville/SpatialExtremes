#include "header.h"

double whittleMatern(double *dist, int nPairs, double sill, double range,
		     double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When ans != 0.0, the whittle-matern parameters are ill-defined.
  
  int i;
  
  //Some preliminary steps: Valid points?
  if (smooth < EPS)
    return R_pow_di(1 - smooth + EPS, 2) * MINF;

  else if (smooth > 150)
    //Required because it could lead to infinite rho values
    return R_pow_di(smooth - 150, 2) * MINF;

  if (range <= 0.0)
    return R_pow_di(1 - range, 2) * MINF;

  if (sill <= 0.0)
    return R_pow_di(1 - sill, 2) * MINF;
  
  else if (sill > 1)
    //the 1.02 factor is here to avoid problem with non
    //feasible region
    return R_pow_di(sill, 2) * MINF;
  
  for (i=0;i<nPairs;i++){

    rho[i] = sill * R_pow(2, 1 - smooth) / gammafn(smooth) *
      R_pow(dist[i] / range, smooth) * 
      bessel_k(dist[i] / range, smooth, 1);
    
  }

  return 0.0;
}

double cauchy(double *dist, int nPairs, double sill, double range,
	      double smooth, double *rho){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When ans != 0.0, the cauchy parameters are ill-defined.

  int i;
  
  //Some preliminary steps: Valid points?
  if (smooth < 0)
    return R_pow_di(1 - smooth, 2) * MINF;

  if (range <= 0.0)
    return R_pow_di(1 - range, 2) * MINF;

  if (sill <= 0.0)
    return R_pow_di(1 - sill, 2) * MINF;
  
  else if (sill > 1)
    return R_pow_di(sill, 2) * MINF;

  for (i=0;i<nPairs;i++)
    rho[i] = sill * R_pow(1 + dist[i] * dist[i] / range / range, -smooth);
    
  return 0.0;
}

double powerExp(double *dist, int nPairs, double sill, double range,
		double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  int i;
    
  //Some preliminary steps: Valid points?
  if (smooth < 0)
    return R_pow_di(1 - smooth, 2) * MINF;

  else if (smooth >= 2.0)
    return R_pow_di(smooth - 1, 2) * MINF;

  if (range <= 0.0)
    return R_pow_di(1 - range, 2) * MINF;

  if (sill <= 0.0)
    return R_pow_di(1 - sill, 2) * MINF;
  
  else if (sill > 1)
    return R_pow_di(sill, 2) * MINF;
  
  for (i=0;i<nPairs;i++)
    rho[i] = sill * exp(-R_pow(dist[i] / range, smooth));
    
  return 0.0;
}

double mahalDistFct(double *distVec, int nPairs, double *cov11,
		    double *cov12, double *cov22, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When ans != 0.0, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  double det;

  det = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (*cov11 <= 0)
    return R_pow_di(1 - *cov11, 2) * MINF;

  if (*cov22 <= 0)
    return R_pow_di(1 - *cov22, 2) * MINF;
  
  if (det <= 1e-10)
    return R_pow_di(1 - det + 1e-10, 2) * MINF;
  
  for (i=0;i<nPairs;i++){

    mahal[i] = (*cov11 * distVec[nPairs + i] * distVec[nPairs + i] -
		2 * *cov12 * distVec[i] * distVec[nPairs + i] +
		*cov22 * distVec[i] * distVec[i]) / det;
    
    mahal[i] = sqrt(mahal[i]);
  }
  
  return 0.0;
}

double mahalDistFct3d(double *distVec, int nPairs, double *cov11,
		      double *cov12, double *cov13, double *cov22, 
		      double *cov23, double *cov33, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 3D
  //When ans != 0.0, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  double det, detMin;

  det = *cov11 * *cov22 * *cov33 - R_pow_di(*cov12, 2) * *cov33 -
    *cov11 * R_pow_di(*cov23, 2) + 2 * *cov12 * *cov13 * *cov23 -
    R_pow_di(*cov13, 2) * *cov22;
  detMin = *cov11 * *cov22 - R_pow_di(*cov12, 2);
  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (det <= 1e-10)
    return R_pow_di(1 - det + 1e-10, 2) * MINF;

  if (*cov11 <= 0)
    return R_pow_di(1 - *cov11, 2) * MINF;

  if (detMin <= 0)
    return R_pow_di(1 - detMin, 2);
  
  for (i=0;i<nPairs;i++){

    mahal[i] = (*cov11 * *cov22 * distVec[2 * nPairs + i] * distVec[2 * nPairs + i] -
		*cov12 * *cov12 * distVec[2 * nPairs + i] * distVec[2 * nPairs + i] -
		2 * *cov11 * *cov23 * distVec[nPairs + i] * distVec[2 * nPairs + i] +
		2 * *cov12 * *cov13 * distVec[nPairs + i] * distVec[2 * nPairs + i] +
		2 * *cov12 * *cov23 * distVec[i] * distVec[2 * nPairs + i] -
		2 * *cov13 * *cov22 * distVec[i] * distVec[2 * nPairs + i] +
		*cov11 * *cov33 * distVec[nPairs + i] * distVec[nPairs + i] - 
		*cov13 * *cov13 * distVec[nPairs + i] * distVec[nPairs + i] -
		2 * *cov12 * *cov33 * distVec[i] * distVec[nPairs + i] +
		2 * *cov13 * *cov23 * distVec[i] * distVec[nPairs + i] +
		*cov22 * *cov33 * distVec[i] * distVec[i] -
		*cov23 * *cov23 * distVec[i] * distVec[i]) / det;

    mahal[i] = sqrt(mahal[i]);
  }
  
  return 0.0;
}

double geomCovariance(double *dist, int nPairs, int covmod,
		      double sigma2, double sill, double range,
		      double smooth, double *rho){

  //This function computes the geometric gaussian covariance function
  //between each pair of locations.
  //When ans != 0.0, the parameters are ill-defined.
  int i;
  double ans = 0.0;

  switch (covmod){
  case 1:
    ans = whittleMatern(dist, nPairs, sill, range, smooth, rho);
    break;
  case 2:
    ans = cauchy(dist, nPairs, sill, range, smooth, rho);
    break;
  case 3:
    ans = powerExp(dist, nPairs, sill, range, smooth, rho);
    break;
  }

  if (ans != 0.0)
    return ans;

  for (i=0;i<nPairs;i++)
    rho[i] = sqrt(2 * sigma2 * (1 - rho[i]));

  return ans;
}

double nsgeomCovariance(double *dist, int nSite, int covmod,
			double *sigma2, double sill, double range,
			double smooth, double *rho){
  
  int i, j, currentPair = 0, nPairs = nSite * (nSite - 1) / 2;
  double ans = 0.0;

  switch (covmod){
  case 1:
    ans = whittleMatern(dist, nPairs, sill, range, smooth, rho);
    break;
  case 2:
    ans = cauchy(dist, nPairs, sill, range, smooth, rho);
    break;
  case 3:
    ans = powerExp(dist, nPairs, sill, range, smooth, rho);
    break;
  }

  if (ans != 0.0)
    return ans;

  for (i=0;i<(nSite-1);i++){
    for (j=i+1;j<nSite;j++){
      rho[currentPair] = sqrt(sigma2[i] - 2 * sqrt(sigma2[i] * sigma2[j]) * 
			      rho[currentPair] + sigma2[j]);
      currentPair++;
    }
  }

  return ans;
}
      

  
  
