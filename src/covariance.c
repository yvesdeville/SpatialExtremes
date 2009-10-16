#include "header.h"

double whittleMatern(double *dist, int nPairs, double sill, double range,
		     double smooth, double *rho){

  //This function computes the whittle-matern covariance function
  //between each pair of locations.
  //When ans != 0.0, the whittle-matern parameters are ill-defined.
  
  int i;
  const double cst = sill * R_pow(2, 1 - smooth) / gammafn(smooth),
    irange = 1 / range;

  //Some preliminary steps: Valid points?
  if (smooth < EPS)
    return (1 - smooth + EPS) * (1 - smooth + EPS) * MINF;

  else if (smooth > 150)
    /* Required because it could lead to infinite rho values - gamma
    function */
    return (smooth - 149) * (smooth - 149) * MINF;

  if (range < 1e-2)
    return (1.01 - range) * (1.01 - range) * MINF;

  if (sill <= 1e-2)
    return (1.01 - sill) * (1.01 - sill) * MINF;
  
  else if (sill > 1)
    return sill *sill * MINF;
  
  for (i=nPairs;i--;){
    double cst2;
    cst2 = dist[i] * irange;
    rho[i] = cst * R_pow(cst2, smooth) * bessel_k(cst2, smooth, 1);
  }

  return 0.0;
}

double cauchy(double *dist, int nPairs, double sill, double range,
	      double smooth, double *rho){

  //This function computes the cauchy covariance function between each
  //pair of locations.
  //When ans != 0.0, the cauchy parameters are ill-defined.

  int i;
  const double irange2 = 1 / (range * range);
  
  //Some preliminary steps: Valid points?
  if (smooth < 0)
    return (1 - smooth) * (1 - smooth) * MINF;
  
  if (range <= 0.0)
    return (1 - range) * (1 - range)* MINF;

  if (sill <= 0.0)
    return (1 - sill) * (1 - sill) * MINF;
  
  else if (sill > 1)
    return sill * sill * MINF;

  for (i=nPairs;i--;)
    rho[i] = sill * R_pow(1 + dist[i] * dist[i] * irange2, -smooth);
    
  return 0.0;
}

double powerExp(double *dist, int nPairs, double sill, double range,
		double smooth, double *rho){

  //This function computes the powered exponential covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  int i;
  const double irange = 1 / range;
    
  //Some preliminary steps: Valid points?
  if ((smooth < 0) || (smooth > 2.0))
    return (1 - smooth) * (1 - smooth) * MINF;

  if (range <= 0.0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0.0)
    return (1 - sill) * (1 - sill) * MINF;
  
  else if (sill > 1)
    return sill * sill * MINF;
  
  for (i=nPairs;i--;)
    rho[i] = sill * exp(-R_pow(dist[i] * irange, smooth));
    
  return 0.0;
}

double bessel(double *dist, int nPairs, int dim, double sill,
	      double range, double smooth, double *rho){
  //This function computes the bessel covariance function
  //between each pair of locations.
  //When ans != 0.0, the powered exponential parameters are ill-defined.

  int i;
  const double irange = 1 / range, cst = sill * R_pow(2, smooth) * gammafn(smooth + 1);

  //Some preliminary steps: Valid points?
  if (smooth < (0.5 * (dim - 2)))
    return (1 + 0.5 * (dim - 2) - smooth) * (1 + 0.5 * (dim - 2) - smooth) * MINF;

  else if (smooth > 30)
    //Require as bessel_j will be numerically undefined
    return (smooth - 29) * (smooth - 29) * MINF;

  if (range <= 0)
    return (1 - range) * (1 - range) * MINF;

  if (sill <= 0)
    return (1 - sill) * (1 - sill) * MINF;

  else if (sill > 1)
    return sill * sill * MINF;

  for (i=nPairs;i--;){
    double cst2;
    cst2 = dist[i] * irange;

    rho[i] = cst * R_pow(cst2, -smooth) * bessel_j(cst2, smooth);
  }

  return 0.0;
}
  
double mahalDistFct(double *distVec, int nPairs, double *cov11,
		    double *cov12, double *cov22, double *mahal){
  //This function computes the mahalanobis distance between each pair
  //of locations. Currently this function is only valid in 2D
  //When ans != 0.0, the covariance matrix and/or the mahalanobis
  //distance is ill-defined.
  
  int i;
  const double det = *cov11 * *cov22 - *cov12 * *cov12,
    idet = 1 / det;

  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (*cov11 <= 0)
    return (1 - *cov11) * (1 - *cov11) * MINF;

  if (*cov22 <= 0)
    return (1 - *cov22) * (1 - *cov22) * MINF;
  
  if (det <= 1e-10)
    return (1 - det + 1e-10) * (1 - det + 1e-10) * MINF;
  
  for (i=nPairs;i--;){
    mahal[i] = (*cov11 * distVec[nPairs + i] * distVec[nPairs + i] -
		2 * *cov12 * distVec[i] * distVec[nPairs + i] +
		*cov22 * distVec[i] * distVec[i]) * idet;
    
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
  const double det = *cov11 * *cov22 * *cov33 - *cov12 * *cov12 * *cov33 -
    *cov11 * *cov23 * *cov23 + 2 * *cov12 * *cov13 * *cov23 -
    *cov13 * *cov13 * *cov22,
    detMin = *cov11 * *cov22 - *cov12 * *cov12,
    idet = 1 / det;

  //We test if the covariance matrix is *not* nonnegative
  //definite e.g. all minor determinant are negative or 0
  if (det <= 1e-10)
    return (1 - det + 1e-10) * (1 - det + 1e-10) * MINF;

  if (*cov11 <= 0)
    return (1 - *cov11) * (1 - *cov11) * MINF;

  if (detMin <= 0)
    return (1 - detMin) * (1 - detMin) * MINF;
  
  for (i=nPairs;i--;){

    mahal[i] = ((*cov22 * *cov33 - *cov23 * *cov23) * distVec[i] * distVec[i] +
		2 * (*cov13 * *cov23 - *cov12 * *cov33) * distVec[i] * distVec[nPairs + i] +
		2 * (*cov12 * *cov23 - *cov13 * *cov22) * distVec[i] * distVec[2 *nPairs + i] +
		(*cov11 * *cov33 - *cov13 * *cov13) * distVec[nPairs + i] * distVec[nPairs + i] +
		2 * (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + i] * distVec[2 * nPairs + i] +
		(*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + i] * distVec[2 * nPairs + i]) *
      idet;

    mahal[i] = sqrt(mahal[i]);
  }
  
  return 0.0;
}

double geomCovariance(double *dist, int nPairs, int dim, int covmod,
		      double sigma2, double sill, double range,
		      double smooth, double *rho){

  //This function computes the geometric gaussian covariance function
  //between each pair of locations.
  //When ans != 0.0, the parameters are ill-defined.
  int i;
  const double twiceSigma2 = 2 * sigma2;
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
  case 4:
    ans = bessel(dist, nPairs, dim, sill, range, smooth, rho);
    break;
  }

  if (ans != 0.0)
    return ans;

  for (i=nPairs;i--;)    
    rho[i] = sqrt(twiceSigma2 * (1 - rho[i]));
  
  return ans;
}

double nsgeomCovariance(double *dist, int nSite, int dim, int covmod,
			double *sigma2, double sill, double range,
			double smooth, double *rho){
  
  int i, j, currentPair = 0;
  const int nPairs = nSite * (nSite - 1) / 2;
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
  case 4:
    ans = bessel(dist, nPairs, dim, sill, range, smooth, rho);
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
      

  
  
