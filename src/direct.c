#include "header.h"

void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
		 double *nugget, double *sill, double *range,
		 double *smooth, double *covMat){

  int i, j, currentPair, nPairs, effnSite = *nSite, zero = 0;
  const double one = 1, dzero = 0;
  double *dist, *rho, flag = 0, *coordGrid;

  if (*grid)
    effnSite = R_pow_di(effnSite, *dim);

  nPairs = effnSite * (effnSite - 1) / 2;

  dist = (double *)R_alloc(nPairs, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  coordGrid = (double *)R_alloc(effnSite * *dim, sizeof(double));

  if (*grid){
    //Coord specify a grid
    for (i=*nSite;i--;)
      for (j=*nSite;j--;){
	coordGrid[j + i * *nSite] = coord[i];
	coordGrid[*nSite * (*nSite + i) + j] = coord[j];
      }

    distance(coordGrid, dim, &effnSite, &zero, dist);
  }

  else
    //Coord don't specify a grid
    distance(coord, dim, nSite, &zero, dist);

  switch (*covmod){
  case 1:
    flag = whittleMatern(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 2:
    flag = cauchy(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 3:
    flag = powerExp(dist, nPairs, dzero, one, *range, *smooth, rho);
    break;
  case 4:
    flag = bessel(dist, nPairs, *dim, dzero, one, *range, *smooth, rho);
    break;
  case 6:
    if (*grid)
      flag = fbm(coordGrid, dist, *dim, effnSite, one, *range, *smooth, rho);

    else
      flag = fbm(coord, dist, *dim, effnSite, one, *range, *smooth, rho);

    break;
  }
  
  if (flag != 0.0)
    error("The covariance parameters seem to be ill-defined. Please check\n");
  
  //Fill the non-diagonal elements of the covariance matrix
  currentPair = -1;
  for (i=0;i<(effnSite-1);i++){
    for (j=i+1;j<effnSite;j++){
      currentPair++;
      covMat[effnSite * i + j] = covMat[effnSite * j + i] = *sill * rho[currentPair];
    }
  }
  
  //Fill the diagonal elements of the covariance matrix
  if (*covmod == 6){
    //Fractional brownian
    double irange = 1 / *range;

    if (*grid){
      for (i=effnSite;i--;){
	covMat[i * (effnSite + 1)] = 0;
	
	for (j=*dim;j--;)
	  covMat[i * (effnSite + 1)] += coordGrid[i + j * effnSite] * coordGrid[i + j * effnSite];
	
	covMat[i * (effnSite + 1)] = 2 * R_pow(sqrt(covMat[i * (effnSite + 1)]) * irange, *smooth);
      }
    }

    else {
      for (i=effnSite;i--;){
	covMat[i * (effnSite + 1)] = 0;
	
	for (j=*dim;j--;)
	  covMat[i * (effnSite + 1)] += coord[i + j * effnSite] * coord[i + j * effnSite];
	
	covMat[i * (effnSite + 1)] = 2 * R_pow(sqrt(covMat[i * (effnSite + 1)]) * irange, *smooth);
      }
    }
  }
  
  else
    for (i=effnSite;i--;)
      covMat[i * (effnSite + 1)] = *sill + *nugget;

  return;
}
  
void direct(int *n, int *nSite, int *grid, int *covmod, double *coord, int *dim,
	    double *nugget, double *sill, double *range, double *smooth,
	    double *ans){

  int neffSite = *nSite, i, j, k;
  double *covmat, one = 1, zero = 0, *d, *u, *v;

  if (*grid)
    neffSite = R_pow_di(neffSite, *dim);

  covmat = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  d = (double *)R_alloc(neffSite, sizeof(double));
  u = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  v = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  
  buildcovmat(nSite, grid, covmod, coord, dim, nugget, sill, range,
	      smooth, covmat);

  /* Compute the singular value decomposition of the covariance
     matrix.

     This piece of code is strongly inspired from Lapack.c */
  double *xvals = (double *) R_alloc(neffSite * neffSite, sizeof(double));
  Memcpy(xvals, covmat, neffSite * neffSite);
  
  {
    int info = 0, lwork = -1;
    int *iwork = (int *) R_alloc(8 * neffSite, sizeof(int));
    double *work, tmp;
      
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgesdd)("A", &neffSite, &neffSite, xvals, &neffSite, d, u,
		     &neffSite, v, &neffSite, &tmp, &lwork, iwork, &info);
    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dgesdd");

    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));

    F77_CALL(dgesdd)("A", &neffSite, &neffSite, xvals, &neffSite, d, u,
		     &neffSite, v, &neffSite, work, &lwork, iwork, &info);
    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dgesdd");
  }

  /* Compute the square root of the covariance matrix */
  // a) First compute diag(sqrt(d)) %*% u
  for (i=0;i<neffSite;i++){
    double dummy = sqrt(d[i]);
    
    for (j=0;j<neffSite;j++)
      u[i + neffSite * j] *= dummy;
  }

  // b) Then compute v^T %*% diag(sqrt(d)) %*% u and put it in covmat
  F77_CALL(dgemm)("T", "N", &neffSite, &neffSite, &neffSite, &one,
		  v, &neffSite, u, &neffSite, &zero, covmat, &neffSite);
  
  /* Simulation part */
  GetRNGstate();

  if (*grid){
    for (i=0;i<*n;i++){
      for (j=0;j<neffSite;j++)
	d[j] = norm_rand();
      
      for (j=0;j<neffSite;j++){
	double sum = 0;
	for (k=0;k<neffSite;k++)
	  sum += d[k] * covmat[j + k * neffSite];
	
	ans[j + i * neffSite] = sum;
      }
    }
  }

  else {
    for (i=0;i<*n;i++){
      for (j=0;j<neffSite;j++)
	d[j] = norm_rand();
      
      for (j=0;j<neffSite;j++){
	double sum = 0;
	for (k=0;k<neffSite;k++)
	  sum += d[k] * covmat[j + k * neffSite];
	
	ans[i + j * *n] = sum;
      }
    }
  }
  
  PutRNGstate();

  return;
}

    
