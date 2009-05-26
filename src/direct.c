#include "header.h"

void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
		 double *nugget, double *sill, double *range,
		 double *smooth, double *covMat){

  int i, j, currentPair, nPairs, effnSite, zero = 0;
  double *dist, *rho, flag = 0;

  if (*grid)
    effnSite = R_pow_di(*nSite, *dim);
    
  else
    effnSite = *nSite;

  nPairs = effnSite * (effnSite - 1) / 2;

  dist = (double *)R_alloc(nPairs, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));

  if (*grid){
    //Coord specify a grid
    double *coordGrid;
    coordGrid = (double *)R_alloc(effnSite * *dim, sizeof(double));

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
    flag = whittleMatern(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 2:
    flag = cauchy(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 3:
    flag = powerExp(dist, nPairs, *sill, *range, *smooth, rho);
    break;
  case 4:
    flag = bessel(dist, nPairs, *dim, *sill, *range, *smooth, rho);
    break;
  }
  
  if (flag != 0.0)
    error("The covariance parameters seem to be ill-defined. Please check\n");
  
  currentPair = -1;
  for (i=0;i<(effnSite-1);i++){
    for (j=i+1;j<effnSite;j++){
      currentPair++;
      covMat[effnSite * i + j] = rho[currentPair];
      covMat[effnSite * j + i] = rho[currentPair];
    }
  }
  
  //Fill the diagonal of sill + nugget
  for (i=effnSite;i--;)
    covMat[i * (effnSite + 1)] = *sill + *nugget;

  return;
}
  
void direct(int *n, int *nSite, int *grid, int *covmod, double *coord, int *dim,
	    double *nugget, double *sill, double *range, double *smooth,
	    double *ans){

  int neffSite, i, j, k, lwork, info = 0;
  double *covmat, one = 1, zero = 0, *work, *xvals,
    tmp, *d, *u, *v, sum, dummy;

  if (*grid)
    neffSite = R_pow_di(*nSite, *dim);

  else
    neffSite = *nSite;

  covmat = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  d = (double *)R_alloc(neffSite, sizeof(double));
  u = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  v = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  xvals = (double *) R_alloc(neffSite * neffSite, sizeof(double));
  
  buildcovmat(nSite, grid, covmod, coord, dim, nugget, sill, range,
	      smooth, covmat);

  /* Compute the singular value decomposition of the covariance
     matrix.

     This piece of code is strongly inspired from Lapack.c */
  
  Memcpy(xvals, covmat, neffSite * neffSite);
  
  {
    int *iwork= (int *) R_alloc(8 * *nSite, sizeof(int));
    
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
    dummy = sqrt(d[i]);
    
    for (j=0;j<neffSite;j++)
      u[i + neffSite * j] *= dummy;
  }

  // b) Then compute v^T %*% diag(sqrt(d)) %*% u and put it in covmat
  F77_CALL(dgemm)("T", "N", &neffSite, &neffSite, &neffSite, &one,
		  v, &neffSite, u, &neffSite, &zero, covmat, &neffSite);
  
  /* Simulation part */
  GetRNGstate();
  
  for (i=0;i<*n;i++){
    for (j=0;j<neffSite;j++)
      d[j] = norm_rand();

    for (j=0;j<neffSite;j++){
      sum = 0;
      for (k=0;k<neffSite;k++)
	sum += d[k] * covmat[j + k * neffSite];
      
      ans[i + j * *n] = sum;
    }
  }

  PutRNGstate();

  return;
}

    
