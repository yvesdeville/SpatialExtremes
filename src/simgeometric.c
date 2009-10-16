#include "header.h"

void rgeomtbm(double *coord, int *nObs, int *nSite, int *dim,
	      int *covmod, int *grid, double *sigma2, double *sill,
	      double *range, double *smooth, double *uBound,
	      int *nlines, double *ans){
  /* This function generates random fields from the geometric model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
    sigma2: the variance of the geometric gaussian process
      sill: the sill parameter - (1 - sill) is a nugget effect
     range: the range parameter
    smooth: the smooth parameter
    uBound: the uniform upper bound for the stoch. proc.
    nlines: the number of lines used in the TBM algo
       ans: the generated random field */

  int i;
  double *lines;

  //rescale the coordinates
  for (i=(*nSite * *dim);i--;){
    const double irange = 1 / *range;
    coord[i] = coord[i] * irange;
  }

  lines = (double *)R_alloc(3 * *nlines, sizeof(double));
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);
  
  GetRNGstate();
  if (*grid){
    //coord defines a grid
    int neffSite = R_pow_di(*nSite, *dim);
    for (i=*nObs;i--;){
      int nKO = neffSite;
      double poisson = 0;

      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	int j;
	double sigma = sqrt(*sigma2),// uBound = exp(3.5 * sigma - 0.5 * *sigma2),
	  nugget = 1 - *sill, ipoisson, u, v, w, angle, norm, thresh,
	  *gp;

	gp = (double *)R_alloc(neffSite, sizeof(double));

	/* ------- Random rotation of the lines ----------*/
	u = unif_rand() - 0.5;
	v = unif_rand() - 0.5;
	w = unif_rand() - 0.5;
	angle = runif(0, M_2PI);
	
	norm = sqrt(u * u + v * v + w * w);
	
	u /= norm;
	v /= norm;
	w /= norm;
	
	rotation(lines, nlines, &u, &v, &w, &angle);
	/* -------------- end of rotation ---------------*/
	
	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	
	for (j=neffSite;j--;)
	  gp[j] = 0;

	tbmcore(nSite, &neffSite, dim, covmod, grid, coord, &nugget,
		sill, range, smooth, nlines, lines, gp);
	
	nKO = neffSite;
	for (j=neffSite;j--;){
	  if (thresh > ans[j + i * neffSite])
	    ans[j + i * neffSite] = fmax2(exp(sigma * gp[j] - 0.5 * *sigma2) *
					  ipoisson, ans[j + i * neffSite]);
	  
	  else
	    nKO--;
	  
	}
      }
    }
  }

  else{
    //coord doesn't define a grid
    int neffSite = *nSite;
    for (i=*nObs;i--;){
      double poisson = 0;
      int nKO = neffSite;
      
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	int j;
	double sigma = sqrt(*sigma2),// uBound = exp(3.5 * sigma - 0.5 * *sigma2),
	  nugget = 1 - *sill, ipoisson, u, v, w, angle, norm, thresh, *gp;

	gp = (double *)R_alloc(neffSite, sizeof(double));

	/* ------- Random rotation of the lines ----------*/
	u = unif_rand() - 0.5;
	v = unif_rand() - 0.5;
	w = unif_rand() - 0.5;
	angle = runif(0, M_2PI);
	
	norm = sqrt(u * u + v * v + w * w);
	
	u /= norm;
	v /= norm;
	w /= norm;
	
	rotation(lines, nlines, &u, &v, &w, &angle);
	/* -------------- end of rotation ---------------*/

	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	for (j=neffSite;j--;)
	  gp[j] = 0;
	
	tbmcore(nSite, &neffSite, dim, covmod, grid, coord, &nugget,
		sill, range, smooth, nlines, lines, gp);
	
	nKO = neffSite;
	for (j=*nSite;j--;){
	  if (thresh > ans[i + j * *nObs])
	    ans[i + j * *nObs] = fmax2(exp(sigma * gp[j] - 0.5 * *sigma2) *
				       ipoisson, ans[i + j * *nObs]);
	  
	  else
	    nKO--;
	}
      }
    }
  }

  PutRNGstate();

  return;
}

void rgeomdirect(double *coord, int *nObs, int *nSite, int *dim,
		 int *covmod, int *grid, double *sigma2, double *sill,
		 double *range, double *smooth, double *uBound,
		 double *ans){
  /* This function generates random fields for the geometric model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
    sigma2: the variance of the geometric gaussian process
      sill: the sill parameter - (1 - sill) is a nugget effect
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field */

  const double sigma = sqrt(*sigma2);//, uBound = exp(3.5 * sigma - 0.5 * *sigma2);
  int i, j, k, lwork, nKO, info = 0, neffSite;
  double poisson, ipoisson, thresh, *gp, nugget = 1 - *sill,
    *covmat, one = 1, zero = 0, *work, *xvals, tmp, *d, *u,
    *v, sum, dummy;

  if (*grid)
    neffSite = R_pow_di(*nSite, *dim);

  else
    neffSite = *nSite;

  covmat = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  d = (double *)R_alloc(neffSite, sizeof(double));
  u = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  v = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  xvals = (double *) R_alloc(neffSite * neffSite, sizeof(double));
  gp = (double *)R_alloc(neffSite, sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, &nugget, sill, range,
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

  /*--------------- end of singular value decomposition ---------------*/

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
  
  GetRNGstate();
  if (*grid){
    //coord defines a grid
    for (i=*nObs;i--;){
      poisson = 0;
      nKO = neffSite;
      
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	for (j=neffSite;j--;)
	  d[j] = norm_rand();
	  
	for (j=neffSite;j--;){
	  sum = 0;
	  for (k=neffSite;k--;)
	    sum += d[k] * covmat[j + k * neffSite];
	  
	  gp[j] = sum;
	}
	
	nKO = neffSite;
	for (j=neffSite;j--;){
	  if (thresh > ans[j + i * neffSite])
	    ans[j + i * neffSite] = fmax2(exp(sigma * gp[j] - 0.5 * *sigma2) *
					  ipoisson, ans[j + i * neffSite]);
	  
	  else
	    nKO--;
	  
	}
      }
    }
  }

  else{
    //coord doesn't define a grid
    for (i=*nObs;i--;){
      poisson = 0;
      nKO = *nSite;
      
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	for (j=neffSite;j--;)
	  d[j] = norm_rand();
	  
	for (j=neffSite;j--;){
	  sum = 0;
	  for (k=neffSite;k--;)
	    sum += d[k] * covmat[j + k * neffSite];
	  
	  gp[j] = sum;
	}
	
	nKO = *nSite;
	for (j=*nSite;j--;){
	  if (thresh > ans[i + j * *nObs])
	    ans[i + j * *nObs] = fmax2(exp(sigma * gp[j] - 0.5 * *sigma2) *
				       ipoisson, ans[i + j * *nObs]);
	  
	  else
	    nKO--;
	}
      }
    }
  }

  PutRNGstate();

  return;
}
