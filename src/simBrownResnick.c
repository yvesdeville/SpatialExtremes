#include "header.h"

void rbrowndirect(double *coord, double *bounds, int *nObs, int *nSite,
		  int *dim, int *grid, double *range, double *smooth,
		  double *ans){
  /* This function generates random fields for the geometric model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
     nSite: the number of locations
       dim: the random field is generated in R^dim
      grid: Does coord specifies a grid?
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field */

  int i, j, k, neffSite, covmod = 6, lagi = 1, lagj = 1, nSimTot = 250;
  const double irange = 1 / *range, inflate = R_pow(4, 1 / *smooth) * *range;
  double one = 1, zero = 0, tmp;

  //Inflate the region
  for (i=*dim;i--;){
    bounds[2 * i] -= inflate;
    bounds[2 * i + 1] += inflate;
  }

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  double *covmat = (double *)R_alloc(neffSite * neffSite, sizeof(double)),
    *d = (double *)R_alloc(neffSite, sizeof(double)),
    *u = (double *)R_alloc(neffSite * neffSite, sizeof(double)),
    *v = (double *)R_alloc(neffSite * neffSite, sizeof(double)),
    *xvals = (double *) R_alloc(neffSite * neffSite, sizeof(double)),
    *gp = (double *)R_alloc(neffSite, sizeof(double));

  GetRNGstate();

  for (i=*nObs;i--;){
    double poisson = 0;
    //int nKO = neffSite;
    
    //while (nKO){
    for (int nSim=nSimTot; nSim--;){
      double *shift = (double *) R_alloc(*dim, sizeof(double));
      double *shiftedCoord = (double *) R_alloc(*dim * *nSite, sizeof(double));
      double *vario = (double *) R_alloc(neffSite, sizeof(double));

      // Shift the locations
      for (j=*dim;j--;)
	shift[j] = runif(bounds[2 * j], bounds[2 * j + 1]);
      
      for (j=*nSite;j--;)
	for (k=*dim;k--;)
	  shiftedCoord[k * *nSite + j] = coord[k * *nSite + j] - shift[k];
      
      // Compute the variogram
      distance2orig(shiftedCoord, *nSite, *dim, vario, *grid);
      
      for (j=neffSite;j--;)
	vario[j] = R_pow(vario[j] * irange, *smooth);
      
      // Compute the covariance matrix for the shifted locations
      buildcovmat(nSite, grid, &covmod, shiftedCoord, dim, &zero, &one, range,
		  smooth, covmat);
      
      /* Compute the singular value decomposition of the covariance
	 matrix.
	 
	 This piece of code is strongly inspired from Lapack.c */
      
      Memcpy(xvals, covmat, neffSite * neffSite);

      {
	int *iwork= (int *) R_alloc(8 * neffSite, sizeof(int));
	
	/* ask for optimal size of work array */
	int lwork = -1, info = 0;
	F77_CALL(dgesdd)("A", &neffSite, &neffSite, xvals, &neffSite, d, u,
			 &neffSite, v, &neffSite, &tmp, &lwork, iwork, &info);
	if (info != 0)
	  error("error code %d from Lapack routine '%s'", info, "dgesdd");
	
	lwork = (int) tmp;
	double *work = (double *) R_alloc(lwork, sizeof(double));
	
	F77_CALL(dgesdd)("A", &neffSite, &neffSite, xvals, &neffSite, d, u,
			   &neffSite, v, &neffSite, work, &lwork, iwork, &info);
	if (info != 0)
	  error("error code %d from Lapack routine '%s'", info, "dgesdd");
      }
      
      /*--------------- end of singular value decomposition ---------------*/
      
      /* Compute the square root of the covariance matrix */
      // a) First compute diag(sqrt(d)) %*% u
      for (j=neffSite;j--;){
	double dummy = sqrt(d[j]);
	
	for (k=neffSite;k--;)
	  u[j + neffSite * k] *= dummy;
      }
      
      // b) Then compute v^T %*% diag(sqrt(d)) %*% u and put it in covmat
      F77_CALL(dgemm)("T", "N", &neffSite, &neffSite, &neffSite, &one,
		      v, &neffSite, u, &neffSite, &zero, covmat, &neffSite);
      
      poisson += exp_rand();
      double ipoisson = -log(poisson);//, thresh = 4 + ipoisson;
	
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      for (j=neffSite;j--;)
	d[j] = norm_rand();
	
      for (j=neffSite;j--;){
	double sum = 0;
	for (k=neffSite;k--;)
	  sum += d[k] * covmat[j + k * neffSite];
	  
	gp[j] = sum;
      }
	
      //nKO = neffSite;
      for (j=neffSite;j--;){
	ans[j * lagj + i * lagi] = fmax2(gp[j] - vario[j] + ipoisson,
					 ans[j * lagj + i * lagi]);
	//nKO -= (thresh <= ans[j * lagj +  i * lagi]);
      }
    }
  }

  PutRNGstate();

  for (i=*nObs * neffSite;i--;)
    ans[i] = M_SQRT2 * exp(ans[i]);

  return;
}
