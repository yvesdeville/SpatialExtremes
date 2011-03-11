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

  int i, neffSite, lagi = 1, lagj = 1;
  const double loguBound = log(*uBound), halfSigma2 = 0.5 * *sigma2;
  double sigma = sqrt(*sigma2), nugget = 1 - *sill,
    *lines = (double *)R_alloc(3 * *nlines, sizeof(double));

  if (*grid){
    neffSite = R_pow_di(*nSite, *dim);
    lagi = neffSite;
  }

  else{
    neffSite = *nSite;
    lagj = *nObs;
  }

  //rescale the coordinates
  for (i=(*nSite * *dim);i--;){
    const double irange = 1 / *range;
    coord[i] = coord[i] * irange;
  }
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);
  
  GetRNGstate();
 
  for (i=*nObs;i--;){
    int nKO = neffSite;
    double poisson = 0;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      int j;
      double *gp = (double *)R_alloc(neffSite, sizeof(double));
      
      /* ------- Random rotation of the lines ----------*/
      double u = unif_rand() - 0.5,
	v = unif_rand() - 0.5,
	w = unif_rand() - 0.5,
	angle = runif(0, M_2PI),	
	inorm = 1 / sqrt(u * u + v * v + w * w);
      
      u *= inorm;
      v *= inorm;
      w *= inorm;
      
      rotation(lines, nlines, &u, &v, &w, &angle);
      /* -------------- end of rotation ---------------*/
      
      poisson += exp_rand();
      double ipoisson = -log(poisson),
	thresh = loguBound + ipoisson;
      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      memset(gp, 0, neffSite * sizeof(double));
      tbmcore(nSite, &neffSite, dim, covmod, grid, coord, &nugget,
	      sill, range, smooth, nlines, lines, gp);
      
      nKO = neffSite;
      double ipoissonMinusHalfSigma2 = ipoisson - halfSigma2;
      for (j=neffSite;j--;){
	ans[j * lagj + i * lagi] = fmax2(sigma * gp[j] + ipoissonMinusHalfSigma2,
					 ans[j * lagj + i * lagi]);
	
	nKO -= (thresh <= ans[j * lagj + i * lagi]);
	
      }
    }
  }

  PutRNGstate();

  /* So fare we generate a max-stable process with standard Gumbel
     margins. Switch to unit Frechet ones */
  for (i=*nObs * neffSite;i--;)
    ans[i] = exp(ans[i]);

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

  int i, j, k, neffSite, lagi = 1, lagj = 1;
  const double loguBound = log(*uBound), halfSigma2 = 0.5 * *sigma2;
  double sigma = sqrt(*sigma2), nugget = 1 - *sill, one = 1, zero = 0;

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

  buildcovmat(nSite, grid, covmod, coord, dim, &nugget, sill, range,
	      smooth, covmat);
  
  /* Compute the singular value decomposition of the covariance
     matrix.

     This piece of code is strongly inspired from Lapack.c */
  
  Memcpy(xvals, covmat, neffSite * neffSite);
  
  {
    int *iwork= (int *) R_alloc(8 * neffSite, sizeof(int));
    double tmp;

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
  for (i=0;i<neffSite;i++){
    double dummy = sqrt(d[i]);
    
    for (j=0;j<neffSite;j++)
      u[i + neffSite * j] *= dummy;
  }

  // b) Then compute v^T %*% diag(sqrt(d)) %*% u and put it in covmat
  F77_CALL(dgemm)("T", "N", &neffSite, &neffSite, &neffSite, &one,
		  v, &neffSite, u, &neffSite, &zero, covmat, &neffSite);
  
  GetRNGstate();
  
  for (i=*nObs;i--;){
    double poisson = 0;
    int nKO = neffSite;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      poisson += exp_rand();
      double ipoisson = -log(poisson), thresh = loguBound + ipoisson;
	
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
      
      nKO = neffSite;
      double ipoissonMinusHalfSigma2 = ipoisson - halfSigma2;
      for (j=neffSite;j--;){
	ans[j * lagj + i * lagi] = fmax2(sigma * gp[j] + ipoissonMinusHalfSigma2,
				      ans[j * lagj + i * lagi]);
	
	nKO -= (thresh <= ans[j * lagj + i * lagi]);
	  
      }
    }
  }
  
  PutRNGstate();

  /* So fare we generate a max-stable process with standard Gumbel
     margins. Switch to unit Frechet ones */
  for (i=*nObs * neffSite;i--;)
    ans[i] = exp(ans[i]);

  return;
}

void rgeomcirc(int *nObs, int *ngrid, double *steps, int *dim,
	       int *covmod, double *sigma2, double *sill, double *range,
	       double *smooth, double *uBound, double *ans){
  /* This function generates random fields from the geometric model

     nObs: the number of observations to be generated
    ngrid: the number of locations along one axis
      dim: the random field is generated in R^dim
   covmod: the covariance model
     sill: the sill parameter - (1 - sill) is a nugget effect
    range: the range parameter
   smooth: the smooth parameter
   uBound: the uniform upper bound for the stoch. proc.
      ans: the generated random field */

  int i, j, k = -1, nbar = R_pow_di(*ngrid, *dim), r, m;
  const double loguBound = log(*uBound), halfSigma2 = 0.5 * *sigma2;
  double sigma = sqrt(*sigma2), nugget = 1 - *sill, *rho, *irho, *dist;

  //Below is a table of highly composite numbers
  int HCN[39] = {1, 2, 4, 6, 12, 24, 36, 48, 60, 120, 180, 240,
		 360, 720, 840, 1260, 1680, 2520, 5040, 7560,
		 10080, 15120, 20160, 25200, 27720, 45360, 50400,
		 55440, 83160, 110880, 166320, 221760, 277200,
		 332640, 498960, 554400, 665280, 720720, 1081080};

    
  /* Find the smallest size m for the circulant embedding matrix */
  {
    int dummy = 2 * (*ngrid - 1);
    do {
      k++;
      m = HCN[k];
    } while (m < dummy);
  }
  
  /* ---------- beginning of the embedding stage ---------- */
  int mbar = m * m, halfM = m / 2, notPosDef = 0;
  do {
    dist = (double *)R_alloc(mbar, sizeof(double));

    notPosDef = 0;
    //Computation of the distance
    for (r=mbar;r--;){
      i = r % m;
      j = r / m;
      
      if (i > halfM)
	i -= m;
      
      if (j > halfM)
	j -= m;
      
      dist[r] = hypot(steps[0] * i, steps[1] * j);
    }

    //Computations of the covariances
    rho = (double *)R_alloc(mbar, sizeof(double));
    irho = (double *)R_alloc(mbar, sizeof(double));
    memset(irho, 0, mbar * sizeof(double));

    switch (*covmod){
    case 1:
      whittleMatern(dist, mbar, *sill, *range, *smooth, rho);
      break;
    case 2:
      cauchy(dist, mbar, *sill, *range, *smooth, rho);
      break;
    case 3:
      powerExp(dist, mbar, *sill, *range, *smooth, rho);
      break;
    case 4:
      bessel(dist, mbar, *dim, *sill, *range, *smooth, rho);
      break;
    }

    /* Compute the eigen values to check if the circulant embbeding
       matrix is positive definite */

    /* Note : The next lines is only valid for 2d random fields. I
       need to change if there are m_1 \neq m_2 as I suppose that m_1
       = m_2 = m */
    int maxf, maxp;

    fft_factor(m, &maxf, &maxp);
    double *work = (double *)R_alloc(4 * maxf, sizeof(double));
    int *iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, m, m, 1, -1, work, iwork);

    fft_factor(m, &maxf, &maxp);
    work = (double *)R_alloc(4 * maxf, sizeof(double));
    iwork = (int *)R_alloc(maxp, sizeof(int));
    fft_work(rho, irho, 1, m, m, -1, work, iwork);

    //Check if the eigenvalues are all positive
    for (i=mbar;i--;){
      notPosDef |= (rho[i] <= 0) || (fabs(irho[i]) > 0.001);
    }

    if (notPosDef){
      k++;
      m = HCN[k];
      halfM = m / 2;
      mbar = m * m;
    }

    if (k > 30)
      error("Impossible to embbed the covariance matrix");
    
  } while (notPosDef);
  /* --------- end of the embedding stage --------- */

  /* Computation of the square root of the eigenvalues */
  for (i=mbar;i--;){
    rho[i] = sqrt(rho[i]);
    irho[i] = 0;//No imaginary part
  }

  int mdag = m / 2 + 1, mdagbar = mdag * mdag;
  double isqrtMbar = 1 / sqrt(mbar);

  double *a = (double *)R_alloc(mbar, sizeof(double));
  double *ia = (double *)R_alloc(mbar, sizeof(double));
  
  GetRNGstate();
  for (i=*nObs;i--;){
    int nKO = nbar;
    double poisson = 0;
    
    while (nKO) {
      /* The stopping rule is reached when nKO = 0 i.e. when each site
	 satisfies the condition in Eq. (8) of Schlather (2002) */
      int j;
      double *gp = (double *)R_alloc(nbar, sizeof(double));
      
      poisson += exp_rand();
      double ipoisson = -log(poisson), thresh = loguBound + ipoisson;
      
      /* We simulate one realisation of a gaussian random field with
	 the required covariance function */
      circcore(rho, a, ia, m, halfM, mdag, mdagbar, *ngrid, nbar, isqrtMbar, nugget, gp);
      
      nKO = nbar;
      double ipoissonMinusHalfSigma2 = ipoisson - halfSigma2;
      for (j=nbar;j--;){
	ans[j + i * nbar] = fmax2(sigma * gp[j] + ipoissonMinusHalfSigma2,
				  ans[j + i * nbar]);
	nKO -= (thresh <= ans[j + i * nbar]);
	
      }
    }
  }
  
  PutRNGstate();

  /* So fare we generate a max-stable process with standard Gumbel
     margins. Switch to unit Frechet ones */
  for (i=*nObs * nbar;i--;)
    ans[i] = exp(ans[i]);
  
  return;
}
