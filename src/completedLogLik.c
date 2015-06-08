#include "header.h"

void completellik(int *nObs, int *nSite, double *data, int *parts, int *partSizes,
		  double *cov, double *weights){
  /* This function computes the log density (Schlather case) */
  *weights=0;
  double *y = malloc(*nSite * sizeof(double));
  int *part = malloc(*nSite * sizeof(int));

  for (int i=0;i<*nObs;i++){

    memcpy(part, parts + i * *nSite, *nSite * sizeof(int));

    for (int j=0;j<*nSite;j++)
      y[j] = data[i + *nObs * j];

    for (int j=0;j<partSizes[i];j++){
      //Loop over all sets of the current partition

      int r = 0;
      for (int k=0;k<*nSite;k++)
	r += (part[k] == j);

      int nr = *nSite - r,
	*tau = malloc(r * sizeof(int)),
	*taubar = malloc(nr * sizeof(int));

      gettau(nSite, part, &j, tau);

      if (r < *nSite) {
	gettaubar(nSite, part, &j, taubar);

	/* Define the parameters of the distribution of the extremal
	   function. For the Schlather case, it is a multivariate
	   Student distribution. Hence below we compute the degree of
	   freedom (DoF), the Cholesky decomposition of the scale
	   matrix (scaleMat) and the location parameter (mu) */
	double DoF = 0,
	    *scaleMat = malloc(nr * nr * sizeof(double)),
	    *mu = malloc(nr * sizeof(double));
	  getParametersSC(tau, taubar, &r, &nr, cov, y, &DoF, mu, scaleMat);

	  /* Compute the proba values, i.e., the probabilities that one
	     extremal function is smaller than the hitting bounds at the
	     remaining conditioning locations. */
	  double prob = 0;
	  computeprobaSC(&DoF, mu, scaleMat, y, &nr, taubar, &prob);
	  *weights += log(prob);

	  free(scaleMat); free(mu);
      }

      /* Compute the intensity function values, i.e.,
	 \lambda_{x_{\tau_j}}(z_{\tau_j}). */
      double f = 0;
      getfvaluesSC(y, nSite, &r, tau, cov, &f);
      *weights += f;
      free(tau); free(taubar);
    }

    /* Add the contribution of the sub-extremal functions */
    double subExtCont=0;
    int nMC = 1000;
    pschlather(y, nSite, cov, &subExtCont, &nMC);
    *weights += subExtCont;
  }

  free(y); free(part);
  return;
}

void llikschlather(int *nObs, int *nSite, double *data, int *nPart,
		   int *parts, int *partSizes, double *cov, double *weights){
  /* This function computes the log likelihood (Schlather case) */
  *weights=1;
  double *y = malloc(*nSite * sizeof(double));
  int *part = malloc(*nSite * sizeof(int));

  for (int i=0;i<*nObs;i++){
    for (int j=0;j<*nSite;j++)
      y[j] = data[i + *nObs * j];

    double contrib = 0;
    for (int j=0;j<*nPart;j++){
      double tmp=1;
      memcpy(part, parts + j * *nSite, *nSite * sizeof(int));

      for (int k=0;k<partSizes[j];k++){
	//Loop over all sets of the current partition

	int r = 0;
	for (int l=0;l<*nSite;l++)
	  r += (part[l] == k);

	int nr = *nSite - r,
	  *tau = malloc(r * sizeof(int)),
	  *taubar = malloc(nr * sizeof(int));

	gettau(nSite, part, &k, tau);

	if (r < *nSite) {
	  gettaubar(nSite, part, &k, taubar);

	  /* Define the parameters of the distribution of the extremal
	     function. For the Schlather case, it is a multivariate
	     Student distribution. Hence below we compute the degree
	     of freedom (DoF), the Cholesky decomposition of the scale
	     matrix (scaleMat) and the location parameter (mu) */
	  double DoF = 0,
	    *scaleMat = malloc(nr * nr * sizeof(double)),
	    *mu = malloc(nr * sizeof(double));
	  getParametersSC(tau, taubar, &r, &nr, cov, y, &DoF, mu, scaleMat);

	  /* Compute the proba values, i.e., the probabilities that one
	     extremal function is smaller than the hitting bounds at the
	     remaining conditioning locations. */
	  double prob = 0;
	  computeprobaSC(&DoF, mu, scaleMat, y, &nr, taubar, &prob);
	  tmp *= prob;

	  free(scaleMat); free(mu);
	}

	/* Compute the intensity function values, i.e.,
	   \lambda_{x_{\tau_j}}(z_{\tau_j}). */
	double f = 0;
	getfvaluesSC(y, nSite, &r, tau, cov, &f);
	tmp *= exp(f);
	free(tau); free(taubar);
      }

      contrib += tmp;
    }

    /* Add the contribution of the sub-extremal functions */
    double subExtCont=0;
    int nMC = 1000;
    pschlather(y, nSite, cov, &subExtCont, &nMC);
    *weights *= exp(subExtCont) * contrib;
  }

  free(y); free(part);
  return;
}


void pschlather(double *q, int *dim, double *cov, double *prob, int *nMC){

  /* Compute log Pr{Z(x_1) < z_1, ..., Z(x_k) < z_k} = -E{max_j
     Y(x_j) / z_j} by Monte Carlo */
  *prob=0;
  int one = 1;
  double *iq = malloc(*dim * sizeof(double));
  for (int i=0;i<*dim;i++)
    iq[i] = 1 / q[i];


  /* Compute the Cholesky decomposition of the covariance matrix */
  double *covChol = malloc(*dim * *dim * sizeof(double));
  memcpy(covChol, cov, *dim * *dim * sizeof(double));

  int info = 0;
  F77_CALL(dpotrf)("U", dim, covChol, dim, &info);

  if (info != 0)
    error("error code %d from Lapack routine '%s'", info, "dpotrf");

  GetRNGstate();
  double max, min,
    *sim = malloc(*dim * sizeof(double));

  /* We will use an antithetic approach since Z = -Z in dist where Z
     \sim N(0, 1). Now since max (covChol^T (-Z)) = - min covChol^T Z
     this doesn't actually require any additional computations---apart
     from the minimum... */

  for (int i=0;i<*nMC;i++){
    for (int j=0;j<*dim;j++)
      sim[j] = norm_rand();

    F77_CALL(dtrmv)("U", "T", "N", dim, covChol, dim, sim, &one);

    max = 0;
    min = 0;
    for (int j=0;j<*dim;j++){
      double tmp = sim[j] * iq[j];
      max = fmax2(max, tmp);
      min = fmin2(min, tmp);
    }

    *prob += max - min;
  }

  free(sim);free(iq);free(covChol);
  PutRNGstate();

  *prob *= -0.5 * M_SQRT2 * M_SQRT_PI / (double) *nMC;
  return;
}
