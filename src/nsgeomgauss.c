#include "header.h"

void nsgeomgaussfull(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, double *locs, double *scales, double *shapes,
		     double *sigma2dsgnmat, double *sigma2coeff, int *nsigma2coeff,
		     double *sill, double *range, double *smooth, int *fitmarge,
		     double *dns){
  //This is the non-stationary geometric gaussian model. It computes
  //the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i;
  double *jac, *rho, *frech, *sigma2;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  sigma2 = (double *)R_alloc(*nSite, sizeof(double));

  *dns = 1.0;

  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (i=0;i<*nSite;i++){
      if (scales[i] <= 0){
	//printf("scales <= 0!!!\n");
	*dns += R_pow_di(1 - scales[i], 2);
	scales[i] = 1e-3;
      }
      
      if (shapes[i] <= -1){
	//printf("shapes <= -1!!!\n");
	*dns += R_pow_di(shapes[i], 2);
	shapes[i] = -0.9;
      }
    }
  }

  //Stage 0: Compute the sigma2 at each location
  dsgnmat2Sigma2(sigma2dsgnmat, sigma2coeff, *nSite,
		 *nsigma2coeff, sigma2);
   
  //Stage 1: Compute the covariance at each location
  *dns += nsgeomCovariance(dist, *nSite, *covmod, sigma2, *sill, *range,
			   *smooth, rho);
    
  //Stage 2: Transformation to unit Frechet
  if (*fitmarge)
    *dns += gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		      jac, frech);
    
  else {
    for (i=0;i<(*nSite * *nObs);i++){
      frech[i] = data[i];
      jac[i] = 0.0;
    }
  }

  *dns *= lpliksmith(frech, rho, jac, *nObs, *nSite);

  return;

}
