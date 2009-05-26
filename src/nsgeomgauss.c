#include "header.h"

void nsgeomgaussfull(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, int *dim, double *locs, double *scales, double *shapes,
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

  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (i=0;i<*nSite;i++){
      if ((scales[i] <= 0) || (shapes[i] <= -1)){
	*dns = MINF;
	return;
      }
    }
  }

  //Stage 0: Compute the sigma2 at each location
  dsgnmat2Sigma2(sigma2dsgnmat, sigma2coeff, *nSite,
		 *nsigma2coeff, sigma2);
   
  //Stage 1: Compute the covariance at each location
  *dns = nsgeomCovariance(dist, *nSite, *dim, *covmod, sigma2, *sill, *range,
			  *smooth, rho);

  if (*dns != 0.0)
    return;
    
  //Stage 2: Transformation to unit Frechet
  if (*fitmarge){
    *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);

    if (*dns != 0.0)
      return;
    
    *dns = lpliksmith(frech, rho, jac, *nObs, *nSite);
  }
    
  else {
    for (i=0;i<(*nSite * *nObs);i++)
      jac[i] = 0.0;
    
    *dns = lpliksmith(data, rho, jac, *nObs, *nSite);
  }  

  if (!R_FINITE(*dns))
    *dns = MINF;

  return;

}

void nsgeomgaussdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			int *dim, double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
			double *locpenalty, double *scaledsgnmat, double *scalepenmat,
			int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
			double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
			double *sigma2dsgnmat, int *nsigma2coeff, double *loccoeff, double *scalecoeff,
			double *shapecoeff, double *sigma2coeff, double *sill, double *range,
			double *smooth, double *dns){
  //This is the geometric gaussian model
  //The GEV parameters are defined using a polynomial response surface
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  double *jac, *rho, *locs, *scales, *shapes, *frech, *sigma2;
    
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  sigma2 = (double *)R_alloc(*nSite, sizeof(double));
  
  //Stage 0: Compute the sigma2 at each location
  dsgnmat2Sigma2(sigma2dsgnmat, sigma2coeff, *nSite,
		 *nsigma2coeff, sigma2);
  
  //Stage 1: Compute the covariance at each location
  *dns = nsgeomCovariance(dist, *nSite, *dim, *covmod, sigma2, *sill, *range,
			  *smooth, rho);

  if (*dns != 0.0)
    return;
    
  //Stage 2: Compute the GEV parameters using the design matrix
  *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
		       loccoeff, scalecoeff, shapecoeff, *nSite,
		       *nloccoeff, *nscalecoeff, *nshapecoeff,
		       locs, scales, shapes);

  if (*dns != 0.0)
    return;

  //Stage 3: Transformation to unit Frechet
  *dns = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		   jac, frech);

  if (*dns != 0.0)
    return;
    
  *dns = lpliksmith(frech, rho, jac, *nObs, *nSite);
  
  //Stage 5: Removing the penalizing terms (if any)
  // 1- For the location parameter
  if (*locpenalty > 0)
    *dns -= penalization(locpenmat, loccoeff, *locpenalty,
			 *nloccoeff, *npparloc);
  
  // 2- For the scale parameter
  if (*scalepenalty > 0)    
    *dns -= penalization(scalepenmat, scalecoeff, *scalepenalty,
			 *nscalecoeff, *npparscale);
  
  // 3- For the shape parameter
  if (*shapepenalty > 0)
    *dns -= penalization(shapepenmat, shapecoeff, *shapepenalty,
			 *nshapecoeff, *npparshape);
  
  if (!R_FINITE(*dns))
    *dns = MINF;

  return;
  
}
