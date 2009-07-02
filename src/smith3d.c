#include "header.h"

void smithfull3d(double *data, double *distVec, int *nSite,
		 int *nObs, double *locs, double *scales, double *shapes,
		 double *cov11, double *cov12, double *cov13, double *cov22,
		 double *cov23, double *cov33, int *fitmarge, double *dns){
  //This is the Smith model - 3d case. It's a wrapper to several
  //sub-functions. It's named xxxfull as it either assume that the
  //margins are unit Frechet, or the GEV parameters are estimated at
  //each locations.

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i;
  double *jac, *mahalDist, *frech;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));

  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (i=0;i<*nSite;i++){
      if ((scales[i] <= 0) | (shapes[i] <= -1)){
	*dns = MINF;
	return;
      }
    }
  }

  //Stage 1: Computing the Mahalanobis distance
  *dns = mahalDistFct3d(distVec, nPairs, cov11, cov12, cov13,
			cov22, cov23, cov33, mahalDist);

  if (*dns != 0.0)
    return;
  
  //Stage 2: Transformation to unit Frechet
  if (*fitmarge){
    *dns = gev2frech(data, *nObs, *nSite, locs, scales,
		     shapes, jac, frech);

    if (*dns != 0.0)
      return;
    
    *dns = lpliksmith(frech, mahalDist, jac, *nObs, *nSite);
  }
  
  else {
    for (i=0;i<(*nSite * *nObs);i++)
      jac[i] = 0.0;    

    *dns = lpliksmith(data, mahalDist, jac, *nObs, *nSite);
  }

  if (!R_FINITE(*dns))
    *dns = MINF;
  
  return;
}

void smithdsgnmat3d(double *data, double *distVec, int *nSite, int *nObs, 
		    double *locdsgnmat, double *locpenmat, int *nloccoeff,
		    int *npparloc, double *locpenalty, double *scaledsgnmat,
		    double *scalepenmat, int *nscalecoeff, int *npparscale,
		    double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		    int *nshapecoeff, int *npparshape, double *shapepenalty,
		    double *loccoeff, double *scalecoeff, double *shapecoeff,
		    double *cov11, double *cov12, double *cov13, double *cov22,
		    double *cov23, double *cov33, double *dns){
  //This is the Smith's model - 3d case. It's named xxxdsgnmat as
  //either linear models or p-splines are used for the gev parameters.
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  double *jac, *mahalDist, *locs, *scales, *shapes, *frech;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  
  //Stage 1: Computing the Mahalanobis distance
  *dns = mahalDistFct3d(distVec, nPairs, cov11, cov12, cov13,
			cov22, cov23, cov33, mahalDist);

  if (*dns != 0.0)
    return;

  //Stage 2: Computing the GEV parameters using the design matrix
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

  *dns = lpliksmith(frech, mahalDist, jac, *nObs, *nSite);
    
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
