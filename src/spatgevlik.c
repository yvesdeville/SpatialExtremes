#include "header.h"

void spatgevlik(double *data, double *covariables, int *nSite, int *nObs,
		double *locdsgnmat, double *locpenmat, int *nloccoeff,
		int *npparloc, double *locpenalty, double *scaledsgnmat,
		double *scalepenmat, int *nscalecoeff, int *npparscale,
		double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		int *nshapecoeff, int *npparshape, double *shapepenalty,
		double *loccoeff, double *scalecoeff, double *shapecoeff,
		double *dns){
  
  /* This function computes the log-likelihood for a spatial GEV model
     i.e. the GEV paramters are defined through a response surface
     that could be either is simple linear model or a p-spline. */

  int i, j;
  double *locs, *scales, *shapes;

  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));

  //Stage 1: Computing the GEV parameters using the design matrix
  *dns = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
		      loccoeff, scalecoeff, shapecoeff, *nSite,
		      *nloccoeff, *nscalecoeff, *nshapecoeff,
		      locs, scales, shapes);

  if (*dns != 0.0){
    *dns = MINF;
    return;
  }

  //Stage 2: Compute the log-likelihood (assuming independence between
  //stations)
  
  for (i=0;i<*nSite;i++){

    if (fabs(shapes[i]) <= 1e-6){
      for (j=0;j<*nObs;j++){
	data[j + i * *nObs] = (data[j + i * *nObs] - locs[i]) /
	  scales[i];
	
	*dns += -log(scales[i]) - data[j + i * *nObs] -
	  exp(-data[j + i * *nObs]);
      }
    }
    
    else{
      for (j=0;j<*nObs;j++){
	data[j + i * *nObs] = 1 + shapes[i] * (data[j + i * *nObs] - locs[i]) /
	  scales[i];
	
	if (data[j + i * *nObs] <= 0){
	  *dns = MINF;
	  return;
	}
	
	*dns += -log(scales[i]) - R_pow(data[j + i * *nObs], -1/shapes[i]) -
	  (1/shapes[i] + 1) * log(data[j + i * *nObs]);
      }
    }
  }

  //Stage 3: Penalizing the llik if p-splines are used
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

  return;
}

	  
