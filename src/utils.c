#include "header.h"

void distance(double *coord, int *nDim, int *nSite,
	      int *vec, double *dist){

  //This function computes either the euclidean distance or the
  //distance vector between each pair of locations
  const int nPair = *nSite * (*nSite - 1) / 2;
  int i, j, k, currentPair = 0;

  if (*vec){
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	for (k=0;k<*nDim;k++)
	  dist[k * nPair + currentPair] = coord[k * *nSite + j] -
	    coord[k * *nSite + i];

	currentPair++;
      }
    }
  }

  else{
    memset(dist, 0, nPair * sizeof(double));

    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	for (k=0;k<*nDim;k++)
	  dist[currentPair] += (coord[i + k * *nSite] -	coord[j + k * *nSite]) *
	    (coord[i + k * *nSite] - coord[j + k * *nSite]);

	dist[currentPair] = sqrt(dist[currentPair]);
	currentPair++;
      }
    }
  }
}

void distance2orig(double *coord, int n, int dim, double *dist, int grid){
  //It computes the l_2 norm of points in R^d i.e. sqrt(x[1]^2 + ... + x[d]^2)
  int i, j;

  if (grid){
    //Only works with two dimensional grids!!!
    int current = -1;
    double dummy;
    memset(dist, 0, R_pow_di(n, dim) * sizeof(double));

    for (i=0;i<n;i++){
      dummy = coord[i] * coord[i];

      for (j=0;j<n;j++){
	current++;
	dist[current] = sqrt(dummy + coord[n + j] * coord[n + j]);
      }
    }
  }

  else {
    memset(dist, 0, n * sizeof(double));

    for (i=n;i--;){
      for (j=dim;j--;)
	dist[i] += coord[i + j * n] * coord[i + j * n];

      dist[i] = sqrt(dist[i]);
    }
  }
}

double gev2frech(double *data, int nObs, int nSite, double *locs,
		 double *scales, double *shapes, double *jac,
		 double *frech){

  //This function transforms the GEV observations to unit Frechet ones
  //and computes the log of the jacobian of each transformation
  //When ans > 0.0, the GEV parameters are invalid.

  int i, j;

  for (i=nSite;i--;){
    double iscale = 1 / scales[i], logScale = log(scales[i]);

    if (shapes[i] == 0.0){
      for (j=nObs;j--;){
	frech[i * nObs + j] = (data[i * nObs + j] - locs[i]) * iscale;
	jac[i * nObs + j] = frech[i * nObs + j] - logScale;
	frech[i * nObs + j] = exp(frech[i * nObs + j]);
      }
    }

    else {
      double ishape = 1 / shapes[i];
      for (j=nObs;j--;){
	frech[i * nObs + j] = 1 + shapes[i] * (data[i * nObs + j] - locs[i]) * iscale;

	if (frech[i * nObs + j] <= 0)
	  return MINF;

	jac[i * nObs + j] = (ishape - 1) * log(frech[i * nObs + j]) - logScale;
	frech[i * nObs + j] = R_pow(frech[i * nObs + j], ishape);
      }
    }
  }
  return 0.0;
}

double gev2frechTrend(double *data, int nObs, int nSite, double *locs, double *scales,
		      double *shapes, double *trendlocs, double *trendscales,
		      double *trendshapes,double *jac, double *frech){

  /* This function transforms the GEV observations to unit Frechet
     ones with a temporal trend and computes the log of the jacobian
     of each transformation

     When ans > 0.0, the GEV parameters are invalid. */

  int i, j;

  for (i=nSite;i--;){
    for (j=nObs;j--;){
      double loc = locs[i] + trendlocs[j], scale = scales[i] + trendscales[j],
	shape = shapes[i] + trendshapes[j], iscale = 1 / scale, logScale = log(scale),
	ishape = 1 / shape;

      if (shape == 0.0){
	frech[i * nObs + j] = (data[i * nObs + j] - loc) * iscale;
	jac[i * nObs + j] = frech[i * nObs + j] - logScale;
	frech[i * nObs + j] = exp(frech[i * nObs + j]);
      }

      else {
	frech[i * nObs + j] = 1 + shape * (data[i * nObs + j] - loc) * iscale;

	if (frech[i * nObs + j] <= 0)
	  return MINF;

	jac[i * nObs + j] = (ishape - 1) * log(frech[i * nObs + j]) - logScale;
	frech[i * nObs + j] = R_pow(frech[i * nObs + j], ishape);
      }
    }
  }

  return 0.0;
}

double dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat, double *shapedsgnmat,
		     double *loccoeff, double *scalecoeff, double *shapecoeff,
		     int nSite, int nloccoeff, int nscalecoeff, int nshapecoeff,
		     double *locs, double *scales, double *shapes){

  //This function computes the GEV parameters from the design matrix
  //when ans > 0.0, the GEV parameters are invalid
  int i, j;

  memset(locs, 0, nSite * sizeof(double));
  memset(scales, 0, nSite * sizeof(double));
  memset(shapes, 0, nSite * sizeof(double));

  for (i=nSite;i--;){

    for (j=nloccoeff;j--;)
      locs[i] += loccoeff[j] * locdsgnmat[i + nSite * j];

    for (j=nscalecoeff;j--;)
      scales[i] += scalecoeff[j] * scaledsgnmat[i + nSite * j];

    if (scales[i] <= 0)
      return MINF;

    for (j=nshapecoeff;j--;)
      shapes[i] += shapecoeff[j] * shapedsgnmat[i + nSite * j];

    if (shapes[i] <= -1)
      return MINF;
  }

  return 0.0;
}

void dsgnmat2temptrend(double *dsgnmatloc, double *dsgnmatscale, double *dsgnmatshape,
		       double *loccoeff, double *scalecoeff, double *shapecoeff, int nSite,
		       int nObs, int *usetempcov, int nloccoeff, int nscalecoeff,
		       int nshapecoeff, double *trendlocs, double *trendscales,
		       double *trendshapes){

  //This function computes the temporal trend for each GEV parameters
  //from the design matrix.
  int i, j;

  memset(trendlocs, 0, nObs * sizeof(double));
  memset(trendscales, 0, nObs * sizeof(double));
  memset(trendshapes, 0, nObs * sizeof(double));

  if (usetempcov[0])
    for (i=nObs;i--;)
      for (j=nloccoeff;j--;)
	trendlocs[i] += loccoeff[j] * dsgnmatloc[i + nObs * j];

  if (usetempcov[1])
    for (i=nObs;i--;)
      for (j=nscalecoeff;j--;)
	trendscales[i] += scalecoeff[j] * dsgnmatscale[i + nObs * j];

  if (usetempcov[2])
    for (i=nObs;i--;)
      for (j=nshapecoeff;j--;)
	trendshapes[i] += shapecoeff[j] * dsgnmatshape[i + nObs * j];

  return;
}

void gev(double *prob, int *n, double *locs, double *scales, double *shapes,
	 double *quant){

  //This function computes the GEV quantiles

  int i;
  const double mlogProb = -log(*prob);

  for (i=*n;i--;){

    if (scales[i] <= 0){
      quant[i] = R_NaReal;

    }

    if (shapes[i] == 0)
      quant[i] = locs[i] - scales[i] * log(mlogProb);

    else
      quant[i] = locs[i] + scales[i] * (R_pow(mlogProb, -shapes[i]) - 1) /
	shapes[i];
  }
}


void dsgnmat2Alpha(double *alphadsgnmat, double *alphacoeff,
		   int nSite, int nalphacoeff, double *alphas){

  //This function computes the 'alpha' values from the design matrix
  //the 'alphas' are used in the schlatherind model
  //We use the expit function to ensure that the alphas always lie in
  //(0,1)
  int i, j;

  memset(alphas, 0, nSite * sizeof(double));

  for (i=nSite;i--;){
    for (j=nalphacoeff;j--;)
      alphas[i] += alphacoeff[j] * alphadsgnmat[i + nSite * j];

    //We use the expit function to ensure that the alphas lie in [0,1]
    alphas[i] = exp(alphas[i]) / (1 + exp(alphas[i]));

  }

  return;
}

void dsgnmat2Sigma2(double *sigma2dsgnmat, double *sigma2coeff,
		    int nSite, int nsigma2coeff, double *sigma2){

  //This function computes the 'sigma2' values from the design matrix
  //the 'sigma2' are used in the non-stationary geometric gaussian model
  int i, j;

  memset(sigma2, 0, nSite * sizeof(double));
  for (i=nSite;i--;){
    for (j=nsigma2coeff;j--;)
      sigma2[i] += sigma2coeff[j] * sigma2dsgnmat[i + nSite * j];

    //We use a log link function to ensure that the sigma2s lie are positive
    sigma2[i] = exp(sigma2[i]);

  }

  return;
}

double gev2unif(double *data, int nObs, int nSite, double *locs,
		double *scales, double *shapes, double *unif,
		double *logdens){

  /* This function transforms the GEV observations to U(0,1) ones.

     When ans > 0.0, the GEV parameters are invalid. */

  int i, j;
  
  for (i=nSite;i--;){
    double iscale = 1 / scales[i], logIscale = log(iscale);

    if (shapes[i] == 0.0)
      for (j=nObs;j--;){
	unif[i * nObs + j] = (data[i * nObs + j] - locs[i]) * iscale;
	logdens[i * nObs + j] = logIscale - unif[i * nObs + j] -
	  exp(-unif[i * nObs + j]);
	unif[i * nObs + j] = exp(-exp(-unif[i * nObs + j]));
      }

    else {
      double ishape = 1 / shapes[i];
      for (j=nObs;j--;){
	unif[i * nObs + j] = 1 + shapes[i] * (data[i * nObs + j] - locs[i]) *
	  iscale;

	if (unif[i * nObs + j] <= 0)
	  return MINF;

	logdens[i * nObs + j] = logIscale - (1 + ishape) * log(unif[i * nObs + j]) - 
	  R_pow(unif[i * nObs + j], -ishape);
	unif[i * nObs + j] = exp(-R_pow(unif[i * nObs + j], -ishape));
      }
    }
  }
  
  return 0.0;
}

double gev2unifTrend(double *data, int nObs, int nSite, double *locs,
		     double *scales, double *shapes, double *trendlocs,
		     double *trendscales, double *trendshapes, double *unif,
		     double *logdens){

  /* This function transforms the GEV observations to U(0,1) ones with
     a temporal trend.

     When ans > 0.0, the GEV parameters are invalid. */

  int i, j;

  for (i=nSite;i--;){
    for (j=nObs;j--;){
      double loc = locs[i] + trendlocs[j], scale = scales[i] + trendscales[j],
	shape = shapes[i] + trendshapes[j], iscale = 1 / scale,
	logIscale = log(iscale), ishape = 1 / shape;

      if (shape == 0.0){
	unif[i * nObs + j] = (data[i * nObs + j] - loc) * iscale;
	logdens[i * nObs + j] = logIscale - unif[i * nObs + j] -
	  exp(-unif[i * nObs + j]);
	unif[i * nObs + j] = exp(-exp(-unif[i * nObs + j]));
      }

      else {
	unif[i * nObs + j] = 1 + shape * (data[i * nObs + j] - loc) * iscale;

	if (unif[i * nObs + j] <= 0)
	  return MINF;

	logdens[i * nObs + j] = logIscale - (1 + ishape) * log(unif[i * nObs + j]) -
	  R_pow(unif[i * nObs + j], -ishape);
	unif[i * nObs + j] = exp(-R_pow(unif[i * nObs + j], -ishape));
      }
    }
  }

  return 0.0;
}
