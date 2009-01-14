#include "header.h"      
  
void smithstderr(double *data, double *distVec, int *nSite,
		 int *nObs, double *locdsgnmat, int *nloccoeff,
		 double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		 int *nshapecoeff, double *loccoeff, double *scalecoeff,
		 double *shapecoeff, double *cov11, double *cov12,
		 double *cov22, int *fitmarge, double *hess, double *grad){

  //This is the Smith model. It computes the hessian of the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, l, currentPair = -1;
  double c1, c2, dAa, dAz1, dAz2, B, dBa, dBz1, dBz2, C,
    dCa, dCz1, dCz2, D, dDa, dDz1, dDz2, *mahalDist, *locs,
    *scales, *shapes, jacCommonSigma, det, *jac, *frech,
    dz1loc, dz2loc, dz1scale, dz2scale, dz1shape, dz2shape, dE, flag;
  //c1, c2 are useful quantities
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //mahalanobis distance
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to cov11, cov12 or cov22
  //det is the determinant of the covariance matrix

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  det = *cov11 * *cov22 - *cov12 * *cov12;

  //Computing the Mahalanobis distance
  flag = mahalDistFct(distVec, nPairs, cov11, cov12,
		      cov22, mahalDist);

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];
  
  //Stage 2: Hessian computations
  // a- Covariance matrix part
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;
	
      for (k=*nObs;k--;){
    	double ifrech1 = 1 / frech[k + i * *nObs],
	  ifrech2 = 1 / frech[k + j * *nObs],
	  ifrech1Square = ifrech1 * ifrech1,
	  ifrech2Square = ifrech2 * ifrech2;
	  
	c1 = log(frech[k + j * *nObs] * ifrech1) * imahal +
	  0.5 * mahalDist[currentPair];
	c2 = mahalDist[currentPair] - c1;
	 
	double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);
	    
	//A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	  dnormc2 * imahal * ifrech2Square;
	C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	  dnormc1 * imahal * ifrech1Square;
	D = c2 * dnormc1 * ifrech2 * imahalSquare * ifrech1Square +
	  c1 * dnormc2 * ifrech1 * imahalSquare * ifrech2Square;

	dAa = - c2 * dnormc1 * ifrech1 * imahal - c1 * dnormc2 * imahal * ifrech2;
	dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	  (1 + c1 * c2 ) * dnormc1 * ifrech1 * ifrech2 * imahalSquare;
	dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square + 
	  (1 + c1 * c2) * dnormc2 * ifrech1 * ifrech2 * imahalSquare;
	dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahalSquare * imahal *
	  ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	  imahalSquare * imahal * ifrech1 * ifrech2Square;

	jacCommonSigma = dAa + (dBa * C + B * dCa + dDa) / (B*C + D);
	 
	hess[k * nPairs + currentPair] = -(*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) *
	  (*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) / 
	  (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;
	hess[(*nObs + k) * nPairs + currentPair] = (*cov11 * distVec[nPairs + currentPair] - 
						     *cov12 * distVec[currentPair]) * 
	  (*cov12 * distVec[nPairs + currentPair] - *cov22 * distVec[currentPair]) /
	  (det * det * mahalDist[currentPair]) * jacCommonSigma;
	hess[(2 * *nObs + k) * nPairs + currentPair] = -(*cov11 * distVec[nPairs + currentPair] - 
							 *cov12 * distVec[currentPair]) *
	  (*cov11 * distVec[nPairs + currentPair] -  *cov12 * distVec[currentPair]) / 
	  (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	grad[k] += hess[k * nPairs + currentPair];
	grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
      }
    }
  }

  if (*fitmarge){
    // b- Marginal part
    currentPair = -1;
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	currentPair++;

	double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;

	for (k=*nObs;k--;){
	 
	  double ifrech1 = 1 / frech[k + i * *nObs],
	    ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1,
	    ifrech2Square = ifrech2 * ifrech2;
	    
	  c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;

	  double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);
	    
	 
	  //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	  B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square;
	  C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square;
	  D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square;

	  double iBCplusD = 1 / (B * C + D);

	  dAz1 = dnormc1 * imahal * ifrech1Square + pnormc1 * ifrech1Square -
	    dnormc2 * imahal * ifrech1 * ifrech2;
	  dAz2 = dnormc2 * imahal * ifrech2Square + pnormc2 * ifrech2Square -
	    dnormc1 * imahal * ifrech1 * ifrech2;
	  dBz1 = D;
	  dBz2 = (mahalDist[currentPair] + c1) * dnormc1 * imahalSquare *
	    ifrech1 * ifrech2Square - 2 * pnormc2 * ifrech2 * ifrech2Square -
	    (2 * mahalDist[currentPair] + c1) * dnormc2 * imahalSquare * ifrech2 *
	    ifrech2Square;
	  dCz1 = (mahalDist[currentPair] + c2) * dnormc2 * imahalSquare *
	    ifrech1Square * ifrech2 - 2 * pnormc1 * ifrech1Square * ifrech1 -
	    (2 * mahalDist[currentPair] + c2) * dnormc1 * imahalSquare * ifrech1 *
	    ifrech1Square;
	  dCz2 = D;
	  dDz1 = (1 - c2 * (mahalDist[currentPair] + c2)) * dnormc1 *
	    imahal * imahalSquare * ifrech1 * ifrech1Square * ifrech2 -
	    (1 + c1 * (mahalDist[currentPair] + c2)) * dnormc2 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;
	  dDz2 = (1 - c1 * (mahalDist[currentPair] + c1)) * dnormc2 *
	    imahal * imahalSquare * ifrech1 * ifrech2 * ifrech2Square -
	    (1 + c2 * (mahalDist[currentPair] + c1)) * dnormc1 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;
	 
	  for (l=0;l<*nloccoeff;l++){
	    dE = (shapes[i] - 1) * locdsgnmat[i + *nSite * l] /
	      scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      (shapes[j] - 1) * locdsgnmat[j + *nSite * l] /
	      scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) /
	      scales[i] * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) /
	      scales[j] * locdsgnmat[j + *nSite * l];

	    hess[((3 + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1loc + dAz2 * dz2loc) +
	      ((dBz1 * dz1loc + dBz2 * dz2loc) * C + B * 
	       (dCz1 * dz1loc + dCz2 * dz2loc) + (dDz1 * dz1loc + dDz2 * dz2loc)) * iBCplusD + dE;
	    grad[(3 + l) * *nObs + k] += hess[((3 + l) * *nObs + k) * nPairs + currentPair];
	    
	  }
	  
	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (locs[i] - scales[i] - data[k + i * *nObs]) /
	      scales[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      scaledsgnmat[j + *nSite * l] * (locs[j] - scales[j] - data[k + j * *nObs]) /
	      scales[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) *
	      (data[k + i * *nObs] - locs[i]) / scales[i] / scales[i] *
	      scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) *
	      (data[k + j * *nObs] - locs[j]) / scales[j] / scales[j] *
	      scaledsgnmat[j + *nSite * l];

	    hess[((3 + *nloccoeff + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1scale + dAz2 * dz2scale) +
	      ((dBz1 * dz1scale + dBz2 * dz2scale) * C + B * 
	       (dCz1 * dz1scale + dCz2 * dz2scale) + (dDz1 * dz1scale + dDz2 * dz2scale)) * iBCplusD + dE;
	    grad[(3 + *nloccoeff + l) * *nObs + k] += hess[((3 + *nloccoeff + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shapes[i] + (1/shapes[i] - 1) * (data[k + i * *nObs] - locs[i]) *
	      shapedsgnmat[i + *nSite * l] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shapes[j] + (1/shapes[j] - 1) * (data[k + j * *nObs] - locs[j]) *
	      shapedsgnmat[j + *nSite * l] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);	      

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] * 
	      (-log(frech[k + i * *nObs]) / shapes[i] + (data[k + i * *nObs] - locs[i]) /
	       shapes[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] * 
	      (-log(frech[k + j * *nObs]) / shapes[j] + (data[k + j * *nObs] - locs[j]) /
	       shapes[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]));

	    hess[((3 + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1shape + dAz2 * dz2shape) +
	      ((dBz1 * dz1shape + dBz2 * dz2shape) * C + B * 
	       (dCz1 * dz1shape + dCz2 * dz2shape) + (dDz1 * dz1shape + dDz2 * dz2shape)) *
	      iBCplusD + dE;
	    grad[(3 + *nloccoeff + *nscalecoeff + l) * *nObs + k] += 
	      hess[((3 + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }
   
  return;
}

void smithstderr3d(double *data, double *distVec, int *nSite,
		   int *nObs, double *locdsgnmat, int *nloccoeff,
		   double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *loccoeff, double *scalecoeff,
		   double *shapecoeff, double *cov11, double *cov12, double *cov13,
		   double *cov22, double *cov23, double *cov33, int *fitmarge, double *hess,
		   double *grad){

  //This is the Smith model. It computes the hessian of the pairwise log-likelihood
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, l, currentPair = -1;
  double c1, c2, dAa, dAz1, dAz2, B, dBa, dBz1, dBz2, C,
    dCa, dCz1, dCz2, D, dDa, dDz1, dDz2, *mahalDist, *locs,
    *scales, *shapes, jacCommonSigma, det, *jac, *frech,
    dz1loc, dz2loc, dz1scale, dz2scale, dz1shape, dz2shape, dE, flag;
  //c1, c2 are useful quantities
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //mahalanobis distance
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to cov11, cov12 or cov22
  //det is the determinant of the covariance matrix

  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  det = *cov11 * *cov22 * *cov33 - *cov12 * *cov12 * *cov33 -
    *cov11 * *cov23 * *cov23 + 2 * *cov12 * *cov13 * *cov23 -
    *cov13 * *cov13 * *cov22;

  //Computing the Mahalanobis distance
  flag = mahalDistFct3d(distVec, nPairs, cov11, cov12, cov13,
			cov22, cov23, cov33, mahalDist);

  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];
  
  //Stage 2: Hessian computations
  // a- Covariance matrix part
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
      currentPair++;
      
      double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;

      for (k=*nObs;k--;){
	double ifrech1 = 1 / frech[k + i * *nObs],
	  ifrech2 = 1 / frech[k + j * *nObs],
	  ifrech1Square = ifrech1 * ifrech1,
	  ifrech2Square = ifrech2 * ifrech2;	  

	c1 = log(frech[k + j * *nObs] * ifrech1) * imahal +
	  0.5 * mahalDist[currentPair];
	c2 = mahalDist[currentPair] - c1;
	 
	double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);

	//A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	  dnormc2 * imahal * ifrech2Square;
	C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	  dnormc1 * imahal * ifrech1Square;
	D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	  c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square;

	dAa = - c2 * dnormc1 * imahal * ifrech1 - c1 * dnormc2 * imahal * ifrech2;
	dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	  (1 + c1 * c2 ) * dnormc1 * imahalSquare * ifrech1 * ifrech2;
	dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square +
	  (1 + c1 * c2) * dnormc2 * imahalSquare * ifrech1 * ifrech2;
	dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahal * imahalSquare *
	  ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	  imahalSquare * imahal * ifrech1 * ifrech2Square;

	jacCommonSigma = dAa + (dBa * C + B * dCa + dDa) / (B*C + D);

	hess[k * nPairs + currentPair] = -((*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
		    (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair] +
		    (*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair]) *
	  ((*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
	   (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair] +
	   (*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair]) /
	  (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	hess[(*nObs + k) * nPairs + currentPair]= ((*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair] +
			    (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair] +
			    (*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair]) *
	  ((*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair] +
	   (*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
	   (*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair]) /
	  (det * det * mahalDist[currentPair]) * jacCommonSigma;

	hess[(2 * *nObs + k) * nPairs + currentPair] = -((*cov22 * *cov33 - *cov23 * *cov23) * distVec[currentPair] +
				(*cov13 * *cov23 - *cov12 * *cov33) * distVec[nPairs + currentPair] +
				(*cov12 * *cov23 - *cov13 * *cov22) * distVec[2 * nPairs + currentPair]) *
	  ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	   (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair] +
	   (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair]) /
	  (det * det * mahalDist[currentPair]) * jacCommonSigma;
	
	hess[(3 * *nObs + k) * nPairs + currentPair] = -
	  ((*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair] +
	   (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair] +
	   (*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair]) *
	  ((*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair] +
	   (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair] +
	   (*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair]) / 
	  (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	hess[(4 * *nObs + k) * nPairs + currentPair] = 
	  ((*cov11 * *cov23 - *cov12 * *cov13) * distVec[2 * nPairs + currentPair] +
	   (*cov12 * *cov33 - *cov13 * *cov23) * distVec[currentPair] +
	   (*cov13 * *cov13 - *cov11 * *cov33) * distVec[nPairs + currentPair]) *
	  ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	   (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair] +
	   (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair]) /
	  (det * det * mahalDist[currentPair]) * jacCommonSigma;

	hess[(5 * *nObs + k) * nPairs + currentPair] = -
	  ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	   (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair] +
	   (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair]) *
	  ((*cov11 * *cov22 - *cov12 * *cov12) * distVec[2 * nPairs + currentPair] +
	   (*cov12 * *cov13 - *cov11 * *cov23) * distVec[nPairs + currentPair] +
	   (*cov12 * *cov23 - *cov13 * *cov22) * distVec[currentPair]) /
	  (2 * det * det * mahalDist[currentPair]) * jacCommonSigma;

	grad[k] += hess[k * nPairs + currentPair];
	grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	grad[4 * *nObs + k] += hess[(4 * *nObs + k) * nPairs + currentPair];
	grad[5 * *nObs + k] += hess[(5 * *nObs + k) * nPairs + currentPair];
      }
    }
  }

  if (*fitmarge){
    // b- Marginal part
    currentPair = -1;
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	 
	currentPair++;

	double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;

	for (k=*nObs;k--;){
	  double ifrech1 = 1 / frech[k + i * *nObs],
	    ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1,
	    ifrech2Square = ifrech2 * ifrech2;	    

	  c1 = log(frech[k + j * *nObs] * ifrech1) * imahal +
	    0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;
	 
	  double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);
	    
	  //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	  B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square;
	  C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square;
	  D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square;

	  double iBCplusD = 1 / (B * C + D);

	  dAz1 = dnormc1 * imahal * ifrech1Square + pnormc1 * ifrech1Square -
	    dnormc2 * imahal * ifrech1 * ifrech2;
	  dAz2 = dnormc2 * imahal * ifrech2Square + pnormc2 * ifrech2Square -
	    dnormc1 * imahal * ifrech1 * ifrech2;
	  dBz1 = D;
	  dBz2 = (mahalDist[currentPair] + c1) * dnormc1 * imahalSquare *
	    ifrech1 * ifrech2Square - 2 * pnormc2 * ifrech2 * ifrech2Square -
	    (2 * mahalDist[currentPair] + c1) * dnormc2 * imahalSquare * ifrech2 *
	    ifrech2Square;
	  dCz1 = (mahalDist[currentPair] + c2) * dnormc2 * imahalSquare *
	    ifrech1Square * ifrech2 - 2 * pnormc1 * ifrech1Square * ifrech1 -
	    (2 * mahalDist[currentPair] + c2) * dnormc1 * imahalSquare * ifrech1 *
	    ifrech1Square;
	  dCz2 = D;
	  dDz1 = (1 - c2 * (mahalDist[currentPair] + c2)) * dnormc1 *
	    imahal * imahalSquare * ifrech1 * ifrech1Square * ifrech2 -
	    (1 + c1 * (mahalDist[currentPair] + c2)) * dnormc2 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;
	  dDz2 = (1 - c1 * (mahalDist[currentPair] + c1)) * dnormc2 *
	    imahal * imahalSquare * ifrech1 * ifrech2 * ifrech2Square -
	    (1 + c2 * (mahalDist[currentPair] + c1)) * dnormc1 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;
	 
	  for (l=0;l<*nloccoeff;l++){
	    dE = (shapes[i] - 1) * locdsgnmat[i + *nSite * l] /
	      scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      (shapes[j] - 1) * locdsgnmat[j + *nSite * l] /
	      scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) /
	      scales[i] * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) /
	      scales[j] * locdsgnmat[j + *nSite * l];	    

	    hess[((6 + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1loc + dAz2 * dz2loc) +
	      ((dBz1 * dz1loc + dBz2 * dz2loc) * C + B * 
	       (dCz1 * dz1loc + dCz2 * dz2loc) + (dDz1 * dz1loc + dDz2 * dz2loc)) * iBCplusD + dE;
	    grad[(6 + l) * *nObs + k] +=  hess[((6 + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (locs[i] - scales[i] - data[k + i * *nObs]) /
	      scales[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      scaledsgnmat[j + *nSite * l] * (locs[j] - scales[j] - data[k + j * *nObs]) /
	      scales[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) *
	      (data[k + i * *nObs] - locs[i]) / scales[i] / scales[i] *
	      scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) *
	      (data[k + j * *nObs] - locs[j]) / scales[j] / scales[j] *
	      scaledsgnmat[j + *nSite * l];

	    hess[((6 + *nloccoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1scale + dAz2 * dz2scale) +
	      ((dBz1 * dz1scale + dBz2 * dz2scale) * C + B * 
	       (dCz1 * dz1scale + dCz2 * dz2scale) + (dDz1 * dz1scale + dDz2 * dz2scale)) *
	      iBCplusD + dE;
	    grad[(6 + *nloccoeff + l) * *nObs + k] +=
	      hess[((6 + *nloccoeff + l) * *nObs + k) * nPairs + currentPair];
	    
	  }

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shapes[i] + (1/shapes[i] - 1) * (data[k + i * *nObs] - locs[i]) *
	      shapedsgnmat[i + *nSite * l] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shapes[j] + (1/shapes[j] - 1) * (data[k + j * *nObs] - locs[j]) *
	      shapedsgnmat[j + *nSite * l] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);	      

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] * 
	      (-log(frech[k + i * *nObs]) / shapes[i] + (data[k + i * *nObs] - locs[i]) /
	       shapes[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] * 
	      (-log(frech[k + j * *nObs]) / shapes[j] + (data[k + j * *nObs] - locs[j]) /
	       shapes[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]));

	    hess[((6 + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1shape + dAz2 * dz2shape) +
	      ((dBz1 * dz1shape + dBz2 * dz2shape) * C + B * 
	       (dCz1 * dz1shape + dCz2 * dz2shape) + (dDz1 * dz1shape + dDz2 * dz2shape)) *
	      iBCplusD + dE;
	    grad[(6 + *nloccoeff + *nscalecoeff + l) * *nObs + k] +=
	      hess[((6 + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }
   
  return;
}


void schlatherstderr(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, double *locdsgnmat, int *nloccoeff,
		     double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		     int *nshapecoeff, double *loccoeff, double *scalecoeff,
		     double *shapecoeff, double *sill, double *range, double *smooth,
		     double *smooth2, int *fitmarge, double *hess, double *grad){

  /* This is the Schlather model. It computes the hessian of the
     pairwise log-likelihood */
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, l, currentPair = -1, nCorPar = 3;
  double c1, dArho, dAz1, dAz2, B, dBrho, dBz1, dBz2, C,
    dCrho, dCz1, dCz2, D, dDrho, dDz1, dDz2, *rho, *locs,
    *scales, *shapes, jacCommonRho, *jac, *frech, dz1loc,
    dz2loc, dz1scale, dz2scale, dz1shape, dz2shape, dE, flag;
  //c1 is a useful quantity
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //covariance function
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to cov11, cov12 or cov22
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  //Stage 0: Compute the covariance at each location
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
    //Here we use 0 for dim as we don't care for the computation of
    //the hessian
    flag = bessel(dist, nPairs, 0, *sill, *range, *smooth, rho);
    break;
  case 5:
    flag = caugen(dist, nPairs, *sill, *range, *smooth, *smooth2, rho);
  }
  
  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];
  
  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
	
      currentPair++;

      for (k=*nObs;k--;){
	c1 = sqrt(frech[k + i * *nObs] * frech[k + i * *nObs] + 
		  frech[k + j * *nObs] * frech[k + j * *nObs] -
		  2 * frech[k + i * *nObs] * frech[k + j * *nObs] *
		  rho[currentPair]);

	B = (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1);
	C = (- rho[currentPair] * frech[k + i * *nObs] + c1 +
	     frech[k + j * *nObs]) / (2 * c1 * frech[k + i * *nObs] *
				      frech[k + i * *nObs]);
	D = (-rho[currentPair] * frech[k + j * *nObs] + c1 +
	     frech[k + i * *nObs]) / (2 * c1 * frech[k + j * *nObs] *
				      frech[k + j * *nObs]);

	dArho =  1 / (2 * c1);
	dBrho = - rho[currentPair] / (c1 * c1 * c1) + 3 * 
	  (1 - rho[currentPair] * rho[currentPair]) * frech[k + i * *nObs] *
	  frech[k + j * *nObs] / (2 * c1 * c1 * c1 * c1 * c1);
	dCrho = (-frech[k + i * *nObs] + frech[k + j * *nObs] *
		 rho[currentPair]) / (2 * c1 * c1 * c1);
	dDrho = (-frech[k + j * *nObs] + frech[k + i * *nObs] *
		 rho[currentPair]) / (2 * c1 * c1 * c1);		   

	jacCommonRho = dArho + (dBrho + dCrho * D + C * dDrho) / (B + C * D);
	 
	hess[k * nPairs + currentPair] = rho[currentPair] / *sill * jacCommonRho;
	grad[k] += hess[k * nPairs + currentPair];

	switch (*covmod){
	case 1:
	  //i.e. Whittle-Matern
	  hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * 
	    (-2 * *smooth / *range + dist[currentPair] * 
	     bessel_k(dist[currentPair] / *range, *smooth + 1, 1) / 
	     (bessel_k(dist[currentPair] / *range, *smooth, 1) * *range * *range)) *
	    jacCommonRho;
	  //The Whittle-Matern covariance function is not
	  //differentiable w.r.t. to the smooth parameter
	  hess[(2 * *nObs + k) * nPairs + currentPair] = R_NaReal;
	  break;
	case 2:
	  //i.e. cauchy
	  hess[(*nObs + k) * nPairs + currentPair] = 2 * dist[currentPair] * dist[currentPair] *
	    *sill * *smooth / (*range * *range * *range) * 
	    R_pow(1 + dist[currentPair] * dist[currentPair] / (*range * *range),
		  - *smooth - 1) * jacCommonRho; 
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] * 
	    log1p(dist[currentPair] *dist[currentPair] / ( *range * *range)) *
	    jacCommonRho;
	  break;
	case 3:
	  //i.e. powered exponential
	  hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range *
	    R_pow(dist[currentPair] / *range, *smooth) * jacCommonRho;
	  hess[(2 * *nObs + k) * nPairs + currentPair] -= rho[currentPair] *
	    R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) *
	    jacCommonRho;
	  break;
	case 4:
	  //i.e. Bessel
	  hess[(*nObs + k) * nPairs + currentPair] = *sill * 
	    R_pow(2 * *range / dist[currentPair], *smooth) *
	    dist[currentPair] / (*range * *range) * gammafn(*smooth + 1) *
	    bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommonRho;
	  hess[(2 * *nObs + k) * nPairs + currentPair] = R_NaReal;
	  break;
	case 5:
	  //i.e. Generalized Cauchy
	  hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range /
	    (1 + R_pow(*range / dist[currentPair], *smooth2)) *
	    jacCommonRho;
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	    rho[currentPair] / *smooth2 *jacCommonRho;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
				  *smooth2 - log(dist[currentPair] / *range) /
				  (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	    *smooth / *smooth2 * rho[currentPair] * jacCommonRho;
	  grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
	  //The caugen has 2 smooth paramaters so add + 1 to nCorPar.
	  nCorPar = 4;
	  break;
	}

	grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];
	grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
      }
    }
  }

  if (*fitmarge){
    // b- Marginal part
    currentPair = -1;
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	
	currentPair++;
	
	for (k=*nObs;k--;){
	  c1 = sqrt(frech[k + i * *nObs] * frech[k + i * *nObs] + 
		    frech[k + j * *nObs] * frech[k + j * *nObs] -
		    2 * frech[k + i * *nObs] * frech[k + j * *nObs] *
		    rho[currentPair]);
	  
	  B = (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1);
	  C = - (rho[currentPair] * frech[k + i * *nObs] - c1 -
		 frech[k + j * *nObs]) / (2 * c1 * frech[k + i * *nObs] *
					  frech[k + i * *nObs]);
	  D = - (rho[currentPair] * frech[k + j * *nObs] - c1 -
		 frech[k + i * *nObs]) / (2 * c1 * frech[k + j * *nObs] *
					  frech[k + j * *nObs]);
	  
	  dAz1 = (frech[k + j * *nObs] + c1 - rho[currentPair] * 
		  frech[k + i * *nObs]) / (2 * c1 * frech[k + i * *nObs] *
					   frech[k + i * *nObs]);
	  dAz2 = (frech[k + i * *nObs] + c1 - rho[currentPair] * 
		  frech[k + j * *nObs]) / (2 * c1 * frech[k + j * *nObs] *
					   frech[k + j * *nObs]);
	  dBz1 = 3 * (rho[currentPair] * rho[currentPair] - 1) * 
	    (frech[k + i * *nObs] - rho[currentPair] * 
	     frech[k + j * *nObs]) / (2 * c1 * c1 * c1 * c1 * c1);
	  dBz2 = 3 * (rho[currentPair] * rho[currentPair] - 1) * 
	    (frech[k + j * *nObs] - rho[currentPair] * 
	     frech[k + i * *nObs]) / (2 * c1 * c1 * c1 * c1 * c1);
	  dCz1 = (2 * rho[currentPair] * frech[k + i * *nObs] * frech[k + i * *nObs] *
		  frech[k + i * *nObs] + 6 * frech[k + i * *nObs] * frech[k + j * *nObs] * 
		  frech[k + j * *nObs] * rho[currentPair] -
		  3 * frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + j * *nObs] *
		  (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 -
		  2 * frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + j * *nObs]) /
	    (2 * c1 * c1 * c1 * frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + i * *nObs]);
	  dCz2 = B;
	  dDz1 = B;
	  dDz2 = (2 * rho[currentPair] * frech[k + j * *nObs] * frech[k + j * *nObs] *
		  frech[k + j * *nObs] + 6 * frech[k + j * *nObs] * frech[k + i * *nObs] * 
		  frech[k + i * *nObs] * rho[currentPair] -
		  3 * frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + i * *nObs] *
		  (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 -
		  2 * frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + i * *nObs]) /
	    (2 * c1 * c1 * c1 * frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + j * *nObs]);
	  	 
	  for (l=0;l<*nloccoeff;l++){
	    dE = (shapes[i] - 1) * locdsgnmat[i + *nSite * l] /
	      scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      (shapes[j] - 1) * locdsgnmat[j + *nSite * l] /
	      scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);
	    
	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) /
	      scales[i] * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) /
	      scales[j] * locdsgnmat[j + *nSite * l];

	    hess[((nCorPar + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1loc + dAz2 * dz2loc) + 
	      ((dBz1 * dz1loc + dBz2 * dz2loc) +
	       (dCz1 * dz1loc + dCz2 * dz2loc) * D +
	       (dDz1 * dz1loc + dDz2 * dz2loc) * C) /
	      (B + C * D) + dE;
	    grad[(nCorPar + l) * *nObs + k] += hess[((nCorPar + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (locs[i] - scales[i] - data[k + i * *nObs]) /
	      scales[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      scaledsgnmat[j + *nSite * l] * (locs[j] - scales[j] - data[k + j * *nObs]) /
	      scales[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) *
	      (data[k + i * *nObs] - locs[i]) / scales[i] / scales[i] *
	      scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) *
	      (data[k + j * *nObs] - locs[j]) / scales[j] / scales[j] *
	      scaledsgnmat[j + *nSite * l];

	    hess[((nCorPar + *nloccoeff + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1scale + dAz2 * dz2scale) +
	      ((dBz1 * dz1scale + dBz2 * dz2scale) + (dCz1 * dz1scale + dCz2 * dz2scale) * D +
	       (dDz1 * dz1scale + dDz2 * dz2scale) * C) / (B + C * D) + dE;
	    grad[(nCorPar + *nloccoeff + l) * *nObs + k] +=
	      hess[((nCorPar + *nloccoeff + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shapes[i] + (1/shapes[i] - 1) * (data[k + i * *nObs] - locs[i]) *
	      shapedsgnmat[i + *nSite * l] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shapes[j] + (1/shapes[j] - 1) * (data[k + j * *nObs] - locs[j]) *
	      shapedsgnmat[j + *nSite * l] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] * 
	      (-log(frech[k + i * *nObs]) / shapes[i] + (data[k + i * *nObs] - locs[i]) /
	       shapes[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] * 
	      (-log(frech[k + j * *nObs]) / shapes[j] + (data[k + j * *nObs] - locs[j]) /
	       shapes[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]));

	    hess[((nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) +
						     (dCz1 * dz1shape + dCz2 * dz2shape) * D +
						     (dDz1 * dz1shape + dDz2 * dz2shape) * C) /	      
	      (B + C * D) + dE;
	    grad[(nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k] += 
	      hess[((nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  return;
}

void schlatherindstderr(int *covmod, double *data, double *dist, int *nSite,
			int *nObs, double *locdsgnmat, int *nloccoeff,
			double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
			int *nshapecoeff, double *loccoeff, double *scalecoeff,
			double *shapecoeff, double *alpha, double *sill, double *range,
			double *smooth, double *smooth2, int *fitmarge, double *hess,
			double *grad){

  /* This is the independent Schlather model. It computes the hessian
    of the pairwise log-likelihood */
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, l, currentPair = -1, nCorPar = 4;
  double c1, dAalpha, dArho, dAz1, dAz2, B, dBalpha, dBrho,
    dBz1, dBz2, C, dCalpha, dCrho, dCz1, dCz2, D, dDalpha,
    dDrho, dDz1, dDz2, *rho, *locs, *scales, *shapes, jacCommonRho,
    *jac, *frech, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape,
    dz2shape, dE, flag;
  //c1 is a useful quantity
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //covariance function
  //jacCommon is the common part of all the jacobians - i.e. the common
  //part when deriving with respect to cov11, cov12 or cov22
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  rho = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));

  //Stage 0: Compute the covariance at each location
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
    //Here we use 0 for dim as we don't care for the computation of
    //the hessian
    flag = bessel(dist, nPairs, 0, *sill, *range, *smooth, rho);
    break;
  case 5:
    flag = caugen(dist, nPairs, *sill, *range, *smooth, *smooth2, rho);
    break;
  }
  
  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
    
  }

  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];
  
  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
	
      currentPair++;
      for (k=*nObs;k--;){
	c1 = sqrt(frech[k + i * *nObs] * frech[k + i * *nObs] + 
		  frech[k + j * *nObs] * frech[k + j * *nObs] -
		  2 * frech[k + i * *nObs] * frech[k + j * *nObs] *
		  rho[currentPair]);

	B = (1 - *alpha) * (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1);
	C = (*alpha - 1) * (rho[currentPair] * frech[k + i * *nObs] - c1 -
			    frech[k + j * *nObs]) / (2 * c1 * frech[k + i * *nObs] *
						     frech[k + i * *nObs]) +
	  *alpha / (frech[k + i * *nObs] * frech[k + i * *nObs]);
	D = (*alpha - 1) * (rho[currentPair] * frech[k + j * *nObs] - c1 -
			    frech[k + i * *nObs]) / (2 * c1 * frech[k + j * *nObs] *
						     frech[k + j * *nObs]) +
	  *alpha / (frech[k + j * *nObs] * frech[k + j * *nObs]);
	
	dAalpha = (c1 - frech[k + i * *nObs] - frech[k + j * *nObs]) /
	  (2 * frech[k + i * *nObs] * frech[k + j * *nObs]);
	dArho =  (1 - *alpha) / (2 * c1);

	dBalpha = (rho[currentPair] * rho[currentPair] - 1) / (2 * c1 * c1 * c1);
	dBrho = (3 * (1 - rho[currentPair] * rho[currentPair]) * frech[k + i * *nObs] *
		 frech[k + j * *nObs] / (2 * c1 * c1 * c1 * c1 * c1) - rho[currentPair] /
		 (c1 * c1 * c1)) * (1 - *alpha);

	dCalpha = (rho[currentPair] * frech[k + i * *nObs] + c1 - frech[k + j * *nObs]) /
	  (2 * c1 * frech[k + i * *nObs] * frech[k + i * *nObs]);
	dCrho = (*alpha - 1) * (frech[k + i * *nObs] - frech[k + j * *nObs] *
				rho[currentPair]) / (2 * c1 * c1 * c1);

	dDalpha = (rho[currentPair] * frech[k + j * *nObs] + c1 - frech[k + i * *nObs]) /
	  (2 * c1 * frech[k + j * *nObs] * frech[k + j * *nObs]);
	dDrho = (*alpha - 1) * (frech[k + j * *nObs] - frech[k + i * *nObs] *
				rho[currentPair]) / (2 * c1 * c1 * c1);		   

	jacCommonRho = dArho + (dBrho + dCrho * D + dDrho * C) / (B + C * D);
	
	hess[k * nPairs + currentPair] = dAalpha + (dBalpha + dCalpha * D + dDalpha * C) / (B + C * D);
	hess[(*nObs + k) * nPairs + currentPair] = rho[currentPair] / *sill * jacCommonRho;

	grad[k] += hess[k * nPairs + currentPair];
	grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];

	switch (*covmod){
	case 1:
	  //i.e. Whittle-Matern
	  hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] * 
	    (-2 * *smooth / *range + dist[currentPair] * 
	     bessel_k(dist[currentPair] / *range, *smooth + 1, 1) / 
	     bessel_k(dist[currentPair] / *range, *smooth, 1) / (*range * *range)) *
	    jacCommonRho;
	  //The Whittle-Matern covariance function is not
	  //differentiable w.r.t. to the smooth parameter
	  hess[(3 * *nObs + k) * nPairs + currentPair] = R_NaReal;
	  break;
	case 2:
	  //i.e. cauchy
	  hess[(2 * *nObs + k) * nPairs + currentPair] = 2 * dist[currentPair] * dist[currentPair] *
	    *sill * *smooth / (*range * *range * *range) * 
	    R_pow(1 + dist[currentPair] * dist[currentPair] / (*range * *range),
		  - *smooth - 1) * jacCommonRho; 
	  hess[(3 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] * 
	    log1p(dist[currentPair] * dist[currentPair] / (*range * *range)) *
	    jacCommonRho;
	  break;
	case 3:
	  //i.e. powered exponential
	  hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth /
	    *range * R_pow(dist[currentPair] / *range, *smooth) * jacCommonRho;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = -rho[currentPair] *
	    R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) *
	    jacCommonRho;
	  break;
	case 4:
	  //i.e. Bessel
	  hess[(2 * *nObs + k) * nPairs + currentPair] = *sill *
	    R_pow(2 * *range / dist[currentPair], *smooth) *
	    dist[currentPair] / (*range * *range) * gammafn(*smooth + 1) *
	    bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommonRho;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = R_NaReal;
	  break;
	case 5:
	  //i.e. Generalized Cauchy
	  hess[(2 * *nObs + k) * nPairs + currentPair] = rho[currentPair] * *smooth / *range /
	    (1 + R_pow(*range / dist[currentPair], *smooth2)) *
	    jacCommonRho;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = -log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	    rho[currentPair] / *smooth2 *jacCommonRho;
	  hess[(4 * *nObs + k) * nPairs + currentPair] = (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
							  *smooth2 - log(dist[currentPair] / *range) /
							  (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	    *smooth / *smooth2 * rho[currentPair] * jacCommonRho;
	  grad[4 * *nObs + k] += hess[(4 * *nObs + k) * nPairs + currentPair];
	  //The caugen has 2 smooth paramaters so add + 1 to nCorPar.
	  nCorPar = 5;
	  break;
	}
	grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
      }
    }
  }

  if (*fitmarge){
    // b- Marginal part
    currentPair = -1;
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	  
	currentPair++;
	for (k=*nObs;k--;){
	  c1 = sqrt(frech[k + i * *nObs] * frech[k + i * *nObs] + 
		    frech[k + j * *nObs] * frech[k + j * *nObs] -
		    2 * frech[k + i * *nObs] * frech[k + j * *nObs] *
		    rho[currentPair]);
	  
	  B = (1 - *alpha) * (1 - rho[currentPair] * rho[currentPair]) / (2 * c1 * c1 * c1);
	  C = (*alpha - 1) * (rho[currentPair] * frech[k + i * *nObs] - c1 -
			      frech[k + j * *nObs]) / (2 * c1 * frech[k + i * *nObs] *
						       frech[k + i * *nObs]) +
	    *alpha / (frech[k + i * *nObs] * frech[k + i * *nObs]);
	  D = (*alpha - 1) * (rho[currentPair] * frech[k + j * *nObs] - c1 -
			      frech[k + i * *nObs]) / (2 * c1 * frech[k + j * *nObs] *
						       frech[k + j * *nObs]) +
	    *alpha / (frech[k + j * *nObs] * frech[k + j * *nObs]);
	  
	  dAz1 = C;
	  dAz2 = D;
	  dBz1 = 3 * (1 - *alpha) * (rho[currentPair] * rho[currentPair] - 1) * 
	    (frech[k + i * *nObs] - rho[currentPair] * frech[k + j * *nObs]) / 
	    (2 * c1 * c1 * c1 * c1 * c1);
	  dBz2 = 3 * (1 - *alpha) * (rho[currentPair] * rho[currentPair] - 1) * 
	    (frech[k + j * *nObs] - rho[currentPair] * frech[k + i * *nObs]) /
	    (2 * c1 * c1 * c1 * c1 * c1);
	  dCz1 = (2 * rho[currentPair] * frech[k + i * *nObs] * frech[k + i * *nObs] *
		  frech[k + i * *nObs] + 6 * frech[k + i * *nObs] * frech[k + j * *nObs] * 
		  frech[k + j * *nObs] * rho[currentPair] -
		  3 * frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + j * *nObs] *
		  (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 -
		  2 * frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + j * *nObs]) /
	    (2 * c1 * c1 * c1 * frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + i * *nObs]) *
	    (1 - *alpha) - 2 * *alpha / (frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + i * *nObs]);
	  dCz2 = B;
	  dDz1 = dCz2;
	  dDz2 = (2 * rho[currentPair] * frech[k + j * *nObs] * frech[k + j * *nObs] *
		  frech[k + j * *nObs] + 6 * frech[k + j * *nObs] * frech[k + i * *nObs] * 
		  frech[k + i * *nObs] * rho[currentPair] -
		  3 * frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + i * *nObs] *
		  (1 + rho[currentPair] * rho[currentPair]) - 2 * c1 * c1 * c1 -
		  2 * frech[k + i * *nObs] * frech[k + i * *nObs] * frech[k + i * *nObs]) /
	    (2 * c1 * c1 * c1 * frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + j * *nObs]) *
	    (1 - *alpha) - 2 * *alpha / (frech[k + j * *nObs] * frech[k + j * *nObs] * frech[k + j * *nObs]);
	  	 
	  for (l=0;l<*nloccoeff;l++){
	    dE = (shapes[i] - 1) * locdsgnmat[i + *nSite * l] /
	      scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      (shapes[j] - 1) * locdsgnmat[j + *nSite * l] /
	      scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);
	    
	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) /
	      scales[i] * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) /
	      scales[j] * locdsgnmat[j + *nSite * l];

	    hess[((nCorPar + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1loc + dAz2 * dz2loc) + 
	      ((dBz1 * dz1loc + dBz2 * dz2loc) + 
	       (dCz1 * dz1loc + dCz2 * dz2loc) * D +
	       C * (dDz1 * dz1loc + dDz2 * dz2loc)) / (B + C * D) + dE;
	    grad[(nCorPar + l) * *nObs + k] += hess[((nCorPar + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (locs[i] - scales[i] - data[k + i * *nObs]) /
	      scales[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      scaledsgnmat[j + *nSite * l] * (locs[j] - scales[j] - data[k + j * *nObs]) /
	      scales[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);
	    
	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) *
	      (data[k + i * *nObs] - locs[i]) / (scales[i] * scales[i]) *
	      scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) *
	      (data[k + j * *nObs] - locs[j]) / (scales[j] * scales[j]) *
	      scaledsgnmat[j + *nSite * l];

	    hess[((nCorPar + *nloccoeff + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1scale + dAz2 * dz2scale) + 
	      ((dBz1 * dz1scale + dBz2 * dz2scale) + 
	       (dCz1 * dz1scale + dCz2 * dz2scale) * D +
	       C * (dDz1 * dz1scale + dDz2 * dz2scale)) / (B + C * D) + dE;
	    grad[(nCorPar + *nloccoeff + l) * *nObs + k] +=
	      hess[((nCorPar + *nloccoeff + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shapes[i] + (1/shapes[i] - 1) * (data[k + i * *nObs] - locs[i]) *
	      shapedsgnmat[i + *nSite * l] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shapes[j] + (1/shapes[j] - 1) * (data[k + j * *nObs] - locs[j]) *
	      shapedsgnmat[j + *nSite * l] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] * 
	      (-log(frech[k + i * *nObs]) / shapes[i] + (data[k + i * *nObs] - locs[i]) /
	       shapes[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] * 
	      (-log(frech[k + j * *nObs]) / shapes[j] + (data[k + j * *nObs] - locs[j]) /
	       shapes[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]));

	    hess[((nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1shape + dAz2 * dz2shape) + ((dBz1 * dz1shape + dBz2 * dz2shape) +
						     (dCz1 * dz1shape + dCz2 * dz2shape) * D +
						     (dDz1 * dz1shape + dDz2 * dz2shape) * C) /	      
	      (B + C * D) + dE;
	    grad[(nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k] += 
	      hess[((nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }

  return;
}

void geomgaussstderr(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, double *locdsgnmat, int *nloccoeff,
		     double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		     int *nshapecoeff, double *loccoeff, double *scalecoeff,
		     double *shapecoeff, double *sigma2, double *sill, double *range,
		     double *smooth, double *smooth2, int *fitmarge, double *hess,
		     double *grad){
  /* This function computes the hessian of the log-pairwise
     likelihood for the geometric gaussian model.
     
     Remember that this model shares the same bivariate density as the
     Smith model except that the Mahalanobis distance is modified. */
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, l, currentPair = -1, nCorPar = 4;
  double c1, c2, dAa, dAz1, dAz2, B, dBa, dBz1, dBz2, C, dCa, dCz1,
    dCz2, D, dDa, dDz1, dDz2, rho, *mahalDist, *locs, *scales, *shapes,
    jacCommon, *jac, *frech, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape,
    dz2shape, dE, flag;
  //c1 is a useful quantity
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //covariance function
  //jacCommon is the common part of all the jacobians
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  
  //Stage 0: Compute the covariance at each location
  flag = geomCovariance(dist, nPairs, 0, *covmod, *sigma2, 1.0e6, *sill, *range,
			*smooth, *smooth2, mahalDist);
  
  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
      
  }
  
  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];
  
  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;

      rho = 1 - mahalDist[currentPair] * mahalDist[currentPair] /
	  (2 * *sigma2);

      for (k=*nObs;k--;){
	double ifrech1 = 1 / frech[k + i * *nObs],
	  ifrech2 = 1 / frech[k + j * *nObs],
	  ifrech1Square = ifrech1 * ifrech1,
	  ifrech2Square = ifrech2 * ifrech2;	  

	c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair];
	c2 = mahalDist[currentPair] - c1;

	double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);
      
	//A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	  dnormc2 * imahal * ifrech2Square;
	C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	  dnormc1 * imahal * ifrech1Square;
	D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	  c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square;

	dAa = - c2 * dnormc1 * imahal * ifrech1 - c1 * dnormc2 * imahal * ifrech2;
	dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	  (1 + c1 * c2 ) * dnormc1 * imahalSquare * ifrech1 * ifrech2;
	dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square +
	  (1 + c1 * c2) * dnormc2 * imahalSquare * ifrech1 * ifrech2;
	dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahalSquare * imahal *
	  ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	  imahalSquare * imahal * ifrech1 * ifrech2Square;

	jacCommon = dAa + (dBa * C + B * dCa + dDa) / (B*C + D);

	hess[k * nPairs + currentPair] = mahalDist[currentPair] / (2 * *sigma2) * jacCommon;
	hess[(*nObs + k) * nPairs + currentPair] = -*sigma2 * rho / (mahalDist[currentPair] * *sill) *
	  jacCommon;

	grad[k] += hess[k * nPairs + currentPair];
	grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];

	switch (*covmod){
	case 1:
	  //i.e. Whittle-Matern
	  hess[(2 * *nObs + k) * nPairs + currentPair] = *sigma2 * rho / mahalDist[currentPair] * 
	    (2 * *smooth / *range - dist[currentPair] * 
	     bessel_k(dist[currentPair] / *range, *smooth + 1, 1) / 
	     (bessel_k(dist[currentPair] / *range, *smooth, 1) * *range * *range)) *
	    jacCommon;
	  //The Whittle-Matern covariance function is not
	  //differentiable w.r.t. to the smooth parameter
	  hess[(3 * *nObs + k) * nPairs + currentPair] = R_NaReal;
	  break;
	case 2:
	  //i.e. cauchy
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -2 * *sigma2 * dist[currentPair] * dist[currentPair] *
	    *sill * *smooth * R_pow(1 + dist[currentPair] * dist[currentPair] /
				    (*range * *range), - *smooth - 1) /
	    (*range * *range * *range * mahalDist[currentPair]) * jacCommon; 
	  hess[(3 * *nObs + k) * nPairs + currentPair] = *sigma2 * rho * 
	    log1p(dist[currentPair] *dist[currentPair] / ( *range * *range)) /
	    mahalDist[currentPair] * jacCommon;
	  break;
	case 3:
	  //i.e. powered exponential
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -*sigma2 * rho * R_pow(dist[currentPair] / *range, *smooth) *
	    *smooth / (*range * mahalDist[currentPair]) * jacCommon;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = *sigma2 * rho *
	    R_pow(dist[currentPair] / *range, *smooth) * log(dist[currentPair] / *range) /
	    mahalDist[currentPair] * jacCommon;
	  break;
	case 4:
	  //i.e. Bessel
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -*sigma2 * *sill * 
	    R_pow(2 * *range / dist[currentPair], *smooth) * gammafn(*smooth + 1) *
	    dist[currentPair] / (*range * *range * mahalDist[currentPair]) * 
	    bessel_j(dist[currentPair] / *range, *smooth + 1) * jacCommon;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = R_NaReal;
	  break;
	case 5:
	  //i.e. Generalized Cauchy
	  hess[(2 * *nObs + k) * nPairs + currentPair] = -*sigma2 * rho * *smooth / 
	    (*range * (1 + R_pow(*range / dist[currentPair], *smooth2)) * mahalDist[currentPair]) *
	    jacCommon;
	  hess[(3 * *nObs + k) * nPairs + currentPair] = *sigma2 * log1p(R_pow(dist[currentPair] / *range, *smooth2)) *
	    rho / (*smooth2 * mahalDist[currentPair])* jacCommon;
	  hess[(4 * *nObs + k) * nPairs + currentPair] = -*sigma2 * (log1p(R_pow(dist[currentPair] / *range, *smooth2)) /
				  *smooth2 - log(dist[currentPair] / *range) /
				  (1 + R_pow(*range / dist[currentPair], *smooth2))) *
	    *smooth / *smooth2 * rho /mahalDist[currentPair] * jacCommon;
	  grad[4 * *nObs + k] += hess[(4 * *nObs + k) * nPairs + currentPair];
	  //The caugen has one more smooth parameter so
	  nCorPar = 5;
	  break;
	}
	grad[2 * *nObs + k] += hess[(2 * *nObs + k) * nPairs + currentPair];
	grad[3 * *nObs + k] += hess[(3 * *nObs + k) * nPairs + currentPair];
      }
    }
  }

  if (*fitmarge){
    // b- Marginal part
    currentPair = -1;
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	 
	currentPair++;
	double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;

	for (k=*nObs;k--;){
	  double ifrech1 = 1 / frech[k + i * *nObs],
	  ifrech2 = 1 / frech[k + j * *nObs],
	  ifrech1Square = ifrech1 * ifrech1,
	    ifrech2Square = ifrech2 * ifrech2;
	  	 
	  c1 = log(frech[k + j * *nObs] * ifrech1) * imahal +
	    0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;
	 
	  double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);

	  //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	  B = - dnormc1 *imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square;
	  C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square;
	  D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square;

	  double iBCplusD = 1 / (B * C + D);
	  
	  dAz1 = dnormc1 * imahal * ifrech1Square + pnormc1 * ifrech1Square -
	    dnormc2 * imahal * ifrech1 * ifrech2;
	  dAz2 = dnormc2 * imahal * ifrech2Square + pnormc2 * ifrech2Square -
	    dnormc1 * imahal * ifrech1 * ifrech2;
	  dBz1 = D;
	  dBz2 = (mahalDist[currentPair] + c1) * dnormc1 * imahalSquare * ifrech1 *
	    ifrech2Square - 2 * pnormc2 * ifrech2 * ifrech2Square -
	    (2 * mahalDist[currentPair] + c1) * dnormc2 * imahalSquare * ifrech2 *
	    ifrech2Square;
	  dCz1 = (mahalDist[currentPair] + c2) * dnormc2 * imahalSquare * ifrech1Square *
	    ifrech2 - 2 * pnormc1 * ifrech1 * ifrech1Square -
	    (2 * mahalDist[currentPair] + c2) * dnormc1 * imahalSquare * ifrech1 *
	    ifrech1Square;
	  dCz2 = D;
	  dDz1 = (1 - c2 * (mahalDist[currentPair] + c2)) * dnormc1 * imahal *
	    imahalSquare * ifrech1 * ifrech1Square * ifrech2 -
	    (1 + c1 * (mahalDist[currentPair] + c2)) * dnormc2 * imahal * imahalSquare *
	    ifrech1Square * ifrech2Square;
	  dDz2 = (1 - c1 * (mahalDist[currentPair] + c1)) * dnormc2 * imahal * imahalSquare *
	    ifrech1 * ifrech2 * ifrech2Square -
	    (1 + c2 * (mahalDist[currentPair] + c1)) * dnormc1 * imahal * imahalSquare * 
	    ifrech1Square * ifrech2Square;
	 
	  for (l=0;l<*nloccoeff;l++){
	    dE = (shapes[i] - 1) * locdsgnmat[i + *nSite * l] /
	      scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      (shapes[j] - 1) * locdsgnmat[j + *nSite * l] /
	      scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) /
	      scales[i] * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) /
	      scales[j] * locdsgnmat[j + *nSite * l];

	    hess[((nCorPar + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1loc + dAz2 * dz2loc) +
	      ((dBz1 * dz1loc + dBz2 * dz2loc) * C + B * 
	       (dCz1 * dz1loc + dCz2 * dz2loc) + (dDz1 * dz1loc + dDz2 * dz2loc)) * iBCplusD + dE;
	    grad[(nCorPar + l) * *nObs + k] += hess[((nCorPar + l) * *nObs + k) * nPairs + currentPair];
	  }
	  
	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (locs[i] - scales[i] - data[k + i * *nObs]) /
	      scales[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      scaledsgnmat[j + *nSite * l] * (locs[j] - scales[j] - data[k + j * *nObs]) /
	      scales[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) *
	      (data[k + i * *nObs] - locs[i]) / scales[i] / scales[i] *
	      scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) *
	      (data[k + j * *nObs] - locs[j]) / scales[j] / scales[j] *
	      scaledsgnmat[j + *nSite * l];

	    hess[((nCorPar + *nloccoeff + l) * *nObs + k) * nPairs + currentPair] =
	      (dAz1 * dz1scale + dAz2 * dz2scale) + 
	      ((dBz1 * dz1scale + dBz2 * dz2scale) * C + B * 
	       (dCz1 * dz1scale + dCz2 * dz2scale) + (dDz1 * dz1scale + dDz2 * dz2scale)) *
	      iBCplusD + dE;
	    grad[(nCorPar + *nloccoeff + l) * *nObs + k] +=
	      hess[((nCorPar + *nloccoeff + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shapes[i] + (1/shapes[i] - 1) * (data[k + i * *nObs] - locs[i]) *
	      shapedsgnmat[i + *nSite * l] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shapes[j] + (1/shapes[j] - 1) * (data[k + j * *nObs] - locs[j]) *
	      shapedsgnmat[j + *nSite * l] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);	      

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] * 
	      (-log(frech[k + i * *nObs]) / shapes[i] + (data[k + i * *nObs] - locs[i]) /
	       shapes[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] * 
	      (-log(frech[k + j * *nObs]) / shapes[j] + (data[k + j * *nObs] - locs[j]) /
	       shapes[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]));

	    hess[((nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1shape + dAz2 * dz2shape) +
	      ((dBz1 * dz1shape + dBz2 * dz2shape) * C + B * 
	       (dCz1 * dz1shape + dCz2 * dz2shape) + (dDz1 * dz1shape + dDz2 * dz2shape)) *
	      iBCplusD + dE;
	    grad[(nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k] +=
	      hess[((nCorPar + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }
   
  return;
}

void brownresnickstderr(double *data, double *dist, int *nSite, int *nObs,
			double *locdsgnmat, int *nloccoeff, double *scaledsgnmat,
			int *nscalecoeff, double *shapedsgnmat, int *nshapecoeff,
			double *loccoeff, double *scalecoeff, double *shapecoeff,
			double *range, double *smooth, int *fitmarge, double *hess,
			double *grad){
  /* This function computes the hessian of the log-pairwise
     likelihood for the Brown-Resnick model.
     
     Remember that this model shares the same bivariate density as the
     Smith model except that the Mahalanobis distance is modified. */
  
  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i, j, k, l, currentPair = -1;
  double c1, c2, dAa, dAz1, dAz2, B, dBa, dBz1, dBz2, C, dCa, dCz1,
    dCz2, D, dDa, dDz1, dDz2, *mahalDist, *locs, *scales, *shapes,
    jacCommon, *jac, *frech, dz1loc, dz2loc, dz1scale, dz2scale, dz1shape,
    dz2shape, dE, flag;
  //c1 is a useful quantity
  //A, B, C, D are part of the log-bivariate density
  //dB, dC and dD are their derivatives with respect to the
  //covariance function
  //jacCommon is the common part of all the jacobians
  
  jac = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  frech = (double *)R_alloc(*nObs * *nSite, sizeof(double));
  
  //Stage 0: Compute the covariance at each location
  flag = brownResnick(dist, nPairs, *range, *smooth, mahalDist);
  
  //Compute the GEV parameters using the design matrix
  if (*fitmarge){
    
    flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
			 loccoeff, scalecoeff, shapecoeff, *nSite,
			 *nloccoeff, *nscalecoeff, *nshapecoeff,
			 locs, scales, shapes);
    
    //Stage 1: Transformation to unit Frechet
    flag = gev2frech(data, *nObs, *nSite, locs, scales, shapes,
		     jac, frech);
      
  }
  
  else
    for (i=(*nSite * *nObs);i--;)
      frech[i] = data[i];
  
  //Stage 2: Hessian computations;
  // a- Covariance part
  for (i=0;i<(*nSite-1);i++){
    for (j=i+1;j<*nSite;j++){
      currentPair++;

      double imahal = 1 / mahalDist[currentPair],
	  imahalSquare = imahal * imahal;

      for (k=*nObs;k--;){
	double ifrech1 = 1 / frech[k + i * *nObs],
	  ifrech2 = 1 / frech[k + j * *nObs],
	  ifrech1Square = ifrech1 * ifrech1,
	  ifrech2Square = ifrech2 * ifrech2;
	  
	c1 = log(frech[k + j * *nObs] * ifrech1) * imahal +
	  0.5 * mahalDist[currentPair];
	c2 = mahalDist[currentPair] - c1;
	 
	double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);
	    
	//A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	  dnormc2 * imahal * ifrech2Square;
	C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	  dnormc1 * imahal * ifrech1Square;
	D = c2 * dnormc1 * ifrech2 * imahalSquare * ifrech1Square +
	  c1 * dnormc2 * ifrech1 * imahalSquare * ifrech2Square;

	dAa = - c2 * dnormc1 * ifrech1 * imahal - c1 * dnormc2 * imahal * ifrech2;
	dBa = (c1 * c1 - 1) * dnormc2 * imahalSquare * ifrech2Square +
	  (1 + c1 * c2 ) * dnormc1 * ifrech1 * ifrech2 * imahalSquare;
	dCa = (c2 * c2 - 1) * dnormc1 * imahalSquare * ifrech1Square + 
	  (1 + c1 * c2) * dnormc2 * ifrech1 * ifrech2 * imahalSquare;
	dDa = (c1 - c1 * c2 * c2 - 2 * c2) * dnormc1 * imahalSquare * imahal *
	  ifrech1Square * ifrech2 + (c2 - c1 * c1 *c2 - 2 * c1) * dnormc2 *
	  imahalSquare * imahal * ifrech1 * ifrech2Square;

	jacCommon = dAa + (dBa * C + B * dCa + dDa) / (B*C + D);

	hess[k * nPairs + currentPair] = -0.5 * *smooth * mahalDist[currentPair] / *range * jacCommon;
	hess[(*nObs + k) * nPairs + currentPair] = 0.5 * log(dist[currentPair] / *range) *
	  mahalDist[currentPair] * jacCommon;

	grad[k] += hess[k * nPairs + currentPair];
	grad[*nObs + k] += hess[(*nObs + k) * nPairs + currentPair];

      }
    }
  }

  if (*fitmarge){
    // b- Marginal part
    currentPair = -1;
    for (i=0;i<(*nSite-1);i++){
      for (j=i+1;j<*nSite;j++){
	 
	currentPair++;
	double imahal = 1 / mahalDist[currentPair],
	    imahalSquare = imahal * imahal;

	for (k=*nObs;k--;){
	  double ifrech1 = 1 / frech[k + i * *nObs],
	    ifrech2 = 1 / frech[k + j * *nObs],
	    ifrech1Square = ifrech1 * ifrech1,
	    ifrech2Square = ifrech2 * ifrech2;	    

	  c1 = log(frech[k + j * *nObs] * ifrech1) * imahal + 0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;

	  double dnormc1 = dnorm(c1, 0., 1., 0),
	    pnormc1 = pnorm(c1, 0., 1., 1, 0),
	    dnormc2 = dnorm(c2, 0., 1., 0),
	    pnormc2 = pnorm(c2, 0., 1., 1, 0);
	    
	 
	  //A = - pnormc1 * ifrech1 - pnormc2 * ifrech2;
	  B = - dnormc1 * imahal * ifrech1 * ifrech2 + pnormc2 * ifrech2Square +
	    dnormc2 * imahal * ifrech2Square;
	  C = - dnormc2 * imahal * ifrech1 * ifrech2 + pnormc1 * ifrech1Square +
	    dnormc1 * imahal * ifrech1Square;
	  D = c2 * dnormc1 * imahalSquare * ifrech1Square * ifrech2 +
	    c1 * dnormc2 * imahalSquare * ifrech1 * ifrech2Square;

	  double iBCplusD = 1 / (B * C + D);

	  dAz1 = dnormc1 * imahal * ifrech1Square + pnormc1 * ifrech1Square -
	    dnormc2 * imahal * ifrech1 * ifrech2;
	  dAz2 = dnormc2 * imahal * ifrech2Square + pnormc2 * ifrech2Square -
	    dnormc1 * imahal * ifrech1 * ifrech2;
	  dBz1 = D;
	  dBz2 = (mahalDist[currentPair] + c1) * dnormc1 * imahalSquare *
	    ifrech1 * ifrech2Square - 2 * pnormc2 * ifrech2 * ifrech2Square -
	    (2 * mahalDist[currentPair] + c1) * dnormc2 * imahalSquare * ifrech2 *
	    ifrech2Square;
	  dCz1 = (mahalDist[currentPair] + c2) * dnormc2 * imahalSquare *
	    ifrech1Square * ifrech2 - 2 * pnormc1 * ifrech1Square * ifrech1 -
	    (2 * mahalDist[currentPair] + c2) * dnormc1 * imahalSquare * ifrech1 *
	    ifrech1Square;
	  dCz2 = D;
	  dDz1 = (1 - c2 * (mahalDist[currentPair] + c2)) * dnormc1 *
	    imahal * imahalSquare * ifrech1 * ifrech1Square * ifrech2 -
	    (1 + c1 * (mahalDist[currentPair] + c2)) * dnormc2 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;
	  dDz2 = (1 - c1 * (mahalDist[currentPair] + c1)) * dnormc2 *
	    imahal * imahalSquare * ifrech1 * ifrech2 * ifrech2Square -
	    (1 + c2 * (mahalDist[currentPair] + c1)) * dnormc1 *
	    imahal * imahalSquare * ifrech1Square * ifrech2Square;
	 
	  for (l=0;l<*nloccoeff;l++){
	    dE = (shapes[i] - 1) * locdsgnmat[i + *nSite * l] /
	      scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      (shapes[j] - 1) * locdsgnmat[j + *nSite * l] /
	      scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1loc = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) /
	      scales[i] * locdsgnmat[i + *nSite * l];
	    dz2loc = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) /
	      scales[j] * locdsgnmat[j + *nSite * l];

	    hess[((2 + l) * *nObs + k) * nPairs + currentPair] = (dAz1 * dz1loc + dAz2 * dz2loc) +
	      ((dBz1 * dz1loc + dBz2 * dz2loc) * C + B * 
	       (dCz1 * dz1loc + dCz2 * dz2loc) + (dDz1 * dz1loc + dDz2 * dz2loc)) * iBCplusD + dE;
	    grad[(2 + l) * *nObs + k] += hess[((2 + l) * *nObs + k) * nPairs + currentPair];
	  }
	  
	  for (l=0;l<*nscalecoeff;l++){
	    dE = scaledsgnmat[i + *nSite * l] * (locs[i] - scales[i] - data[k + i * *nObs]) /
	      scales[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) +
	      scaledsgnmat[j + *nSite * l] * (locs[j] - scales[j] - data[k + j * *nObs]) /
	      scales[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);

	    dz1scale = - R_pow(frech[k + i * *nObs], 1 - shapes[i]) *
	      (data[k + i * *nObs] - locs[i]) / scales[i] / scales[i] *
	      scaledsgnmat[i + *nSite * l];
	    dz2scale = - R_pow(frech[k + j * *nObs], 1 - shapes[j]) *
	      (data[k + j * *nObs] - locs[j]) / scales[j] / scales[j] *
	      scaledsgnmat[j + *nSite * l];

	    hess[((2 + *nloccoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1scale + dAz2 * dz2scale) +
	      ((dBz1 * dz1scale + dBz2 * dz2scale) * C + B * 
	       (dCz1 * dz1scale + dCz2 * dz2scale) + (dDz1 * dz1scale + dDz2 * dz2scale)) *
	      iBCplusD + dE;
	    grad[(2 + *nloccoeff + l) * *nObs + k] +=
	      hess[((2 + *nloccoeff + l) * *nObs + k) * nPairs + currentPair];
	  }

	  for (l=0;l<*nshapecoeff;l++){
	    dE = -shapedsgnmat[i + *nSite * l] * log(frech[k + i * *nObs]) /
	      shapes[i] + (1/shapes[i] - 1) * (data[k + i * *nObs] - locs[i]) *
	      shapedsgnmat[i + *nSite * l] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]) -
	      shapedsgnmat[j + *nSite * l] * log(frech[k + j * *nObs]) /
	      shapes[j] + (1/shapes[j] - 1) * (data[k + j * *nObs] - locs[j]) *
	      shapedsgnmat[j + *nSite * l] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]);	      

	    dz1shape = frech[k + i * *nObs] * shapedsgnmat[i + *nSite * l] * 
	      (-log(frech[k + i * *nObs]) / shapes[i] + (data[k + i * *nObs] - locs[i]) /
	       shapes[i] / scales[i] / R_pow(frech[k + i * *nObs], shapes[i]));
	    dz2shape = frech[k + j * *nObs] * shapedsgnmat[j + *nSite * l] * 
	      (-log(frech[k + j * *nObs]) / shapes[j] + (data[k + j * *nObs] - locs[j]) /
	       shapes[j] / scales[j] / R_pow(frech[k + j * *nObs], shapes[j]));

	    hess[((2 + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair] = 
	      (dAz1 * dz1shape + dAz2 * dz2shape) +
	      ((dBz1 * dz1shape + dBz2 * dz2shape) * C + B * 
	       (dCz1 * dz1shape + dCz2 * dz2shape) + (dDz1 * dz1shape + dDz2 * dz2shape)) *
	      iBCplusD + dE;
	    grad[(2 + *nloccoeff + *nscalecoeff + l) * *nObs + k] +=
	      hess[((2 + *nloccoeff + *nscalecoeff + l) * *nObs + k) * nPairs + currentPair];
	  }
	}
      }
    }
  }
   
  return;
}


void spatgevstderr(double *data, int *nSite, int *nObs, double *locdsgnmat,
		   int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		   double *shapedsgnmat, int *nshapecoeff, double *loccoeff,
		   double *scalecoeff, double *shapecoeff, double *hess, 
		   double *grad){

  int i, j, k;
  double *locs, *scales, *shapes, flag;
    
  locs = (double *)R_alloc(*nSite, sizeof(double));
  scales = (double *)R_alloc(*nSite, sizeof(double));
  shapes = (double *)R_alloc(*nSite, sizeof(double));
  
  //Stage 0. Compute the GEV parameters using the design matrix
  flag = dsgnmat2Param(locdsgnmat, scaledsgnmat, shapedsgnmat,
		       loccoeff, scalecoeff, shapecoeff, *nSite,
		       *nloccoeff, *nscalecoeff, *nshapecoeff,
		       locs, scales, shapes);
    
  //Stage 1. Compute the gradient
  for (i=0;i<*nObs;i++){
    for (j=0;j<*nSite;j++){

      double dataTrans = 1 + shapes[j] * (data[j * *nObs + i] - locs[j]) / scales[j];

      for (k=0;k<*nloccoeff;k++){
	hess[(k * *nObs + i) * *nSite + j] = 
	  ((1 + shapes[j]) / dataTrans - R_pow(dataTrans, - 1 / shapes[j] - 1)) *
	  locdsgnmat[k * *nSite + j] / scales[j];
	grad[k * *nObs + i] += hess[(k * *nObs + i) * *nSite + j];
      }

      for (k=0;k<*nscalecoeff;k++){
	hess[((*nloccoeff + k) * *nObs + i) * *nSite + j] = 
	  (-1 + (1 + shapes[j]) * (data[j * *nObs + i] - locs[j]) / (scales[j] * dataTrans) -
	   R_pow(dataTrans, - 1 / shapes[j] - 1) * (data[j * *nObs + i] - locs[j]) / scales[j]) *
	  scaledsgnmat[k * *nSite + j] / scales[j];
	grad[(*nloccoeff + k) * *nObs + i] += hess[((*nloccoeff + k) * *nObs + i) * *nSite + j];
      }

      for (k=0;k<*nshapecoeff;k++){
	hess[((*nloccoeff + *nscalecoeff + k) * *nObs + i) * *nSite + j] = 
	  (log1p(dataTrans - 1) / shapes[j] - (1 + shapes[j]) * (data[j * *nObs + i] - locs[j]) / 
	   (scales[j] * dataTrans) - R_pow(dataTrans, - 1 / shapes[j]) *
	   (log1p(dataTrans - 1) / shapes[j] - (data[j * *nObs + i] - locs[j]) / 
	    (scales[j] * dataTrans))) * shapedsgnmat[k * *nSite + j] / shapes[j];
	grad[(*nloccoeff + *nscalecoeff + k) * *nObs + i] +=
	  hess[((*nloccoeff + *nscalecoeff + k) * *nObs + i) * *nSite + j];
      }
    }
  }

  return;
}
