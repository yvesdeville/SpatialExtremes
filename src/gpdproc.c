#include "header.h"

void gpdprocfull(double *data, double *distVec, int *nSite,
		 int *nObs, double *excRates, double *threshs, double *scales,
		 double *shapes, double *cov11, double *cov12,
		 double *cov22, int *fitmarge, double *dns){
  //This is the Smith gpd process model. It's a wrapper to several
  //sub-functions. It's named xxxfull as it either assume that the
  //margins are unit Frechet, or the GPD parameters are estimated at
  //each locations.

  const int nPairs = *nSite * (*nSite - 1) / 2;
  int i;
  double *jac, *mahalDist, *ugpd;
  
  jac = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  mahalDist = (double *)R_alloc(nPairs, sizeof(double));
  ugpd = (double *)R_alloc(*nSite * *nObs, sizeof(double));
  
  //Some preliminary steps: Valid points?
  if (*fitmarge){
    for (i=*nSite;i--;){
      if ((scales[i] <= 0) || (shapes[i] <= -1)){
	*dns = MINF;
	return;
      }
    }
  }

  
  //Stage 1: Computing the Mahalanobis distance
  *dns = mahalDistFct(distVec, nPairs, cov11, cov12, cov22, mahalDist);

  if (*dns != 0.0)
      return;

  //Check if the parameters define an appropriate upper bound
  double uBound = 0.5 * M_1_PI / sqrt(*cov11 * *cov22 - *cov12 * *cov12);
  for (i=*nSite;i--;){
    if ((uBound * excRates[i]) > 1){
      *dns = MINF;
      return;
    }
  }

  //Stage 2: Transformation to unit GPD if any
  if (*fitmarge){
    *dns = gpd2ugpd(data, *nObs, *nSite, excRates, threshs, scales,
		    shapes, jac, ugpd);

    if (*dns != 0.0)
      return;

    *dns = lpliksmithgpd(ugpd, mahalDist, jac, excRates, *nObs, *nSite);
  }
  
  else {
    memset(jac, 0, *nSite * *nObs * sizeof(double));
    *dns = lpliksmithgpd(data, mahalDist, jac, excRates, *nObs, *nSite);
  }
  
  return;
}

double gpd2ugpd(double *data, int nObs, int nSite, double *excRates,
		double *threshs, double *scales, double *shapes,
		double *jac, double *ugpd){

  //This function transforms the GPD observations to unit GPD ones
  //and computes the log of the jacobian of each transformation
  //When ans > 0.0, the GPD parameters are invalid.
  
  int i, j;
      
  for (i=nSite;i--;){
    double iscale = 1 / scales[i], iexcRate = 1 / excRates[i];

    if (shapes[i] == 0.0){
            
      for (j=nObs;j--;){

	if (data[i * nObs + j] > threshs[i]){
	  ugpd[i * nObs + j] = (data[i * nObs + j] - threshs[i]) * iscale;
	  jac[i * nObs + j] = ugpd[i * nObs + j] - log(scales[i] * excRates[i]);
	  ugpd[i * nObs + j] = exp(ugpd[i * nObs + j]) * iexcRate;
	  
	}

	else {
	  ugpd[i * nObs + j] = iexcRate;
	  jac[i * nObs + j] = 0;
	}
      }
    }
      
    else {
      double ishape = 1 / shapes[i];

      for (j=nObs;j--;){

	if (data[i * nObs + j] > threshs[i]){
	  ugpd[i * nObs + j] = 1 + shapes[i] * (data[i * nObs + j] - threshs[i]) *
	    iscale;
	  
	  if (ugpd[i * nObs + j] <= 0)
	    return MINF;
	  
	  jac[i * nObs + j] = (ishape - 1) * log(ugpd[i * nObs + j]) -
	    log(scales[i] * excRates[i]);
	  ugpd[i * nObs + j] = R_pow(ugpd[i * nObs + j], ishape) * iexcRate;
	  
	}

	else {
	  ugpd[i * nObs + j] = iexcRate;
	  jac[i * nObs + j] = 0;
	}
      }
    }
  }
  
  return 0.0;
}

double lpliksmithgpd(double *data, double *mahalDist, double *jac,
		     double *excRates, int nObs, int nSite){
  /*This function computes the log-pairwise likelihood for the smith
    gpd process model.

    Rmq: 1 / excRates corresponds to the threshold at the unit GPD
    scale */
  
  int i, j, k, currentPair = -1;
  double c1, c2, dns = 0.0, idata1, idata2, idata1Square, idata2Square,
    imahal;

  for (i=0;i<(nSite - 1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      imahal = 1 / mahalDist[currentPair];
            
      for (k=nObs;k--;){

	/* Now we use the censored likelihood so we need to check if
	   both components are extremes or not. whichCase is a way to
	   know in which region (y_1, y_2) lies */
	int whichCase = ((data[k + i * nObs] * excRates[i]) > 1) +
	  2 * ((data[k + j * nObs] * excRates[j]) > 1);

	switch (whichCase){
	case 0:
	  //Both components are non extremes
	  c1 = log(excRates[i] / excRates[j]) * imahal + 0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;
	  dns += log(1 - excRates[i] * pnorm(c1, 0, 1, 1, 0) - excRates[j] *
		     pnorm(c2, 0, 1, 1, 0)); 
	  break;
	case 1:
	  //Component 1 is extreme, 2 isn't
	  idata1 = 1 / data[k + i * nObs];
	  idata1Square = idata1 * idata1;
	  
	  c1 = -log(data[k + i * nObs] * excRates[j]) * imahal + 0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;
	  
	  if ((fabs(c1) > 38) && (fabs(c2) > 38))
	    /* This means that only data1 or data2 contributes to the
	       log-likelihood. Hence we consider these parameters as
	       unfeasible as the bivariate distribution in this case
	       degenerates
	       
		 Rmq: 38 is the limiting accuracy for dnorm */
	    return MINF;
	  
	  else
	    dns += log(pnorm(c1, 0, 1, 1, 0) * idata1Square + dnorm(c1, 0, 1, 0) *
		       imahal * idata1Square - excRates[j] * dnorm(c2, 0, 1, 0) * imahal *
		       idata1) + jac[k + i * nObs];
	  break;
	case 2:
	  //Component 2 is extreme, 1 isn't
	  idata2 = 1 / data[k + j * nObs];
	  idata2Square = idata2 * idata2;
	  
	  c1 = log(data[k + j * nObs] * excRates[i]) * imahal + 0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;
	  
	  if ((fabs(c1) > 38) && (fabs(c2) > 38))
	    /* This means that only data1 or data2 contributes to the
	       log-likelihood. Hence we consider these parameters as
	       unfeasible as the bivariate distribution in this case
	       degenerates
	       
	       Rmq: 38 is the limiting accuracy for dnorm */
	    return MINF;
	  
	  else
	    dns += log(pnorm(c2, 0, 1, 1, 0) * idata2Square + dnorm(c2, 0, 1, 0) *
		       imahal * idata2Square - excRates[i] * dnorm(c1, 0, 1, 0) * imahal * 
		       idata2) + jac[k + j * nObs];
	  break;
	case 3:
	  //Both are extremes.
	  c1 = log(data[k + j * nObs] / data[k + i * nObs]) * imahal + 0.5 * mahalDist[currentPair];
	  c2 = mahalDist[currentPair] - c1;
	  
	  if ((fabs(c1) > 38) && (fabs(c2) > 38))
	    /* This means that only data1 or data2 contributes to the
	       log-likelihood. Hence we consider these parameters as
	       unfeasible as the bivariate distribution in this case
	       degenerates
	       
	       Rmq: 38 is the limiting accuracy for dnorm */
	    return MINF;
	  
	  else{
	    dns += log(c1 * data[k + i * nObs] * dnorm(c2, 0, 1, 0) +
		       c2 * data[k + j * nObs] * dnorm(c1, 0, 1, 0)) -
	      2 * log(mahalDist[currentPair] * data[k + i * nObs] * data[k + j * nObs]) +
	      jac[k + i * nObs] + jac[k + j * nObs];
	  }	
	  break;
	
	}
      }
    }
  }
  
  return dns;
}
