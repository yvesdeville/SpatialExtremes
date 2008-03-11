#include "header.h"

void extCoeffSmith(double *frech, int *nObs, int *nSite,
		   double *extCoeff){
  int i, j, k, currentPair = 0;

  for (i=0;i<(*nSite - 1);i++){
    for (j=i+1;j<*nSite;j++){
      for (k=0;k<*nObs;k++){
	extCoeff[currentPair] = extCoeff[currentPair] + 
	  fmin2(frech[i * *nObs + k], frech[j * *nObs + k]);
      }
      extCoeff[currentPair] = *nObs / extCoeff[currentPair];

      if (extCoeff[currentPair] > 2)
	extCoeff[currentPair] = 2;

      if (extCoeff[currentPair] < 1)
	extCoeff[currentPair] = 1;

      currentPair++;
    }
  }

  return;
}

void extCoeffST(double *frech, double *xbar, double *z, double *theta,
		int *nObs, double *dns){

  int i;
  double frechScaled;

  for (i=0;i<*nObs;i++){
    frechScaled = fmax2(frech[i] * xbar[0], 
			frech[i + *nObs] * xbar[1]);

    if (frechScaled > *z)
      *dns =  *dns - log(*theta) + *theta / frechScaled;

    else
      *dns = *dns + *theta / *z;

  }

  return;
}

  

  
  
  

