#include "header.h"

void gevlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns){

  //It computes the log-likelihood for the GEV
  int i;
  double *dvec;
  
  dvec = (double *)R_alloc(*n, sizeof(double));

  if( (*scale <= 0) & (*shape < -1)) {
    *dns = -1e6;
    return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - *loc) / *scale;
    
    if(fabs(*shape) <= 1e-6){
      *shape = 0.0;
      dvec[i] = -log(*scale) - data[i] - exp(-data[i]);
    }

    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
	*dns = -1e6;
	return;
      }
      dvec[i] = -log(*scale) - R_pow(data[i], -1 / *shape) -
	(1 / *shape + 1) * log(data[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];

  return;
}

void gpdlik(double *exceed, int *n, double *thresh, double *scale,
	    double *shape, double *dns){
  //It computes the log-likelihood for the GPD
  int i;
  double *dvec;
  
  dvec = (double *)R_alloc(*n, sizeof(double));

  if ((*scale <= 0) && (*shape < -1)) {
    *dns = -1e6;
    return;
  }

  for (i=0;i<*n;i++) {
    exceed[i] = (exceed[i] - *thresh) / *scale;
    
    if (exceed[i] <= 0) {
      *dns = -1e6;
      return;
    }

    if(fabs(*shape) <= 1e-6){
      *shape = 0.0; 
      dvec[i] = -log(*scale) - exceed[i];
    }

    else {
      exceed[i] = 1 + *shape * exceed[i];
      
      if (exceed[i] <= 0) {
	*dns = -1e6;
	return;
      }
      
      dvec[i] = -log(*scale) - (1 / *shape + 1) * log(exceed[i]);
    }
  }
  
  for(i=0;i<*n;i++) 
    *dns = *dns + dvec[i];

  return;
}
