#include "header.h"

void gevlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns){

  //It computes the log-likelihood for the GEV
  int i;
  
  if( (*scale <= 0) || (*shape < -1)) {
    *dns = -1e6;
    return;
  }

  if (fabs(*shape) <= 1e-16){
    for (i=0;i<*n;i++){
      data[i] = (data[i] - *loc) / *scale;
      *dns += -log(*scale) - data[i] - exp(-data[i]);
    }
  }

  else{
    for(i=0;i<*n;i++){
      
      data[i] = 1 + *shape * (data[i] - *loc) / *scale;
      
      if (data[i] <= 0) {
	*dns = -1e6;
	return;
      }
      
      *dns += -log(*scale) - R_pow(data[i], -1 / *shape) -
	(1 / *shape + 1) * log(data[i]);
    }
  }

  return;
}

void gpdlik(double *exceed, int *n, double *thresh, double *scale,
	    double *shape, double *dns){
  //It computes the log-likelihood for the GPD
  int i;
  
  if ((*scale <= 0) || (*shape < -1)) {
    *dns = -1e6;
    return;
  }

  if (fabs(*shape) <= 1e-16){
    for (i=0;i<*n;i++){
      exceed[i] = (exceed[i] - *thresh) / *scale;

      if (exceed[i] <= 0){
	*dns = -1e6;
	return;
      }
      
      *dns += -log(*scale) - exceed[i];
    }
  }

  else{
    for (i=0;i<*n;i++) {
      exceed[i] = (exceed[i] - *thresh) / *scale;
      
      if (exceed[i] <= 0) {
	*dns = -1e6;
	return;
      }
      
      exceed[i] = 1 + *shape * exceed[i];
      
      if (exceed[i] <= 0) {
	*dns = -1e6;
	return;
      }
      
      *dns += -log(*scale) - (1 / *shape + 1) * log(exceed[i]);
    }
  }
  
  return;
}
