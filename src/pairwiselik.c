#include "header.h"

double lplikschlather(double *data, double *rho, double *jac,
		      int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schlather model.

  int i, j, k, currentPair = -1;
  double c1, dns, lFvec, dvecM1, dvecM2, dvecMixed, pbm,
    data1Square, data2Square;
  //c1 is a useful quantity - see documentation

  dns = 0.0;
  for (i=0;i<(nSite-1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      
      for (k=0;k<nObs;k++){
	data1Square = R_pow_di(data[k + i * nObs], 2);
	data2Square = R_pow_di(data[k + j * nObs], 2);

	c1 = sqrt(data1Square + data2Square - 2 * data[k + i * nObs] *
		  data[k + j * nObs] * rho[currentPair]);
	
	//It's the log of the joint CDF
	lFvec = - (1/data[k + i * nObs] + 1/data[k + j * nObs]) *
	  (1 + sqrt(1 - 2 * (rho[currentPair] + 1) * data[k + i * nObs] *
		    data[k + j * nObs] / 
		    R_pow_di(data[k + i * nObs] + data[k + j * nObs], 2))) / 2;
	
	//It's the partial derivative for marge 1
	dvecM1 = -(rho[currentPair] * data[k + i * nObs] - c1 - 
		   data[k + j * nObs]) / 2 / c1 / data1Square;

	//It's the partial derivative for marge 2
	dvecM2 = -(rho[currentPair] * data[k + j * nObs] - c1 - 
		   data[k + i * nObs]) / 2 / c1 / data2Square;

	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.

	//It's the mixed partial derivative
	dvecMixed = (1 - R_pow_di(rho[currentPair], 2)) / 2 / R_pow_di(c1, 3) + 
	  dvecM1 * dvecM2;

	pbm = (1 + fabs(dvecMixed)) * MINF;
	dvecMixed = log(dvecMixed);

	if (!R_FINITE(dvecMixed)){
	  //printf("dvecMixed is erradic\n");
	  dns += pbm;
	}

	else{
	  //Now the final step, multiplying by Fvec and the gradient
	  dns += dvecMixed + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	}
      }
    }
  }

  return dns;

}

double lpliksmith(double *data, double *mahalDist, double *jac,
		  int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //smith model.
  
  int i, j, k, currentPair = -1;
  double c1, c2, dns, lFvec, dvecM1, dvecM2, dvecMixed,
    dnormc1, dnormc2, pnormc1, pnormc2, data1Square,
    data2Square;

  dns = 0.0;
  for (i=0;i<(nSite - 1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      
      for (k=0;k<nObs;k++){
	c1 = log(data[k + j * nObs] / data[k + i * nObs]) /
	  mahalDist[currentPair] + mahalDist[currentPair] / 2;
	c2 = mahalDist[currentPair] - c1;
	data1Square = R_pow_di(data[k + i * nObs], 2);
	data2Square = R_pow_di(data[k + j * nObs], 2);
	
	dnormc1 = dnorm(c1, 0., 1., 0);
	dnormc2 = dnorm(c2, 0., 1., 0);
	pnormc1 = pnorm(c1, 0., 1., 1, 0);
	pnormc2 = pnorm(c2, 0., 1., 1, 0);

	//It's the log of the joint CDF
	lFvec = -pnormc1 / data[k + i * nObs] - pnormc2 / data[k + j * nObs];

	//It's the partial derivative for marge 1
	dvecM1 = (dnormc1 / mahalDist[currentPair]  + pnormc1) /
	  data1Square - dnormc2 / data[k + i * nObs] / data[k + j * nObs] /
	  mahalDist[currentPair];
	
	//It's the partial derivative for marge 2
	dvecM2 = (dnormc2 / mahalDist[currentPair] + pnormc2) / data2Square -
	  dnormc1 / data[k + i * nObs] / data[k + j * nObs] /
	  mahalDist[currentPair];
	
	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.
	
	//It's the mixed partial derivative
	dvecMixed = dvecM1 * dvecM2 + 
	  (c2 * dnormc1 / data1Square / data[k + j * nObs] +
	   c1 * dnormc2 / data2Square / data[k + i * nObs]) / 
	  R_pow_di(mahalDist[currentPair], 2);

	dvecMixed = log(dvecMixed);

	if (!R_FINITE(dvecMixed)){
	  //printf("dvecMixed = %f\n", dvecMixed);
	  //The likelihood is ill-defined if the mahalanobis distance
	  //is either too small or to large thus
	  
	  dns += R_pow_di(1 + 1/mahalDist[currentPair] +
			  mahalDist[currentPair], 2) * MINF;	  
	}

	else{
	  //Now the final step, multiplying by Fvec and the gradient
	  dns += dvecMixed + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	}
      }
    }
  }
  
  return dns;
}

double lplikschlatherind(double *data, double alpha, double *rho,
			 double *jac, int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schlather model allowing for independence.

  int i, j, k, currentPair = -1;
  double c1, dns, lFvec, dvecM1, dvecM2, dvecMixed, pbm,
    data1Square, data2Square;
  //c1 is a useful quantity - see documentation

  dns = 0.0;
  for (i=0;i<(nSite-1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      
      for (k=0;k<nObs;k++){
	data1Square = R_pow_di(data[k + i * nObs], 2);
	data2Square = R_pow_di(data[k + j * nObs], 2);

	c1 = sqrt(data1Square + data2Square - 2 * data[k + i * nObs] *
		  data[k + j * nObs] * rho[currentPair]);
	
	//It's the log of the joint CDF
	lFvec = (alpha - 1) * (1/data[k + i * nObs] + 1/data[k + j * nObs]) *
	  (1 + sqrt(1 - 2 * (rho[currentPair] + 1) * data[k + i * nObs] *
		    data[k + j * nObs] / 
		    R_pow_di(data[k + i * nObs] + data[k + j * nObs], 2))) / 2 -
	  alpha * (1/data[k + i * nObs] + 1/data[k + j * nObs]);
	
	//It's the partial derivative for marge 1
	dvecM1 = (alpha - 1) *(rho[currentPair] * data[k + i * nObs] - c1 - 
		   data[k + j * nObs]) / 2 / c1 / data1Square + alpha /
	  data1Square;

	//It's the partial derivative for marge 2
	dvecM2 = (alpha - 1) * (rho[currentPair] * data[k + j * nObs] - c1 - 
		   data[k + i * nObs]) / 2 / c1 / data2Square + alpha /
	  data2Square;

	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.

	//It's the mixed partial derivative
	dvecMixed = (1 - alpha) * (1 - R_pow_di(rho[currentPair], 2)) / 2 /
	  R_pow_di(c1, 3) + dvecM1 * dvecM2;

	pbm = (1 + fabs(dvecMixed)) * MINF;
	dvecMixed = log(dvecMixed);

	if (!R_FINITE(dvecMixed)){
	  //printf("dvecMixed is erradic\n");
	  dns += pbm;
	}

	else{
	  //Now the final step, multiplying by Fvec and the gradient
	  dns += dvecMixed + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	}
      }
    }
  }

  return dns;

}
