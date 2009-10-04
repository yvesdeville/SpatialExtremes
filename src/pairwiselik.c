#include "header.h"

double lplikschlather(double *data, double *rho, double *jac,
		      int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schlather model.

  int i, j, k, currentPair = -1;
  double c1, dns = 0.0, lFvec, dvecM1, dvecM2, dvecMixed,
    data1Square, data2Square;
  //c1 is a useful quantity - see documentation

  for (i=0;i<(nSite - 1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      
      for (k=nObs;k--;){
	
	data1Square = data[k + i * nObs] * data[k + i * nObs];
	data2Square = data[k + j * nObs] * data[k + j * nObs];

	c1 = sqrt(data1Square + data2Square - 2 * data[k + i * nObs] *
		  data[k + j * nObs] * rho[currentPair]);
	
	//It's the log of the joint CDF
	lFvec = - (data[k + i * nObs] + data[k + j * nObs] + c1) / 
	  (2 * data[k + i * nObs] * data[k + j * nObs]);
	
	//It's the partial derivative for marge 1
	dvecM1 = -(rho[currentPair] * data[k + i * nObs] - c1 - 
		   data[k + j * nObs]) / (2 * c1 * data1Square);

	//It's the partial derivative for marge 2
	dvecM2 = -(rho[currentPair] * data[k + j * nObs] - c1 - 
		   data[k + i * nObs]) / (2 * c1 * data2Square);

	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.

	//It's the mixed partial derivative
	dvecMixed = (1 - rho[currentPair] * rho[currentPair]) / 
	  (2 * c1 * c1 * c1) + dvecM1 * dvecM2;

	//Now the final step, multiplying by Fvec and the gradient
	dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];

      }
    }
  }

  if (!R_FINITE(dns))
    dns = MINF;
  
  return dns;

}

double lpliksmith(double *data, double *mahalDist, double *jac,
		  int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //smith model.
  
  int i, j, k, currentPair = -1;
  double c1, c2, dns = 0.0, lFvec, dvecM1, dvecM2, dvecMixed,
    dnormc1, dnormc2, pnormc1, pnormc2, idata1, idata2, idata1Square,
    idata2Square, imahal, imahalSquare;

  for (i=0;i<(nSite - 1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      imahal = 1 / mahalDist[currentPair];
      imahalSquare = imahal * imahal;
      
      for (k=nObs;k--;){
	idata1 = 1 / data[k + i * nObs];
	idata2 = 1 / data[k + j * nObs];
	idata1Square = idata1 * idata1;
	idata2Square = idata2 * idata2;

	c1 = log(data[k + j * nObs] * idata1) * imahal + mahalDist[currentPair] / 2;
	c2 = mahalDist[currentPair] - c1;
	
	dnormc1 = dnorm(c1, 0., 1., 0);
	dnormc2 = dnorm(c2, 0., 1., 0);
	pnormc1 = pnorm(c1, 0., 1., 1, 0);
	pnormc2 = pnorm(c2, 0., 1., 1, 0);

	//It's the log of the joint CDF
	lFvec = -pnormc1 * idata1 - pnormc2 * idata2;

	//It's the partial derivative for marge 1
	dvecM1 = (dnormc1 * imahal + pnormc1) * idata1Square - dnormc2 *
	  idata1 * idata2 * imahal;
	
	//It's the partial derivative for marge 2
	dvecM2 = (dnormc2 * imahal + pnormc2) * idata2Square - dnormc1 *
	  idata1 * idata2 * imahal;
	
	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.
	
	//It's the mixed partial derivative
	dvecMixed = dvecM1 * dvecM2 + 
	  (data[k + j * nObs] * c2 * dnormc1 + data[k + i * nObs] * c1 * dnormc2) *
	  idata1Square * idata2Square * imahalSquare;

	//Now the final step, multiplying by Fvec and the gradient
	dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	
      }
    }
  }

  if (!R_FINITE(dns))
    dns = MINF;
  
  return dns;
}

double lplikschlatherind(double *data, double alpha, double *rho,
			 double *jac, int nObs, int nSite){
  //This function computes the log-pairwise likelihood for the
  //schlather model allowing for independence.

  int i, j, k, currentPair = -1;
  double c1, dns = 0.0, lFvec, dvecM1, dvecM2, dvecMixed,
    data1Square, data2Square;
  //c1 is a useful quantity - see documentation

  for (i=0;i<(nSite-1);i++){
    for (j=i+1;j<nSite;j++){
      currentPair++;
      
      for (k=nObs;k--;){
	data1Square = data[k + i * nObs] * data[k + i * nObs];
	data2Square = data[k + j * nObs] * data[k + j * nObs];

	c1 = sqrt(data1Square + data2Square - 2 * data[k + i * nObs] *
		  data[k + j * nObs] * rho[currentPair]);
	
	//It's the log of the joint CDF
	lFvec = ((-alpha - 1) * (data[k + i * nObs] + data[k + j * nObs]) +
		 (alpha - 1) * c1) / (2 * data[k + i * nObs] * data[k + j * nObs]);
	
	//It's the partial derivative for marge 1
	dvecM1 = (alpha - 1) *(rho[currentPair] * data[k + i * nObs] - c1 - 
			       data[k + j * nObs]) / (2 * c1 * data1Square) +
	  alpha / data1Square;

	//It's the partial derivative for marge 2
	dvecM2 = (alpha - 1) * (rho[currentPair] * data[k + j * nObs] - c1 - 
				data[k + i * nObs]) / (2 * c1 * data2Square) +
	  alpha / data2Square;

	//Rmq: to have dvecM1 and dvecM2 we have to multiply
	//them by Fvec[i]. It's not done yet as dvecMixed has to be
	//computed first.

	//It's the mixed partial derivative
	dvecMixed = (1 - alpha) * (1 - rho[currentPair] * rho[currentPair]) / 
	  (2 * c1 * c1 * c1) + dvecM1 * dvecM2;

	//Now the final step, multiplying by Fvec and the gradient
	dns += log(dvecMixed) + lFvec + jac[k + i * nObs] + jac[k + j * nObs];
	
      }
    }
  }

  if (!R_FINITE(dns))
    dns = MINF;

  return dns;

}
