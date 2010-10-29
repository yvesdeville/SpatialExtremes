#include "header.h"

void latentgev(int *n, double *data, int *nSite, int *nObs, int *covmod, 
	       int *dim, double *distMat, double *dsgnMat, int *nBeta, double *beta,
	       double *sills, double *ranges, double *smooths, double *gevParams,
	       double *hyperSill, double *hyperRange, double *hyperSmooth,
	       double *hyperBetaMean, double *hyperBetaIcov, double *propGev,
	       double *propRanges, double *propSmooths, double *mcLoc,
	       double *mcScale, double *mcShape, double *accRates,
	       double *extRates, int *thin, int *burnin){
  
  
  int iter = 0, iterThin = 0, idxSite, idxSite2, idxMarge, idxBeta, info = 0,
    oneInt = 1, zeroInt = 0, nSite2 = *nSite * *nSite,
    nPairs = *nSite * (*nSite + 1) / 2,
    *cumBeta = (int *) R_alloc(3, sizeof(int)),
    *cumBeta2 = (int *) R_alloc(3, sizeof(int)),
    *nBeta2 = (int *) R_alloc(3, sizeof(int)),
    lagLoc = nBeta[0] + 3 + *nSite, lagScale = nBeta[1] + 3 + *nSite,
    lagShape = nBeta[2] + 3 + *nSite;

  cumBeta[0] = 0;
  cumBeta[1] = nBeta[0];
  cumBeta[2] = nBeta[0] + nBeta[1];
  cumBeta2[0] = 0;
  cumBeta2[1] = nBeta[0] * nBeta[0];
  cumBeta2[2] = nBeta[0] * nBeta[0] + nBeta[1] * nBeta[1];
  nBeta2[0] = nBeta[0] * nBeta[0];
  nBeta2[1] = nBeta[1] * nBeta[1];
  nBeta2[2] = nBeta[2] * nBeta[2];  
  
  double one = 1.0, zero = 0.0, flag, logDetProp,
    *logDet = (double *) R_alloc(3, sizeof(double)),
    *icovMat = (double *) R_alloc(3 * nSite2, sizeof(double)),
    *icovMatChol = (double *) R_alloc(3 * nSite2, sizeof(double)),
    *GPmean = (double *) R_alloc(3 * *nSite, sizeof(double)),
    *resTop = (double *) R_alloc(*nSite, sizeof(double)),
    *resBottom = (double *) R_alloc(*nSite, sizeof(double)),
    *covariances = (double *) R_alloc(nPairs, sizeof(double)),
    *proposalGEV = (double *) R_alloc(3, sizeof(double)),
    *dataCopy = (double *) R_alloc(*nObs, sizeof(double)),
    *icovMatProp = (double *) R_alloc(nSite2, sizeof(double)),
    *icovMatPropChol = (double *) R_alloc(nSite2, sizeof(double));

  
  memset(logDet, 0, 3 * sizeof(double));
  memset(icovMat, 0, 3 * nSite2 * sizeof(double));
  memset(GPmean, 0, 3 * *nSite * sizeof(double));
  memset(icovMatProp, 0, nSite2 * sizeof(double));
  memset(icovMatPropChol, 0, nSite2 * sizeof(double));

  /*----------------------------------------------------*/
  //                                                    \\
  //           Compute some initial objects             \\
  //                                                    \\
  /*----------------------------------------------------*/

  // a. The inverse of the covariance matrices
  for (idxMarge=0;idxMarge<3;idxMarge++){

    switch(covmod[idxMarge]){
    case 1:
      flag = whittleMatern(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
			   smooths[idxMarge], covariances);
      break;
    case 2:
      flag = cauchy(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
		    smooths[idxMarge], covariances);
      break;
    case 3:
      flag = powerExp(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
		      smooths[idxMarge], covariances);
      break;
    case 4:
      flag = bessel(distMat, nPairs, *dim, sills[idxMarge], ranges[idxMarge],
		    smooths[idxMarge], covariances);
      break;
    }

    if (flag != 0)
      error("The starting values (covariance parameter) are ill-defined. Please check\n");

    /* We need to fill in the upper triangular part of icovMat with
       covariances */
    {
      int current=-1;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
	  current++;
	  icovMat[idxSite + idxSite2 * *nSite + idxMarge * nSite2] = covariances[current];
	}
    }

    // Compute the Cholesky decomposition of the covariance matrices
    F77_CALL(dpotrf)("U", nSite, icovMat + idxMarge * nSite2, nSite, &info);
    
    if (flag != 0)
      error("Impossible to get the Cholesky decomp. from the starting values\n");

    /* Compute the log of the determinant of the proposal
       cov. mat. using the sum of the square of the diagonal elements of
       the Cholesky decomposition */
    for (idxSite2=0;idxSite2<*nSite;idxSite2++)
      logDet[idxMarge] += 2 * log(icovMat[idxSite2 * (*nSite + 1) + idxMarge *
					  nSite2]);
          
    /* Compute the inverse of this Cholesky decomposition. 
    
       WARNING: At this stage icovMat is the inverse of the covariance
       matrices but only the upper diagonal elements are stored */
    F77_CALL(dpotri)("U", nSite, icovMat + idxMarge * nSite2, nSite, &info);

    /* Now we compute its Cholesky decomposition as we will use it a
       lot */
    memcpy(icovMatChol + idxMarge * nSite2, icovMat + idxMarge * nSite2,
	   nSite2 * sizeof(double));
    F77_CALL(dpotrf)("U", nSite, icovMatChol + idxMarge * nSite2, nSite, &info);
  }

  // b. The mean of the Gaussian processes
  for (idxMarge=0;idxMarge<3;idxMarge++)
    for (idxSite=0;idxSite<*nSite;idxSite++)
      for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	GPmean[idxSite + idxMarge * *nSite] += 
	  dsgnMat[idxBeta * *nSite + idxSite + cumBeta[idxMarge] * *nSite] *
	  beta[cumBeta[idxMarge] + idxBeta];
  

  /*----------------------------------------------------*/
  //                                                    \\
  //               Starting the MCMC algo               \\
  //                                                    \\
  /*----------------------------------------------------*/

  GetRNGstate();
  while (iterThin<*n){
    
    /*----------------------------------------------------*/
    //                                                    \\
    //           Updating the GEV parameters              \\
    //                                                    \\
    /*----------------------------------------------------*/
    
    for (idxSite=0;idxSite<*nSite;idxSite++){
      double topGEV = 0, bottomGEV = 0, topGP = 0, bottomGP = 0;
      
      for (idxMarge=0;idxMarge<3;idxMarge++)
	proposalGEV[idxMarge] = gevParams[idxMarge * *nSite + idxSite] +
	  propGev[idxMarge] * (unif_rand() - 0.5);

      gevlik(data + idxSite * *nObs, nObs, proposalGEV, proposalGEV + 1,
	     proposalGEV + 2, &topGEV);

      if (topGEV == -1e6){
	extRates[0]++;
	continue;
      }

      gevlik(data + idxSite * *nObs, nObs, gevParams + idxSite, gevParams +
	     *nSite + idxSite, gevParams + 2 * *nSite + idxSite, &bottomGEV);

      for (idxMarge=0;idxMarge<3;idxMarge++){
	for (idxSite2=0;idxSite2<*nSite;idxSite2++)
	  resBottom[idxSite2] = gevParams[idxSite2 + idxMarge * *nSite] -
	    GPmean[idxSite2 + idxMarge * *nSite];
	
	memcpy(resTop, resBottom, *nSite * sizeof(double));
	resTop[idxSite] = proposalGEV[idxMarge] - GPmean[idxSite + idxMarge *
							 *nSite];

	F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatChol + idxMarge * nSite2,
			nSite, resTop, &oneInt);
	F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatChol + idxMarge * nSite2,
			nSite, resBottom, &oneInt);

	for (idxSite2=0;idxSite2<*nSite;idxSite2++){
	  topGP += resTop[idxSite2] * resTop[idxSite2];
	  bottomGP += resBottom[idxSite2] * resBottom[idxSite2];
	}	
      }

      topGP *= -0.5;
      bottomGP *= -0.5;

      if (unif_rand() < exp(topGEV - bottomGEV + topGP - bottomGP)){
	for (idxMarge=0;idxMarge<3;idxMarge++)
	  gevParams[idxSite + idxMarge * *nSite] = proposalGEV[idxMarge];
	  
	accRates[0]++;
      }
    }
  	    

    /*----------------------------------------------------*/
    //                                                    \\
    //        Updating the regression parameters          \\
    //                (conjugate prior)                   \\ 
    //                                                    \\
    /*----------------------------------------------------*/

    for (idxMarge=0;idxMarge<3;idxMarge++){
      
      /* conjCovMat is the covariance matrix for the conjugate
	 distribution i.e. MVN

	 conjCovMatChol is its Cholesky decomposition */
      double *dummy = malloc(*nSite * nBeta[idxMarge] * sizeof(double)),
	*conjCovMat = malloc(nBeta2[idxMarge] * sizeof(double)),
	*conjCovMatChol = malloc(nBeta2[idxMarge] * sizeof(double));
      
      memcpy(conjCovMat, hyperBetaIcov + cumBeta2[idxMarge],
	     nBeta2[idxMarge] * sizeof(double));
      memcpy(dummy, dsgnMat + *nSite * cumBeta[idxMarge],
	     *nSite * nBeta[idxMarge] * sizeof(double));
      
      // Compute dummy = icovMatChol %*% dsgnMat
      F77_CALL(dtrmm)("L", "U", "N", "N", nSite, nBeta + idxMarge, &one,
		      icovMatChol + idxMarge * nSite2, nSite, dummy, nSite);
      
      /* Compute conjCovMat = dummy^T %*% dummy + conjCovMat
	 
	 WARNING: Only the upper diagonal elements will be stored */
      F77_CALL(dsyrk)("U", "T", nBeta + idxMarge, nSite, &one, dummy, nSite,
		      &one, conjCovMat, nBeta + idxMarge);

      /* Compute the inverse of the conjCovMat
	 
	 WARNING: This will be  the inverse but only the upper diagonal
	 elements are stored. */
      F77_CALL(dpotrf)("U", nBeta + idxMarge, conjCovMat, nBeta + idxMarge,
		       &info);
      F77_CALL(dpotri)("U", nBeta + idxMarge, conjCovMat, nBeta + idxMarge,
		       &info);

      // Compute its Cholesky decomposition
      memcpy(conjCovMatChol, conjCovMat, nBeta2[idxMarge] * sizeof(double));
      
      F77_CALL(dpotrf)("U", nBeta + idxMarge, conjCovMatChol, nBeta + idxMarge,
		       &info);
      
      // Compute dummy2 = icovMatChol %*% (locs or scales or shapes)
        double *dummy2 = malloc(*nSite * sizeof(double));
      memcpy(dummy2, gevParams + idxMarge * *nSite, *nSite * sizeof(double));
      
      F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatChol + idxMarge * nSite2,
		      nSite, dummy2, &oneInt);

      //Compute dummy3 = hyperBetaIcov %*% hyperBetaMean
      double *dummy3 = malloc(nBeta[idxMarge] * sizeof(double));
      memset(dummy3, 0, nBeta[idxMarge] * sizeof(double));
      F77_CALL(dsymv)("U", nBeta + idxMarge, &one, hyperBetaIcov +
		      cumBeta2[idxMarge], nBeta + idxMarge, hyperBetaMean +
		      cumBeta[idxMarge], &oneInt, &zero, dummy3, &oneInt);
      
      // Compute dummy3 = dummy3 + dummy^T %*% dummy2 (dummy2 is a vector)
      F77_CALL(dgemv)("T", nSite, nBeta + idxMarge, &one, dummy, nSite, dummy2,
		      &oneInt, &one, dummy3, &oneInt);
      
      // conjMean is the mean for the conjugate distribution i.e. MVN
      double *conjMean = malloc(nBeta[idxMarge] * sizeof(double));
      memset(conjMean, 0, nBeta[idxMarge] * sizeof(double));
      // Compute conjMean = conjCovMat %*% dummy3
      F77_CALL(dsymv)("U", nBeta  + idxMarge, &one, conjCovMat, nBeta +
		      idxMarge, dummy3, &oneInt, &zero, conjMean, &oneInt);

      /* The new state is a realisation from the MVN(conjMean,
	 conjCovMat) so we simulate it from the Cholesky
	 decomposition */
      
      double *stdNormal = malloc(nBeta[idxMarge] * sizeof(double));
      for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	stdNormal[idxBeta] = norm_rand();
      
      F77_CALL(dtrmv)("U", "T", "N", nBeta + idxMarge, conjCovMatChol,
		      nBeta + idxMarge, stdNormal, &oneInt);
      
      for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	beta[cumBeta[idxMarge] + idxBeta] = stdNormal[idxBeta] +
	  conjMean[idxBeta];
      
      //The last step is to update the mean of the GP
      for (idxSite=0;idxSite<*nSite;idxSite++){
	GPmean[idxSite + idxMarge * *nSite] = 0;
	
	for (idxBeta=0;idxBeta<nBeta[idxMarge];idxBeta++)
	  GPmean[idxSite + idxMarge * *nSite] += dsgnMat[idxBeta * *nSite + idxSite +
							 cumBeta[idxMarge] * *nSite] *
	    beta[cumBeta[idxMarge] + idxBeta];
      }

      free(dummy);
      free(conjCovMat);
      free(conjCovMatChol);
      free(dummy2);
      free(dummy3);
      free(conjMean);
      free(stdNormal);
    }


    /*----------------------------------------------------*/
    //                                                    \\
    //        Updating the sills (conjugate prior)        \\
    //                                                    \\
    /*----------------------------------------------------*/

    for (idxMarge=0;idxMarge<3;idxMarge++){
      for (idxSite=0;idxSite<*nSite;idxSite++)
	resTop[idxSite] = gevParams[idxSite + idxMarge * *nSite] -
	  GPmean[idxSite + idxMarge * *nSite];

      // Compute resTop = icovMatChol %*% resTop
      F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatChol + idxMarge * nSite2,
		      nSite, resTop, &oneInt); 

      double shape = 0.5 * *nSite + hyperSill[2 * idxMarge];
      double scale = hyperSill[1 + 2 * idxMarge];
      for (idxSite=0;idxSite<*nSite;idxSite++)
	scale += 0.5 * sills[idxMarge] * resTop[idxSite] * resTop[idxSite];

      /* Rmk: If Y ~ Gamma(shape = shape, rate = 1 / scale) then X :=
	 1 / Y \sim IGamma(shape = shape, scale = scale) */
      sills[idxMarge] = 1 / rgamma(shape,  1 / scale);

      // Now we need to update the covariance matrix and its inverse
      switch(covmod[idxMarge]){
      case 1:
	flag = whittleMatern(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
			     smooths[idxMarge], covariances);
	break;
      case 2:
	flag = cauchy(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
		      smooths[idxMarge], covariances);
	break;
      case 3:
	flag = powerExp(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
			smooths[idxMarge], covariances);
	break;
      case 4:
	flag = bessel(distMat, nPairs, *dim, sills[idxMarge], ranges[idxMarge],
		      smooths[idxMarge], covariances);
	break;
      }

      /* We need to fill in the upper triangular part of icovMat with
	 covariances */
      {
	int current=-1;
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
	    current++;
	    icovMat[idxSite + idxSite2 * *nSite + idxMarge * nSite2] = covariances[current];
	  }
      }

      // Cholesky decomposition of the covariance matrices
      F77_CALL(dpotrf)("U", nSite, icovMat + idxMarge * nSite2, nSite, &info);

      // Compute the log of the determinant of the proposal cov. mat.
      logDet[idxMarge] = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	logDet[idxMarge] += 2 * log(icovMat[idxSite * (1 + *nSite) + idxMarge *
					    nSite2]);

      // Inverse of the proposal cov. mat 
      F77_CALL(dpotri)("U", nSite, icovMat + idxMarge * nSite2, nSite, &info);

      // And its Cholesky decomposition
      memcpy(icovMatChol + idxMarge * nSite2, icovMat + idxMarge * nSite2,
	     nSite2 * sizeof(double));
      F77_CALL(dpotrf)("U", nSite, icovMatChol + idxMarge * nSite2, nSite,
		       &info);
    }


    /*----------------------------------------------------*/
    //                                                    \\
    //          Updating the ranges (M.-H. step)          \\
    //                                                    \\
    /*----------------------------------------------------*/
    
    for (idxMarge=0;idxMarge<3;idxMarge++){
      if (propRanges[idxMarge] == 0)
	continue;

      //double rangeProp = ranges[idxMarge] + propRanges[idxMarge] *
      //(unif_rand() - 0.5);
      double rangeProp = rlnorm(log(ranges[idxMarge]), propRanges[idxMarge]),
	logpropRatio = log(rangeProp / ranges[idxMarge]);

      switch(covmod[idxMarge]){
      case 1:
	flag = whittleMatern(distMat, nPairs, sills[idxMarge], rangeProp,
			     smooths[idxMarge], covariances);
	break;
      case 2:
	flag = cauchy(distMat, nPairs, sills[idxMarge], rangeProp,
		      smooths[idxMarge], covariances);
	break;
      case 3:
      flag = powerExp(distMat, nPairs, sills[idxMarge], rangeProp,
		      smooths[idxMarge], covariances);
      break;
      case 4:
	flag = bessel(distMat, nPairs, *dim, sills[idxMarge], rangeProp,
		      smooths[idxMarge], covariances);
	break;
      }
      
      if (flag != 0){
	extRates[1 + idxMarge]++;
	continue;
      }
      
      /* We need to fill in the upper triangular part of icovMatProp
	 with covariances */
      {
	int current=-1;
	for (idxSite=0;idxSite<*nSite;idxSite++)
	  for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
	    current++;
	    icovMatProp[idxSite + idxSite2 * *nSite] = covariances[current];
	  }
      }

      // Cholesky decomposition of the proposal cov. mat.
      F77_CALL(dpotrf)("U", nSite, icovMatProp, nSite, &info);

      if (info != 0){
	extRates[1 + idxMarge]++;
	continue;
      }

      // Log of the determinant of the proposal cov. mat.
      logDetProp = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
	logDetProp += 2 * log(icovMatProp[idxSite * (1 + *nSite)]);

      // Inverse of the proposal cov. mat.
      F77_CALL(dpotri)("U", nSite, icovMatProp, nSite, &info);

      // And its Cholesky decomposition
      memcpy(icovMatPropChol, icovMatProp, nSite2 * sizeof(double));
      F77_CALL(dpotrf)("U", nSite, icovMatPropChol, nSite, &info);

      if (info != 0){
	extRates[1 + idxMarge]++;
	continue;
      }

      for (idxSite=0;idxSite<*nSite;idxSite++)
	resBottom[idxSite] = gevParams[idxSite + idxMarge * *nSite] -
	  GPmean[idxSite + idxMarge * *nSite];

      memcpy(resTop, resBottom, *nSite * sizeof(double));
      
      F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatChol + idxMarge * nSite2,
		      nSite, resBottom, &oneInt);
      F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatPropChol, nSite, resTop,
		      &oneInt);

      double top = logDetProp, bottom = logDet[idxMarge],
	logpriorRatio = (hyperRange[2 * idxMarge] - 1) * 
	log(rangeProp / ranges[idxMarge]) + (ranges[idxMarge] - rangeProp) /
	hyperRange[2 * idxMarge + 1];

      for (idxSite=0;idxSite<*nSite;idxSite++){
	top += resTop[idxSite] * resTop[idxSite];
	bottom += resBottom[idxSite] * resBottom[idxSite];
      }

      top *= -0.5;
      bottom *= -0.5;

      if (unif_rand() < exp(top - bottom + logpriorRatio + logpropRatio)){
	ranges[idxMarge] = rangeProp;
	logDet[idxMarge] = logDetProp;
	memcpy(icovMat + idxMarge * nSite2, icovMatProp, nSite2 *
	       sizeof(double));
	memcpy(icovMatChol + idxMarge * nSite2, icovMatPropChol, nSite2 *
	       sizeof(double));
	accRates[1 + idxMarge]++;
      }    
    }

    /*----------------------------------------------------*/
    //                                                    \\
    //         Updating the smooths (M.-H. step)          \\
    //                                                    \\
    /*----------------------------------------------------*/
    
    for (idxMarge=0;idxMarge<3;idxMarge++){
      if (propSmooths[idxMarge] == 0)
	continue;

      //double smoothProp = smooths[idxMarge] + propSmooths[idxMarge] *
      //(unif_rand() - 0.5);
      double smoothProp = rlnorm(log(smooths[idxMarge]), propSmooths[idxMarge]),
	logpropRatio = log(smoothProp / smooths[idxMarge]);
    
      switch(covmod[idxMarge]){
      case 1:
	flag = whittleMatern(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
			     smoothProp, covariances);
	break;
      case 2:
	flag = cauchy(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
		      smoothProp, covariances);
	break;
      case 3:
	flag = powerExp(distMat, nPairs, sills[idxMarge], ranges[idxMarge],
			smoothProp, covariances);
	break;
      case 4:
	flag = bessel(distMat, nPairs, *dim, sills[idxMarge], ranges[idxMarge],
		      smoothProp, covariances);
	break;
      }
            
      if (flag != 0){
    	extRates[4 + idxMarge]++;
    	continue;
      }
      
      /* We need to fill in the upper triangular part of icovMatProp
    	 with covariances */
      {
    	int current=-1;
    	for (idxSite=0;idxSite<*nSite;idxSite++)
    	  for (idxSite2=idxSite;idxSite2<*nSite;idxSite2++){
    	    current++;
    	    icovMatProp[idxSite + idxSite2 * *nSite] = covariances[current];
    	  }
      }
    
      // Cholesky decomposition of the proposal cov. mat.
      F77_CALL(dpotrf)("U", nSite, icovMatProp, nSite, &info);
    
      if (info != 0){
    	extRates[4 + idxMarge]++;
    	continue;
      }
    
      // Log of the determinant of the proposal cov. mat.
      logDetProp = 0;
      for (idxSite=0;idxSite<*nSite;idxSite++)
    	logDetProp += 2 * log(icovMatProp[idxSite * (1 + *nSite)]);
    
      // Inverse of the proposal cov. mat.
      F77_CALL(dpotri)("U", nSite, icovMatProp, nSite, &info);
    
      // And its Cholesky decomposition
      memcpy(icovMatPropChol, icovMatProp, nSite2 * sizeof(double));
      F77_CALL(dpotrf)("U", nSite, icovMatPropChol, nSite, &info);
    
      if (info != 0){
    	extRates[4 + idxMarge]++;
    	continue;
      }
    
      for (idxSite=0;idxSite<*nSite;idxSite++)
    	resBottom[idxSite] = gevParams[idxSite + idxMarge * *nSite] -
    	  GPmean[idxSite + idxMarge * *nSite];
    
      memcpy(resTop, resBottom, *nSite * sizeof(double));
      
      F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatChol + idxMarge * nSite2,
    		      nSite, resBottom, &oneInt);
      F77_CALL(dtrmv)("U", "N", "N", nSite, icovMatPropChol, nSite, resTop,
    		      &oneInt);
    
      double top = logDetProp, bottom = logDet[idxMarge],
    	logpriorRatio = (hyperSmooth[2 * idxMarge] - 1) * 
    	log(smoothProp / smooths[idxMarge]) + (smooths[idxMarge] - smoothProp) /
    	hyperSmooth[2 * idxMarge + 1];
    
      for (idxSite=0;idxSite<*nSite;idxSite++){
    	top += resTop[idxSite] * resTop[idxSite];
    	bottom += resBottom[idxSite] * resBottom[idxSite];
      }
    
      top *= -0.5;
      bottom *= -0.5;
    
      if (unif_rand() < exp(top - bottom + logpriorRatio + logpropRatio)){
    	smooths[idxMarge] = smoothProp;
    	logDet[idxMarge] = logDetProp;
    	memcpy(icovMat + idxMarge * nSite2, icovMatProp, nSite2 *
    	       sizeof(double));
    	memcpy(icovMatChol + idxMarge * nSite2, icovMatPropChol, nSite2 *
    	       sizeof(double));
    	accRates[4 + idxMarge]++;
      }    
    }

    iter++;

    //Need to store the new state into the mc object.
    if ((iter > *burnin) & ((iter % *thin) == 0)){
      mcLoc[nBeta[0] + iterThin * lagLoc] = sills[0];
      mcLoc[nBeta[0] + 1 + iterThin * lagLoc] = ranges[0];
      mcLoc[nBeta[0] + 2 + iterThin * lagLoc] = smooths[0];
      
      mcScale[nBeta[1] + iterThin * lagScale] = sills[1];
      mcScale[nBeta[1] + 1 + iterThin * lagScale] = ranges[1];
      mcScale[nBeta[1] + 2 + iterThin * lagScale] = smooths[1];
      
      mcShape[nBeta[2] + iterThin * lagShape] = sills[2];
      mcShape[nBeta[2] + 1 + iterThin * lagShape] = ranges[2];
      mcShape[nBeta[2] + 2 + iterThin * lagShape] = smooths[2];
      
      for (idxBeta=0;idxBeta<nBeta[0];idxBeta++)
	mcLoc[idxBeta + iterThin * lagLoc] = beta[idxBeta];
      
      for (idxBeta=0;idxBeta<nBeta[1];idxBeta++)
	mcScale[idxBeta + iterThin * lagScale] = beta[cumBeta[1] + idxBeta];
      
      for (idxBeta=0;idxBeta<nBeta[2];idxBeta++)
	mcShape[idxBeta + iterThin * lagShape] = beta[cumBeta[2] + idxBeta];
      
      for (idxSite=0;idxSite<*nSite;idxSite++){
	mcLoc[nBeta[0] + 3 + idxSite + iterThin * lagLoc] = gevParams[idxSite];
	mcScale[nBeta[1] + 3 + idxSite + iterThin * lagScale] = gevParams[*nSite + idxSite];
	mcShape[nBeta[2] + 3 + idxSite + iterThin * lagShape] = gevParams[2 * *nSite + idxSite];
      }
      iterThin++;
    }
  }
  GetRNGstate();
  
  for (int i=7;i--;){
    accRates[i] /= (double) iter;
    extRates[i] /= (double) iter;
  }

  return;
}
  
void DIC(int *nChain, int *nSite, int *nObs, double *data, double *chainLoc,
	 double *chainScale, double *chainShape, double *postLoc,
	 double *postScale, double *postShape, double *dic, double *effNpar,
	 double *dbar){

  int i,j;
  for (i=*nChain;i--;)
    for (j=*nSite;j--;){
      double dummy = 0;
      gevlik(data + j * *nObs, nObs, chainLoc + j * *nChain + i,
	     chainScale + j * *nChain + i, chainShape + j * *nChain + i,
	     &dummy);
      *dbar += dummy;
    }

  *dbar *= -2 / ((double) *nChain);
  
  for (i=*nSite;i--;){
    double dummy = 0;
    gevlik(data + i * *nObs, nObs, postLoc + i, postScale + i, postShape + i,
	   &dummy);
    *effNpar += dummy;
  }

  *effNpar = *dbar + 2 * *effNpar;
  *dic = *effNpar + *dbar;

  return;
}
  
