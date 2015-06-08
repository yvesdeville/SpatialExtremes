#include "header.h"



void empiricalConcProb(double *data, int *nSite, int *nObs, int *blockSize,
		       int *nBlock, double *concProb){
  /* This function compputes the pairwise concurence probability
     assuming the data are in the max domain of attraction */
  
  double *dataBlock = malloc(*blockSize * *nSite * sizeof(double));
  // will store the data corresponding to a given block

  for (int k=0;k<*nBlock;k++){//k is the block number
    double blockMax[2];//will store the pairwise block maxima
    int currentPair = 0,
      blockMaxIdx[2];//will store the indices for the pairwise block maxima

    for (int i=0;i<*blockSize;i++)
      for (int j=0;j<*nSite;j++)
	dataBlock[i + j * *blockSize] = data[k * *blockSize + i + j * *nObs];

    for (int i=0;i<(*nSite-1);i++){
      // Compute the block maxima for site i (and its index) and block number k
      blockMax[0] = dataBlock[i * *blockSize];
      blockMaxIdx[0] = 0;

      for (int l=1;l<*blockSize;l++){
	if (blockMax[0] < dataBlock[l + i * *blockSize]){
	  blockMax[0] = dataBlock[l + i * *blockSize];
	  blockMaxIdx[0] = l;
	}
      }
      
      for (int j=i+1;j<*nSite;j++){
	// Compute the block maxima for site j (and its index) and block number k
	blockMax[1] = dataBlock[j * *blockSize];
	blockMaxIdx[1] = 0;

	for (int l=1;l<*blockSize;l++){
	  if (blockMax[1] < dataBlock[l + j * *blockSize]){
	    blockMax[1] = dataBlock[l + j * *blockSize];
	    blockMaxIdx[1] = l;
	  }
	}

	if (blockMaxIdx[0] == blockMaxIdx[1])
	  // Concurent extremes
	  concProb[currentPair]++;

	currentPair++;
      }
    }
  }

  for (int i=0;i<(*nSite * (*nSite - 1) / 2);i++)
    concProb[i] /= (double) *nBlock;      

  free(dataBlock);
  return;
}

void concProbKendall(double *data, int *nSite, int *nObs, double *concProb){
  
  int currentPair=0;
  for (int i=0;i<(*nSite-1);i++){
    for (int j=i+1;j<*nSite;j++){

      concProb[currentPair]=0;
      int nEffObs = 0;//the number of effective contribution (because of possible missing values)
      for (int k=0;k<(*nObs-1);k++){

	if (ISNA(data[k + i * *nObs]) || ISNA(data[k + j * *nObs]))
	  continue;
	
	for (int l=k+1;l<*nObs;l++){

	  if (ISNA(data[l + i * *nObs]) || ISNA(data[l + j * *nObs]))
	    continue;
	  
	  nEffObs++;
	  concProb[currentPair] += sign(data[k + i * *nObs] - data[l + i * *nObs]) *
	    sign(data[k + j * *nObs] - data[l + j * *nObs]);
	}
      }

      concProb[currentPair] *= nEffObs == 0 ? NA_REAL : 1.0 / ((double) nEffObs);
      currentPair++;
    }
  }

  return;
}


void empiricalBootConcProb(double *data, int *nSite, int *nObs, int *blockSize,
			   double *concProb){

  const double normCst = lchoose(*nObs, *blockSize);
  int currentPair=0;
  for (int i=0;i<(*nSite-1);i++){
    for (int j=i+1; j<*nSite; j++){
      // For each pair compute the estimator
      concProb[currentPair] = 0;

      for (int k=0;k<*nObs;k++){
	int d = 0;
	
	for (int l=0;l<*nObs;l++)
	  d += (data[l + i * *nObs] < data[k + i * *nObs]) && (data[l + j * *nObs] < data[k + j * *nObs]);

	concProb[currentPair] += exp(lchoose(d, *blockSize - 1) - normCst);
      }

      currentPair++;
    }
  }
}
