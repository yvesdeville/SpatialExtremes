#include "header.h"

void rschlathertbm(double *coord, int *nObs, int *nSite, int *dim,
		   int *covmod, int *grid, double *sill, double *range,
		   double *smooth, double *uBound, int *nlines,
		   double *ans){
  /* This function generates random fields from the Schlather model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      sill: the sill parameter - (1 - sill) is a nugget effect
     range: the range parameter
    smooth: the smooth parameter
    uBound: the uniform upper bound for the stoch. proc.
    nlines: the number of lines used for the TBM algo
       ans: the generated random field */

  int i;
  double *lines;

  //rescale the coordinates
  for (i=(*nSite * *dim);i--;){
    const double irange = 1 / *range;
    coord[i] = coord[i] * irange;
  }

  lines = (double *)R_alloc(3 * *nlines, sizeof(double));
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);
  
  GetRNGstate();
  if (*grid){
    //coord defines a grid
    int neffSite = R_pow_di(*nSite, *dim);
    for (i=*nObs;i--;){
      int nKO = neffSite;
      double poisson = 0;

      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	int j;
	double nugget = 1 - *sill, ipoisson, u, v, w,
	  angle, inorm, thresh, *gp;

	gp = (double *)R_alloc(neffSite, sizeof(double));

	/* ------- Random rotation of the lines ----------*/
	u = unif_rand() - 0.5;
	v = unif_rand() - 0.5;
	w = unif_rand() - 0.5;
	angle = runif(0, M_2PI);
	
	inorm = 1 / sqrt(u * u + v * v + w * w);
	
	u *= inorm;
	v *= inorm;
	w *= inorm;
	
	rotation(lines, nlines, &u, &v, &w, &angle);
	/* -------------- end of rotation ---------------*/
	
	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	
	for (j=neffSite;j--;)
	  gp[j] = 0;

	tbmcore(nSite, &neffSite, dim, covmod, grid, coord, &nugget,
		sill, range, smooth, nlines, lines, gp);
	
	nKO = neffSite;
	for (j=neffSite;j--;){
	  if (thresh > ans[j + i * neffSite])
	    ans[j + i * neffSite] = fmax2(gp[j] * ipoisson, ans[j + i * neffSite]);
	  
	  else
	    nKO--;
	  
	}
      }
    }

    //Lastly we multiply by 1 / E[max(0, Y)]
    for (i=(neffSite * *nObs);i--;){
      const double imean = M_SQRT2 * M_SQRT_PI;
      ans[i] *= imean;
    }
  }

  else{
    //coord doesn't define a grid
    int neffSite = *nSite;
    for (i=*nObs;i--;){
      double poisson = 0;
      int nKO = neffSite;
      
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	int j;
	double nugget = 1 - *sill, ipoisson, u, v, w,
	  angle, inorm, thresh, *gp;

	gp = (double *)R_alloc(neffSite, sizeof(double));

	/* ------- Random rotation of the lines ----------*/
	u = unif_rand() - 0.5;
	v = unif_rand() - 0.5;
	w = unif_rand() - 0.5;
	angle = runif(0, M_2PI);
	
	inorm = 1 / sqrt(u * u + v * v + w * w);
	
	u *= inorm;
	v *= inorm;
	w *= inorm;
	
	rotation(lines, nlines, &u, &v, &w, &angle);
	/* -------------- end of rotation ---------------*/

	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	for (j=neffSite;j--;)
	  gp[j] = 0;
	
	tbmcore(nSite, &neffSite, dim, covmod, grid, coord, &nugget,
		sill, range, smooth, nlines, lines, gp);
	
	nKO = neffSite;
	for (j=*nSite;j--;){
	  if (thresh > ans[i + j * *nObs])
	    ans[i + j * *nObs] = fmax2(gp[j] * ipoisson, ans[i + j * *nObs]);
	  
	  else
	    nKO--;
	}
      }
    }

    //Lastly we multiply by 1 / E[max(0, Y)]
    for (i=(neffSite * *nObs);i--;){
      const double imean = M_SQRT2 * M_SQRT_PI;
      ans[i] *= imean;
    }
  }

  PutRNGstate();

  return;
}

void rschlatherdirect(double *coord, int *nObs, int *nSite, int *dim,
		      int *covmod, int *grid, double *sill, double *range,
		      double *smooth, double *uBound, double *ans){
  /* This function generates random fields for the Schlather model

     coord: the coordinates of the locations
      nObs: the number of observations to be generated
    nSite: the number of locations
       dim: the random field is generated in R^dim
    covmod: the covariance model
      grid: Does coord specifies a grid?
      sill: the sill parameter - (1 - sill) is a nugget effect
     range: the range parameter
    smooth: the smooth parameter
       ans: the generated random field */

  int i, j, k, lwork, nKO, info = 0, neffSite;
  double poisson, ipoisson, thresh, *gp, nugget = 1 - *sill,
    *covmat, one = 1, zero = 0, *work, *xvals, tmp, *d, *u,
    *v, sum, dummy;

  if (*grid)
    neffSite = R_pow_di(*nSite, *dim);

  else
    neffSite = *nSite;

  covmat = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  d = (double *)R_alloc(neffSite, sizeof(double));
  u = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  v = (double *)R_alloc(neffSite * neffSite, sizeof(double));
  xvals = (double *) R_alloc(neffSite * neffSite, sizeof(double));
  gp = (double *)R_alloc(neffSite, sizeof(double));

  buildcovmat(nSite, grid, covmod, coord, dim, &nugget, sill, range,
	      smooth, covmat);
  
  /* Compute the singular value decomposition of the covariance
     matrix.

     This piece of code is strongly inspired from Lapack.c */
  
  Memcpy(xvals, covmat, neffSite * neffSite);
  
  {
    int *iwork= (int *) R_alloc(8 * *nSite, sizeof(int));
    
    /* ask for optimal size of work array */
    lwork = -1;
    F77_CALL(dgesdd)("A", &neffSite, &neffSite, xvals, &neffSite, d, u,
		     &neffSite, v, &neffSite, &tmp, &lwork, iwork, &info);
    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dgesdd");

    lwork = (int) tmp;
    work = (double *) R_alloc(lwork, sizeof(double));

    F77_CALL(dgesdd)("A", &neffSite, &neffSite, xvals, &neffSite, d, u,
		     &neffSite, v, &neffSite, work, &lwork, iwork, &info);
    if (info != 0)
      error("error code %d from Lapack routine '%s'", info, "dgesdd");
  }

  /*--------------- end of singular value decomposition ---------------*/

  /* Compute the square root of the covariance matrix */
  // a) First compute diag(sqrt(d)) %*% u
  for (i=0;i<neffSite;i++){
    dummy = sqrt(d[i]);
    
    for (j=0;j<neffSite;j++)
      u[i + neffSite * j] *= dummy;
  }

  // b) Then compute v^T %*% diag(sqrt(d)) %*% u and put it in covmat
  F77_CALL(dgemm)("T", "N", &neffSite, &neffSite, &neffSite, &one,
		  v, &neffSite, u, &neffSite, &zero, covmat, &neffSite);
  
  GetRNGstate();
  if (*grid){
    //coord defines a grid
    for (i=*nObs;i--;){
      poisson = 0;
      nKO = neffSite;
      
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	for (j=neffSite;j--;)
	  d[j] = norm_rand();
	  
	for (j=neffSite;j--;){
	  sum = 0;
	  for (k=neffSite;k--;)
	    sum += d[k] * covmat[j + k * neffSite];
	  
	  gp[j] = sum;
	}
	
	nKO = neffSite;
	for (j=neffSite;j--;){
	  if (thresh > ans[j + i * neffSite])
	    ans[j + i * neffSite] = fmax2(gp[j] * ipoisson, ans[j + i * neffSite]);
	  
	  else
	    nKO--;
	  
	}
      }
    }
  }

  else{
    //coord doesn't define a grid
    for (i=*nObs;i--;){
      poisson = 0;
      nKO = *nSite;
      
      while (nKO) {
	/* The stopping rule is reached when nKO = 0 i.e. when each site
	   satisfies the condition in Eq. (8) of Schlather (2002) */
	poisson += exp_rand();
	ipoisson = 1 / poisson;
	thresh = *uBound * ipoisson;
	
	/* We simulate one realisation of a gaussian random field with
	   the required covariance function */
	for (j=neffSite;j--;)
	  d[j] = norm_rand();
	  
	for (j=neffSite;j--;){
	  sum = 0;
	  for (k=neffSite;k--;)
	    sum += d[k] * covmat[j + k * neffSite];
	  
	  gp[j] = sum;
	}
	
	nKO = *nSite;
	for (j=*nSite;j--;){
	  if (thresh > ans[i + j * *nObs])
	    ans[i + j * *nObs] = fmax2(gp[j] * ipoisson, ans[i + j * *nObs]);
	  
	  else
	    nKO--;
	}
      }
    }
  }

  PutRNGstate();
  //Lastly we multiply by 1 / E[max(0, Y)]
  for (i=(neffSite * *nObs);i--;){
    const double imean = M_SQRT2 * M_SQRT_PI;
    ans[i] *= imean;
  }
  
  return;
}

void tbmcore(int *nsite, int *neffSite, int *dim, int *covmod,
	     int *grid, double *coord, double *nugget, double *sill,
	     double *range, double *smooth, int *nlines, double *lines,
	     double *ans){
  /* This function is the same as the tbm function except that it
     generates only one realisation of the random field. This is only
     useful for the call by rschlathertbm - for CPU reasons. */

  int j;
  
  if (*grid){
    //coord defines a grid
    int k, l, m;
    double freq, eucProd, u1, u2, G, phase;
    double cl, lj, cl1, cl2;
    double halfSmooth = 0.5 * *smooth, ismooth = 1 / *smooth;
    double u3;

    switch (*dim){
    case 2:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd =  cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 2:
	//Cauchy
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3) * rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 3:
	//Powered exponential
	for (j=*nlines;j--;){
	  u1 = rexp(1);
	  u2 = runif(-M_PI_2, M_PI_2);
	  u3 = halfSmooth * (u2 - M_PI_2);

	  G = fabs(sin(u3) * R_pow(cos(u2), -ismooth) *
		   R_pow(cos(u2 - u3) / u1, (2 - *smooth) * ismooth));
	  
	  freq = sqrt(rchisq(3) * M_SQRT_3 / G);
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 4:
	//Bessel
	for (j=*nlines;j--;){
	  freq = sqrt(beta(1.5, *smooth - 0.5) * *range);
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      case 5:
	//Gaussian
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3));
	  phase = M_2PI * unif_rand();
	  
	  lj = lines[2 * *nlines + j];
	  
	  for (k=*nsite;k--;){
	    cl = coord[k] * lines[*nlines + j];
	    for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lj;
	      ans[k  + l * *nsite] += cos(freq * eucProd + phase);
	    }
	  }
	}
	break;
      }
      break;
      
    case 3:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 2:
	//Cauchy
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3) * rgamma(*smooth ,1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 3:
	//Powered exponential
	for (j=*nlines;j--;){
	  u1 = rexp(1);
	  u2 = runif(-M_PI_2, M_PI_2);
	  G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -1 / *smooth) *
		   R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			 (2 - *smooth) / *smooth));
	
	  freq = sqrt(rchisq(3) * M_SQRT_3 / G);
	  phase = M_2PI * unif_rand();
	
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 4:
	//Bessel
	for (j=*nlines;j--;){
	  freq = sqrt(beta(1.5, *smooth - 0.5) * *range);
	  phase = M_2PI * unif_rand();
	
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}
	break;
      case 5:
	//Gaussian
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3));
	  phase = M_2PI * unif_rand();
	
	  for (k=*nsite;k--;){
	    cl1 = coord[k] * lines[j];
	    for (l=*nsite;l--;){
	      cl2 = coord[*nsite + l] * lines[*nlines + j];
	      for (m=*nsite;m--;){
		eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		ans[k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	}    
      }
    }
  }  

  else{
    //coord doesn't define a grid
    int k;
    double freq, eucProd, u1, u2, G, phase;
    switch (*dim){
    case 2:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 2:
	 //Cauchy
	 for (j=*nlines;j--;){
	   freq = sqrt(2 * rchisq(3) * rgamma(*smooth, 1));
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 3:
	 //Powered exponential
	 for (j=*nlines;j--;){
	   u1 = rexp(1);
	   u2 = runif(-M_PI_2, M_PI_2);
	   G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -1 / *smooth) *
		    R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			  (2 - *smooth) / *smooth));
	
	   freq = sqrt(rchisq(3) * M_SQRT_3 / G);
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 4:
	 //Bessel
	 for (j=*nlines;j--;){
	   freq = sqrt(beta(1.5, *smooth - 0.5) * *range);
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
       case 5:
	 //Gaussian
	 for (j=*nlines;j--;){
	   freq = sqrt(2 * rchisq(3));
	   phase = M_2PI * unif_rand();
	
	   for (k=*nsite;k--;){
	     eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	     ans[k] += cos(freq * eucProd + phase);
	   }
	 }
	 break;
      }
      break;    
    case 3:
      switch (*covmod){
      case 1:
	//Whittle-Matern
	for (j=*nlines;j--;){
	  freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 2:
	//Cauchy
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3) * rgamma(*smooth ,1));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 3:
	//Powered exponential
	for (j=*nlines;j--;){
	  u1 = rexp(1);
	  u2 = runif(-M_PI_2, M_PI_2);
	  G = fabs(sin(0.5 * *smooth * (u2 - M_PI_2)) * R_pow(cos(u2), -1 / *smooth) *
		   R_pow(cos(u2 - 0.5 * *smooth * (u2 - M_PI_2)) / u1, 
			 (2 - *smooth) / *smooth));
	  
	  freq = sqrt(rchisq(3) * M_SQRT_3 / G);
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 4:
	//Bessel
	for (j=*nlines;j--;){
	  freq = sqrt(beta(1.5, *smooth - 0.5) * *range);
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}
	break;
      case 5:
	//Gaussian
	for (j=*nlines;j--;){
	  freq = sqrt(2 * rchisq(3));
	  phase = M_2PI * unif_rand();
	  
	  for (k=*nsite;k--;){
	    eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
	      coord[2 * *nsite + k] * lines[2 * *nlines + j];
	    ans[k] += cos(freq * eucProd + phase);
	  }
	}    
      }
    }
  }  

  for (j=*neffSite;j--;){
    const double normConst = sqrt(2 * *sill / *nlines);
    ans[j] *= normConst;
  }

  if (*nugget != 0){
    const double sdnugget = sqrt(*nugget);
    for (j=*neffSite;j--;)
      ans[j] += sdnugget * norm_rand();
  }
  
  return;
}
