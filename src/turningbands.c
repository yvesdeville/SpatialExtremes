#include "header.h"

void tbm(int *nobs, int *nsite, int *dim, int *covmod, int *grid, 
	 double *coord, double *nugget, double *sill, double *range,
	 double *smooth, int *nlines, double *ans){
  
  const double normConst = sqrt(2 * *sill / *nlines),
    sdnugget = sqrt(*nugget);
  int i, j, k, l, m, ngrid = *nsite;
  double *lines, freq, eucProd, u1, u2, G, phase,
    irange = 1 / *range, u, v, w, inorm, angle;

  //rescale the coordinates
  for (i=(*nsite * *dim);i--;)
    coord[i] = coord[i] * irange;


  if (*grid)
    for (i=(*nobs * R_pow_di(*nsite, *dim));i--;)
      ans[i] = 0;

  else
    for (i=(*nobs * *nsite);i--;)
      ans[i] = 0;

  lines = (double *)R_alloc(3 * *nlines, sizeof(double));
  
  if ((*covmod == 3) && (*smooth == 2))
    //This is the gaussian case
    *covmod = 5;

  //Generate lines
  vandercorput(nlines, lines);

  GetRNGstate();
  if (*grid){
    //coord defines a grid
    ngrid = R_pow_di(*nsite, *dim);

    for (i=*nobs;i--;){
      //Random rotation of the lines
      u = unif_rand() - 0.5;
      v = unif_rand() - 0.5;
      w = unif_rand() - 0.5;
      angle = runif(0, M_2PI);

      inorm = 1 / sqrt(u * u + v * v + w * w);
    
      u *= inorm;
      v *= inorm;
      w *= inorm;
    
      rotation(lines, nlines, &u, &v, &w, &angle);

      //Turning bands part
      if (*dim == 2){
	double cl;
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd =  cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * ngrid + k  + l * *nsite] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	  break;
	case 2:
	  //Cauchy
	  for (j=*nlines;j--;){
	    freq = sqrt(2 * rchisq(3) * rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
	      eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
	      ans[i * ngrid + k  + l * *nsite] += cos(freq * eucProd + phase);
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
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * ngrid + k  + l * *nsite] += cos(freq * eucProd + phase);
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
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * ngrid + k  + l * *nsite] += cos(freq * eucProd + phase);
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
	      cl = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		eucProd = cl + coord[*nsite + l] * lines[2 * *nlines + j];
		ans[i * ngrid + k  + l * *nsite] += cos(freq * eucProd + phase);
	      }
	    }
	  }
	  break;
	}
      }

      else if (*dim == 3){
	double cl1, cl2;
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * ngrid + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * ngrid + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * ngrid + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * ngrid + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
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
	      cl1 = coord[k] * lines[*nlines + j];
	      for (l=*nsite;l--;){
		cl2 = coord[*nsite + l] * lines[*nlines + j];
		for (m=*nsite;m--;){
		  eucProd = cl1 + cl2 + coord[2 * *nsite + m] * lines[2 * *nlines + j];
		  ans[i * ngrid + k + *nsite * ( l + m * *nsite)] += cos(freq * eucProd + phase);
		}
	      }
	    }
	  }    
	}
      }
    }  
  }

  else{
    //coord doesn't define a grid
    for (i=*nobs;i--;){
      //Random rotation of the lines
      u = unif_rand() - 0.5;
      v = unif_rand() - 0.5;
      w = unif_rand() - 0.5;
      angle = runif(0, M_2PI);

      inorm = 1 / sqrt(u * u + v * v + w * w);
    
      u *= inorm;
      v *= inorm;
      w *= inorm;
    
      rotation(lines, nlines, &u, &v, &w, &angle);

      //Turning bands part
      if (*dim == 2){
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[*nlines + j] + coord[*nsite + k] * lines[2 * *nlines + j];
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
	    }
	  }
	  break;
	}
      }

      else if (*dim == 3){
	switch (*covmod){
	case 1:
	  //Whittle-Matern
	  for (j=*nlines;j--;){
	    freq = sqrt(0.5 * rchisq(3) / rgamma(*smooth, 1));
	    phase = M_2PI * unif_rand();
	
	    for (k=*nsite;k--;){
	      eucProd = coord[k] * lines[j] + coord[*nsite + k] * lines[*nlines + j] +
		coord[2 * *nsite + k] * lines[2 * *nlines + j];
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
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
	      ans[i + k * *nobs] += cos(freq * eucProd + phase);
	    }
	  }    
	}
      }
    }  
  }

  for (i=(ngrid * *nobs);i--;)
    ans[i] *= normConst;

  if (*nugget != 0)
    for (i=(ngrid * *nobs);i--;)
      ans[i] += sdnugget * norm_rand();
  
  PutRNGstate();
  return;
}
