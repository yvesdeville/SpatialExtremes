#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

#define MINF -1.0e15
#define EPS DBL_EPSILON

///////////////////////////////////
//  From schlather.c
//
void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   int *dim, int *weighted, double *weights, double *locs,
		   double *scales, double *shapes, double *sill, double *range,
		   double *smooth, double *smooth2, int *fitmarge, double *dns);
void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   int *dim, int *weighted, double *weights, double *locs,
		   double *scales, double *shapes, double *sill, double *range,
		   double *smooth, double *smooth2, int *fitmarge,double *dns);

///////////////////////////////////
//  From schlatherind.c
//
void schlatherindfull(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, int *dim, int *weighted, double *weights,
		      double *locs, double *scales, double *shapes, 
		      double *alpha, double *sill, double *range, double *smooth,
		      double *smooth2, int *fitmarge,double *dns);
void schlatherinddsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			 int *dim, int *weighted, double *weights, double *locdsgnmat,
			 double *locpenmat, int *nloccoeff, int *npparloc,
			 double *locpenalty, double *scaledsgnmat, double *scalepenmat,
			 int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
			 double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
			 double *loccoeff, double *scalecoeff, double *shapecoeff, double *alpha,
			 double *sill, double *range, double *smooth, double *smooth2, double *dns);

///////////////////////////////////
//  From geomgauss.c
//
void geomgaussfull(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
		   int *weighted, double *weights, double *locs, double *scales, double *shapes,
		   double *sigma2, double *sigma2Bound, double *sill, double *range,
		   double *smooth, double *smooth2, int *fitmarge,double *dns);
void geomgaussdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		      int *dim, int *weighted, double *weights, double *locdsgnmat,
		      double *locpenmat, int *nloccoeff, int *npparloc, double *locpenalty,
		      double *scaledsgnmat, double *scalepenmat, int *nscalecoeff, int *npparscale,
		      double *scalepenalty, double *shapedsgnmat, double *shapepenmat, int *nshapecoeff,
		      int *npparshape, double *shapepenalty, double *loccoeff, double *scalecoeff,
		      double *shapecoeff, double *sigma2, double *sigma2Bound, double *sill,
		      double *range, double *smooth, double *smooth2, double *dns);

///////////////////////////////////
//  From nsgeomgauss.c
//
void nsgeomgaussfull(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, int *dim, double *locs, double *scales, double *shapes,
		     double *sigma2dsgnmat, double *sigma2coeff, int *nsigma2coeff,
		     double *sill, double *range, double *smooth, double *smooth2, int *fitmarge,
		     double *dns);
void nsgeomgaussdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			int *dim, double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
			double *locpenalty, double *scaledsgnmat, double *scalepenmat,
			int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
			double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
			double *sigma2dsgnmat, int *nsigma2coeff, double *loccoeff, double *scalecoeff,
			double *shapecoeff, double *sigma2coeff, double *sill, double *range,
			double *smooth, double *smooth2, double *dns);

///////////////////////////////////
//  From brownResnick.c
//
void brownresnickfull(double *data, double *dist, int *nSite, int *nObs, int *weighted,
		      double *weights, double *locs, double *scales, double *shapes,
		      double *range, double *smooth, int *fitmarge, double *dns);
void brownresnickdsgnmat(double *data, double *dist, int *nSite, int *nObs, int *weighted,
			 double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
			 int *npparloc, double *locpenalty, double *scaledsgnmat, double *scalepenmat,
			 int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
			 double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
			 double *loccoeff, double *scalecoeff, double *shapecoeff, double *range,
			 double *smooth, double *dns);

///////////////////////////////////
//  From smith.c
//
void smithfull(double *data, double *distVec, int *nSite, int *nObs, int *weighted, double *weights,
	       double *locs, double *scales, double *shapes, double *cov11, double *cov12,
	       double *cov22, int *fitmarge, double *dns);
void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		  double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat, double *scalepenmat,
		  int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
		  double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
		  double *loccoeff, double *scalecoeff, double *shapecoeff, double *cov11,
		  double *cov12, double *cov22, double *dns);

///////////////////////////////////
//  From smith3d.c
//
void smithfull3d(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		 double *weights, double *locs, double *scales, double *shapes,
		 double *cov11, double *cov12, double *cov13, double *cov22,
		 double *cov23, double *cov33, int *fitmarge, double *dns);
void smithdsgnmat3d(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		    double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
		    int *npparloc, double *locpenalty, double *scaledsgnmat,
		    double *scalepenmat, int *nscalecoeff, int *npparscale,
		    double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		    int *nshapecoeff, int *npparshape, double *shapepenalty,
		    double *loccoeff, double *scalecoeff, double *shapecoeff,
		    double *cov11, double *cov12, double *cov13, double *cov22,
		    double *cov23, double *cov33, double *dns);

///////////////////////////////////
//  From utils.c
//
void distance(double *coord, int *nDim, int *nSite,
	      int *vec, double *dist);
double gev2frech(double *data, int nObs, int nSite, double *locs,
		 double *scales, double *shapes, double *jac, double *frech);
double dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
		     double *shapedsgnmat, double *loccoeff, 
		     double *scalecoeff, double *shapecoeff,
		     int nSite, int nloccoeff, int nscalecoeff,
		     int nshapecoeff, double *locs, double *scales,
		     double *shapes);
void dsgnmat2Alpha(double *alphadsgnmat, double *alphacoeff, 
		   int nSite, int nalphacoeff, double *alphas);
void dsgnmat2Sigma2(double *sigma2dsgnmat, double *sigma2coeff, 
		    int nSite, int nsigma2coeff, double *sigma2);
void gev(double *prob, int *n, double *locs, double *scales, double *shapes,
	 double *quant);

///////////////////////////////////
//  From univllik.c
//
void gevlik(double *data, int *n, double *loc, double *scale,
	    double *shape, double *dns);
void gpdlik(double *exceed, int *n, double *thresh, double *scale,
	    double *shape, double *dns);

///////////////////////////////////
//  From covariance.c
//
double whittleMatern(double *dist, int nPairs, double sill, double range,
		     double smooth, double *rho);
double cauchy(double *dist, int nPairs, double sill, double range,
	      double smooth, double *rho);
double caugen(double *dist, int nPairs, double sill, double range,
	      double smooth, double smooth2, double *rho);
double powerExp(double *dist, int nPairs, double sill, double range,
		double smooth, double *rho);
double bessel(double *dist, int nPairs, int dim, double sill,
	      double range, double smooth, double *rho);
double mahalDistFct(double *distVec, int nPairs, double *cov11,
		    double *cov12, double *cov22, double *mahal);
double mahalDistFct3d(double *distVec, int nPairs, double *cov11,
		      double *cov12, double *cov13, double *cov22, 
		      double *cov23, double *cov33, double *mahal);
double geomCovariance(double *dist, int nPairs, int dim, int covmod,
		      double sigma2, double sigma2Bound, double sill,
		      double range, double smooth, double smooth2,
		      double *rho);
double nsgeomCovariance(double *dist, int nSite, int dim, int covmod,
			double *sigma2, double sill, double range,
			double smooth, double smooth2, double *rho);
double brownResnick(double *dist, int nPairs, double range, double smooth,
		    double *rho);

///////////////////////////////////
//  From mcmc.c
//
SEXP gibbs(SEXP n, SEXP np, SEXP thin, SEXP init,
	   SEXP psd, SEXP f, SEXP rho);

///////////////////////////////////
//  From pairwiselik.c
//
double lpliksmith(double *data, double *rho, double *jac,
		  int nObs, int nSite);
double lplikschlather(double *data, double *rho, double *jac,
		      int nObs, int nSite);
double lplikschlatherind(double *data, double alpha, double *rho,
			 double *jac, int nObs, int nSite);

///////////////////////////////////
//  From weightedPairwiselik.c
//
double wlplikschlather(double *data, double *rho, double *jac,
		       int nObs, int nSite, double *weights);
double wlpliksmith(double *data, double *mahalDist, double *jac,
		   int nObs, int nSite, double *weights);
double wlplikschlatherind(double *data, double alpha, double *rho,
			  double *jac, int nObs, int nSite, double *weights);

///////////////////////////////////
//  From penalizations.c
//
double penalization(double *penmat, double *beta, double pencoeff, int n,
		    int nppar);
double penalization2(double *penmat, double *beta, double pencoeff, int n,
		     int nppar);


///////////////////////////////////
//  From extcoeff.c
//
void extCoeffSmith(double *frech, int *nObs, int *nSite,
		   double *extCoeff);
void extCoeffST(double *frech, double *xBar, double *z, double *theta,
		int *nObs, double *dns);

///////////////////////////////////
//  From fitcovmat.c
//
void fitcovmat2d(double *cov11, double *cov12, double *cov22,
		 int *nPairs, double *dist, double *extcoeff,
		 double *weights, double *ans);
void fitcovmat3d(double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33,
		 int *nPairs, double *dist, double *extcoeff,
		 double *weights, double *ans);
void fitcovariance(int *covmod, double *sill, double *range, double *smooth,
		   double *smooth2, int *nPairs, int *dim, double *distVec,
		   double *extcoeff, double *weights, double *ans);
void fiticovariance(int *covmod, double *alpha, double *sill, double *range,
		    double *smooth, double *smooth2, int *nPairs, int *dim,
		    double *dist, double *extcoeff, double *weights, double *ans);
void fitgcovariance(int *covmod, double *sigma2, double *sigma2Bound, double *sill,
		    double *range, double *smooth, double *smooth2, int *nPairs,
		    int *dim, double *dist, double *extcoeff, double *weights,
		    double *ans);
void fitbrcovariance(double *range, double *smooth, int *nPairs,
		     double *dist, double *extcoeff, double *weights,
		     double *ans);

///////////////////////////////////
//  From spatgevlik.c
//
void spatgevlik(double *data, double *covariables, int *nSite, int *nObs,
		double *locdsgnmat, double *locpenmat, int *nloccoeff,
		int *npparloc, double *locpenalty, double *scaledsgnmat,
		double *scalepenmat, int *nscalecoeff, int *npparscale,
		double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		int *nshapecoeff, int *npparshape, double *shapepenalty,
		double *loccoeff, double *scalecoeff, double *shapecoeff,
		double *dns);

///////////////////////////////////
//  From madogram.c
//
void madogram(double *data, int *nObs, int *nSite, double *mado);
void lmadogram(double *data, int *nObs, int *nSite, double *lambda,
	       int *nLambda, double *lmado);

///////////////////////////////////
//  From simsmith.c
//
void rsmith1d(double *coord, double *center, double *edge, int *nObs,
	      int *nSites, double *var, double *ans);
void rsmith2d(double *coord, double *center, double *edge, int *nObs,
	      int *nSites, int *grid, double *cov11, double *cov12,
	      double *cov22, double *ans);

///////////////////////////////////
//  From direct.c
//
void buildcovmat(int *nSite, int *grid, int *covmod, double *coord, int *dim,
		 double *nugget, double *sill, double *range, double *smooth,
		 double *covMat);
void direct(int *n, int *nSite, int *grid, int *covmod, double *coord, int *dim,
	    double *nugget, double *sill, double *range, double *smooth,
	    double *ans);

///////////////////////////////////
//  From randomlines.c
//
void vandercorput(int *n, double *coord);
void rotation(double *coord, int *n, double *u, double *v, double *w,
	      double *angle);

///////////////////////////////////
//  From turningbands.c
//
void tbm(int *nobs, int *nsite, int *dim, int *covmod, int *grid, 
	 double *coord, double *nugget, double *sill, double *range,
	 double *smooth, int *nlines, double *ans);

///////////////////////////////////
//  From simschlather.c
//
void rschlathertbm(double *coord, int *nObs, int *nSites, int *dim,
		   int *covmod, int *grid, double *sill, double *range,
		   double *smooth, double *uBound, int *nlines,
		   double *ans);
void rschlatherdirect(double *coord, int *nObs, int *nSites, int *dim,
		      int *covmod, int *grid, double *sill, double *range,
		      double *smooth, double *uBound, double *ans);
void tbmcore(int *nsite, int *neffSite, int *dim, int *covmod,
	     int *grid, double *coord, double *nugget, double *sill,
	     double *range, double *smooth, int *nlines, double *lines,
	     double *ans);

///////////////////////////////////
//  From simgeometric.c
//
void rgeomtbm(double *coord, int *nObs, int *nSite, int *dim,
	      int *covmod, int *grid, double *sigma2, double *sill,
	      double *range, double *smooth, double *uBound, int *nlines,
	      double *ans);
void rgeomdirect(double *coord, int *nObs, int *nSite, int *dim,
		 int *covmod, int *grid, double *sigma2, double *sill,
		 double *range, double *smooth, double *uBound,
		 double *ans);

///////////////////////////////////
//  From gpdproc.c
//
void gpdprocfull(double *data, double *distVec, int *nSite,
		 int *nObs, double *excRates, double *threshs, double *scales,
		 double *shapes, double *cov11, double *cov12,
		 double *cov22, int *fitmarge, double *dns);
double gpd2ugpd(double *data, int nObs, int nSite, double *excRates,
		double *threshs, double *scales, double *shapes,
		double *jac, double *ugpd);
double lpliksmithgpd(double *data, double *mahalDist, double *jac,
		     double *excRates, int nObs, int nSite);

///////////////////////////////////
//  From standardErrors.c
//
void smithstderr(double *data, double *distVec, int *nSite,
		 int *nObs, double *locdsgnmat, int *nloccoeff,
		 double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		 int *nshapecoeff, double *loccoeff, double *scalecoeff,
		 double *shapecoeff, double *cov11, double *cov12,
		 double *cov22, int *fitmarge, double *hess, double *grad);
void smithstderr3d(double *data, double *distVec, int *nSite,
		   int *nObs, double *locdsgnmat, int *nloccoeff,
		   double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *loccoeff, double *scalecoeff,
		   double *shapecoeff, double *cov11, double *cov12, double *cov13,
		   double *cov22, double *cov23, double *cov33, int *fitmarge, double *hess,
		   double *grad);
void schlatherstderr(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, double *locdsgnmat, int *nloccoeff,
		     double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		     int *nshapecoeff, double *loccoeff, double *scalecoeff,
		     double *shapecoeff, double *sill, double *range, double *smooth,
		     double *smooth2, int *fitmarge, double *hess, double *grad);
void schlatherindstderr(int *covmod, double *data, double *dist, int *nSite,
			int *nObs, double *locdsgnmat, int *nloccoeff,
			double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
			int *nshapecoeff, double *loccoeff, double *scalecoeff,
			double *shapecoeff, double *alpha, double *sill, double *range,
			double *smooth, double *smooth2, int *fitmarge, double *hess,
			double *grad);
void geomgaussstderr(int *covmod, double *data, double *dist, int *nSite,
		     int *nObs, double *locdsgnmat, int *nloccoeff,
		     double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		     int *nshapecoeff, double *loccoeff, double *scalecoeff,
		     double *shapecoeff, double *sigma2, double *sill, double *range,
		     double *smooth, double *smooth2, int *fitmarge, double *hess,
		     double *grad);
void brownresnickstderr(double *data, double *dist, int *nSite, int *nObs,
			double *locdsgnmat, int *nloccoeff, double *scaledsgnmat,
			int *nscalecoeff, double *shapedsgnmat, int *nshapecoeff,
			double *loccoeff, double *scalecoeff, double *shapecoeff,
			double *range, double *smooth, int *fitmarge, double *hess,
			double *grad);
void spatgevstderr(double *data, int *nSite, int *nObs, double *locdsgnmat,
		   int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		   double *shapedsgnmat, int *nshapecoeff, double *loccoeff,
		   double *scalecoeff, double *shapecoeff, double *hess, 
		   double *grad);

///////////////////////////////////
//  From weightedStandardErrors.c
//
void wsmithstderr(double *data, double *distVec, int *nSite,
		  int *nObs, double *locdsgnmat, int *nloccoeff,
		  double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		  int *nshapecoeff, double *loccoeff, double *scalecoeff,
		  double *shapecoeff, double *cov11, double *cov12,
		  double *cov22, int *fitmarge, double *weights, double *hess,
		  double *grad);
void wsmithstderr3d(double *data, double *distVec, int *nSite,
		    int *nObs, double *locdsgnmat, int *nloccoeff,
		    double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		    int *nshapecoeff, double *loccoeff, double *scalecoeff,
		    double *shapecoeff, double *cov11, double *cov12, double *cov13,
		    double *cov22, double *cov23, double *cov33, int *fitmarge,
		    double *weights, double *hess, double *grad);
void wschlatherstderr(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, double *locdsgnmat, int *nloccoeff,
		      double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		      int *nshapecoeff, double *loccoeff, double *scalecoeff,
		      double *shapecoeff, double *sill, double *range, double *smooth,
		      double *smooth2, int *fitmarge, double *weights, double *hess,
		      double *grad);
void wschlatherindstderr(int *covmod, double *data, double *dist, int *nSite,
			 int *nObs, double *locdsgnmat, int *nloccoeff,
			 double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
			 int *nshapecoeff, double *loccoeff, double *scalecoeff,
			 double *shapecoeff, double *alpha, double *sill, double *range,
			 double *smooth, double *smooth2, int *fitmarge, double *weights,
			 double *hess, double *grad);
void wgeomgaussstderr(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, double *locdsgnmat, int *nloccoeff,
		      double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		      int *nshapecoeff, double *loccoeff, double *scalecoeff,
		      double *shapecoeff, double *sigma2, double *sill, double *range,
		      double *smooth, double *smooth2, int *fitmarge, double *weights,
		      double *hess, double *grad);
void wbrownresnickstderr(double *data, double *dist, int *nSite, int *nObs,
			 double *locdsgnmat, int *nloccoeff, double *scaledsgnmat,
			 int *nscalecoeff, double *shapedsgnmat, int *nshapecoeff,
			 double *loccoeff, double *scalecoeff, double *shapecoeff,
			 double *range, double *smooth, int *fitmarge, double *weights,
			 double *hess, double *grad);
