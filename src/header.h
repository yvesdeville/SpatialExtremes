#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#define MINF -1.0e15
#define EPS DBL_EPSILON

///////////////////////////////////
//  From schlather.c
//
void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   double *locs, double *scales, double *shapes, double *sill,
		   double *range, double *smooth, int *fitmarge, double *dns);
void schlatherdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		      double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
		      double *locpenalty, double *scaledsgnmat, double *scalepenmat,
		      int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
		      double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
		      double *loccoeff, double *scalecoeff, double *shapecoeff, double *sill,
		      double *range, double *smooth, double *dns);

///////////////////////////////////
//  From schlatherind.c
//
void schlatherindfull(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, double *locs, double *scales, double *shapes,
		      double *alpha, double *sill, double *range, double *smooth,
		      int *fitmarge,double *dns);
void schlatherinddsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			 double *locdsgnmat, double *locpenmat, int *nloccoeff, int *npparloc,
			 double *locpenalty, double *scaledsgnmat, double *scalepenmat,
			 int *nscalecoeff, int *npparscale, double *scalepenalty, double *shapedsgnmat,
			 double *shapepenmat, int *nshapecoeff, int *npparshape, double *shapepenalty,
			 double *loccoeff, double *scalecoeff, double *shapecoeff, double *alpha,
			 double *sill, double *range, double *smooth, double *dns);

///////////////////////////////////
//  From smith.c
//
void smithfull(double *data, double *distVec, int *nSite,
	       int *nObs, double *locs, double *scales, double *shapes,
	       double *cov11, double *cov12, double *cov22, int *fitmarge,
	       double *dns);
void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, 
		  double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat,
		  double *scalepenmat, int *nscalecoeff, int *npparscale,
		  double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		  int *nshapecoeff, int *npparshape, double *shapepenalty,
		  double *loccoeff, double *scalecoeff, double *shapecoeff,
		  double *cov11, double *cov12, double *cov22, double *dns);

///////////////////////////////////
//  From smith3d.c
//
void smithfull3d(double *data, double *distVec, int *nSite,
		 int *nObs, double *locs, double *scales, double *shapes,
		 double *cov11, double *cov12, double *cov13, double *cov22,
		 double *cov23, double *cov33, int *fitmarge, double *dns);
void smithdsgnmat3d(double *data, double *distVec, int *nSite, int *nObs, 
		    double *locdsgnmat, double *locpenmat, int *nloccoeff,
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
double powerExp(double *dist, int nPairs, double sill, double range,
		double smooth, double *rho);
double mahalDistFct(double *distVec, int nPairs, double *cov11,
		    double *cov12, double *cov22, double *mahal);
double mahalDistFct3d(double *distVec, int nPairs, double *cov11,
		      double *cov12, double *cov13, double *cov22, 
		      double *cov23, double *cov33, double *mahal);
double geomCovariance(double *dist, int nPairs, int covmod,
		      double sigma2, double sill, double range,
		      double smooth, double *rho);

///////////////////////////////////
//  From mcmc.c
//
SEXP gibbs(SEXP n, SEXP np, SEXP thin, SEXP init,
	   SEXP psd, SEXP f, SEXP rho);

///////////////////////////////////
//  From gradients.c
//
void smithgrad(double *data, double *distVec, int *nSite,
	       int *nObs, double *locdsgnmat, int *nloccoeff,
	       double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
	       int *nshapecoeff, double *loccoeff, double *scalecoeff,
	       double *shapecoeff, double *cov11, double *cov12,
	       double *cov22, int *fitmarge, double *grad);
void smithgrad3d(double *data, double *distVec, int *nSite,
		 int *nObs, double *locdsgnmat, int *nloccoeff,
		 double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		 int *nshapecoeff, double *loccoeff, double *scalecoeff,
		 double *shapecoeff, double *cov11, double *cov12, double *cov13,
		 double *cov22, double *cov23, double *cov33, int *fitmarge, double *grad);
void schlathergrad(int *covmod, double *data, double *dist, int *nSite,
		   int *nObs, double *locdsgnmat, int *nloccoeff,
		   double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *loccoeff, double *scalecoeff,
		   double *shapecoeff, double *sill, double *range, double *smooth,
		   int *fitmarge, double *grad);
void schlatherindgrad(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, double *locdsgnmat, int *nloccoeff,
		      double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		      int *nshapecoeff, double *loccoeff, double *scalecoeff,
		      double *shapecoeff, double *alpha, double *sill, double *range,
		      double *smooth, int *fitmarge, double *grad);
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
		   int *nPairs, double *distVec, double *extcoeff,
		   double *weights, double *ans);
