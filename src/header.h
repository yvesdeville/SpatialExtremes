#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>
#include <R_ext/Applic.h>

#define MINF -1.0e15
#define EPS DBL_EPSILON

///////////////////////////////////
//  From schlather.c
//
void schlatherfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   int *dim, int *weighted, double *weights, double *locs,
		   double *scales, double *shapes, double *sill, double *range,
		   double *smooth, double *smooth2, int *fitmarge, double *dns);
void schlatherdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
		      int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
		      int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
		      double *scalepenmat, int *nscalecoeff, int *npparscale,
		      double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		      int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
		      double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
		      int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
		      double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
		      double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
		      int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
		      double *loccoeff, double *scalecoeff, double *shapecoeff,
		      double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
		      double *sill, double *range, double *smooth, double *smooth2, double *dns);

///////////////////////////////////
//  From schlatherind.c
//
void schlatherindfull(int *covmod, double *data, double *dist, int *nSite,
		      int *nObs, int *dim, int *weighted, double *weights,
		      double *locs, double *scales, double *shapes,
		      double *alpha, double *sill, double *range, double *smooth,
		      double *smooth2, int *fitmarge,double *dns);
void schlatherinddsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
			 int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
			 int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
			 double *scalepenmat, int *nscalecoeff, int *npparscale,
			 double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
			 int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
			 double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
			 int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
			 double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
			 double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
			 int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
			 double *loccoeff, double *scalecoeff, double *shapecoeff,
			 double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
			 double *alpha, double *sill, double *range, double *smooth,
			 double *smooth2, double *dns);

///////////////////////////////////
//  From geomgauss.c
//
void geomgaussfull(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
		   int *weighted, double *weights, double *locs, double *scales, double *shapes,
		   double *sigma2, double *sigma2Bound, double *sill, double *range,
		   double *smooth, double *smooth2, int *fitmarge,double *dns);
void geomgaussdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
		      int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
		      int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
		      double *scalepenmat, int *nscalecoeff, int *npparscale,
		      double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		      int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
		      double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
		      int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
		      double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
		      double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
		      int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
		      double *loccoeff, double *scalecoeff, double *shapecoeff,
		      double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
		      double *sigma2, double *sigma2Bound, double *sill, double *range,
		      double *smooth, double *smooth2, double *dns);

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
			 int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
			 int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
			 double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
			 int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
			 double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
			 double *temppenaltyshape, double *loccoeff, double *scalecoeff,
			 double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
			 double *tempcoeffshape,double *range, double *smooth, double *dns);

///////////////////////////////////
//  From smith.c
//
void smithfull(double *data, double *distVec, int *nSite, int *nObs, int *weighted, double *weights,
	       double *locs, double *scales, double *shapes, double *cov11, double *cov12,
	       double *cov22, int *fitmarge, double *dns);
void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		  double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat,
		  double *scalepenmat, int *nscalecoeff, int *npparscale,
		  double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		  int *nshapecoeff, int *npparshape, double *shapepenalty,
		  int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
		  int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
		  double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
		  int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
		  double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
		  double *temppenaltyshape, double *loccoeff, double *scalecoeff,
		  double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
		  double *tempcoeffshape, double *cov11, double *cov12, double *cov22,
		  double *dns);

///////////////////////////////////
//  From smith3d.c
//
void smithfull3d(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		 double *weights, double *locs, double *scales, double *shapes,
		 double *cov11, double *cov12, double *cov13, double *cov22,
		 double *cov23, double *cov33, int *fitmarge, double *dns);
void smithdsgnmat(double *data, double *distVec, int *nSite, int *nObs, int *weighted,
		  double *weights, double *locdsgnmat, double *locpenmat, int *nloccoeff,
		  int *npparloc, double *locpenalty, double *scaledsgnmat,
		  double *scalepenmat, int *nscalecoeff, int *npparscale,
		  double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		  int *nshapecoeff, int *npparshape, double *shapepenalty,
		  int *usetempcov, double *tempdsgnmatloc, double *temppenmatloc,
		  int *ntempcoeffloc, int *nppartempcoeffloc, double *temppenaltyloc,
		  double *tempdsgnmatscale, double *temppenmatscale, int *ntempcoeffscale,
		  int *nppartempcoeffscale, double *temppenaltyscale, double *tempdsgnmatshape,
		  double *temppenmatshape, int *ntempcoeffshape, int *nppartempcoeffshape,
		  double *temppenaltyshape, double *loccoeff, double *scalecoeff,
		  double *shapecoeff, double *tempcoeffloc, double *tempcoeffscale,
		  double *tempcoeffshape, double *cov11, double *cov12, double *cov22,
		  double *dns);

///////////////////////////////////
//  From extremalt.c
//
void extremaltfull(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		   int *dim, int *weighted, double *weights, double *locs, double *scales,
		   double *shapes, double *sill, double *range, double *smooth, double *smooth2,
		   double *df, int *fitmarge, double *dns);
void extremaltdsgnmat(int *covmod, double *data, double *dist, int *nSite, int *nObs, int *dim,
		      int *weighted, double *weights, double *locdsgnmat, double *locpenmat,
		      int *nloccoeff, int *npparloc, double *locpenalty, double *scaledsgnmat,
		      double *scalepenmat, int *nscalecoeff, int *npparscale,
		      double *scalepenalty, double *shapedsgnmat, double *shapepenmat,
		      int *nshapecoeff, int *npparshape, double *shapepenalty, int *usetempcov,
		      double *tempdsgnmatloc, double *temppenmatloc, int *ntempcoeffloc,
		      int *nppartempcoeffloc, double *temppenaltyloc, double *tempdsgnmatscale,
		      double *temppenmatscale, int *ntempcoeffscale, int *nppartempcoeffscale,
		      double *temppenaltyscale, double *tempdsgnmatshape, double *temppenmatshape,
		      int *ntempcoeffshape, int *nppartempcoeffshape, double *temppenaltyshape,
		      double *loccoeff, double *scalecoeff, double *shapecoeff,
		      double *tempcoeffloc, double *tempcoeffscale, double *tempcoeffshape,
		      double *sill, double *range, double *smooth, double *smooth2, double *df,
		      double *dns);

///////////////////////////////////
//  From utils.c
//
void distance(double *coord, int *nDim, int *nSite, int *vec, double *dist);
void distance2orig(double *coord, int n, int dim, double *dist, int grid);
double gev2frech(double *data, int nObs, int nSite, double *locs,
		 double *scales, double *shapes, double *jac, double *frech);
double gev2frechTrend(double *data, int nObs, int nSite, double *locs, double *scales,
		      double *shapes, double *trendlocs, double *trendscales,
		      double *trendshapes,double *jac, double *frech);
double dsgnmat2Param(double *locdsgnmat, double *scaledsgnmat,
		     double *shapedsgnmat, double *loccoeff,
		     double *scalecoeff, double *shapecoeff,
		     int nSite, int nloccoeff, int nscalecoeff,
		     int nshapecoeff, double *locs, double *scales,
		     double *shapes);
double dsgnmat2Param2(double *locdsgnmat, double *scaledsgnmat, double *shapedsgnmat,
		      double *tempdsgnmatLoc, double *tempdsgnmatScale, double *tempdsgnmatShape,
		      double *loccoeff, double *scalecoeff, double *shapecoeff,
		      double *tempcoeffLoc, double *tempcoeffScale, double *tempcoeffShape,
		      int nSite, int nObs, int *usetempcov, int nloccoeff, int nscalecoeff,
		      int nshapecoeff, int ntempcoeffLoc, int ntempcoeffScale,
		      int ntempcoeffShape, double *locs, double *scales, double *shapes);
void dsgnmat2temptrend(double *dsgnmatloc, double *dsgnmatscale, double *dsgnmatshape,
		       double *loccoeff, double *scalecoeff, double *shapecoeff, int nSite,
		       int nObs, int *usetempcov, int nloccoeff, int nscalecoeff,
		       int nshapecoeff, double *trendlocs, double *trendscales,
		       double *trendshapes);
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
double whittleMatern(double *dist, int n, double sill, double range,
		     double smooth, double *rho);
double cauchy(double *dist, int n, double sill, double range,
	      double smooth, double *rho);
double caugen(double *dist, int n, double sill, double range,
	      double smooth, double smooth2, double *rho);
double powerExp(double *dist, int n, double sill, double range,
		double smooth, double *rho);
double bessel(double *dist, int n, int dim, double sill,
	      double range, double smooth, double *rho);
double mahalDistFct(double *distVec, int n, double *cov11,
		    double *cov12, double *cov22, double *mahal);
double mahalDistFct3d(double *distVec, int n, double *cov11,
		      double *cov12, double *cov13, double *cov22,
		      double *cov23, double *cov33, double *mahal);
double geomCovariance(double *dist, int n, int dim, int covmod,
		      double sigma2, double sigma2Bound, double sill,
		      double range, double smooth, double smooth2,
		      double *rho);
double nsgeomCovariance(double *dist, int nSite, int dim, int covmod,
			double *sigma2, double sill, double range,
			double smooth, double smooth2, double *rho);
double brownResnick(double *dist, int n, double range, double smooth,
		    double *rho);
double fbm(double *coord, double *dist, int dim, int nSite, double sill, double range,
	   double smooth, double *rho);

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
double lplikextremalt(double *data, double *rho, double df, double *jac,
		      int nObs, int nSite);

///////////////////////////////////
//  From weightedPairwiselik.c
//
double wlplikschlather(double *data, double *rho, double *jac,
		       int nObs, int nSite, double *weights);
double wlpliksmith(double *data, double *mahalDist, double *jac,
		   int nObs, int nSite, double *weights);
double wlplikschlatherind(double *data, double alpha, double *rho,
			  double *jac, int nObs, int nSite, double *weights);
double wlplikextremalt(double *data, double *rho, double df, double *jac,
		       int nObs, int nSite, double *weights);

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
void fittcovariance(int *covmod, double *sill, double *range, double *smooth,
		    double *smooth2, double *DoF, int *nPairs, int *dim, double *dist,
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
		int *usetempcov, double *tempdsgnmatLoc, double *temppenmatLoc,
		int *ntempcoeffLoc, int *nppartempcoeffLoc, double *temppenaltyLoc,
		double *tempdsgnmatScale, double *temppenmatScale, int *ntempcoeffScale,
		int *nppartempcoeffScale, double *temppenaltyScale, double *tempdsgnmatShape,
		double *temppenmatShape, int *ntempcoeffShape, int *nppartempcoeffShape,
		double *temppenaltyShape, double *loccoeff, double *scalecoeff,
		double *shapecoeff, double *tempcoeffLoc, double *tempcoeffScale,
		double *tempcoeffShape, double *dns);

///////////////////////////////////
//  From madogram.c
//
void madogram(double *data, int *nObs, int *nSite, double *mado);
void variogram(double *data, int *nObs, int *nSite, double *vario);
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
void circcore(double *rho, double *a, double *ia, int m, int halfM, int mdag,
	      int mdagbar, int ngrid, int nbar, double isqrtMbar, double nugget,
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
void rgeomcirc(int *nObs, int *ngrid, double *steps, int *dim,
	       int *covmod, double *sigma2, double *sill, double *range,
	       double *smooth, double *uBound, double *ans);

///////////////////////////////////
//  From simBrownResnick.c
//
void rbrowndirect(double *coord, double *bounds, int *nObs, int *nSite,
		  int *dim, int *grid, double *range, double *smooth,
		  double *ans);

///////////////////////////////////
//  From simextremalt.c
//
void rextremalttbm(double *coord, int *nObs, int *nSite, int *dim,
		   int *covmod, int *grid, double *sill, double *range,
		   double *smooth, double *DoF, int *blockSize, int *nlines,
		   double *ans);
void rextremaltdirect(double *coord, int *nObs, int *nSite, int *dim,
		      int *covmod, int *grid, double *sill, double *range,
		      double *smooth, double *DoF, int *blockSize, double *ans);
void rextremaltcirc(int *nObs, int *ngrid, double *steps, int *dim,
		    int *covmod, double *sill, double *range,
		    double *smooth, double *DoF, int *blockSize, double *ans);

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
void smithstderr(double *data, double *distVec, int *nSite, int *nObs, double *locdsgnmat,
		 int *nloccoeff, double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		 int *nshapecoeff, double *tempdsgnmatloc, int *ntemploccoeff,
		 double *tempdsgnmatscale, int *ntempscalecoeff, double *tempdsgnmatshape,
		 int *ntempshapecoeff,double *loccoeff, double *scalecoeff, double *shapecoeff,
		 double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff,
		 double *cov11, double *cov12, double *cov22, int *fitmarge, int *usetempcov,
		 double *weights, double *hess, double *grad);
void smithstderr3d(double *data, double *distVec, int *nSite, int *nObs, double *locdsgnmat,
		   int *nloccoeff, double *scaledsgnmat, int *nscalecoeff, double *shapedsgnmat,
		   int *nshapecoeff, double *tempdsgnmatloc, int *ntemploccoeff,
		   double *tempdsgnmatscale, int *ntempscalecoeff, double *tempdsgnmatshape,
		   int *ntempshapecoeff, double *loccoeff, double *scalecoeff, double *shapecoeff,
		   double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff, double *cov11,
		   double *cov12, double *cov13, double *cov22, double *cov23, double *cov33, int *fitmarge,
		   int *usetempcov, double *weights, double *hess, double *grad);
void schlatherstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		     double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		     double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		     int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		     double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
		     double *scalecoeff, double *shapecoeff, double *temploccoeff,
		     double *tempscalecoeff, double *tempshapecoeff, double *sill, double *range,
		     double *smooth, double *smooth2, int *fitmarge, int *usetempcov, double *weights,
		     double *hess, double *grad);
void schlatherindstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
			double *locdsgnmat, int *nloccoeff, double *scaledsgnmat,
			int *nscalecoeff, double *shapedsgnmat,	int *nshapecoeff,
			double *tempdsgnmatloc, int *ntemploccoeff, double *tempdsgnmatscale,
			int *ntempscalecoeff, double *tempdsgnmatshape, int *ntempshapecoeff,
			double *loccoeff, double *scalecoeff, double *shapecoeff,
			double *temploccoeff, double *tempscalecoeff, double *tempshapecoeff,
			double *alpha, double *sill, double *range, double *smooth,
			double *smooth2, int *fitmarge, int *usetempcov, double *weights, double *hess,
			double *grad);
void geomgaussstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		     double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		     double *shapedsgnmat, int *nshapecoeff,  double *tempdsgnmatloc,
		     int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		     double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
		     double *scalecoeff, double *shapecoeff, double *temploccoeff,
		     double *tempscalecoeff, double *tempshapecoeff, double *sigma2,
		     double *sill, double *range, double *smooth, double *smooth2,
		     int *fitmarge, int *usetempcov, double *weights, double *hess, double *grad);
void brownresnickstderr(double *data, double *dist, int *nSite, int *nObs, double *locdsgnmat,
			int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
			double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
			int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
			double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
			double *scalecoeff, double *shapecoeff, double *temploccoeff,
			double *tempscalecoeff, double *tempshapecoeff, double *range,
			double *smooth, int *fitmarge, int *usetempcov, double *weights,
			double *hess, double *grad);
void spatgevstderr(double *data, int *nSite, int *nObs, double *locdsgnmat,
		   int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		   double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		   int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		   double *tempdsgnmatshape, int *ntempshapecoeff,  double *loccoeff,
		   double *scalecoeff, double *shapecoeff, double *temploccoeff,
		   double *tempscalecoeff, double *tempshapecoeff, int *usetempcov,
		   double *hess, double *grad);
void extremaltstderr(int *covmod, double *data, double *dist, int *nSite, int *nObs,
		     double *locdsgnmat, int *nloccoeff, double *scaledsgnmat, int *nscalecoeff,
		     double *shapedsgnmat, int *nshapecoeff, double *tempdsgnmatloc,
		     int *ntemploccoeff, double *tempdsgnmatscale, int *ntempscalecoeff,
		     double *tempdsgnmatshape, int *ntempshapecoeff, double *loccoeff,
		     double *scalecoeff, double *shapecoeff, double *temploccoeff,
		     double *tempscalecoeff, double *tempshapecoeff, double *sill, double *range,
		     double *smooth, double *smooth2, double *df, int *fitmarge, int *usetempcov,
		     double *weights, double *hess, double *grad);

///////////////////////////////////
//  From standardErrorsCommonPart.c
//
void marginalPartSmith(int *start, int *nObs, int *nSite, double *data, double *frech,
		       double *mahalDist, double *locs, double *scales, double *shapes,
		       double *trendlocs, double *trendscales, double *trendshapes,
		       int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
		       int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
		       double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
		       double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
		       double *hess, double *grad);
void marginalPartSchlat(int *start, int *nObs, int *nSite, double *data, double *frech,
			double *rho, double *locs, double *scales, double *shapes,
			double *trendlocs, double *trendscales, double *trendshapes,
			int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
			int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
			double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
			double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
			double *hess, double *grad);
void marginalPartiSchlat(int *start, int *nObs, int *nSite, double *data, double *frech,
			 double *alpha, double *rho, double *locs, double *scales, double *shapes,
			 double *trendlocs, double *trendscales, double *trendshapes,
			 int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
			 int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
			 double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
			 double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
			 double *hess, double *grad);
void marginalPartExtremalt(int *start, int *nObs, int *nSite, double *data, double *frech,
			   double *df, double *rho, double *locs, double *scales, double *shapes,
			   double *trendlocs, double *trendscales, double *trendshapes,
			   int *nloccoeff, int *nscalecoeff, int *nshapecoeff, int *ntemploccoeff,
			   int *ntempscalecoeff, int *ntempshapecoeff, double *locdsgnmat,
			   double *scaledsgnmat, double *shapedsgnmat, double *tempdsgnmatloc,
			   double *tempdsgnmatscale, double *tempdsgnmatshape, double *weights,
			   double *hess, double *grad);

///////////////////////////////////
//  From circulant.c
//
void circemb(int *nsim, int *ngrid, double *steps, int *dim, int *covmod,
	     double *nugget, double *sill, double *range, double *smooth,
	     double *ans);

///////////////////////////////////
//  From latentVariable.c
//
void latentgev(int *n, double *data, int *nSite, int *nObs, int *covmod, 
	       int *dim, double *distMat, double *dsgnMat, int *nBeta, double *beta,
	       double *sills, double *ranges, double *smooths, double *gevParams,
	       double *hyperSill, double *hyperRange, double *hyperSmooth,
	       double *hyperBetaMean, double *hyperBetaIcov, double *propGev,
	       double *propRanges, double *propSmooths, double *mcLoc,
	       double *mcScale, double *mcShape, double *accRates,
	       double *extRates, int *thin, int *burnin);
void DIC(int *n, int *nSite, int *nObs, double *data, double *chainLoc,
	 double *chainScale, double *chainShape, double *postLoc,
	 double *postScale, double *postShape, double *dic, double *effNpar,
	 double *dbar);
