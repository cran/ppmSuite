/****************************************************************************************
 * Copyright (c) 2014 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit PPMx models when covariate values are missing
 * The trick is to carry a label vector indicating when a covariate
 * is missing.
 *
 ****************************************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <Rmath.h>

#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*****************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
* nobs = number of response values
* ncon = number of continuous covariates
* ncat = number of categorical covariates
* Mcon = nobs x ncon logical vector indicating missing of continuous covariates.
* Mcat = nobs x ncat logical vector indicating missing of categorical covariates.
* Cvec = ncat x 1 vector indicating the number of categories for each categorical covariate.
* PPM = logical indicating if PPM or PPMx should be used.
* M = double indicating value of M associated with cohesion.
* cohesion = integer indicating which cohesion function to use
*   1 - DP stype
*   2 - uniform
* similarity_function = int indicating similarity function to be used
		1 - auxiliary model
		2 - double dipper
		3 - alpha*exp(-variance)
* consim = integer indicating which similarity to use if 1 or 2 is selected above.
    1 - NN
    2 - NNIG
* y = nobs x 1 vector that contains response values
* Xcon = nobs x ncon contiguous double vector that contains continuous covariate values
* Xcat = nobs x ncat contiguous int vector that contains categorical covariate values
* npred = integer indicating number of out of sample predictions
* Xconp = npred x ncon matrix of continuous covariates to make predictions
* Xcatp = npred x ncat matrix of categorical covariates to make predictions
* Mconp = npred x ncon matrix of logicals indicating if continuous covariate is missing
* Mcatp = npred x ncat matrix of logicals indicating if categorical covariate is missing
* simParms = vector containing similarity functions that are ordered in the following way;
* dissimtn = dissimilarity matrix of training covariates for Gower
* dissimtt = dissimilarity matrix of testing covariates for Gower
* calibrate = integer indicating how similarity function should be calibrated
    0 - no calibration
    1 - "calibration" similarity function
    2 - "coarsened" similarity function
* modelPriors = priors from likelihood model
*
* Output:
* mu = nout x ncov matrix of MCMC iterates
* sig2 = nout vector of MCMC iterates
* mu0 = nout vector of MCMC iterates
* sig20 = nout vector of MCMC iterates
* Si = nout vector of MCMC iterates
* nclus = nout vector of MCMC iterates
* like = scalar with lpml value
* lpml = scalar with lpml value
* waic = scalar containing waic value
* ispred
* ppred
* ppredclass

*****************************************************************************************/

void gaussian_ppmx_missing(
          double *y, int *nobs,
          double *Xcon, int *Mcon, int *ncon,
          int *Xcat, int *Mcat,  int *ncat, int *Cvec,
          int *npred,
          double *Xconp, int *Mconp,
          int *Xcatp, int *Mcatp,
          double *M, int *meanModel, double *modelPriors, double *simParms,
          int *PPM, int *cohesion, int *similarity_function, int *consim,
          double *dissimtn, double *dissimtt, int *calibrate,  double *mh,
          int *verbose, int *draws, int *burn, int *thin,
          int *Si, int *nclus, double *mu, double *sig2, double *beta,
          double *mu0, double *sig20,
          double *like, double *waic, double *lpml,
          double *ispred, double *ppred, int *predclass,
          double *rbpred, double *predclass_prob){



  // i - MCMC index
  // ii - MCMC save index
  // j - individual index
  // jj - second player index (for double for loops)
  // c - categorical variable index
  // p - number of covariates index
  // pp - second covariate index
  // k - cluster index
  // t - subset of covariates index

  int i, ii, j, jj, jjj, c, p, pp, k, b, bb;
  int nout = (*draws-*burn)/(*thin);
  int ncov = *ncon + *ncat;

  // determine the categorical variable with max number of categories
  int max_C;
  max_C = Cvec[0];
  for(p = 0; p < (*ncat); p++){
    if(max_C < Cvec[p]) max_C = Cvec[p];
  }
  if(*ncat == 0) max_C = 1;

  // MLEs needed if calibration 1 option is selected.
  // Also, to facilitate selecting values for the similarity
  // function, I standardize the continuous covariates to have
  // zero mean and unit standard deviation.
  double *mnmle = R_VectorInit(*ncon, 0.0);
  double *s2mle = R_VectorInit(*ncon, 1.0);
  double sum, sum2, nobserved;

  if(!(*PPM)){

    for(p = 0; p < *ncon; p++){
      sum = 0.0, sum2=0.0, nobserved=0.0;
      for(j = 0; j < *nobs; j ++){
        if(Mcon[j*(*ncon)+p] == 0){
          sum = sum + Xcon[j*(*ncon) + p];
          sum2 = sum2 + Xcon[j*(*ncon) + p]*Xcon[j*(*ncon) + p];
          nobserved = nobserved + 1;
        }
      }
      for(pp = 0; pp < *npred; pp++){
        if(Mconp[pp*(*ncon)+p] == 0){
          sum = sum + Xconp[pp*(*ncon) + p];
          sum2 = sum2 + Xconp[pp*(*ncon) + p]*Xconp[pp*(*ncon) + p];
          nobserved = nobserved + 1;
        }
      }
      mnmle[p] = sum/nobserved;
      s2mle[p] = sum2/nobserved - mnmle[p]*mnmle[p];
    }
  }

  // Create Xmat to update the covariates from the
  // global regression that is included in the data model
  // Only used if *meanModel == 2;  Note, these are not standardized
  double *fullXmat = R_VectorInit((*nobs)*(ncov), 0.0);
  double *fullXmatp = R_VectorInit((*npred)*(ncov), 0.0);
  if(!(*PPM)){
    for(j = 0; j < *nobs; j++){
	  for(b = 0; b < ncov; b++){
	    if(b < *ncon){
          fullXmat[j*ncov+b] = Xcon[j*(*ncon)+b];
          if(Mcon[j*(*ncon) + b] == 0){
//            Xcon[j*(*ncon) + b] = (Xcon[j*(*ncon)+b] - mnmle[b])/sqrt(s2mle[b]);
	      }
	      if(Mcon[j*(*ncon) + b] == 1){// if missing update auxiliary variable
	        fullXmat[j*(*ncon) + b] = rnorm(mnmle[b], sqrt(s2mle[b]));
	      }
        }
	    if(b >= *ncon){
	      fullXmat[j*ncov+b] = ((double) Xcat[j*(*ncat) + (b-*ncon)]);
	    }
      }
    }


    for(pp = 0; pp < *npred; pp++){
      for(b = 0; b < ncov; b++){
        if(b < *ncon){
          fullXmatp[pp*ncov+b] = Xconp[pp*(*ncon)+b];
          if(Mconp[pp*(*ncon) + b] == 0){
//            Xconp[pp*(*ncon) + b] = (Xconp[pp*(*ncon)+b] - mnmle[b])/sqrt(s2mle[b]);
	      }
	    }
	    if(b >= *ncon){
	      fullXmatp[pp*ncov+b] = ((double) Xcatp[pp*(*ncat)+ (b-*ncon)]);
	    }
      }
    }
  }

  // =====================================================================================
  //
  // Memory vectors to hold a single MCMC iterate
  //
  // =====================================================================================
  double _mu0=0.0, _sig20=1.0;
  double *_beta = R_VectorInit(ncov,0.0);;
  double *_muh = R_VectorInit(*nobs,0.0);;
  double *_sig2h = R_VectorInit(*nobs, 0.5*modelPriors[2]);;
  double *_like = R_VectorInit(*nobs,0.0);;

  int _nclus=0;
  int _Si[*nobs], nh[*nobs];

  // Initialize cluster labels so that two groups exist
  for(j = 0; j < *nobs; j++){
    _Si[j] = rbinom(2, 0.25) + 1;
    nh[j] = 0; // initialize the number of units in each cluster
  }

  // Create vector of cluster sizes
  for(j = 0; j < *nobs; j++){
    nh[_Si[j]-1] = nh[_Si[j]-1] + 1;
  }
  // Count the number of clusters
  for(j = 0; j < *nobs; j++){
    if(nh[j] > 0) _nclus = _nclus + 1;
  }

  // Stuff to compute the posterior predictive
  double *_ispred = R_VectorInit(*nobs, 0.0);

  double *_ppred = R_Vector((*npred));
  double *_rbpred = R_Vector((*npred));
  int _predclass[*npred];
  double *_predclass_prob = R_VectorInit((*npred)*(*nobs+1), 0.0);

  // ===================================================================================
  //
  // scratch vectors of memory needed to update parameters
  //
  // ===================================================================================

  // stuff that I need to update Si (cluster labels);
  int iaux, auxint, nhxtmp;
  int nhctmp[max_C];
  double auxreal, sumxtmp, sumx2tmp, npdN,npdY,npd, mn, s2k, xcontmp, uu, xb=0.0, vb=0.0;
  double mudraw, sdraw, maxph, denph, cprobh;
  double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw;
  double lgcont,lgcatt;

  double *ph = R_VectorInit(*nobs, 0.0);
  double *probh = R_VectorInit(*nobs, 0.0);

  double *gtilN = R_VectorInit((*nobs+1),0.0);
  double *gtilY = R_VectorInit((*nobs+1),0.0);
  double *lgtilN = R_VectorInit((*nobs+1),0.0);
  double *lgtilY = R_VectorInit((*nobs+1),0.0);
  double sgY, sgN,  lgtilNk, lgtilYk, maxgtilY, maxgtilN;

  double *sumx = R_VectorInit((*nobs)*(*ncon),0.0);
  double *sumx2 = R_VectorInit((*nobs)*(*ncon),0.0);
  int nhx[(*nobs)*(*ncon)];
  int nhc[(*nobs)*(*ncat)*max_C];


  // stuff I need to update muh
  double *sumy = R_VectorInit((*nobs),0.0);
  double *sumy2 = R_VectorInit((*nobs),0.0);
  double mstar, s2star;

  // stuff I need to update sig2h and sig20
  double os0, ns0, lln, llo, llr, nsig, osig;

  // stuff I need to update mu0
  double summu, summu2;

  //stuff I need to update beta;  Only used with meanModel = 2
  double ld;
  double *scr1 = R_Vector(*nobs);
  double *scr2 = R_Vector(*nobs);
  double *sumyxbt = R_VectorInit(ncov,0.0);
  double *sumXXp = R_VectorInit(ncov*ncov, 0.0);
  double *Sstar = R_VectorInit(ncov*ncov, 0.0);
  double *Mstar = R_VectorInit(ncov, 0.0);


  // Stuff to compute lpml, likelihood, WAIC, and Rao-Blackwellized density values
  double _lpml, elppdWAIC;
  double *CPOinv = R_VectorInit(*nobs, 0.0);
  double *mnlike = R_VectorInit(*nobs, 0.0);
  double *mnllike = R_VectorInit(*nobs, 0.0);


  // ===================================================================================
  //
  // Prior distribution parameter values
  //
  // ===================================================================================

  // priors for mu0
  double m = modelPriors[0]; double s2 = modelPriors[1];

  // prior values for sig2h and sig20;
  double smin=0, smax=modelPriors[2];
  double s0min=0, s0max=modelPriors[3];

  // priors for beta.  Only used if meanModel==2
  double mb=0; double s2b = 100^2;

  // DP weight parameter
  double Mdp = *M;

  // Similarity function parameters
  // dirichlet denominator parameter
  double *dirweights = R_VectorInit(max_C, simParms[5]);

  double m0=simParms[0];
  double s20=simParms[1];
  double v=simParms[2];
  double k0=simParms[3];
  double nu0=simParms[4];

  // For the variance similarity function
  double alpha = simParms[6];

  // ===================================================================================
  //
  // Initialize the cluster-specific sufficient statistics for continuous covariates
  // and categorical covariates.
  //
  // ===================================================================================
  if(!(*PPM)){
    for(j=0;j<*nobs;j++){
      mnlike[j] = 0.0;
      mnllike[j] = 0.0;
      for(p=0;p<*ncon;p++){
        sumx[j*(*ncon) + p] = 0.0;
        sumx2[j*(*ncon) + p] = 0.0;
        nhx[j*(*ncon) + p] = 0;
      }
      for(p=0;p<*ncat;p++){
        for(c=0; c<max_C; c++){
          nhc[(j*(*ncat) + p)*max_C + c] = 0;
        }
      }
    }

    // Fill in cluster-specific sufficient statistics based on first partition
    for(j = 0; j < *nobs; j++){
      for(p=0; p<*ncon; p++){
        if(Mcon[j*(*ncon)+p] == 0){
          sumx[(_Si[j]-1)*(*ncon) + p] = sumx[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p];
          sumx2[(_Si[j]-1)*(*ncon) + p] = sumx2[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          nhx[(_Si[j]-1)*(*ncon) + p] = nhx[(_Si[j]-1)*(*ncon) + p] + 1;
        }
      }
      for(p=0; p<*ncat; p++){
        if(Mcat[j*(*ncat)+p] == 0){
          nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
            nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] + 1;
        }
      }
    }
  }


  // M-H tuning parameters
  double csigSIG0=mh[0], csigSIG=mh[1];

  if(*verbose){
    Rprintf("nobs = %d\n", *nobs);
    Rprintf("ncon = %d\n", *ncon);
    Rprintf("ncat = %d\n", *ncat);
    Rprintf("ncov = %d\n", ncov);
    Rprintf("npred = %d\n", *npred);
    if(*ncat > 0 ) RprintIVecAsMat("Cvec",Cvec, 1, *ncat);
    Rprintf("Prior values: m = %.2f, s2 = %.2f smax = %.2f, s0max = %.2f\n", m, s2, smax, s0max);
    Rprintf("M = %f\n", Mdp);
    Rprintf("Simlarity values: m0 = %.2f, s20 = %.2f v = %.2f, k0 = %.2f, nu0 = %.2f, a0 = %.2f, alpha = %.2f\n", m0, s20, v, k0, nu0, simParms[5], alpha);
    if(*nobs > 50) RprintVecAsMat("First fifty response values",y, 1, 50);
  }


  ii = 0;

  // =============================================================================================
  //
  // start of the mcmc algorithm;
  //
  // =============================================================================================

  double calc_time = 0.0;
  clock_t  begin = clock();


  for(i = 0; i < *draws; i++){

//    if(*verbose & ((i+1) % 1000 == 0)){
    if(*verbose){
//	      time_t now;
//	      time(&now);

//	      Rprintf("mcmc iter = %d ===================================================== \n", i+1);
//        Rprintf("%s", ctime(&now));

      clock_t ith_iterate = clock();
      calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

      Rprintf("Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*draws))*100.0, calc_time);
//      fflush(stdout);
//      Rprintf("Progress:%.1f%\r", ((double) (i+1) / (double) (*draws))*100.0);


    }

    //////////////////////////////////////////////////////////////////////////////////
    //
    // update the cluster labels using the polya urn scheme of
    // algorithm 8 found in  Radford Neal's
    //	"Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
    //	paper.
    //
    // To accomodate missing here, if a covariate is missing for individual i, then
    // I am simply not evaluating the similarity function for that individual.
    // This is done by dragging around a logical vector that indicates missing or not.
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(j = 0; j < *nobs; j++){

      if(nh[_Si[j]-1] > 1){

        // Observation belongs to a non-singleton ...
        nh[_Si[j]-1] = nh[_Si[j]-1] - 1;


        if(!(*PPM)){
          // need to reduce the sumx sumx2 to by removing the $i$th individual
          for(p = 0; p < *ncon; p++){
            if(Mcon[j*(*ncon)+p] == 0){
              sumx[(_Si[j]-1)*(*ncon) + p] = sumx[(_Si[j]-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
              sumx2[(_Si[j]-1)*(*ncon) + p] = sumx2[(_Si[j]-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
              nhx[(_Si[j]-1)*(*ncon) + p] = nhx[(_Si[j]-1)*(*ncon) + p] - 1;
            }
          }
          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){
            if(Mcat[j*(*ncat)+p] == 0){
              nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
                 nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;
            }
          }
        }

      }else{

        // Observation is a member of a singleton cluster ...
        iaux = _Si[j];

        if(iaux < _nclus){

          // Need to relabel clusters.  I will do this by swapping cluster labels
          // _Si[j] and _nclus along with cluster specific parameters and
          // covariate values;

          // All members of last cluster will be assigned subject i's cluster label
          for(jj = 0; jj < *nobs; jj++){
            if(_Si[jj] == _nclus){
              _Si[jj] = iaux;
            }
          }

          _Si[j] = _nclus;

          // The following steps swaps order of cluster specific parameters
          // so that the newly labeled subjects from previous step retain
          // their correct cluster specific parameters

          auxreal = _muh[iaux-1];
          _muh[iaux-1] = _muh[_nclus-1];
          _muh[_nclus-1] = auxreal;

          auxreal = _sig2h[iaux-1];
          _sig2h[iaux-1] = _sig2h[_nclus-1];
          _sig2h[_nclus-1] = auxreal;


          // the number of members in cluster is also swapped with the last
          nh[iaux-1] = nh[_nclus-1];
          nh[_nclus-1] = 1;

          if(!(*PPM)){
            // need to swap sumx and sumx2
            // here missingness only comes into play with how many subjects are in each
            // group.  Hence need to include nhx.
            for(p = 0; p < *ncon; p++){
              auxreal = sumx[(iaux-1)*(*ncon) + p];
              sumx[(iaux-1)*(*ncon) + p] = sumx[(_nclus-1)*(*ncon) + p];
              sumx[(_nclus-1)*(*ncon) + p] = auxreal;

              auxreal = sumx2[(iaux-1)*(*ncon) + p];
              sumx2[(iaux-1)*(*ncon) + p] = sumx2[(_nclus-1)*(*ncon) + p];
              sumx2[(_nclus-1)*(*ncon) + p] = auxreal;

              auxint = nhx[(iaux-1)*(*ncon) + p];
              nhx[(iaux-1)*(*ncon) + p] = nhx[(_nclus-1)*(*ncon) + p];
              nhx[(_nclus-1)*(*ncon) + p] = auxint;

            }

            // need to swap nhc as well
            for(p = 0; p < *ncat; p++){
              for(c=0; c<max_C; c++){
                auxint = nhc[((iaux-1)*(*ncat) + p)*(max_C) + c];
                nhc[((iaux-1)*(*ncat) + p)*(max_C) + c] = nhc[((_nclus-1)*(*ncat) + p)*(max_C) + c];
                nhc[((_nclus-1)*(*ncat) + p)*(max_C) + c] = auxint;
              }
            }
          }

        }

      	// Now remove the ith obs
      	nh[_nclus-1] = nh[_nclus-1] - 1;

        // need to reduce the sumx sumx2 for the $i$th subject

        if(!(*PPM)){
          for(p = 0; p < *ncon; p++){
            if(Mcon[j*(*ncon)+p] == 0){
              sumx[(_nclus-1)*(*ncon) + p] = sumx[(_nclus-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
              sumx2[(_nclus-1)*(*ncon) + p] = sumx2[(_nclus-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
              nhx[(_nclus-1)*(*ncon) + p] = nhx[(_nclus-1)*(*ncon) + p] - 1;
            }
          }
          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){
            if(Mcat[j*(*ncat)+p] == 0){
               nhc[((_nclus-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
                  nhc[((_nclus-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;
             }
          }
        }

  		// finally reduce the number of clusters
      	_nclus = _nclus - 1;

      }

      // The atoms have been relabeled if necessary and now we need to
      // update Si.


      if(*meanModel==2){
        // create x'beta but only include observed x's when doing this
        xb = 0.0, vb = 0.0;
        for(b = 0; b < ncov; b++){
          if(b < *ncon){
            if(Mcon[j*(*ncon)+b] == 0){// not missing
              xb = xb + (fullXmat[j*ncov+b] - mnmle[b])*_beta[b];
            }
            if(Mcon[j*(*ncon)+b] == 1){// missing covariate
              vb = vb + s2mle[b]*_beta[b]*_beta[b];
            }
          }
          if(b >= *ncon){ // I don't think categorical covariates are possible here.
            if(Mcat[j*(*ncat) + (b-*ncon)] == 0){
              xb = xb + fullXmat[j*ncov+b]*_beta[b];
            }
          }
        }
      }

      // Begin the cluster probabilities

      for(k=0; k<_nclus; k++){
        lgconY = 0.0;
      	lgconN = 0.0;
      	lgcatY = 0.0;
      	lgcatN = 0.0;

      	if(!(*PPM)){

      	  // start by calculating similarity for continuous covariates
      	  for(p=0; p<(*ncon); p++){

            sumxtmp = sumx[k*(*ncon) + p];
            sumx2tmp = sumx2[k*(*ncon) + p];
            nhxtmp = nhx[k*(*ncon) + p];

      	    if(*similarity_function==1){ // Auxilliary
      	      if(*consim==1){ // NN
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 0, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
      	      if(*consim==2){// NNIG
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 0, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }

      	    }
      	    if(*similarity_function==2){ //Double Dipper
      	      if(*consim==1){
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 1, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
      	      if(*consim==2){
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 1, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
      	    }
      	    if(*similarity_function==3){ // variance
      	      lgcont = gsimconEV(sumxtmp, sumx2tmp, nhxtmp, alpha,1);
      	      lgconN = lgconN + lgcont;
      	    }


      	    // now add jth individual back;
      	    if(Mcon[j*(*ncon)+p] == 0){
      	      sumxtmp = sumxtmp + Xcon[j*(*ncon)+p];
      	      sumx2tmp = sumx2tmp + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
      	      nhxtmp = nhxtmp + 1;
            }


      	    if(*similarity_function==1){ // Auxilliary
      	      if(*consim==1){
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 0, 0, 1);
      	        lgconY = lgconY + lgcont;
      	      }
      	      if(*consim==2){
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 0, 0, 1);
      	        lgconY = lgconY + lgcont;
      	      }
      	    }
      	    if(*similarity_function==2){ //Double Dipper
      	      if(*consim==1){
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 1, 0, 1);
      	        lgconY = lgconY + lgcont;
      	      }
      	      if(*consim==2){
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 1, 0, 1);
      	        lgconY = lgconY + lgcont;
      	      }
      	    }
      	    if(*similarity_function==3){ // variance
      	      lgcont = gsimconEV(sumxtmp, sumx2tmp, nhxtmp, alpha,1);
      	      lgconY = lgconY + lgcont;
      	    }

      	  }


      	  // Now calculate similarity for the categorical covariates
      	  for(p=0; p<(*ncat); p++){

            for(c = 0; c < max_C; c++){
              nhctmp[c] = nhc[(k*(*ncat) + p)*(max_C) + c];
            }


      	    if(*similarity_function==1){ // Auxiliary
      	      lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
      	      lgcatN = lgcatN + lgcatt;
      	  	}
      	  	if(*similarity_function==2){// Double dipper
      	      lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
      	      lgcatN = lgcatN + lgcatt;
      	  	}
      	  	if(*similarity_function==3){// Using Entropy instead of variance here
      	  	  lgcatt = 0.0;
      	  	  for(c=0;c<Cvec[p];c++){
          	    if(nhctmp[c]==0){
          	      lgcatt = lgcatt + 0;
          	    }else{
          	      lgcatt = lgcatt + -((double) nhctmp[c]/(double) nh[k])*
                         			  (log((double) nhctmp[c]/(double) nh[k])/log(2));
          	    }
      	  	  }
      	  	  lgcatN = lgcatN + -(alpha)*lgcatt;
      	  	}

            // include the categorical covariate in the kth cluster
      	    if(Mcat[j*(*ncat)+p] == 0){
      	  	  nhctmp[Xcat[j*(*ncat)+p]] = nhctmp[Xcat[j*(*ncat)+p]] + 1;
            }

      	  	if(*similarity_function==1){
      	  	  lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
      	  	  lgcatY = lgcatY + lgcatt;
      	  	}
      	  	if(*similarity_function==2){
      	  	  lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
      	  	  lgcatY = lgcatY + lgcatt;
      	  	}
      	  	if(*similarity_function==3){// Using Entropy instead of variance here
      	  	  lgcatt = 0.0;
      	  	  for(c=0;c<Cvec[p];c++){
          	    if(nhctmp[c]==0){
          		  lgcatt = lgcatt + 0;
          		}else{
          		  lgcatt = lgcatt + -((double) nhctmp[c]/(double) nh[k]+1)*
          							  (log((double) nhctmp[c]/(double) nh[k]+1)/log(2));
          		}
      	  	  }
      	  	  lgcatY = lgcatY + -(alpha)*lgcatt;
      	  	}
      	  }


      	  gtilY[k] = lgconY + lgcatY;
      	  gtilN[k] = lgconN + lgcatN;

      	  //////////////////////////////////////////////////////////
      	  // Gower Compute similarity values for gower dissimilarity
      	  //////////////////////////////////////////////////////////
      	  if(*similarity_function==4){ // I don't think this one needs to be callibrated.
      	    npd=0.0;
      	  	lgconY = 0.0;
      	  	lgconN = 0.0;
      	    for(jj = 0; jj < *nobs; jj++){

      	  	  if((_Si[jj] == k+1) & (jj != j)){
      	  	    lgconY = lgconY + dissimtn[jj*(*nobs) + j];
      	  	    for(jjj = 0; jjj < jj; jjj++){
      	  		  if((_Si[jjj] == k+1) & (jjj != j)){
      	  		    lgconN = lgconN + dissimtn[jj*(*nobs) + jjj];
      	  		    lgconY = lgconY + dissimtn[jj*(*nobs) + jjj];
      	  		    npd=npd+1;
      	  		  }
      	  		}
      	  	  }
      	  	}

      	  	npdN = nh[k]*(nh[k]-1)/2;
      	  	if(npdN == 0) npdN = 1.0;
      	  	npdY = (nh[k]+1)*(nh[k])/2;

      	  	// This is cluster-mean Gower dissimilarity
      	  	lgconN = -(alpha)*lgconN/(npdN);
      	  	lgconY = -(alpha)*lgconY/(npdY);

      	  	// I will use cluster-total gower dissimilarity
      	  	lgconN = -(alpha)*lgconN;
      	  	lgconY = -(alpha)*lgconY;

      	  }
      	  //////////////////////////////////////////////////////////
      	  // End of Gower similarity
      	  //////////////////////////////////////////////////////////
      	} // THIS ENDS THE PPMX PART.

        mn = _muh[k];
        s2k = _sig2h[k];
        if(*meanModel==2){
          mn = _muh[k] + xb;
          s2k = _sig2h[k] + vb;
        }

      	// Compute the unnormalized cluster probabilities
      	// Note that if PPMx = FALSE then
      	// lgcatY = lgcatN = lgconY = lgconN = 0;
      	ph[k] = dnorm(y[j], mn, sqrt(s2k), 1) +
      		    	log((double) nh[k]) + // cohesion part
      		        lgcatY - lgcatN + // Categorical part only nonzero if PPMx=TRUE
      				lgconY - lgconN;  // Continuous part only nonzero if PPMx=TRUE



      	if(*calibrate == 2){
      		ph[k] = dnorm(y[j], mn, sqrt(_sig2h[k]), 1) +
      	            	log((double) nh[k]) +
      	            	(1/((double)*ncon + (double)*ncat))* //PPMx part only nonzero if PPMx=TRUE
      	            	(lgcatY + lgconY - lgcatN - lgconN);
      	}
      	// Uniform cohesion
      	if(*cohesion==2){
      		ph[k] = ph[k] - log((double) nh[k]);
      	}


      }


      // Need to consider allocating subject to new cluster
      mudraw = rnorm(_mu0, sqrt(_sig20));
      sdraw = runif(smin, smax);

      lgcondraw = 0.0;
      lgcatdraw = 0.0;
      if(!(*PPM)){

        // similarity for continuous covariate
        for(p=0;p<(*ncon);p++){

          if(Mcon[j*(*ncon)+p] == 0){
            xcontmp = Xcon[j*(*ncon)+p];
            if(*similarity_function==1){ // Auxilliary
              if(*consim==1){
                lgcont = gsimconNN(m0,v,s20,xcontmp,xcontmp*xcontmp, mnmle[p],1,0,0, 1);
                lgcondraw = lgcondraw + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp, xcontmp*xcontmp,mnmle[p],s2mle[p], 1, 0,0, 1);
                lgcondraw = lgcondraw + lgcont;
              }
            }
            if(*similarity_function==2){ // Double Dipper
              if(*consim==1){
          	    lgcont = gsimconNN(m0,v,s20,xcontmp,xcontmp*xcontmp, mnmle[p], 1, 1, 0, 1);
          	    lgcondraw = lgcondraw + lgcont;
          	  }
          	  if(*consim==2){
          	    lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp, xcontmp*xcontmp,mnmle[p],s2mle[p], 1, 1, 0, 1);
          	    lgcondraw = lgcondraw + lgcont;
          	  }
            }
            if(*similarity_function==3){ // Variance
              lgcont = gsimconEV(xcontmp, xcontmp*xcontmp, 1,alpha,1);
          	  lgcondraw = lgcondraw + lgcont;
            }
          }
          if(*similarity_function==4){ // Dissimilarity
            lgcondraw = -(alpha)*0;
          }

        }
        // similarity for categorical covariate
        for(p=0;p<(*ncat);p++){
          for(c=0;c<Cvec[p];c++){nhctmp[c] = 0;}

          if(Mcat[j*(*ncat)+p] == 0){ // only consider this if observation is not missing.
            nhctmp[Xcat[j*(*ncat)+p]] = 1;


            if(*similarity_function==1){
              lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
          	  lgcatdraw = lgcatdraw + lgcatt;
            }
            if(*similarity_function==2){
          	  lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
          	  lgcatdraw = lgcatdraw + lgcatt;
            }
            if(*similarity_function==3){
          	  lgcatdraw = lgcatdraw + -(alpha)*0;
            }

          }

          if(*similarity_function==4){
            lgcatdraw = -(alpha)*0;
          }
        }
        gtilY[_nclus] = lgcondraw + lgcatdraw;
        gtilN[_nclus] = lgcondraw + lgcatdraw;
      } // THIS ENDS THE PPMX PART.

      mn = mudraw;
      s2k = sdraw*sdraw;
      if(*meanModel == 2){
        mn = mudraw + xb;
        s2k = sdraw*sdraw + vb;
      }

      // Note that if PPMx = FALSE, then
      // lgcondraw = lgcondraw =  0;


      ph[_nclus] = dnorm(y[j],mn,sqrt(s2k),1) +
                       	log(Mdp) +
                       	lgcondraw +
                       	lgcatdraw;

      if(*calibrate==2){
      	ph[_nclus] = dnorm(y[j],mn,sqrt(s2k),1) +
      	                 	log(Mdp) +
      	                 	(1/((double)*ncon + (double)*ncat))*(lgcondraw + lgcatdraw);
      }

      if(*cohesion==2){
      	ph[_nclus] = ph[_nclus] - log(Mdp);
      }

      /////////////////////////////////////////////////////////////////////////////
      // This is the calibration used when the similarity is normalized
      /////////////////////////////////////////////////////////////////////////////
      if((*calibrate==1) & (*PPM != 1)){
        maxgtilN = gtilN[0];
      	maxgtilY = gtilY[0];
      	for(k=1; k < _nclus+1; k++){

      	  if(maxgtilN < gtilN[k]) maxgtilN = gtilN[k];

      	  if(k < _nclus){
      	    if(maxgtilY < gtilY[k]) maxgtilY = gtilY[k];
      	  }
      	}

      	sgY=0.0;
      	sgN=0.0;
      	for(k=0; k<_nclus+1; k++){

      	  lgtilN[k] = gtilN[k] - maxgtilN;
      	  sgN = sgN + exp(lgtilN[k]);

      	  if(k < _nclus){// If x is included in an existing cluster in cannot be a singleton
      	    lgtilY[k] = gtilY[k] - maxgtilY;
      	    sgY = sgY + exp(lgtilY[k]);
      	  }
      	}

      	// Calibrate the unnormalized cluster probabilities
      	for(k=0; k<_nclus; k++){
      		lgtilNk = lgtilN[k] - log(sgN);
      		lgtilYk = lgtilY[k] - log(sgY);

            mn = _muh[k];
            s2k = _sig2h[k];
            if(*meanModel == 2){
              mn = _muh[k] + xb;
              s2k = _sig2h[k] + vb;
            }

      		ph[k] = dnorm(y[j], mn, sqrt(s2k), 1) +
              	    	log((double) nh[k]) +  // Cohesion part
      				    lgtilYk - lgtilNk; //This takes into account both cont and cat vars

      		if(*cohesion==2){
      	    	ph[k] = ph[k] - log((double) nh[k]);
          	}

      	}

      	// calibration for a singleton
        mn = mudraw;
        s2k = sdraw*sdraw;
        if(*meanModel == 2){
          mn = mudraw + xb;
          s2k = sdraw*sdraw + vb;
        }

      	ph[_nclus] = dnorm(y[j],mn,sqrt(s2k),1) +
                       	    log(Mdp) +
      					    lgtilN[_nclus] - log(sgN);

      	if(*cohesion==2){// Note with a uniform cohesion, for a new cluster
      	                 // the value of log(c({nclus_iter}}) = log(1) = 0;
      		ph[_nclus] = ph[_nclus] - log(Mdp);
      	}

      }
      /////////////////////////////////////////////////////////////////////////////
      // End of calibration used when the similarity is normalized
      /////////////////////////////////////////////////////////////////////////////


      maxph = ph[0];
      for(k = 1; k < _nclus+1; k++){
      	if(maxph < ph[k]) maxph = ph[k];
      }

      denph = 0.0;
      for(k = 0; k < _nclus+1; k++){
      	ph[k] = exp(ph[k] - maxph);
      	denph = denph + ph[k];
      }

      for(k = 0; k < _nclus+1; k++){
      	probh[k] = ph[k]/denph;
      }


      uu = runif(0.0,1.0);

      cprobh= 0.0;
      iaux = _nclus+1;
      for(k = 0; k < _nclus+1; k++){
      	cprobh = cprobh + probh[k];
      	if (uu < cprobh){
      		iaux = k+1;
      		break;
      	}
      }


      if(iaux <= _nclus){

      	_Si[j] = iaux;
      	nh[_Si[j]-1] = nh[_Si[j]-1] + 1;

      }else{

      	_nclus = _nclus + 1;
      	_Si[j] = _nclus;
      	nh[_Si[j]-1] = 1;

      	_muh[_Si[j]-1] = mudraw;
      	_sig2h[_Si[j]-1] = sdraw*sdraw;
      }

      // need to now add the xcon to the cluster to which it was assigned;
      if(!(*PPM)){
        for(p = 0; p < *ncon; p++){
          if(Mcon[j*(*ncon)+p] == 0){
            sumx[(_Si[j]-1)*(*ncon) + p] = sumx[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p];
            sumx2[(_Si[j]-1)*(*ncon) + p] = sumx2[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
            nhx[(_Si[j]-1)*(*ncon) + p] = nhx[(_Si[j]-1)*(*ncon) + p] + 1;
          }
        }
        // need to now add the xcat to the cluster to which it was assigned;
        for(p = 0; p < *ncat; p++){
          if(Mcat[j*(*ncat)+p] == 0){
            nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
               nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] + 1;
          }
        }
      }

    }


    //////////////////////////////////////////////////////////////////////////////////
    //
    // update sumy and sumy2
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < _nclus; k++){
      sumy[k] = 0.0;
      sumy2[k] = 0.0;
    }


    for(j = 0; j < *nobs; j++){


      // sum with-in this for loop to save some computing time;
      // This are used to update muh which is next
      if(*meanModel==1){
        sumy[_Si[j]-1] = sumy[_Si[j]-1] + y[j];
        sumy2[_Si[j]-1] = sumy2[_Si[j]-1] + y[j]*y[j];
      }
      if(*meanModel==2){ // This needs to be updated using Matt's stuff
        xb = 0.0;
        for(b = 0; b < ncov; b++){
          if(b < *ncon){
            xb = xb + (fullXmat[j*ncov+b] - mnmle[b])*_beta[b]; // this has auxiliary values
          }
          if(b >= *ncon){
            xb = xb + fullXmat[j*ncov+b]*_beta[b]; // not sure yet how to ``center'' categorical covariate
          }
        }

        sumy[_Si[j]-1] = sumy[_Si[j]-1] + (y[j] - xb);
        sumy2[_Si[j]-1] = sumy2[_Si[j]-1] + (y[j] - xb)*(y[j] - xb);

      }
    }



    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update sig2h  cluster specific variance parameters with metropolis and
    // Uniform prior.
    //
    //////////////////////////////////////////////////////////////////////////////////
    for(k = 0; k < _nclus; k++){
      osig = sqrt(_sig2h[k]);
      nsig = rnorm(osig,csigSIG);

      if((nsig > 0) & (nsig < smax)){

        llo = -(nh[k])*log(osig) - 1/(2*osig*osig)*(sumy2[k] - 2*_muh[k]*sumy[k] + nh[k]*_muh[k]*_muh[k]);
        lln = -(nh[k])*log(nsig) - 1/(2*nsig*nsig)*(sumy2[k] - 2*_muh[k]*sumy[k] + nh[k]*_muh[k]*_muh[k]);

        llo = llo + dunif(osig, smin, smax, 1);
        lln = lln + dunif(nsig, smin, smax, 1);

        llr = lln - llo;
        uu = runif(0,1);

        if(log(uu) < llr){
        	_sig2h[k] = nsig*nsig;
        }
      }
    }



    //////////////////////////////////////////////////////////////////////////////////
    //
    // update muh's cluster specific means
    //
    //////////////////////////////////////////////////////////////////////////////////
    summu = 0.0, summu2=0.0;
    for(k = 0; k < _nclus; k++){

      s2star = 1/((double) nh[k]/_sig2h[k] + 1/_sig20);
      mstar = s2star*( (1/_sig2h[k])*sumy[k] + (1/_sig20)*_mu0);

      _muh[k] = rnorm(mstar, sqrt(s2star));

      summu = summu + _muh[k];  // This is used in the updating of mu0 and sig20
      summu2 = summu2 + _muh[k]*_muh[k];  // This is used in the updating of sig20
    }


    //////////////////////////////////////////////////////////////////////////////////////
    //
    // Update mu0  prior mean of muh
    //
    //////////////////////////////////////////////////////////////////////////////////////

    s2star = 1/(((double) _nclus/_sig20) + (1/s2));
    mstar = s2star*((1/_sig20)*summu + (1/s2)*m);

    _mu0 = rnorm(mstar, sqrt(s2star));


    //////////////////////////////////////////////////////////////////////////////////////
    //
    // Update sig20  prior variance of muh
    //
    //////////////////////////////////////////////////////////////////////////////////////
    os0 = sqrt(_sig20);
    ns0 = rnorm(os0,csigSIG0);
    if((ns0 > 0) & (ns0 < s0max)){

      llo = -(_nclus)*log(os0) - 1/(2*os0*os0)*(summu2 - 2*_mu0*summu + _nclus*_mu0*_mu0);
      lln = -(_nclus)*log(ns0) - 1/(2*ns0*ns0)*(summu2 - 2*_mu0*summu + _nclus*_mu0*_mu0);

      llo = llo + dunif(os0, s0min, s0max, 1);
      lln = lln + dunif(ns0, s0min, s0max, 1);

      llr = lln - llo;
      uu = runif(0,1);

      if(log(uu) < llr){
      	_sig20 = ns0*ns0;
      }

    }



    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update beta for each cluster configuration.
    // they will need to be saved in beta_iter.
    //
    //////////////////////////////////////////////////////////////////////////////////

    if(*meanModel==2){


//      Rprintf("updating beta \n");

      for(b = 0; b < ncov; b++){
        sumyxbt[b] = 0.0;
      	for(bb = 0; bb < ncov; bb++){
      		sumXXp[b*ncov+bb] = 0.0;
      	}
      }


      for(j = 0; j < *nobs; j++){
      	for(b = 0; b < ncov; b++){

      	  sumyxbt[b] = sumyxbt[b] +  (1/_sig2h[_Si[j]-1])*(fullXmat[j*ncov+b] - mnmle[b])*
      					                                  (y[j] - _muh[_Si[j]-1]);

      	  for(bb = 0; bb < ncov; bb++){

            sumXXp[b*ncov+bb] = sumXXp[b*ncov+bb] + (1/_sig2h[_Si[j]-1])*
      						                         (fullXmat[j*ncov+b] - mnmle[b])
      						                        *(fullXmat[j*ncov+bb] - mnmle[bb]);

      	  }
        }
      }


      for(b = 0; b < ncov; b++){
      	sumyxbt[b] = sumyxbt[b] + 1/s2b*mb;
      	// If this is commented out it implies that I am using a p(beta) propto 1 prior
      	// and updating mubeta and s2beta does nothing as they are not included in the
      	// prior.
      	for(bb = 0; bb < ncov; bb++){
          Sstar[b*ncov + bb] = sumXXp[b*ncov+bb];

      	  // If this is commented out it implies that I am using a p(beta) propto 1 prior
          if(b==bb) Sstar[b*ncov + bb] = sumXXp[b*ncov+bb] + 1/s2b;

      	}
      }


     cholesky(Sstar, ncov, &ld);
     inverse_from_cholesky(Sstar, scr1, scr2, ncov); //Sstar is now an inverse;


      matrix_product(Sstar, sumyxbt, Mstar, ncov, 1, ncov);

      cholesky(Sstar, ncov , &ld);

      ran_mvnorm(Mstar, Sstar, ncov, scr1, scr2);


      for(b = 0; b < ncov; b++){
        _beta[b]=scr2[b];
      }
//      RprintVecAsMat("beta", _beta, 1, ncov);


      // Update the auxiliar x values
      for(j = 0; j < *nobs; j++){
	    for(b = 0; b < ncov; b++){
	      if(b < *ncon){
	        if(Mcon[j*(*ncon) + b] == 1){// if missing update auxiliary variable
	          fullXmat[j*(*ncon) + b] = rnorm(mnmle[b], sqrt(s2mle[b]));
	        }
          }
	      if(b >= *ncon){ // not sure what to do here;
	        fullXmat[j*ncov+b] = ((double) Xcat[j*(*ncat) + (b-*ncon)]);
	      }
        }
      }
    }


    ////////////////////////////////////////////////////////////////////////////////////////////
    //
    // in sample prediction to assess model fit
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    if((i >= (*burn)) & ((i) % *thin ==0)){

      for(j = 0; j < *nobs; j++){

        mn = _muh[_Si[j]-1];
        s2k = _sig2h[_Si[j]-1];

        if(*meanModel==2){
          // create x'beta but only include observed x's when doing this
          xb = 0.0, vb = 0.0;
          for(b = 0; b < ncov; b++){
            if(b < *ncon){
              if(Mcon[j*(*ncon)+b] == 0){// not missing
                xb = xb + (fullXmat[j*ncov+b] - mnmle[b])*_beta[b];
              }
              if(Mcon[j*(*ncon)+b] == 1){// missing covariate
                vb = vb + s2mle[b]*_beta[b]*_beta[b];
              }
            }
            if(b >= *ncon){ // I don't think categorical covariates are possible here.
              if(Mcat[j*(*ncat) + (b-*ncon)] == 0){
                xb = xb + fullXmat[j*ncov+b]*_beta[b];
              }
            }
          }
          mn = _muh[_Si[j]-1] + xb;
          s2k = _sig2h[_Si[j]-1] + vb;
        }

        _ispred[j] = rnorm(mn, sqrt(s2k));


        /////////////////////////////////////////////
        //
        // Compute the CPO and lpml using the mixture
        //
        /////////////////////////////////////////////

        _like[j] = dnorm(y[j], mn, sqrt(s2k), 0);

        // These are needed for WAIC
        mnlike[j] = mnlike[j] + (_like[j])/(double) nout;
        mnllike[j] = mnllike[j] + log(_like[j])/(double) nout;

        CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/_like[j]);

      }

    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    //
    // out of sample prediction using posterior predictive?
    //
    ////////////////////////////////////////////////////////////////////////////////////////////

    if((i >= (*burn)) & ((i) % *thin ==0)){
      for(pp = 0; pp < *npred; pp++){

        for(k = 0; k < _nclus; k++){

      	  lgconN=0.0, lgconY=0.0;
      	  lgcatN=0.0, lgcatY=0.0;

      	  if(!(*PPM)){

      	    for(p=0; p<(*ncon); p++){
              sumxtmp = sumx[k*(*ncon) + p];
              sumx2tmp = sumx2[k*(*ncon) + p];
              nhxtmp = nhx[k*(*ncon) + p];

      	  	  if(*similarity_function==1){
      	  	    if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 0, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx2tmp, sumx2[k], mnmle[p], s2mle[p], nhxtmp, 0, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==2){
      	  	  	if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumx2tmp, sumx2tmp, mnmle[p], nhxtmp, 1, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 1, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==3){
      	  	    lgcont = gsimconEV(sumxtmp, sumx2tmp, nhxtmp, alpha, 1);
      	  	  	lgconN = lgconN + lgcont;
      	  	  }

      	  	  // now add ppth prediction to cluster;
      	  	  if(Mconp[pp*(*ncon)+p] == 0){
      	  	    sumxtmp = sumxtmp + Xconp[pp*(*ncon)+p];
      	  	    sumx2tmp = sumx2tmp + Xconp[pp*(*ncon)+p]*Xconp[pp*(*ncon)+p];
      	  	    nhxtmp = nhxtmp + 1;
              }

      	  	  if(*similarity_function==1){ // Auxilliary
      	  	  	if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 0, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 0, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==2){ // Double Dipper
      	  	  	if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nhxtmp, 1, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nhxtmp, 1, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==3){ // Variance
      	  	    lgcont = gsimconEV(sumxtmp, sumx2tmp, nhxtmp,alpha, 1);
      	  	    lgconY = lgconY + lgcont;
      	  	  }
      	  	} // This ends the loop through ncon continuous covariates


      	  	for(p=0; p<(*ncat); p++){

              for(c = 0; c < max_C; c++){
                nhctmp[c] = nhc[(k*(*ncat) + p)*(max_C) + c];
              }


      	  	  if(*similarity_function==1){
      	  	    lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
      	  	  	lgcatN = lgcatN + lgcatt;
      	  	  }
      	  	  if(*similarity_function==2){
      	  	  	lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
      	  	  	lgcatN = lgcatN + lgcatt;
      	  	  }
      	  	  if(*similarity_function==3){
      	  	  	lgcatt = 0.0;
                for(c=0;c<Cvec[p];c++){
                  if(nhctmp[c]==0){
                    lgcatt = lgcatt + 0;
                  }else{
                    lgcatt = lgcatt + -((double) nhctmp[c]/(double) nh[k])*
                                       (log((double) nhctmp[c]/(double)  nh[k])/log(2));
                  }
                }
                lgcatN = lgcatN + -(alpha)*lgcatt;
      	  	  }

              // include the categorical covariate in the kth cluster
              if(Mcatp[pp*(*ncat)+p] == 0){
      	  	    nhctmp[Xcatp[pp*(*ncat)+p]] = nhctmp[Xcatp[pp*(*ncat)+p]] + 1;
              }

      	  	  if(*similarity_function==1){
      	  	    lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
      	  	    lgcatY = lgcatY + lgcatt;
      	  	  }
      	  	  if(*similarity_function==2){
      	  	  	lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
      	  	  	lgcatY = lgcatY + lgcatt;
      	  	  }
      	  	  if(*similarity_function==3){// Use entropy
      	  	  	lgcatt = 0.0;
                for(c=0;c<Cvec[p];c++){
                  if(nhctmp[c]==0){
                    lgcatt = lgcatt + 0;
                  }else{
                  	lgcatt = lgcatt + -((double) nhctmp[c]/(double) nh[k]+1)*
                  	                   (log((double) nhctmp[c]/(double) nh[k]+1)/log(2));
                  }
                }
                lgcatY = lgcatY + -(alpha)*lgcatt;
      	  	  }

      	  	} // This ends the loop through ncat categorical covariates

      	  	// These are for calibration 1
      	    gtilY[k] = lgconY + lgcatY;
      	    gtilN[k] = lgconN + lgcatN;

           	//////////////////////////////////////////////////////////
      	   	// Gower Compute similarity values for gower dissimilarity
      	   	//////////////////////////////////////////////////////////
            if(*similarity_function==4){ //.
              npd=0.0;
              lgconY = 0.0;
              lgconN = 0.0;
              for(j = 0; j < *nobs; j++){
                if(_Si[j] == k+1){

                  lgconY = lgconY + dissimtt[pp*(*nobs) + j];

                  for(jj = 0; jj < j; jj++){
                    if(_Si[jj] == k+1){
                      lgconN = lgconN + dissimtn[j*(*nobs) + jj];
                      lgconY = lgconY + dissimtn[j*(*nobs) + jj];
                    }
                  }
                }
              }

              npdN = nh[k]*(nh[k]-1)/2;
              if(npdN==0)npdN=1;
              npdY = (nh[k]+1)*(nh[k])/2;

              // This is cluster-mean Gower dissimilarity, but the next is used
              lgconN = -(alpha)*lgconN/(npdN);
              lgconY = -(alpha)*lgconY/(npdY);

              // Just Use the cluster-total Gower dissimilarity
              lgconN = -(alpha)*lgconN;
              lgconY = -(alpha)*lgconY;
            }
            //////////////////////////////////////////////////////////
      	   	// End of Gower similarity values for gower dissimilarity
      	   	//////////////////////////////////////////////////////////

      	  } // THIS ENDS THE PPMX PART

      	  // Note that if PPMx = FALSE, then
      	  // lgcatY = lgcatN = lgconY = lgconN = 0;

      	  ph[k] = log((double) nh[k]) +
      	       	 	lgcatY - lgcatN +
      	  			lgconY - lgconN;

      	  if(*calibrate == 2){

      	  	ph[k] = log((double) nh[k]) +
      	  	            (1/((double)*ncon + (double)*ncat))*
      	  	            (lgcatY + lgconY - lgcatN - lgconN);
      	  }

          if(*cohesion==2) ph[k] =  ph[k] - log((double) nh[k]);
        } // This ends loop through existing clusters.



      	// Now evaluate probability of new subject being assigned to own cluster
        lgcondraw = 0.0;
        lgcatdraw = 0.0;
        if(!(*PPM)){
          for(p=0;p<*ncon;p++){
            if(Mconp[pp*(*ncon)+p] == 0){
              xcontmp = Xconp[pp*(*ncon)+p];
              if(*similarity_function==1){
                if(*consim==1){
                  lgcondraw = lgcondraw + gsimconNN(m0,v,s20,xcontmp,xcontmp*xcontmp, mnmle[p],1,0,0,1);
                }
                if(*consim==2){
              	  lgcondraw = lgcondraw + gsimconNNIG(m0, k0, nu0, s20,xcontmp,xcontmp*xcontmp, mnmle[p],s2mle[p],1,0,0,1);
                }
              }
              if(*similarity_function==2){
                if(*consim==1){
              	  lgcondraw = lgcondraw + gsimconNN(m0,v,s20,xcontmp,xcontmp*xcontmp, mnmle[p],1,1,0,1);
                }
                if(*consim==2){
              	  lgcondraw = lgcondraw + gsimconNNIG(m0, k0, nu0, s20,xcontmp,xcontmp*xcontmp,mnmle[p],s2mle[p],1, 1, 0,1);
                }
              }
              if(*similarity_function==3){
                lgcondraw = lgcondraw + gsimconEV(xcontmp,xcontmp*xcontmp,1,alpha,1);
              }
            }

          }

          for(p=0;p<(*ncat);p++){
            for(c=0;c<Cvec[p];c++){nhctmp[c] = 0;}

            if(Mcatp[pp*(*ncat)+p] == 0){
              nhctmp[Xcatp[pp*(*ncat)+p]] = 1;

          	  if(*similarity_function==1){
          	    lgcatdraw = lgcatdraw + gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
          	  }
          	  if(*similarity_function==2){
          	    lgcatdraw = lgcatdraw + gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
          	  }
              if(*similarity_function==3){
                lgcatdraw = lgcatdraw + -(alpha)*0;
              }
            }
          }
          gtilY[_nclus] = lgcatdraw + lgcondraw;
          gtilN[_nclus] = lgcatdraw + lgcondraw;

        }

      	// Note if PPMx = FALSE, then logcon0 = lgcat0 = 0;
      	ph[_nclus] = log((double) Mdp) + lgcondraw + lgcatdraw;

      	// Gower similarity
      	if(*similarity_function==4) ph[_nclus] = log((double) Mdp) + log(1);

      	// Calibration through coarsening
      	if(*calibrate==2){
      		ph[_nclus] = log(Mdp) +
      							(1/((double)*ncon + (double)*ncat))*(lgcondraw + lgcatdraw);
      	}

      	// The uniform cohesion
      	if(*cohesion==2) ph[_nclus] = ph[_nclus] - log((double) Mdp);


      	//////////////////////////////////////////////////////////////////////////
      	// This is the calibration used when the similarity is standardized by
      	//////////////////////////////////////////////////////////////////////////

        if((*calibrate==1) & (*PPM != 1)){

          maxgtilN = gtilN[0];
          maxgtilY = gtilY[0];
          for(k=1; k < _nclus+1; k++){
            if(maxgtilN < gtilN[k]) maxgtilN = gtilN[k];

            if(k < _nclus){
              if(maxgtilY < gtilY[k]) maxgtilY = gtilY[k];
            }
          }


          sgY=0.0;
          sgN=0.0;

          for(k=0; k<_nclus+1; k++){

            lgtilN[k] = gtilN[k] - maxgtilN;
            sgN = sgN + exp(lgtilN[k]);

            if(k < _nclus){
              lgtilY[k] = gtilY[k] - maxgtilY;
              sgY = sgY + exp(lgtilY[k]);
            }
          }


          for(k=0; k<_nclus; k++){
            lgtilNk = lgtilN[k] - log(sgN);
            lgtilYk = lgtilY[k] - log(sgY);

            ph[k] = log((double) nh[k]) + lgtilYk - lgtilNk; //This takes into account both cont and cat vars

            if(*cohesion==2){
              ph[k] = ph[k] - log((double) nh[k]);
            }
          }
          // calibration for a singleton
          ph[_nclus] =  log(Mdp) + lgtilN[_nclus] - log(sgN);

          if(*cohesion==2){// Note with a uniform cohesion, for a new cluster
                           // the value of log(c({_nclus}}) = log(1) = 0;
           	ph[_nclus] = ph[_nclus] - log(Mdp);
          }
        }
      	//////////////////////////////////////////////////////////////////////////
      	// End of calibration used when the similarity is standardized by
      	//////////////////////////////////////////////////////////////////////////


        maxph = ph[0];
        for(k = 1; k < _nclus+1; k++){
          if(ph[k] > maxph) maxph=ph[k];
        }

        denph = 0.0;
        for(k = 0; k < _nclus+1; k++){
          ph[k] = exp(ph[k] - maxph);
          denph = denph + ph[k];
        }

        for(k = 0; k < _nclus+1; k++){
          probh[k] = ph[k]/denph;
        }

        uu = runif(0.0,1.0);

        cprobh= 0.0;

        for(k = 0; k < _nclus+1; k++){
          cprobh = cprobh + probh[k];
          if (uu < cprobh){
            iaux = k+1;
            break;
          }
        }

        mudraw = rnorm(_mu0, sqrt(_sig20));
        sdraw = runif(smin, smax);

        if(iaux <= _nclus){
          mudraw = _muh[(iaux-1)];
          sdraw = sqrt(_sig2h[(iaux-1)]);
        }else{
          mudraw = rnorm(_mu0,sqrt(_sig20));
          sdraw = runif(smin, smax);
        }

        mn = mudraw;
        s2k = sdraw*sdraw;
        if(*meanModel==2){
          // create x'beta but only include observed x's when doing this
          xb = 0.0, vb = 0.0;
          for(b = 0; b < ncov; b++){
            if(b < *ncon){
              if(Mconp[pp*(*ncon)+b] == 0){// not missing
                xb = xb + (fullXmatp[pp*ncov+b] - mnmle[b])*_beta[b];
              }
              if(Mconp[pp*(*ncon)+b] == 1){// missing covariate
                vb = vb + s2mle[b]*_beta[b]*_beta[b];
              }
            }
            if(b >= *ncon){ // I don't think categorical covariates are possible here.
              if(Mcat[j*(*ncat) + (b-*ncon)] == 0){
                xb = xb + fullXmat[j*ncov+b]*_beta[b];
              }
            }
          }
          mn = mudraw + xb;
          s2k = s2k + vb;
        }

        _ppred[pp] = rnorm(mn, sqrt(s2k));
        _predclass[pp] = iaux;

        mn = 0.0;
        for(k = 0; k < _nclus; k++){
          if(*meanModel == 1) mn = mn +  _muh[k]*probh[k];
          if(*meanModel == 2) mn = mn +  (_muh[k] + xb)*probh[k];
          _predclass_prob[pp*(*nobs) + k] = probh[k];
        }

        if(*meanModel == 1) mn = mn + rnorm(_mu0,sqrt(_sig20))*probh[_nclus];
        if(*meanModel == 2) mn = mn + rnorm(_mu0+xb,sqrt(_sig20))*probh[_nclus];

        _rbpred[pp] = mn;

        _predclass_prob[pp*(*nobs) + _nclus] = probh[_nclus];
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////
    //
    // Store MCMC iterates
    //
    //////////////////////////////////////////////////////////////////////////////////////
    if((i >= (*burn)) & ((i) % *thin ==0)){

      mu0[ii] = _mu0;
      sig20[ii] = _sig20;
      nclus[ii] = _nclus;

      for(j=0; j<*nobs; j++){
        Si[ii + nout*j] = _Si[j];
        like[ii + nout*j] = _like[j];
        ispred[ii + nout*j] = _ispred[j];

        mu[ii + nout*j] = _muh[_Si[j]-1];
        sig2[ii + nout*j] = _sig2h[_Si[j]-1];
      }

      for(b = 0; b < ncov; b++){
        beta[ii + nout*b] = _beta[b];
      }

      for(pp = 0; pp < *npred; pp++){
        ppred[ii + nout*pp] = _ppred[pp];
        predclass[ii + nout*pp] = _predclass[pp];

        rbpred[ii + nout*pp] = _rbpred[pp];
      }
      for(pp = 0; pp < (*nobs)*(*npred); pp++){
        predclass_prob[ii + nout*pp] = _predclass_prob[pp];
        _predclass_prob[pp] = 0.0;
      }

      ii = ii + 1;
    }

  }


  //////////////////////////////////////////////////////////////////////////////////
  // calculate LPML
  //////////////////////////////////////////////////////////////////////////////////
  _lpml=0.0;
  for(j = 0; j < *nobs; j++){
    _lpml = _lpml + log(1/CPOinv[j]);
  }
  lpml[0] = _lpml;

  ////////////////////////////////////////////////////////////////////////////////////////////
  // Computing WAIC  (see Gelman article in lit review folder)
  ////////////////////////////////////////////////////////////////////////////////////////////
  elppdWAIC = 0.0;
  for(j = 0; j < *nobs; j++){
    elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
  }
  waic[0] = -2*elppdWAIC;




}



SEXP GAUSSIAN_PPMX_MISSING(SEXP y, SEXP nobs,
                           SEXP Xcon, SEXP Mcon, SEXP ncon,
                           SEXP Xcat, SEXP Mcat, SEXP ncat, SEXP Cvec,
                           SEXP npred,
                           SEXP Xconp, SEXP Mconp,
                           SEXP Xcatp, SEXP Mcatp,
                           SEXP M, SEXP meanModel, SEXP modelPriors, SEXP simParms,
                           SEXP PPM, SEXP cohesion, SEXP similarity_function, SEXP consim,
                           SEXP dissimtn, SEXP dissimtt, SEXP calibrate, SEXP mh,
                           SEXP verbose, SEXP draws, SEXP burn, SEXP thin){
  int nprot = 0;
  int _nobs = asInteger(nobs);
  int _ncon = asInteger(ncon);
  int _ncat = asInteger(ncat);
  int _npred = asInteger(npred);

  int _meanModel = asInteger(meanModel);
  int _PPM = asInteger(PPM);
  int _cohesion = asInteger(cohesion);
  int _similarity_function = asInteger(similarity_function);
  int _consim = asInteger(consim);
  int _calibrate = asInteger(calibrate);
  int _verbose = asInteger(verbose);

  int _niter = asInteger(draws);
  int _nburn = asInteger(burn);
  int _nthin = asInteger(thin);

  double _M = asReal(M);


  double nout = (_niter-_nburn)/_nthin;
  double ncov = _ncon + _ncat;

  y = PROTECT(coerceVector(y, REALSXP)); nprot++;
  Xcon =  PROTECT(coerceVector(Xcon, REALSXP)); nprot++;
  Mcon =  PROTECT(coerceVector(Mcon, INTSXP)); nprot++;
  Xcat =  PROTECT(coerceVector(Xcat, INTSXP)); nprot++;
  Mcat =  PROTECT(coerceVector(Mcat, INTSXP)); nprot++;
  Xconp =  PROTECT(coerceVector(Xconp, REALSXP)); nprot++;
  Mconp =  PROTECT(coerceVector(Mconp, INTSXP)); nprot++;
  Xcatp =  PROTECT(coerceVector(Xcatp, INTSXP)); nprot++;
  Mcatp =  PROTECT(coerceVector(Mcatp, INTSXP)); nprot++;
  Cvec =  PROTECT(coerceVector(Cvec, INTSXP)); nprot++;
  modelPriors = PROTECT(coerceVector(modelPriors, REALSXP)); nprot++;
  simParms = PROTECT(coerceVector(simParms, REALSXP)); nprot++;
  mh = PROTECT(coerceVector(mh, REALSXP)); nprot++;

  dissimtn  = PROTECT(coerceVector(dissimtn, REALSXP)); nprot++;
  dissimtt  = PROTECT(coerceVector(dissimtt, REALSXP)); nprot++;

  SEXP Si = PROTECT(allocMatrix(INTSXP, nout, _nobs)); nprot++;
  SEXP MU = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++;
  SEXP SIG2 = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++;
  SEXP BETA = PROTECT(allocMatrix(REALSXP, nout, ncov)); nprot++;
  SEXP NCLUS = PROTECT(allocMatrix(INTSXP, nout, 1)); nprot++;
  SEXP MU0 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP SIG20 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;

  SEXP ISPRED = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++;
  SEXP PPRED = PROTECT(allocMatrix(REALSXP, nout, _npred)); nprot++;
  SEXP RBPRED = PROTECT(allocMatrix(REALSXP, nout, _npred)); nprot++;
  SEXP PREDCLASS = PROTECT(allocMatrix(INTSXP, nout, _npred)); nprot++;
  SEXP PREDCLASS_PROB = PROTECT(allocMatrix(REALSXP, nout, (_nobs)*(_npred))); nprot++;

  SEXP LIKE = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++;
  SEXP WAIC = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
  SEXP LPML = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;

  double *MUout, *SIG2out, *MU0out, *SIG20out, *BETAout, *LIKEout, *WAICout, *LPMLout;
  double *ISPREDout, *PPREDout, *RBPREDout, *PREDCLASS_PROBout;
  int *Siout, *NCLUSout, *PREDCLASSout;

  MUout = REAL(MU);
  SIG2out = REAL(SIG2);
  MU0out = REAL(MU0);
  SIG20out = REAL(SIG20);
  BETAout = REAL(BETA);
  Siout = INTEGER(Si);
  NCLUSout = INTEGER(NCLUS);

  LIKEout = REAL(LIKE);
  WAICout = REAL(WAIC);
  LPMLout = REAL(LPML);

  ISPREDout = REAL(ISPRED);
  PPREDout = REAL(PPRED);
  RBPREDout = REAL(RBPRED);
  PREDCLASS_PROBout = REAL(PREDCLASS_PROB);
  PREDCLASSout = INTEGER(PREDCLASS);



  GetRNGstate();

  gaussian_ppmx_missing(REAL(y), &_nobs,
                REAL(Xcon), INTEGER(Mcon), &_ncon,
                INTEGER(Xcat), INTEGER(Mcat), &_ncat, INTEGER(Cvec),
                &_npred,
                REAL(Xconp), INTEGER(Mconp),
                INTEGER(Xcatp), INTEGER(Mcatp),
                &_M, &_meanModel, REAL(modelPriors), REAL(simParms),
                &_PPM, &_cohesion, &_similarity_function, &_consim,
                REAL(dissimtn), REAL(dissimtt), &_calibrate, REAL(mh),
                &_verbose, &_niter, &_nburn, &_nthin,
                Siout, NCLUSout, MUout, SIG2out, BETAout,
                MU0out, SIG20out,
                LIKEout, WAICout, LPMLout,
                ISPREDout, PPREDout, PREDCLASSout,
                RBPREDout, PREDCLASS_PROBout);

  PutRNGstate();


  SEXP ans = PROTECT(allocVector(VECSXP, 15)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIG2);
  SET_VECTOR_ELT(ans, 2, BETA);
  SET_VECTOR_ELT(ans, 3, Si);
  SET_VECTOR_ELT(ans, 4, MU0);
  SET_VECTOR_ELT(ans, 5, SIG20);
  SET_VECTOR_ELT(ans, 6, NCLUS);
  SET_VECTOR_ELT(ans, 7, LIKE);
  SET_VECTOR_ELT(ans, 8, WAIC);
  SET_VECTOR_ELT(ans, 9, LPML);
  SET_VECTOR_ELT(ans, 10, ISPRED);
  SET_VECTOR_ELT(ans, 11, PPRED);
  SET_VECTOR_ELT(ans, 12, PREDCLASS);
  SET_VECTOR_ELT(ans, 13, RBPRED);
  SET_VECTOR_ELT(ans, 14, PREDCLASS_PROB);


  SEXP nm = allocVector(STRSXP, 15);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("sig2"));
  SET_STRING_ELT(nm, 2, mkChar("beta"));
  SET_STRING_ELT(nm, 3, mkChar("Si"));
  SET_STRING_ELT(nm, 4, mkChar("mu0"));
  SET_STRING_ELT(nm, 5, mkChar("sig20"));
  SET_STRING_ELT(nm, 6, mkChar("nclus"));
  SET_STRING_ELT(nm, 7, mkChar("like"));
  SET_STRING_ELT(nm, 8, mkChar("WAIC"));
  SET_STRING_ELT(nm, 9, mkChar("lpml"));
  SET_STRING_ELT(nm, 10, mkChar("fitted.values"));
  SET_STRING_ELT(nm, 11, mkChar("ppred"));
  SET_STRING_ELT(nm, 12, mkChar("predclass"));
  SET_STRING_ELT(nm, 13, mkChar("rbpred"));
  SET_STRING_ELT(nm, 14, mkChar("predclass_prob"));

  UNPROTECT(nprot);
  return(ans);
}



