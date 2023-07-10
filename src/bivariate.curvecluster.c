/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit a bivariate hierarchical q-degree (b)-spline model with
 * a fixed number of knots penalized using P-splines.
 *
 * I need to put a brief description of the model here.
 *
 * begin with a fairly flat Gaussian prior on
 *************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*****************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned

* nsubject = number of subjects in data set
* nobs = vector whose entries indicate number of observations per subject
* y1 = vector containing nobs saggital values for each subject
* y2 = vector containing nobs frontal values for each subject
* z = vector containing nobs time (time is incremented by one) values for each player

* K = smoothing matrix for penalized B-splines.
* dimK = dimension of K.

* ncon = number of continuous covariates
* ncat = number of categorical covariates
* Cvec = ncat x 1 vector indicating the number of categories for each categorical covariate.

* Xcon = nsubject x ncon contiguous double vector that contains continuous covariate values
* Xcat = nsubject x ncat contiguous int vector that contains categorical covariate values


* PPM = logical indicating if PPM or PPMx should be used.
* M = double indicating value of M associated with cohesion (scale parameter of DP).

* simularity_function = int indicating similarity function for both 
                        continuous and categorical variable
		1 - auxiliary model
		2 - double dipper
* consim = 1 or 2.  1 implies sim for con var is N-N.  2 implies sim is N-NIG

* npred = integer indicating number of out of sample predictions
* npredobs = integer indicating number of "time" points for which to predict
* Xconp = npred x ncon matrix of continuous covariates to make predictions
* Xcatp = npred x ncat matrix of categorical covariates to make predictions
*
* simparms = vector containing parameter values for similarity functions
* alpha = parameter value associated with the cluster variance similarity
* calibration = integer that determines if the calibrated similarity should be used
*	0 - do not use the calibrated similarity
*	1 - use the calibrated similarity
*****************************************************************************************/


void mcmc_bivariate_curvecluster(int *draws, int *burn, int *thin, //3
                                 int *nsubject, int *nobs, //2
                                 double *y1, double *y2, double *z, //3
                                 double *K, int *dimK, //2
                                 int *ncon, int *ncat, int *Cvec, //3
                                 double *Xcon, int *Xcat, //2
                                 int *PPM, double *M,//2
                                 int *similarity_function, int *consim, //2
                                 int *npred, int *npredobs, double *Xconp, int *Xcatp, //4
                                 double *simParms, double *Aparm, //2
                                 int *calibrate, double *modelPriors, double *Hmat, //3
                                 int *balanced, double *mh, //2
                                 double *beta1, double *beta2, double *beta01, double *beta02, //4
                                 double *sig21, double *sig22, double *mub01, double *mub02, //4
                                 double *sig2b01,double *sig2b02,double *lam1, double *lam2, //4
                                 double *tau21, double *tau22, double *theta1, double *theta2, //4
                                 double *mu1,double *mu2, int *Si, int *nclus, //4
                                 double *ppred1,  double *ppred2, int *predclass, //3
                                 double *llike, double *lpml,double *WAIC){ //3


  // i - MCMC iterate
  // ii - MCMC save iterate
  // j - player iterate
  // jj - second player iterate
  // t - time (cumulative minutes played) per player iterate
  // c - category iteration
  // b - beta and thetah iterate
  // bb - second beta and thetah iterate
  // k - cluster iterate
  // kk - knot iterate
  // p - covariate iterate (for both continuous and categorical);

  int c, i, ii, j, t, p, jj, k, kk, b, bb;
  int pp;
  int csobs = 0, csHobs=0;
  int nb= (*dimK), N ;
  int nout = (*draws - *burn)/(*thin);

  Rprintf("nsubject = %d\n", *nsubject);
  Rprintf("nb = %d\n", nb);
  Rprintf("ncon = %d\n", *ncon);
  Rprintf("ncat = %d\n", *ncat);


  // RprintIVecAsMat("nobs", nobs, 1, *nsubject);
  if(*ncat > 0) RprintIVecAsMat("Cvec", Cvec, 1, *ncat);


  double *zpred = R_Vector(*npredobs);
  for(pp = 0; pp < *npredobs; pp++) zpred[pp] = (double) (pp+1);
  //	RprintVecAsMat("zpred", zpred, 1, *npredobs);
  //	nstrtkntsindx = ((int) floor((double) *nstartknots));

  //	Rprintf("nstrtkntsindx = %d\n", nstrtkntsindx);
  int  max_nobs;


  max_nobs = nobs[0];
  N=0;
  for(j = 0; j < (*nsubject); j++){
    //		Rprintf("j = %d\n", j);
    if(max_nobs < nobs[j]) max_nobs = nobs[j];
    N = N + nobs[j];

  }

  N = csobs;

  Rprintf("max_nobs = %d\n", max_nobs);
  Rprintf("N = %d\n", N);

  int max_C;
  max_C = Cvec[0];
  for(p = 0; p < (*ncat); p++){
    if(max_C < Cvec[p]) max_C = Cvec[p];
  }

  if(*ncat == 0) max_C = 1.0;

  // Rprintf("max_C = %d\n", max_C);



  // ===================================================================================
  //
  // Memory vectors to hold MCMC iterates for non cluster specific parameters
  //
  // ===================================================================================

  double *sig21_iter = R_VectorInit(*nsubject, ((modelPriors[0]-0.0)/2.0)*((modelPriors[0]-0.0)/2.0));
  double *sig22_iter = R_VectorInit(*nsubject, ((modelPriors[0]-0.0)/2.0)*((modelPriors[0]-0.0)/2.0));

  double *beta01_iter = R_VectorInit(*nsubject, 0.0);
  double *beta02_iter = R_VectorInit(*nsubject, 0.0);

  double *beta1_iter = R_VectorInit(nb*(*nsubject),0.0);
  double *beta2_iter = R_VectorInit(nb*(*nsubject),0.0);


  double *mu1_iter = R_VectorInit(nb,0.0);
  double *mu2_iter = R_VectorInit(nb,0.0);

  double mub01_iter = 0.0;
  double mub02_iter = 0.0;
  double sig2b01_iter = 1.0;
  double sig2b02_iter = 1.0;


  int Si_iter[*nsubject];
  int nclus_iter = 0;


  // ===================================================================================
  //
  // Memory vectors to hold MCMC iterates for cluster specific parameters
  //
  // ===================================================================================

  double *tau2h1 = R_VectorInit(*nsubject, 1.0);
  double *lamh1 = R_VectorInit(*nsubject, 0.5);
  double *thetah1 = R_VectorInit(nb*(*nsubject), 0.0);

  double *tau2h2 = R_VectorInit(*nsubject, 1.0);
  double *lamh2 = R_VectorInit(*nsubject, 0.5);
  double *thetah2 = R_VectorInit(nb*(*nsubject), 0.0);


  //	tau2h[0] = 0.01; tau2h[1] = 1.0; tau2h[2] = 100.0;

  int nh[*nsubject];


  // ===================================================================================
  //
  // Initialize a few parameter vectors
  //
  // ===================================================================================

  // Initialize Si according to covariates
  for(j = 0; j < *nsubject; j++){
    //		Si_iter[j] = x[j]*(*ncat2) + x2[j]+1;
    //Si_iter[j] = 1;
    Si_iter[j] = rbinom(3,0.4)+1;
    nh[j] = 0;
  }
  // Initial enumeration of number of players per cluster;
  for(j = 0; j < *nsubject; j++){
    for(k = 0; k < *nsubject; k++){
      if(Si_iter[j] == k+1) nh[k] = nh[k] + 1;
    }
  }
  // Initialize the number of clusters
  for(j = 0; j < *nsubject; j++){
    if(nh[j] > 0) nclus_iter = nclus_iter + 1;
  }



  // ===================================================================================
  //
  // scratch vectors of memory needed to update parameters
  //
  // ===================================================================================
  int  big = max_nobs;
  if(nb > max_nobs) big = nb;
  Rprintf("big = %d\n", big);
  // These are made particularly big to make sure there is enough memory
  double *scr1 = R_Vector((big)*(big));
  double *scr2 = R_Vector((big)*(big));
  double *scr3 = R_Vector((big)*(big));
  double *scr4 = R_Vector((big)*(big));


  // stuff that I need to update Si (cluster labels);
  int iaux=1, auxint, nhctmp[max_C];

  double auxreal, uu, sumxtmp, sumx2tmp, xcontmp;
  double sumsq2, sumsq1, maxph, denph, cprobh;
  double lamdraw1, tau2draw1, lamdraw2, tau2draw2;

  double *thetadraw1 = R_VectorInit(nb, 0.0);
  double *thetadraw2 = R_VectorInit(nb, 0.0);

  double *ph = R_VectorInit(*nsubject, 0.0);
  double *probh = R_VectorInit(*nsubject, 0.0);

  double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw;
  double lgcont,lgcatt;
  //	double mgconN, mgconY, mgcatN, mgcatY, sgconN, sgconY, sgcatY, sgcatN;

  double *gtilN = R_VectorInit((*nsubject+1),0.0);
  double *gtilY = R_VectorInit((*nsubject+1),0.0);
  double sgY, sgN, lgtilN, lgtilY, maxgtilY, maxgtilN;

  double *sumx = R_VectorInit((*nsubject)*(*ncon),0.0);
  double *sumx2 = R_VectorInit((*nsubject)*(*ncon),0.0);
  int  nhc[(*nsubject)*(*ncat)*max_C];


  // stuff I need to update sig2 (subject specific), mub0, sig2b0;
  double astar, bstar1, bstar2, os, ns, sumsq;
  double sumb01,sumb02;

  // stuff that I need for player specific beta's;
  //	double *H = R_VectorInit((*nsubject)*(nb)*(max_nobs),0.0);
  double *H = R_VectorInit((*nsubject)*(nb)*(big),0.0);
  double *tH = R_VectorInit((*nsubject)*(nb)*(big),0.0);
  double *HtH = R_Vector(nb*nb);
  double *Hty1 = R_Vector(nb);
  double *Hty2 = R_Vector(nb);

  double *y1_tmp = R_Vector(big);
  double *y2_tmp = R_Vector(big);

  double sumy1_Hb, s2star1, mstar1;
  double sumy2_Hb, s2star2, mstar2;
  double *y1_b0 = R_Vector(big);
  double *y2_b0 = R_Vector(big);
  double *Hb1 = R_Vector(big);
  double *Hb2 = R_Vector(big);


  // Create the inverse of the K penalty matrix
  double ld, ld1, ld2;
  double *Kinv = R_Vector(nb*nb);
  for(b = 0; b < nb; b++){for(bb = 0; bb < nb; bb++){Kinv[b*nb+bb] = K[b*nb+bb];}}
  cholesky(Kinv, nb, &ld);
  inverse_from_cholesky(Kinv, scr1, scr2, nb);

  // stuff I need for tau2h;
  double *thtmp1 = R_Vector(nb);
  double *thtmp2 = R_Vector(nb);

  // stuff I need for thetah
  double *sumbeta1 = R_Vector(nb);
  double *sumbeta2 = R_Vector(nb);
  double *Mstar1 = R_Vector(nb);
  double *Mstar2 = R_Vector(nb);
  double *Sstar1 = R_Vector(nb*nb);
  double *Sstar2 = R_Vector(nb*nb);
  double *outrmvnorm1 = R_Vector(nb);
  double *outrmvnorm2 = R_Vector(nb);

  //stuff I need for lamh
  double olam1, nlam1, lln, llo, llr, lln1, lln2, llo1, llo2, llr1, llr2, ldo1, ldn1;
  double olam2, nlam2, ldo2, ldn2;
  double *btmp1 = R_Vector(nb);
  double *btmp2 = R_Vector(nb);
  double *nV1 = R_Vector(nb*nb);
  double *nV2 = R_Vector(nb*nb);
  double *oV1 = R_Vector(nb*nb);
  double *oV2 = R_Vector(nb*nb);


  // stuff I need for mu
  double *sumtheta1 = R_Vector(nb);
  double *sumtheta2 = R_Vector(nb);
  double sumtau21;
  double sumtau22;

  // stuff that I need to perform the predictions
//  double *thetatmp1 = R_VectorInit(nb, 0.0);
//  double *thetatmp2 = R_VectorInit(nb, 0.0);
//  double *bpred1 = R_VectorInit(nb, 0.0);
//  double *bpred2 = R_VectorInit(nb, 0.0);
//  double *ppredtmp1 = R_VectorInit(100,0.0);
//  double *ppredtmp2 = R_VectorInit(100,0.0);
//  double lgcon0, lgcat0=0.0;


  // stuff that I need to compute the lpml;
  double llikeval, lpml_iter = 0.0, elppdWAIC;;
  double *CPO = R_VectorInit((*nsubject), 0.0);
  double *llike_iter = R_VectorInit(*nsubject,0.0);
  double *mnlike = R_VectorInit((*nsubject), 0.0);
  double *mnllike = R_VectorInit((*nsubject), 0.0);
  int like0,nout_0=0;;





  if(*balanced==1){
    for(t=0; t < nobs[0]; t++){
      for(kk=0; kk<nb; kk++){
        H[t*(nb) + kk] = Hmat[t*(nb) + kk];
      }
    }
    
    mat_transpose(H, tH, nobs[0], nb);
    matrix_product(tH, H, HtH, nb, nb, nobs[0]);

    // RprintVecAsMat("H", H, nobs[0], nb);
  }



  double *mnmle = R_Vector(*ncon);
  double *s2mle = R_Vector(*ncon);
  double sum, sum2;
  for(p = 0; p < *ncon; p++){
    sum = 0.0, sum2=0.0;
    for(j = 0; j < *nsubject; j ++){
      sum = sum + Xcon[j*(*ncon) + p];
      sum2 = sum2 + Xcon[j*(*ncon) + p]*Xcon[j*(*ncon) + p];
    }

    mnmle[p] = sum/((double) *nsubject);
    s2mle[p] = sum2/((double) *nsubject) - mnmle[p]*mnmle[p];
  }


  // ===================================================================================
  //
  // Prior parameter values
  //
  // ===================================================================================
  // upper bound on sigma21 and sigma22
  double Asig = modelPriors[0]; // this is for Metropolis step;

  Rprintf("Asig = %f\n", Asig);


  // priors for mu
  double s2mu = modelPriors[1];


  // priors for mub0 and sig2b0;
  double mb0 = modelPriors[2], s2b0 = modelPriors[3];
  double ab0 = modelPriors[4], bb0 = modelPriors[5];

  // IG parameters for tau2 Less smoothing (more local behavior) make denominator of
  // bt smaller (as this makes tau2 smaller)
  //	double at = 0.0; double bt = 1.0;
  double at = modelPriors[6]; double bt = modelPriors[7];

  Rprintf("at = %f\n", at);
  Rprintf("bt = %f\n", bt);

  // Uniform prior for lam
  double A1 = Aparm[0];
  double A2 = Aparm[1];
  Rprintf("A1 = %f\n", A1);
  Rprintf("A2 = %f\n", A2);

  // DP weight parameter
  double Mdp = *M;

  // dirichlet denominator parameter
  double *dirweights = R_VectorInit(max_C, simParms[5]);
  //	double m0=0.0, s20=0.5, v=1.0, k0=1.0, nu0=1.0;
  double m0=simParms[0];
  double s20=simParms[1];
  double v2=simParms[2];
  double k0=simParms[3];
  double nu0=simParms[4];
  double alpha=simParms[6]; // parameter for gower and variance similarity


  Rprintf("Mdp = %f\n", Mdp);
  RprintVecAsMat("dirweights", dirweights, 1, max_C);




  // ===================================================================================
  //
  // Initialize the cluster-specific sufficient statistics for continuous covariates
  // and categorical covariates.
  //
  // ===================================================================================
  if(!(*PPM)){
    for(j=0;j<*nsubject;j++){
      mnlike[j] = 0.0;
      mnllike[j] = 0.0;
      for(p=0;p<*ncon;p++){
        sumx[j*(*ncon) + p] = 0.0;
        sumx2[j*(*ncon) + p] = 0.0;
      }
      for(p=0;p<*ncat;p++){
        for(c=0; c<max_C; c++){
          nhc[(j*(*ncat) + p)*max_C + c] = 0;
        }
      }
    }

    // Fill in cluster-specific sufficient statistics based on first partition
    for(j = 0; j < *nsubject; j++){
      for(p=0; p<*ncon; p++){
        sumx[(Si_iter[j]-1)*(*ncon) + p] = sumx[(Si_iter[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p];
        sumx2[(Si_iter[j]-1)*(*ncon) + p] = sumx2[(Si_iter[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
      }
      for(p=0; p<*ncat; p++){
        nhc[((Si_iter[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
          nhc[((Si_iter[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] + 1;
      }
    }
  }

//  RprintIVecAsMat("nhc", nhc, *ncat, max_C);
//  RprintIVecAsMat("Xcat", Xcat, *nobs, *ncat);

  // Candidate density sd's for M-H stuff
  double csigLAM1 = mh[0], csigLAM2 = mh[1], csigSIG=mh[2];

  //	RprintVecAsMat("H", H, nobs[0], nb);

  ii = 0;

  GetRNGstate();


  // ===================================================================================
  //
  // start of the mcmc algorithm;
  //
  // ===================================================================================

  for(i = 0; i < *draws; i++){


    if((i+1) % 10000 == 0){
      time_t now;
      time(&now);

      Rprintf("mcmc iter = %d =========================================== \n", i+1);
      Rprintf("%s", ctime(&now));
      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
      //			RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);
    }



    //////////////////////////////////////////////////////////////////////////////////
    //
    // update the cluster labels using the polya urn scheme of
    // algorithm 8 found in  Radford Neal's
    //	"Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
    //	paper.
    //
    //////////////////////////////////////////////////////////////////////////////////

    //		RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
    //		Rprintf("nclus_iter = %d\n", nclus_iter);
    for(j = 0; j < *nsubject; j++){
   
//      Rprintf("j = %d\n", j);

      if(nh[Si_iter[j]-1] > 1){

        // Observation belongs to a non-singleton ...
        nh[Si_iter[j]-1] = nh[Si_iter[j]-1] - 1;

        if(!(*PPM)){
          // need to reduce the sumx sumx2 to
          for(p = 0; p < *ncon; p++){
            sumx[(Si_iter[j]-1)*(*ncon) + p] = sumx[(Si_iter[j]-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
            sumx2[(Si_iter[j]-1)*(*ncon) + p] = sumx2[(Si_iter[j]-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          }

          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){

            nhc[((Si_iter[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
               nhc[((Si_iter[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;

          }
        }



      }else{

        // Observation is a member of a singleton cluster ...

        iaux = Si_iter[j];
        //				Rprintf("iaux = %d\n", iaux);
        if(iaux < nclus_iter){

          // Need to relabel clusters. I will do this by swapping cluster labels
          // Si_iter[j] and nclus_iter along with cluster specific parameters;


          // All members of last cluster will be assigned subject j's label
          for(jj = 0; jj < *nsubject; jj++){

            if(Si_iter[jj] == nclus_iter){
              Si_iter[jj] = iaux;
            }
          }


          Si_iter[j] = nclus_iter;

          // The following steps swaps order of cluster specific parameters
          // so that the newly labeled subjects from previous step retain
          // their correct cluster specific parameters
          auxreal = tau2h1[iaux-1];
          tau2h1[iaux-1] = tau2h1[nclus_iter-1];
          tau2h1[nclus_iter-1] = auxreal;

          auxreal = tau2h2[iaux-1];
          tau2h2[iaux-1] = tau2h2[nclus_iter-1];
          tau2h2[nclus_iter-1] = auxreal;

          auxreal = lamh1[iaux-1];
          lamh1[iaux-1] = lamh1[nclus_iter-1];
          lamh1[nclus_iter-1] = auxreal;

          auxreal = lamh2[iaux-1];
          lamh2[iaux-1] = lamh2[nclus_iter-1];
          lamh2[nclus_iter-1] = auxreal;

          for(b = 0; b < nb; b++){

            auxreal = thetah1[b*(nclus_iter) + iaux-1];
            thetah1[b*(nclus_iter) + iaux-1] = thetah1[b*(nclus_iter) + nclus_iter-1];
            thetah1[b*(nclus_iter) + nclus_iter-1] = auxreal;

            auxreal = thetah2[b*(nclus_iter) + iaux-1];
            thetah2[b*(nclus_iter) + iaux-1] = thetah2[b*(nclus_iter) + nclus_iter-1];
            thetah2[b*(nclus_iter) + nclus_iter-1] = auxreal;
          }


          // the number of members in cluster is also swapped with the last
          nh[iaux-1] = nh[nclus_iter-1];
          nh[nclus_iter-1] = 1;
          
          if(!(*PPM)){
            // need to swap sumx and sumx2
            for(p = 0; p < *ncon; p++){
              auxreal = sumx[(iaux-1)*(*ncon) + p];
              sumx[(iaux-1)*(*ncon) + p] = sumx[(nclus_iter-1)*(*ncon) + p];
              sumx[(nclus_iter-1)*(*ncon) + p] = auxreal;

              auxreal = sumx2[(iaux-1)*(*ncon) + p];
              sumx2[(iaux-1)*(*ncon) + p] = sumx2[(nclus_iter-1)*(*ncon) + p];
              sumx2[(nclus_iter-1)*(*ncon) + p] = auxreal;

            }

            // need to swap nhc as well
            for(p = 0; p < *ncat; p++){
              for(c=0; c<max_C; c++){
                auxint = nhc[((iaux-1)*(*ncat) + p)*(max_C) + c];
                nhc[((iaux-1)*(*ncat) + p)*(max_C) + c] = nhc[((nclus_iter-1)*(*ncat) + p)*(max_C) + c];
                nhc[((nclus_iter-1)*(*ncat) + p)*(max_C) + c] = auxint;
              }
            }
          }

        }


        // Now remove the ith obs and last cluster;
        nh[nclus_iter-1] = nh[nclus_iter-1] - 1;


        // need to reduce the sumx sumx2
        if(!(*PPM)){
          for(p = 0; p < *ncon; p++){
            sumx[(nclus_iter-1)*(*ncon) + p] = sumx[(nclus_iter-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
            sumx2[(nclus_iter-1)*(*ncon) + p] = sumx2[(nclus_iter-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          }

          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){

             nhc[((nclus_iter-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
                nhc[((nclus_iter-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;

           }
        }
        // Finally reduce the number of clusters
        nclus_iter = nclus_iter - 1;
      }



      // The atoms have been relabeled if necessary and now we need to
      // update Si.

      for(b = 0; b < nb; b++){

        btmp1[b] = beta1_iter[b*(*nsubject) + j];
        btmp2[b] = beta2_iter[b*(*nsubject) + j];

      }


      ///////////////////////////////////
      //
      // Begin the cluster probabilities
      //
      //////////////////////////////////
      for(k = 0 ; k < nclus_iter; k++){

//        Rprintf("k = %d ==================== \n", k);

        for(b = 0; b < nb; b++){

          thtmp1[b] = thetah1[b*(nclus_iter) + k];
          thtmp2[b] = thetah2[b*(nclus_iter) + k];

          for(bb = 0; bb < nb; bb++){

            oV1[b*nb+bb] = 0.0;
            oV2[b*nb+bb] = 0.0;

            if(b == bb){

              oV1[b*nb+bb] = 1/(lamh1[k]*lamh1[k]);
              oV2[b*nb+bb] = 1/(lamh2[k]*lamh2[k]);
              //							oV[b*nb+bb] = 1/(lamh[k]);
            }
          }
        }


        ldo1 = 2.0*nb*log(lamh1[k]);
        ldo2 = 2.0*nb*log(lamh2[k]);


        lgconY = 0.0;
        lgconN = 0.0;
        lgcatY=0.0;
        lgcatN=0.0;
        
        if(!(*PPM)){
        
          for(p=0; p<(*ncon); p++){
  
            sumxtmp = sumx[k*(*ncon) + p];
            sumx2tmp = sumx2[k*(*ncon) + p];
  
            if(*similarity_function==1){ // Auxilliary
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 0, 0, 1);
                lgconN = lgconN + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k], 0, 0, 1);
                lgconN = lgconN + lgcont;
              }
  
            }
  
            if(*similarity_function==2){ //Double Dipper
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 1, 0, 1);
                lgconN = lgconN + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k], 1, 0, 1);
                lgconN = lgconN + lgcont;
              }
            }
  
            if(*similarity_function==3){ // Cluster Variance
              lgcont = gsimconEV(sumxtmp, sumx2tmp, nh[k], alpha,1);
              lgconN = lgconN + lgcont;
            }
  
            // now add jth individual back;
      	    sumxtmp = sumxtmp + Xcon[j*(*ncon)+p];
      	    sumx2tmp = sumx2tmp + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
  
            if(*similarity_function==1){ // Auxilliary
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k]+1, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k]+1, 0, 0, 1);
                lgconY = lgconY + lgcont;
              }
            }
            if(*similarity_function==2){ //Double Dipper
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k]+1, 1, 0, 1);
                lgconY = lgconY + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k]+1, 1, 0, 1);
                lgconY = lgconY + lgcont;
              }
            }
            if(*similarity_function==3){ // variance
              //						Rprintf("gsimconEV = %f\n", gsimconEV(sumx, sumx2, nhtmp, alpha, 0));
              lgcont = gsimconEV(sumxtmp, sumx2tmp, nh[k]+1, alpha,1);
              lgconY = lgconY + lgcont;
            }
  
          }
  
  
     	  // Now calculate similarity for the categorical covariates  
          for(p=0; p<(*ncat); p++){
  
             for(c=0;c<Cvec[p];c++){
              nhctmp[c] = nhc[(k*(*ncat) + p)*(max_C) + c];
            }
  
//            RprintIVecAsMat("nhctmp", nhctmp, 1, Cvec[p]);
//            Rprintf("Cvec[p] = %d\n", Cvec[p]);
            
            if(*similarity_function==1){
              lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 0, 1);
              lgcatN = lgcatN + lgcatt;
            }
            if(*similarity_function==2){
              lgcatt = gsimcatDM(nhctmp, dirweights, Cvec[p], 1, 1);
              lgcatN = lgcatN + lgcatt;
            }
            if(*similarity_function==3){// Using Entropy instead of variance here
              lgcatt = 0.0;
              for(c=0;c<Cvec[p];c++){
                if(nhctmp[c]==0){
                  lgcatt = lgcatt + 0;
                }else{
                  lgcatt = lgcatt + -((double) nhc[c]/(double) nh[k])*
                                      (log((double) nhc[c]/(double) nh[k])/log(2));
                }
  
              }
              lgcatN = lgcatN + -(alpha)*lgcatt;
            }
//            Rprintf("lgactN = %f\n", lgcatN);
  
            // include the categorical covariate in the kth cluster
//            Rprintf("Xcat[j*(*ncat)+p] = %d\n", Xcat[j*(*ncat)+p]);
      	  	nhctmp[Xcat[j*(*ncat)+p]] = nhctmp[Xcat[j*(*ncat)+p]] + 1;
//            RprintIVecAsMat("nhctmp", nhctmp, 1, Cvec[p]);
  
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
                  lgcatt = lgcatt + -((double) nhc[c]/(double) nh[k]+1)*
                                      (log((double) nhc[c]/(double) nh[k]+1)/log(2));
                }
              }
              lgcatY = lgcatY + -(alpha)*lgcatt;
            }
  
//            Rprintf("lgcatY = %f\n", lgcatY);
          }
  
          gtilY[k] = lgconY + lgcatY;
          gtilN[k] = lgconN + lgcatN;
        } // THIS ENDS THE PPMX PART.  


      	// Compute the unnormalized cluster probabilities
      	// Note that if PPMx = FALSE then
      	// lgcatY = lgcatN = lgconY = lgconN = 0;
//      	Rprintf("dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1) = %f\n", dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1));
//      	Rprintf("dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1) = %f\n", dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1));
//      	Rprintf("log((double) nh[k]) = %f\n", log((double) nh[k]));
//      	Rprintf("lgcatY - lgcatN = %f\n", lgcatY - lgcatN);
//      	Rprintf("lgconY - lgconN = %f\n", lgconY - lgconN);
        ph[k] = dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1) +
                dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1) +
                   log((double) nh[k]) + // Scale parameter from DP
                   lgcatY - lgcatN +  // Ratio of categorical variables
                   lgconY - lgconN;   // Ratio of continuous variabels


      }


      lamdraw1 = runif(0,A1);
      lamdraw2 = runif(0,A2);
      tau2draw1 = 1/rgamma(at,bt); // shape and scale  E(tau2) = atbt
      tau2draw2 = 1/rgamma(at,bt); // shape and scale  E(tau2) = atbt

      for(b = 0; b < nb; b++){
        for(bb = 0; bb < nb; bb++){

          oV1[b*nb+bb] = 0.0;
          oV2[b*nb+bb] = 0.0;

          if(b == bb){

            oV1[b*nb+bb] = 1/(lamdraw1*lamdraw1);
            oV2[b*nb+bb] = 1/(lamdraw2*lamdraw2);
          }

          nV1[b*nb+bb] = tau2draw1*Kinv[b*nb+bb];
          nV2[b*nb+bb] = tau2draw2*Kinv[b*nb+bb];

        }
      }

      cholesky(nV1, nb , &ld1);
      cholesky(nV2, nb , &ld2);

      ran_mvnorm(mu1_iter, nV1, nb, scr1, thetadraw1);
      ran_mvnorm(mu2_iter, nV2, nb, scr1, thetadraw2);

      ldo1 = 2.0*nb*log(lamdraw1);
      ldo2 = 2.0*nb*log(lamdraw2);


      lgcondraw = 0.0;
      lgcatdraw = 0.0;
      
      if(!(*PPM)){
        for(p=0;p<(*ncon);p++){
  
          xcontmp = Xcon[j*(*ncon)+p];
  
          if(*similarity_function==1){ // Auxilliary
            if(*consim==1){
              lgcont = gsimconNN(m0,v2,s20,xcontmp,xcontmp*xcontmp, mnmle[p],1,0,0, 1);
              lgcondraw = lgcondraw + lgcont;
            }
            if(*consim==2){
              lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp, xcontmp*xcontmp,mnmle[p],s2mle[p], 1, 0,0, 1);
              lgcondraw = lgcondraw + lgcont;
            }
          }
          if(*similarity_function==2){ // Double Dipper
            if(*consim==1){
              lgcont = gsimconNN(m0,v2,s20,xcontmp,xcontmp*xcontmp, mnmle[p], 1, 1, 0, 1);
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
  
          // similarity for categorical covariate
        for(p=0;p<(*ncat);p++){
  
            for(c=0;c<Cvec[p];c++){nhctmp[c] = 0;}
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
  
        gtilY[nclus_iter] = lgcondraw + lgcatdraw;
        gtilN[nclus_iter] = lgcondraw + lgcatdraw;
      }


      ph[nclus_iter] = dmvnorm(btmp1,thetadraw1,oV1,nb,ldo1,scr1,1) +
                       dmvnorm(btmp2,thetadraw2,oV2,nb,ldo2,scr2,1) +
                         log(Mdp) +
                         lgcondraw +
                         lgcatdraw;



      
      /////////////////////////////////////////////////////////////////////////////
      // This is the calibration used when the similarity is standardized by
      // the sum of all cluster similarity values.
      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////
      if(*calibrate==1){
  
        for(b = 0; b < nb; b++){
  
          thtmp1[b] = thetah1[b*(nclus_iter) + k];
          thtmp2[b] = thetah2[b*(nclus_iter) + k];
  
          for(bb = 0; bb < nb; bb++){
  
            oV1[b*nb+bb] = 0.0;
            oV2[b*nb+bb] = 0.0;
  
            if(b == bb){
  
              oV1[b*nb+bb] = 1/(lamh1[k]*lamh1[k]);
              oV2[b*nb+bb] = 1/(lamh2[k]*lamh2[k]);
            }
          }
        }
  
        ldo1 = 2.0*nb*log(lamh1[k]);
        ldo2 = 2.0*nb*log(lamh2[k]);
  
  
        maxgtilY = gtilY[0];
        maxgtilN = gtilN[0];
        for(k=1; k < nclus_iter+1; k++){
          if(maxgtilN < gtilN[k]) maxgtilN = gtilY[k];
          if(maxgtilY < gtilY[k]) maxgtilY = gtilY[k];
        }

        sgY=0.0;
        sgN=0.0;
        for(k=0; k<nclus_iter+1; k++){
          gtilN[k] = exp(gtilN[k] - maxgtilY);
          sgN = sgN + gtilN[k];
  
          gtilY[k] = exp(gtilY[k] - maxgtilY);
          sgY = sgY + gtilY[k];
        }
  
      	// Calibrate the unnormalized cluster probabilities
        for(k=0; k<nclus_iter; k++){
          lgtilY = log(gtilY[k]) - log(sgY);
          lgtilN = log(gtilN[k]) - log(sgN);
  
          ph[k] = dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1) +
                  dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1) +
                      log((double) nh[k]) +
                      lgtilY - lgtilN; //This takes into account both cont and cat vars
  
        }
  
      	// calibration for a singleton
        ph[nclus_iter] = dmvnorm(btmp1,thetadraw1,oV1,nb,ldo1,scr1,1) +
                         dmvnorm(btmp2,thetadraw2,oV2,nb,ldo2,scr2,1) +
                             log(Mdp) +
                             log(gtilN[nclus_iter]) - log(sgN);
  
      }

      maxph = ph[0];
      for(k = 1; k < nclus_iter+1; k++){
      	if(maxph < ph[k]) maxph = ph[k];
      }
      
      
//      RprintVecAsMat("ph", ph, 1, nclus_iter+1);      
      denph = 0.0;
      for(k = 0; k < nclus_iter+1; k++){
        ph[k] = exp(ph[k] - maxph);
        denph = denph + ph[k];
      }

      for(k = 0; k < nclus_iter+1; k++){
        probh[k] = ph[k]/denph;
      }
//      RprintVecAsMat("probh", probh, 1, nclus_iter+1);      

      uu = runif(0.0,1.0);

      cprobh= 0.0;;
      iaux = nclus_iter+1;
      for(k = 0; k < nclus_iter+1; k++){
        cprobh = cprobh + probh[k];
        if (uu < cprobh){
          iaux = k+1;
          break;
        }
      }

      if(iaux <= nclus_iter){

        Si_iter[j] = iaux;
        nh[Si_iter[j]-1] = nh[Si_iter[j]-1] + 1;

      }else{

        nclus_iter = nclus_iter + 1;
        Si_iter[j] = nclus_iter;
        nh[Si_iter[j]-1] = 1;

        tau2h1[Si_iter[j]-1] = tau2draw1;
        tau2h2[Si_iter[j]-1] = tau2draw2;

        lamh1[Si_iter[j]-1] = lamdraw1;
        lamh2[Si_iter[j]-1] = lamdraw2;

        for(b = 0; b < nb; b++){
          thetah1[b*(nclus_iter) + Si_iter[j]-1] = thetadraw1[b];
          thetah2[b*(nclus_iter) + Si_iter[j]-1] = thetadraw2[b];
        }
      }
//      RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);      
//      RprintIVecAsMat("nh", nh, 1, nclus_iter);      
      
      // need to now add the xcon to the cluster to which it was assigned;
      if(!(*PPM)){
        for(p = 0; p < *ncon; p++){
          sumx[(Si_iter[j]-1)*(*ncon) + p] = sumx[(Si_iter[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p];
          sumx2[(Si_iter[j]-1)*(*ncon) + p] = sumx2[(Si_iter[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
        }

        // need to now add the xcat to the cluster to which it was assigned;
        for(p = 0; p < *ncat; p++){

          nhc[((Si_iter[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
             nhc[((Si_iter[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] + 1;

        }
      }
      //RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
    }


    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update lam (using MH-step) this is the standard deviation not the variance
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      olam1 = lamh1[k];
      nlam1 = rnorm(olam1,csigLAM1);

      olam2 = lamh2[k];
      nlam2 = rnorm(olam2,csigLAM2);

      if((nlam1 > 0) & (nlam2 > 0) & (nlam1 < A1) & (nlam2 < A2)){
        for(b = 0; b < nb; b++){
          thtmp1[b] = thetah1[b*(nclus_iter) + k];
          thtmp2[b] = thetah2[b*(nclus_iter) + k];
          for(bb = 0; bb < nb; bb++){
            oV1[b*nb+bb] = 0.0;
            nV1[b*nb+bb] = 0.0;

            oV2[b*nb+bb] = 0.0;
            nV2[b*nb+bb] = 0.0;
            if(b == bb){

              oV1[b*nb+bb] = 1/(olam1*olam1);
              nV1[b*nb+bb] = 1/(nlam1*nlam1);

              oV2[b*nb+bb] = 1/(olam2*olam2);
              nV2[b*nb+bb] = 1/(nlam2*nlam2);
            }
          }
        }

        ldo1 = 2.0*nb*log(olam1);
        ldn1 = 2.0*nb*log(nlam1);

        ldo2 = 2.0*nb*log(olam2);
        ldn2 = 2.0*nb*log(nlam2);

        lln1 = 0.0, lln2=0.0;
        llo1 = 0.0, llo2=0.0;
        for(j = 0; j < *nsubject; j++){
          if(Si_iter[j] == k+1){
            for(b = 0; b < nb; b++){
              btmp1[b] = beta1_iter[b*(*nsubject) + j];
              btmp2[b] = beta2_iter[b*(*nsubject) + j];
            }

            llo1 = llo1 + dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1);
            llo2 = llo2 + dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1);

            lln1 = lln1 + dmvnorm(btmp1, thtmp1, nV1, nb, ldn1, scr1, 1);
            lln2 = lln2 + dmvnorm(btmp2, thtmp2, nV2, nb, ldn2, scr2, 1);
          }
        }

        llo1 = llo1 + dunif(olam1, 0.0, A1, 1);
        llo2 = llo2 + dunif(olam2, 0.0, A2, 1);

        lln1 = lln1 + dunif(nlam1, 0.0, A1, 1);
        lln2 = lln2 + dunif(nlam2, 0.0, A2, 1);


        llr1 = lln1 - llo1;
        llr2 = lln2 - llo2;

        uu = runif(0.0,1.0);
        if(log(uu) < llr1){
          lamh1[k] = nlam1;
        }
        uu = runif(0.0,1.0);
        if(log(uu) < llr2){
          lamh2[k] = nlam2;
        }
      }
    }


    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // update player specific betas, sigma2, and beta0;								//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    csobs=0;
    csHobs=0;
    for(j = 0; j < *nsubject; j++){

      for(t = 0; t < nobs[j]; t++){

        y1_tmp[t] = y1[csobs];
        y2_tmp[t] = y2[csobs];

        y1_b0[t] = y1_tmp[t]-beta01_iter[j];
        y2_b0[t] = y2_tmp[t]-beta02_iter[j];

        csobs = csobs+1;

      }

      if(*balanced!=1){
        for(t=0; t < nobs[j]*(nb); t++){
		  H[t] = Hmat[csHobs + t];
		}
		csHobs = csHobs + nobs[j]*(nb);

        mat_transpose(H, tH, nobs[j], nb);
        matrix_product(tH, H, HtH, nb, nb, nobs[j]);
	  }



      matrix_product(tH, y1_tmp, Hty1, nb, 1, nobs[j]);
      matrix_product(tH, y2_tmp, Hty2, nb, 1, nobs[j]);

      for(b = 0; b < nb; b++){
        for(bb = 0; bb < nb; bb++){

          Sstar1[b*nb+bb] = (1/sig21_iter[j])*HtH[b*nb+bb];
          Sstar2[b*nb+bb] = (1/sig22_iter[j])*HtH[b*nb+bb];

          if(b == bb){
            Sstar1[b*nb+bb] = (1/sig21_iter[j])*HtH[b*nb+bb] +
                              (1/(lamh1[Si_iter[j]-1]*lamh1[Si_iter[j]-1]));

            Sstar2[b*nb+bb] = (1/sig22_iter[j])*HtH[b*nb+bb] +
                              (1/(lamh2[Si_iter[j]-1]*lamh2[Si_iter[j]-1]));
          }
        }
      }

      cholesky(Sstar1, nb, &ld1);
      inverse_from_cholesky(Sstar1, scr1, scr2, nb);

      cholesky(Sstar2, nb, &ld2);
      inverse_from_cholesky(Sstar2, scr1, scr2, nb);

      matrix_product(tH, y1_b0, scr1,  nb, 1, nobs[j]);
      matrix_product(tH, y2_b0, scr2,  nb, 1, nobs[j]);

      for(b = 0; b < nb; b++){
        scr3[b] = (1/sig21_iter[j])*scr1[b] +
                  (1/(lamh1[Si_iter[j]-1]*lamh1[Si_iter[j]-1]))*thetah1[b*(nclus_iter) + Si_iter[j]-1];

        scr4[b] = (1/sig22_iter[j])*scr2[b] +
                  (1/(lamh2[Si_iter[j]-1]*lamh2[Si_iter[j]-1]))*thetah2[b*(nclus_iter) + Si_iter[j]-1];

      }

      matrix_product(Sstar1, scr3, Mstar1, nb, 1, nb);
      matrix_product(Sstar2, scr4, Mstar2, nb, 1, nb);

      cholesky(Sstar1, nb , &ld1);
      cholesky(Sstar2, nb , &ld2);

      ran_mvnorm(Mstar1, Sstar1, nb, scr1, outrmvnorm1);
      ran_mvnorm(Mstar2, Sstar2, nb, scr2, outrmvnorm2);

      for(b = 0; b < nb; b++){

        beta1_iter[b*(*nsubject) + j] = outrmvnorm1[b];
        btmp1[b] = beta1_iter[b*(*nsubject) + j];

        beta2_iter[b*(*nsubject) + j] = outrmvnorm2[b];
        btmp2[b] = beta2_iter[b*(*nsubject) + j];


      }


      ///////////////////////////////////////////
      //									   //
      // udate sigma2 within the same loop.    //
      //									   //
      ///////////////////////////////////////////
      // Start with sig21
      matrix_product(H, btmp1, Hb1, nobs[j], 1, nb);

      os = sqrt(sig21_iter[j]);
      ns = rnorm(os, csigSIG);

      if(ns > 0.0){

        for(jj = 0; jj < nobs[j]; jj++){
          scr3[jj] = y1_b0[jj] - Hb1[jj];
        }

        sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

        llo = -0.5*nobs[j]*log(os*os) - 0.5*(1/(os*os))*sumsq;
        lln = -0.5*nobs[j]*log(ns*ns) - 0.5*(1/(ns*ns))*sumsq;

        llo = llo + dunif(os, 0.0, Asig, 1);
        lln = lln + dunif(ns, 0.0, Asig, 1);

        llr = lln - llo;

        uu = runif(0.0,1.0);
        if(llr > log(uu)) sig21_iter[j] = ns*ns;

      }


      // Now for sig22
      matrix_product(H, btmp2, Hb2, nobs[j], 1, nb);

      os = sqrt(sig22_iter[j]);
      ns = rnorm(os, csigSIG);

      if(ns > 0.0){

        for(jj = 0; jj < nobs[j]; jj++){
          scr3[jj] = y2_b0[jj] - Hb2[jj];
        }

        sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

        llo = -0.5*nobs[j]*log(os*os) - 0.5*(1/(os*os))*sumsq;
        lln = -0.5*nobs[j]*log(ns*ns) - 0.5*(1/(ns*ns))*sumsq;

        llo = llo + dunif(os, 0.0, Asig, 1);
        lln = lln + dunif(ns, 0.0, Asig, 1);

        llr = lln - llo;

        uu = runif(0.0,1.0);
        if(llr > log(uu)) sig22_iter[j] = ns*ns;

      }


/*
      matrix_product(H, btmp1, Hb1, nobs[j], 1, nb);
      matrix_product(H, btmp2, Hb2, nobs[j], 1, nb);

      sumy1_Hb = 0.0;
      sumy2_Hb = 0.0;

      for(jj = 0; jj < nobs[j]; jj++){

        scr3[jj] = y1_b0[jj] - Hb1[jj];
        scr4[jj] = y2_b0[jj] - Hb2[jj];
        sumy1_Hb = sumy1_Hb + (y1_tmp[jj] - Hb1[jj]);
        sumy2_Hb = sumy2_Hb + (y2_tmp[jj] - Hb2[jj]);
      }

      sumsq1 = inner_product(scr3, 1, scr3, 1, nobs[j]);
      sumsq2 = inner_product(scr4, 1, scr4, 1, nobs[j]);

      astar = 0.5*nobs[j] + asig;
      bstar1 = 0.5*sumsq1 + 1/bsig;
      bstar2 = 0.5*sumsq2 + 1/bsig;



      // these are for the hierarchical variance structure
      //			astar = 0.5*nobs[j] + 0.5*nuh[Si_iter[j]];
      //			bstar = 0.5*sumsq + 0.5;


      //			Rprintf("astar = %f\n", astar);
      //			Rprintf("bstar = %f\n", bstar);
      //bstar is rate and rgamma requires scale hence inverse
      sig21_iter[j] = 1/rgamma(astar, 1/bstar1);
      sig22_iter[j] = 1/rgamma(astar, 1/bstar2);

      //			Rprintf("sig21 = %f\n", sig21_iter[j]);
      //			Rprintf("sig22 = %f\n", sig22_iter[j]);

*/
      ////////////////////////////////////////////
      //									      //
      // update beta0 within in the same loop;  //
      //									      //
      ////////////////////////////////////////////
      sumy1_Hb = 0.0;
      sumy2_Hb = 0.0;

      for(jj = 0; jj < nobs[j]; jj++){
        sumy1_Hb = sumy1_Hb + (y1_tmp[jj] - Hb1[jj]);
        sumy2_Hb = sumy2_Hb + (y2_tmp[jj] - Hb2[jj]);
      }



      s2star1 = 1.0/((nobs[j]/sig21_iter[j]) + (1/sig2b01_iter));
      s2star2 = 1.0/((nobs[j]/sig22_iter[j]) + (1/sig2b02_iter));

      mstar1 = s2star1*((1/sig21_iter[j])*sumy1_Hb + (1/sig2b01_iter)*mub01_iter);
      mstar2 = s2star2*((1/sig22_iter[j])*sumy2_Hb + (1/sig2b02_iter)*mub02_iter);

      beta01_iter[j] = rnorm(mstar1, sqrt(s2star1));
      beta02_iter[j] = rnorm(mstar2, sqrt(s2star2));

      //////////////////////////////////////////////////////////////////////////
      //									      								//
      // Evaluate the likelihood for each observation used to calculate lpml  //
      //									      								//
      //////////////////////////////////////////////////////////////////////////
      if((i > (*burn-1)) & (i % (*thin) == 0)){
        like0=0;
        llikeval = 0.0;
        for(t = 0; t < nobs[j]; t++){
          llikeval = llikeval + dnorm(y1_tmp[t], beta01_iter[j] + Hb1[t], sqrt(sig21_iter[j]), 1)*
                                dnorm(y2_tmp[t], beta02_iter[j] + Hb2[t], sqrt(sig22_iter[j]), 1);
        }

        llike_iter[j] = llikeval;

        // These are needed for WAIC
        mnlike[j] = mnlike[j] + exp(llike_iter[j])/(double) nout;
        mnllike[j] = mnllike[j] + (llike_iter[j])/(double) nout;


        if(exp(llike_iter[j]) < 1e-320) like0=1;

        if(like0==1) nout_0 = nout_0+1;

        if(like0==0){
          CPO[j] = CPO[j] + (1/exp(llike_iter[j]));
          //				Rprintf("like = %f\n", likeval);
        }
      }
    }



    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate mub0 a global intercept of mean of beta0i     						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    sumb01 = 0.0;
    sumb02 = 0.0;
    for(j = 0; j < *nsubject; j++){

      sumb01 = sumb01 + beta01_iter[j];
      sumb02 = sumb02 + beta02_iter[j];

    }

    s2star1 = 1.0/((*nsubject)/sig2b01_iter + (1/s2b0));
    s2star2 = 1.0/((*nsubject)/sig2b02_iter + (1/s2b0));

    mstar1 = s2star1*((1/sig2b01_iter)*(sumb01) + (1/s2b0)*mb0);
    mstar2 = s2star2*((1/sig2b02_iter)*(sumb02) + (1/s2b0)*mb0);

    //		Rprintf("mstar = %f\n", mstar);
    //		Rprintf("s2star = %f\n", s2star);

    mub01_iter = rnorm(mstar1, sqrt(s2star1));
    mub02_iter = rnorm(mstar2, sqrt(s2star2));

    //		Rprintf("mub01 = %f\n", mub01_iter);
    //		Rprintf("mub02 = %f\n", mub02_iter);
    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate sig2b0 a global intercept of mean of beta0i     						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    sumsq1 = 0.0;
    sumsq2 = 0.0;
    for(j = 0; j < *nsubject; j++){

      sumsq1 = sumsq1 + (beta1_iter[j] - mub01_iter)*(beta1_iter[j] - mub01_iter);
      sumsq2 = sumsq2 + (beta2_iter[j] - mub02_iter)*(beta2_iter[j] - mub02_iter);

    }

    astar = 0.5*(*nsubject) + ab0;
    bstar1 = 0.5*sumsq1 + 1/bb0;
    bstar2 = 0.5*sumsq2 + 1/bb0;

    sig2b01_iter = 1/rgamma(astar, 1/bstar1);
    sig2b02_iter = 1/rgamma(astar, 1/bstar2);


    //		Rprintf("sig2b01 = %f\n", sig2b01_iter);
    //		Rprintf("sig2b02 = %f\n", sig2b02_iter);

    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate thetah each of the cluster specific coefficients;						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    //		Rprintf("nb = %d\n", nb);
    //		RprintIVecAsMat("nc", nc, 1, nclus_iter);
    //		RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);
    //		RprintVecAsMat("lamh", lamh, 1, nclus_iter);

    for(k = 0; k < nclus_iter; k++){

      //			Rprintf("k = %d =================== \n", k);
      //			Rprintf("tau2h =%f\n", tau2h[k]);
      //			Rprintf("lamh =%f\n", lamh[k]);
      //			RprintVecAsMat("K", K, nb, nb);
      //			Rprintf("nk = %d\n", nh[k]);
      for(b = 0; b < nb; b++){
        for(bb = 0; bb < nb; bb++){

          Sstar1[b*nb+bb] = (1/tau2h1[k])*K[b*nb+bb];
          Sstar2[b*nb+bb] = (1/tau2h2[k])*K[b*nb+bb];

          if(b == bb){
            Sstar1[b*nb+bb] = ((double) nh[k]/(lamh1[k]*lamh1[k])) +
              (1/tau2h1[k])*K[b*nb+bb];

            Sstar2[b*nb+bb] = ((double) nh[k]/(lamh2[k]*lamh2[k])) +
              (1/tau2h2[k])*K[b*nb+bb];

          }
        }

        sumbeta1[b] = 0.0;
        sumbeta2[b] = 0.0;

      }

      //			RprintVecAsMat("Sstar",Sstar, nb, nb);

      cholesky(Sstar1, nb, &ld1);
      inverse_from_cholesky(Sstar1, scr1, scr2, nb);

      cholesky(Sstar2, nb, &ld2);
      inverse_from_cholesky(Sstar2, scr1, scr2, nb);

      //			RprintVecAsMat("Sstar",Sstar, nb, nb);
      //			RprintVecAsMat("mu_iter = ", mu_iter, 1, nb);

      matrix_product(K, mu1_iter, scr1, nb, 1, nb);
      matrix_product(K, mu2_iter, scr2, nb, 1, nb);

      //			RprintVecAsMat("Kmu",scr1, 1, nb);
      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);

      for(j = 0; j < *nsubject; j++){

        if(Si_iter[j] == k+1){

          for(b = 0; b < nb; b++){

            sumbeta1[b] = sumbeta1[b] + (1/(lamh1[k]*lamh1[k]))*
              beta1_iter[b*(*nsubject) + j];

            sumbeta2[b] = sumbeta2[b] + (1/(lamh2[k]*lamh2[k]))*
              beta2_iter[b*(*nsubject) + j];
          }
        }
      }

      for(b=0; b<nb; b++){
        sumbeta1[b] = sumbeta1[b] + (1/tau2h1[k])*scr1[b];
        sumbeta2[b] = sumbeta2[b] + (1/tau2h2[k])*scr2[b];
      }

      //			RprintVecAsMat("sumbeta1", sumbeta1, 1, nb);


      matrix_product(Sstar1, sumbeta1, Mstar1, nb, 1, nb);
      matrix_product(Sstar2, sumbeta2, Mstar2, nb, 1, nb);

      //			RprintVecAsMat("Mstar1", Mstar1, 1, nb);
      //			RprintVecAsMat("Sstar", Sstar, nb, nb);
      //			RprintVecAsMat("Mstar2", Mstar2, 1, nb);
      //			RprintVecAsMat("Sstar", Sstar, nb, nb);


      cholesky(Sstar1, nb , &ld1);
      cholesky(Sstar2, nb , &ld2);


      ran_mvnorm(Mstar1, Sstar1, nb, scr1, outrmvnorm1);
      ran_mvnorm(Mstar2, Sstar2, nb, scr2, outrmvnorm2);


      //RprintVecAsMat("thetah1", outrmvnorm1, 1, nb);

      for(b = 0; b < nb; b++){

        thetah1[b*(nclus_iter) + k] = outrmvnorm1[b];
        thetah2[b*(nclus_iter) + k] = outrmvnorm2[b];


      }

      //			RprintVecAsMat("thetah1", thetah1, nb, nclus_iter);

    }

    // Rprintf("nclus_iter = %d\n", nclus_iter);
    //		RprintVecAsMat("thetah", thetah, nb, *nsubject);
    //		RprintVecAsMat("thetah1", thetah1, nb, nclus_iter);
    //		RprintVecAsMat("thetah2", thetah2, nb, nclus_iter);
    //		RprintVecAsMat("thetah1", thetah1, 1, nb*(*nsubject));
    //		RprintVecAsMat("thetah2", thetah2, nb, nclus_iter);




    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update tau2 for each of the clusters (P-spline smoothing parameter)
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      //			Rprintf("k = %d =================== \n", k);

      for(b = 0; b < nb; b++){

        thtmp1[b] = thetah1[b*(nclus_iter) + k] - mu1_iter[b];
        thtmp2[b] = thetah2[b*(nclus_iter) + k] - mu2_iter[b];

      }

      //			RprintVecAsMat("thtmp1", thtmp1, 1, nb);
      //			RprintVecAsMat("mu1_iter", mu1_iter, 1, nb);
      //			RprintVecAsMat("thtmp2", thtmp2, 1, nb);
      //			RprintVecAsMat("mu2_iter", mu2_iter, 1, nb);

      sumsq1 = quform(thtmp1,K,nb);
      sumsq2 = quform(thtmp2,K,nb);

      //			Rprintf("sumsq1 = %f\n", sumsq1);
      //			Rprintf("sumsq2 = %f\n", sumsq2);

      astar = 0.5*(nb) + at;
      bstar1 = 1/bt + 0.5*sumsq1;
      bstar2 = 1/bt + 0.5*sumsq2;

      //			Rprintf("astar = %f\n", astar);
      //			Rprintf("bstar1 = %f\n", bstar1);
      //			Rprintf("bstar2 = %f\n", bstar2);

      tau2h1[k] = 1/rgamma(astar, 1/bstar1);// E(tau2) = astarbstar for gamma.  bstar is scale
      tau2h2[k] = 1/rgamma(astar, 1/bstar2);// E(tau2) = astarbstar for gamma.  bstar is scale
      //			tau2h[k] = 0.01;

    }





    //		RprintVecAsMat("tau2h1", tau2h1, 1, nclus_iter);
    //		RprintVecAsMat("tau2h2", tau2h2, 1, nclus_iter);
    //		Rprintf("nclus_iter = %d\n", nclus_iter);







    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate mu center of cluster specific thetah;									//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////

    //		Rprintf("nb = %d\n", nb);
    //		RprintIVecAsMat("nc", nc, 1, nclus_iter);
    //		RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);

    for(b = 0; b < nb; b++){
      sumtheta1[b] = 0.0;
      sumtheta2[b] = 0.0;
    }

    sumtau21 = 0.0;
    sumtau22 = 0.0;

    for(k = 0; k < nclus_iter; k++){
      sumtau21 = sumtau21 + (1/tau2h1[k]);
      sumtau22 = sumtau22 + (1/tau2h2[k]);
      for(b = 0; b < nb; b++){
        sumtheta1[b] = sumtheta1[b] + (1/tau2h1[k])*thetah1[b*(nclus_iter) + k];
        sumtheta2[b] = sumtheta2[b] + (1/tau2h2[k])*thetah2[b*(nclus_iter) + k];
      }
    }

    //		Rprintf("sumtau2 = %f\n", sumtau2);
    //		RprintVecAsMat("sumtheta",sumtheta, 1, nb);


    for(b = 0; b < nb; b++){
      for(bb = 0; bb < nb; bb++){

        Sstar1[b*nb+bb] = sumtau21*K[b*nb+bb];
        Sstar2[b*nb+bb] = sumtau22*K[b*nb+bb];

        if(b == bb){
          Sstar1[b*nb+bb] = sumtau21*K[b*nb+bb]+(1/s2mu) ;
          Sstar2[b*nb+bb] = sumtau22*K[b*nb+bb]+(1/s2mu) ;
        }
      }

    }

    //		RprintVecAsMat("sumtheta",sumtheta, 1, nb);
    //		RprintVecAsMat("Sstar",Sstar, nb, nb);

    matrix_product(K, sumtheta1, scr3, nb, 1, nb);
    matrix_product(K, sumtheta2, scr4, nb, 1, nb);

    //		RprintVecAsMat("sumthetaK",scr3, 1, nb);

    cholesky(Sstar1, nb, &ld1);
    inverse_from_cholesky(Sstar1, scr1, scr2, nb);

    cholesky(Sstar2, nb, &ld2);
    inverse_from_cholesky(Sstar2, scr1, scr2, nb);

    //		RprintVecAsMat("Sstar",Sstar, nb, nb);


    matrix_product(Sstar1, scr3, Mstar1, nb, 1, nb);
    matrix_product(Sstar2, scr4, Mstar2, nb, 1, nb);

    //		RprintVecAsMat("Mstar", Mstar, 1, nb);
    //		RprintVecAsMat("Sstar", Sstar, nb, nb);

    cholesky(Sstar1, nb , &ld1);
    cholesky(Sstar2, nb , &ld2);

    ran_mvnorm(Mstar1, Sstar1, nb, scr1, outrmvnorm1);
    ran_mvnorm(Mstar2, Sstar2, nb, scr2, outrmvnorm2);


    //		RprintVecAsMat("mu", outrmvnorm, 1, nb);

    for(b = 0; b < nb; b++){
      mu1_iter[b] = outrmvnorm1[b];
      mu2_iter[b] = outrmvnorm2[b];
//      mu1_iter[b] = 0;
//      mu2_iter[b] = 0;
    }

    //		RprintVecAsMat("thetah", thetah, nb, nclus_iter);


    //		RprintVecAsMat("mu1_iter", mu1_iter, 1, nb);
    //		RprintVecAsMat("mu2_iter", mu2_iter, 1, nb);



    //////////////////////////////////////////////////////////////////////////////////
    //
    // Posterior predictives.  I'll use ages 19, 20 ,21 and the three role variables;
    //
    //////////////////////////////////////////////////////////////////////////////////
/*
    if((i > (*burn-1)) & i % ((*thin) == 0)){

      for(pp = 0; pp < *npred; pp++){

        //				Rprintf("pp = %d ================================================ \n", pp+1);



        //				Rprintf("nclus_iter = %d\n", nclus_iter);
        for(k = 0; k < nclus_iter; k++){
          //					Rprintf("k = %d  ========== \n", k);


          lgconN=0.0, lgconY=0.0;
          for(p=0; p<(*ncon); p++){
            //						Rprintf("p = %d ====== \n", p) ;
            nhtmp = 0;
            for(j = 0; j < *nobs; j++){
              if(Si_iter[j] == k+1){
                xcontmp[nhtmp] = Xcon[j*(*ncon)+p]; //create cluster specific x-vector
                nhtmp = nhtmp+1;
                //								Rprintf("nhtmp = %d\n", nhtmp);
              }
            }

            //						Rprintf("nhtmp = %d\n", nhtmp);
            //						Rprintf("nh[k] = %d\n", nh[k]);
            //						RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp);

            sumx = 0.0;
            sumx2 = 0.0;
            for(t = 0; t < nhtmp; t++){

              sumx = sumx + xcontmp[t];
              sumx2 = sumx2 + xcontmp[t]*xcontmp[t];

            }

            //						Rprintf("sumx = %f\n", sumx);
            //						Rprintf("sumx2 = %f\n", sumx2);
            if(*similarity_function==1){
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
                //								lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 1, 1);
                lgconN = lgconN + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
                //								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 1, 1);
                lgconN = lgconN + lgcont;
              }
            }
            if(*similarity_function==2){
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
                //								lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 1, 1);
                lgconN = lgconN + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
                //								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 1, 1);
                lgconN = lgconN + lgcont;
              }
            }
            if(*similarity_function==3){
              lgcont = gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
              lgconN = lgconN + lgcont;
            }
            //						Rprintf("lgconN = %f\n", lgconN);


            // now add ppth prediction to cluster;
            //						Rprintf("xconpred[pp] = %f\n", xconpred[pp]);
            xcontmp[nhtmp] = Xconp[pp*(*ncon)+p];
            sumx = sumx + Xconp[pp*(*ncon)+p];
            sumx2 = sumx2 + Xconp[pp*(*ncon)+p]*Xconp[pp*(*ncon)+p];
            nhtmp = nhtmp + 1;

            //						Rprintf("nhtmp = %d\n", nhtmp);
            //						RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp);

            //						Rprintf("sumx = %f\n", sumx);
            //						Rprintf("sumx2 = %f\n", sumx2);
            if(*similarity_function==1){ // Auxilliary
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
                //								lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 1, 1);
                lgconY = lgconY + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
                //								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 1, 1);
                lgconY = lgconY + lgcont;
              }
            }
            if(*similarity_function==2){ // Double Dipper
              if(*consim==1){
                lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
                //								lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 1, 1);
                lgconY = lgconY + lgcont;
              }
              if(*consim==2){
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
                //								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 1, 1);
                lgconY = lgconY + lgcont;
              }
            }
            if(*similarity_function==3){ // Variance
              lgcont = gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
              lgconY = lgconY + lgcont;
            }

            //						Rprintf("lgconY = %f\n", lgconY);

            //						RprintVecAsMat("calmatconpY", calmatconpY, nclus_iter, *ncon);

          }
          //					Rprintf("lgconY - lgconN = %f\n", lgconY - lgconN);




          //					Rprintf("*ncat = %d", *ncat);
          lgcatY=0.0, lgcatN=0.0;

          for(p=0; p<(*ncat); p++){
            //						Rprintf("p = %d ====== \n", p) ;
            for(c=0;c<Cvec[p];c++){nhc[c]=0;}
            nhtmp=0;
            for(j = 0; j < *nobs; j++){
              //							Rprintf("j = %d\n", j);
              //							Rprintf("Si_iter[j] = %d\n", Si_iter[j]);
              //							Rprintf("Xcat[j*(*ncat)+p] = %d\n", Xcat[j*(*ncat)+p]);

              if(Si_iter[j]==k+1){
                nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1; // this needs to be a vector
                nhtmp = nhtmp+1;

                //                              	RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
              }
            }


            //						RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
            //						Rprintf("nhtmp =%d\n", nhtmp);

            if(*similarity_function==1){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
              lgcatN = lgcatN + lgcatt;
            }
            if(*similarity_function==2){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
              lgcatN = lgcatN + lgcatt;
            }
            if(*similarity_function==3){
              //							lgcatt = gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
              lgcatt = 0.0;
              for(c=0;c<Cvec[p];c++){
                if(nhc[c]==0){
                  lgcatt = lgcatt + 0;
                }else{
                  lgcatt = lgcatt + -((double) nhc[c]/(double) nhtmp)*(
                    log((double) nhc[c]/(double) nhtmp)/log(2));
                }
              }
              lgcatN = lgcatN + -(alpha)*lgcatt;
            }
            //						Rprintf("lgcatN = %f\n", lgcatN);


            //						RprintVecAsMat("calmatcatpN", calmatcatpN, nclus_iter, *ncat);

            //						Rprintf("Xcatp[pp*(*ncat)+p] = %d\n", Xcatp[pp*(*ncat)+p]);

            nhc[Xcatp[pp*(*ncat)+p]] = nhc[Xcatp[pp*(*ncat)+p]] + 1;
            nhtmp=nhtmp + 1;

            //						RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);

            if(*similarity_function==1){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
              lgcatY = lgcatY + lgcatt;
            }
            if(*similarity_function==2){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
              lgcatY = lgcatY + lgcatt;
            }
            if(*similarity_function==3){// Use entropy
              //							lgcatt = gsimconEV(sumx, sumx2, nhtmp, alpha,  1);
              lgcatt = 0.0;
              for(c=0;c<Cvec[p];c++){
                if(nhc[c]==0){
                  lgcatt = lgcatt + 0;
                }else{
                  lgcatt = lgcatt + -((double) nhc[c]/(double) nhtmp)*(
                    log((double) nhc[c]/(double) nhtmp)/log(2));
                }
              }
              lgcatY = lgcatY + -(alpha)*lgcatt;
            }
            //						Rprintf("lgcatY = %f\n", lgcatY);

          }


          gtilY[k] = lgconY + lgcatY;
          gtilN[k] = lgconN + lgcatN;


        }

        //				RprintVecAsMat("ph", ph, 1, nclus_iter);
        //				RprintVecAsMat("ph1", ph1, 1, nclus_iter);

        lgcon0=0.0;
        for(p=0;p<*ncon;p++){
          xcontmp[0] = Xconp[pp*(*ncon)+p];
          if(*similarity_function==1){
            if(*consim==1){
              lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0,1);
              //							lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,1,1);
            }
            if(*consim==2){
              lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,0,1);
              //							lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,1,1);
            }
          }
          if(*similarity_function==2){
            if(*consim==1){
              lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,0,1);
              //							lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,1,1);
            }
            if(*consim==2){
              lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 0,1);
              //							lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 1,1);
            }
          }
          if(*similarity_function==3){
            lgcon0 = lgcon0 + gsimconEV(xcontmp[0],xcontmp[0]*xcontmp[0],1,alpha,1);
          }
        }



        lgcat0 = 0.0;

        for(p=0;p<(*ncat);p++){
          for(c=0;c<Cvec[p];c++) nhc[c] = 0;

          nhc[Xcatp[pp*(*ncat)+p]] = 1;
          //					RprintIVecAsMat("nhc =", nhc, 1, Cvec[p]);

          if(*similarity_function==1){
            lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
          }
          if(*similarity_function==2){
            lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
          }
          if(*similarity_function==3){
            //						lgcat0 = lgcat0 + gsimconEV(xcattmp[0], xcattmp[0]*xcattmp[0],  1,alpha, 1);
            lgcat0 = lgcat0 + -(alpha)*0;

          }


        }


        gtilY[nclus_iter] = lgcat0 + lgcon0;
        gtilN[nclus_iter] = lgcat0 + lgcon0;

        ph[nclus_iter] = log((double) Mdp) + lgcon0 + lgcat0;

        if(*similarity_function==4) ph[nclus_iter] = log((double) Mdp) + log(1);


        //				RprintVecAsMat("ph", ph, 1, nclus_iter+1);

        if(*PPM) ph[nclus_iter] = log((double) Mdp);

        denph = 0.0;
        for(k = 0; k < nclus_iter+1; k++){

          //					ph[k] = exp(ph[k] - maxph);
          //					ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
          denph = denph + ph[k];

        }

        //				RprintVecAsMat("ph", ph, 1, nclus_iter+1);

        for(k = 0; k < nclus_iter+1; k++){

          probh[k] = ph[k]/denph;

        }
        //				Rprintf("denph = %f\n", denph);

        //				RprintVecAsMat("probh", probh, 1, nclus_iter+1);

        uu = runif(0.0,1.0);
        //				Rprintf("uu = %f\n", uu);

        cprobh= 0.0;

        for(k = 0; k < nclus_iter+1; k++){

          cprobh = cprobh + probh[k];

          if (uu < cprobh){

            iaux = k+1;
            break;
          }

        }


        //				Rprintf("iaux = %d\n", iaux);
        //				Rprintf("nb = %d\n", nb);
        if(iaux < nclus_iter){

          for(b = 0; b < nb; b++){

            thetatmp1[b] = thetah1[b*(nclus_iter) + iaux-1];
            thetatmp2[b] = thetah2[b*(nclus_iter) + iaux-1];

          }


          for(b = 0; b < nb; b++){
            for(bb = 0; bb < nb; bb++){

              Sstar1[b*nb+bb] = 0.0;
              Sstar2[b*nb+bb] = 0.0;

              // THis is the Cholesky Decomposition as needed in ran_mvnorm
              if(b == bb){
                Sstar1[b*nb+bb] = sqrt(lamh1[iaux-1]*lamh1[iaux-1]);
                Sstar2[b*nb+bb] = sqrt(lamh2[iaux-1]*lamh2[iaux-1]);
              }
            }
          }

          //					RprintVecAsMat("thetatmp = ", thetatmp, 1, nb);
        }else{

          tau2draw1 = 1/rgamma(at,bt);
          tau2draw2 = 1/rgamma(at,bt);
          //					Rprintf("tau2draw = %f\n", tau2draw);
          lamdraw1 = runif(0,A1);
          lamdraw2 = runif(0,A2);
          //					Rprintf("lamdraw = %f\n", lamdraw);

          for(b = 0; b < nb; b++){
            for(bb = 0; bb < nb; bb++){

              nV1[b*nb+bb] = tau2draw1*Kinv[b*nb+bb];
              nV2[b*nb+bb] = tau2draw2*Kinv[b*nb+bb];

              Sstar1[b*nb+bb] = 0.0;
              Sstar2[b*nb+bb] = 0.0;
              // THis is the Cholesky Decomposition as needed in ran_mvnorm
              if(b == bb){
                Sstar1[b*nb+bb] = sqrt(lamdraw1*lamdraw1);
                Sstar2[b*nb+bb] = sqrt(lamdraw2*lamdraw2);
              }
            }


          }


          //					RprintVecAsMat("nV = ", nV, nb, nb);
          cholesky(nV1, nb , &ld1);
          cholesky(nV2, nb , &ld2);

          //					RprintVecAsMat("mu_iter =", mu_iter, 1, nb);
          ran_mvnorm(mu1_iter, nV1, nb, scr1, thetatmp1);
          ran_mvnorm(mu2_iter, nV2, nb, scr2, thetatmp2);




        }


        //				RprintVecAsMat("thetatmp = ", thetatmp, 1, nb);
        //				RprintVecAsMat("Sstar = ", Sstar, nb, nb);


        ran_mvnorm(thetatmp1, Sstar1, nb, scr1, bpred1);
        ran_mvnorm(thetatmp2, Sstar2, nb, scr2, bpred2);


        //				RprintVecAsMat("bpred = ", bpred, 1, nb);

        matrix_product(Hpred, bpred1, ppredtmp1, *npredobs, 1, nb);
        matrix_product(Hpred, bpred2, ppredtmp2, *npredobs, 1, nb);

        //				RprintVecAsMat("ppredtmp", ppredtmp, 1, *npredobs);

      }
    }
*/




    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // Save MCMC iterates															//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////

    if((i > (*burn-1)) & (i % (*thin) == 0)){


      nclus[ii] = nclus_iter;


      mub01[ii] = mub01_iter;
      mub02[ii] = mub02_iter;

      sig2b01[ii] = sig2b01_iter;
      sig2b02[ii] = sig2b02_iter;


      for(b = 0; b < nb; b++){
        mu1[ii*(nb) + b] = mu1_iter[b];
        mu2[ii*(nb) + b] = mu2_iter[b];
      }


      for(j = 0; j < *nsubject; j++){

        Si[ii*(*nsubject) + j] = Si_iter[j];


        sig21[ii*(*nsubject) + j] = sig21_iter[j];
        sig22[ii*(*nsubject) + j] = sig22_iter[j];

        beta01[ii*(*nsubject) + j] = beta01_iter[j];
        beta02[ii*(*nsubject) + j] = beta02_iter[j];

        llike[ii*(*nsubject) + j] = llike_iter[j];

        lam1[ii*(*nsubject) + j] = lamh1[Si_iter[j]-1];
        lam2[ii*(*nsubject) + j] = lamh2[Si_iter[j]-1];

        tau21[ii*(*nsubject) + j] = tau2h1[Si_iter[j]-1];
        tau22[ii*(*nsubject) + j] = tau2h2[Si_iter[j]-1];


        for(b = 0; b < nb; b++){

          beta1[(ii*(nb) + b)*(*nsubject) + j] = beta1_iter[b*(*nsubject) + j];
          beta2[(ii*(nb) + b)*(*nsubject) + j] = beta2_iter[b*(*nsubject) + j];

          theta1[(ii*(nb) + b)*(*nsubject) + j] = thetah1[b*(nclus_iter) + Si_iter[j]-1];
          theta2[(ii*(nb) + b)*(*nsubject) + j] = thetah2[b*(nclus_iter) + Si_iter[j]-1];


        }

      }

      ii = ii+1;

    }


  }

  lpml_iter=0.0;
  for(j = 0; j < *nsubject; j++){
    //		Rprintf("j = %d\n", j);

    //		Rprintf("CPO = %f\n", CPO[j]);

    lpml_iter = lpml_iter - log((1/(double) nout-nout_0)*CPO[j]);

  }

  lpml[0] = lpml_iter;



  elppdWAIC = 0.0;

  for(j = 0; j < *nsubject; j++){
    elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
  }
  WAIC[0] = -2*elppdWAIC;

  PutRNGstate();

}
