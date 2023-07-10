/*************************************************************
 * Copyright (c) 2012 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit a hierarchical q-degree (b)-spline model with
 * a fixed number of knots penalized using P-splines.
 * We introduce covariates by using Product Partition model PPMx
 * that is developed by Fernando and co-authors.
 * The ideas of Fernando were used to implement the PPMx model
 *
 * Recall that for a single subject
 * the model is linear in a nonlinear Mls basis set
 * (see chapters 3,4 in book) :
 *
 *	Y_i = B_k * beta_i  + epsilon
 *	epsilon ~ N(0, sig2_i I)
 *
 *	Priors:
 *
 *	beta_i ~ N(theta_{S_i}, lam2 * I)
 *	theta_j \propto exp{-0.5/tau_jtheta_j K theta_j}
 *
 * B is a matrix whose entries corresond to a particular
 * basis dependent on covariates and K is a smoothing matrix
 * that is defined by the degree of random walk that is
 * employed.
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
* y = vector containing nobs production values for each player
* z = vector containing nobs time (time is incremented by one) values for each player

* K = smoothing matrix for penalized B-splines.
* dimK = dimension of K
* ncon = number of continuous covariates
* ncat = number of categorical covariates
* Cvec = ncat x 1 vector indicating the number of categories for each categorical covariate.

* Xcon = nsubject x ncon contiguous double vector that contains continuous covariate values
* Xcat = nsubject x ncat contiguous int vector that contains categorical covariate values


* PPM = logical indicating if PPM or PPMx should be used.
* M = double indicating value of M associated with cohesion (scale parameter of DP).

* similarity_function = int indicating similarity function of PPMx variable
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
* Hmat = matrix of b-splines basis functions.  Use this only if all subjects are measured
*			same time points.
* balanced = scalar indicating if Hmat above should be used as the design matrix.
* modelPriors = vector containing prior values associated with data hierarchical model
*
* Output:  The following objects contain MCMC iterates associated with named parmaeter
* beta =
* beta0 =
* sig2 =
* mub0 =
* sig2b0 =
* lam =
* tau2 =
* theta =
* mu =
* Si =
* nclus =
* like =
* WAIC =
* lpml =
* ispred =
* ppred =
* ppredclass =
*
*****************************************************************************************/


void mcmc_curvecluster(	int *draws, int *burn, int *thin, //3
                        int *nsubject, int *nobs, //2
                        double *y, double *z, //2
                        double *K, int *dimK, //2
                        int *ncon, int *ncat, int *Cvec, //3
                        double *Xcon, int *Xcat, //2
                        int *PPM, double *M, //2
                        int *similarity_function, int *consim, //3
                        int *npred, int *npredobs, double *Xconp, int *Xcatp, //4
                        double *simParms, double *Aparm, //2
                        int *calibrate, double *modelPriors, //2
                        double *Hmat, int *balanced, double *mh, //3
                        double *Hmat_pred, // 1
                        double *beta, double *beta0, double *sig2, double *mub0, //4
                        double *sig2b0,double *lam, double *tau2, double *theta, //4
                        double *mu, int *Si, int *nclus, //3
                        double *ppred, int *predclass, 
                        double *predlines, double *llike, //2
                        double *lpml, double *WAIC){ //2

  // i - MCMC iterate
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
//  int pp;
  int csobs = 0, csHobs=0, csobs_pred=0, csHobs_pred=0;
  int nb= (*dimK), N, N_pred=0;

  int nout = (*draws - *burn)/(*thin);

  Rprintf("nsubject = %d\n", *nsubject);
  Rprintf("nb = %d\n", nb);
  Rprintf("ncon = %d\n", *ncon);
  Rprintf("ncat = %d\n", *ncat);
  Rprintf("balanced = %d\n", *balanced);

  RprintIVecAsMat("nobs", nobs, 1, *nsubject);

  int  max_nobs;
  double *sumy = R_VectorInit(*nsubject,0.0);

  max_nobs = nobs[0];
  csobs = 0;
  for(j = 0; j < (*nsubject); j++){
    for(t = 0; t < nobs[j]; t++){
      sumy[j] = sumy[j] + y[csobs];
      csobs = csobs + 1;
    }
    if(max_nobs < nobs[j]) max_nobs = nobs[j];
    N_pred = N_pred + nobs[j] + (*npredobs);
  }


  N = csobs;
  Rprintf("N = %d\n", N);
  Rprintf("N_pred = %d\n", N_pred);
  Rprintf("max_nobs = %d\n", max_nobs);
  // RprintVecAsMat("K",K, nb, nb);
  // RprintVecAsMat("y", y, 1, N);
  // RprintVecAsMat("sumy", sumy, 1, *nsubject);

  int max_C;
  max_C = Cvec[0];
  for(p = 0; p < (*ncat); p++){
    if(max_C < Cvec[p]) max_C = Cvec[p];
  }

  if(*ncat == 0) max_C = 1.0;

  Rprintf("max_C = %d\n", max_C);


  // ===================================================================================
  //
  // Memory vectors to hold MCMC iterates for non cluster specific parameters
  //
  // ===================================================================================

  double *sig2_iter = R_VectorInit(*nsubject, ((modelPriors[0]-0.0)/2.0)*((modelPriors[0]-0.0)/2.0));
  double *beta0_iter = R_VectorInit(*nsubject, 0.0);

  double *beta_iter = R_VectorInit(nb*(*nsubject), 0.0);

  double *mu_iter = R_VectorInit(nb,0.0);

  double mub0_iter = 0.0;
  double sig2b0_iter = 1.0;


  int Si_iter[*nsubject];
  int nclus_iter = 0;

  // ===================================================================================
  //
  // Memory vectors to hold MCMC iterates for cluster specific parameters
  //
  // ===================================================================================

  double *tau2h = R_VectorInit(*nsubject, 1.0);
  double *lamh = R_VectorInit(*nsubject, 0.99*(*Aparm));
  double *thetah = R_VectorInit(nb*(*nsubject), 0.0);

  int nh[*nsubject];
//  RprintVecAsMat("sig2_iter", sig2_iter, 1, *nsubject);


  // ===================================================================================
  //
  // Initialize a few parameter vectors
  //
  // ===================================================================================

  // Initialize Si according to covariates
  for(j = 0; j < *nsubject; j++){
    //		Si_iter[j] = x[j]*(*ncat2) + x2[j]+1;
    //Si_iter[j] = 1;
    //Si_iter[j] = j+1;
    Si_iter[j] = rbinom(3,0.4)+1;
    nh[j] = 0;
  }
  // Initial enumeration of number of players per cluster;
  for(j = 0; j < *nsubject; j++){
    nh[Si_iter[j]-1] = nh[Si_iter[j]-1]+1;
  }
  // Initialize the number of clusters
  for(j = 0; j < *nsubject; j++){
    if(nh[j] > 0) nclus_iter = nclus_iter + 1;
  }
  RprintIVecAsMat("Si", Si_iter, 1, *nsubject);
  RprintIVecAsMat("nh", nh, 1, nclus_iter);
  Rprintf("nclus = %d\n", nclus_iter);
  
  // ===================================================================================
  //
  // scratch vectors of memory needed to update parameters
  //
  // ===================================================================================
  int  big = max_nobs;
  if(nb > max_nobs) big = nb;
//  Rprintf("big = %d\n", big);
  // These are made particularly big to make sure there is enough memory
  double *scr1 = R_Vector(big*big);
  double *scr2 = R_Vector(big*big);
  double *scr3 = R_Vector(big*big);


  // stuff that I need to update Si (cluster labels);
  int iaux=1, auxint;
  int nhctmp[max_C];
  double auxreal, uu, sumxtmp, sumx2tmp, xcontmp;
  double sumsq, maxph, denph, cprobh;
  double lamdraw, tau2draw;

  double *thetadraw = R_VectorInit(nb, 0.0);

  double *ph = R_VectorInit(*nsubject, 0.0);
  double *probh = R_VectorInit(*nsubject, 0.0);

  double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw;
  double lgcont,lgcatt;
  //	double mgconN, mgconY, mgcatN, mgcatY, sgconN, sgconY, sgcatY, sgcatN;

  double *gtilN = R_VectorInit(*nsubject+1,0.0);
  double *gtilY = R_VectorInit(*nsubject+1,0.0);
  double sgY, sgN,  lgtilN, lgtilY, maxgtilY, maxgtilN;

  double *sumx = R_VectorInit((*nsubject)*(*ncon),0.0);
  double *sumx2 = R_VectorInit((*nsubject)*(*ncon),0.0);
  int  nhc[(*nsubject)*(*ncat)*max_C];
  
  // stuff I need to update sig2 (player specific), mub0, sig2b0;
  double astar, bstar, os, ns, sumb0;

  // stuff that I need for player specific beta's;
  double *H = R_VectorInit((*nsubject)*(nb)*(big),0.0);
  double *tH = R_VectorInit((*nsubject)*(nb)*(big),0.0);
  double *HtH = R_Vector(nb*nb);
  double *Hty = R_Vector(nb);
//  double *bhat = R_Vector(nb);

  double *y_tmp = R_Vector(big);

  double sumy_Hb, s2star, mstar;
  double *y_b0 = R_Vector(big);
  double *Hb = R_Vector(big);

  double *H_pred = R_VectorInit((*nsubject)*(nb)*(big+(*npredobs)),0.0);
  double *Hb_pred = R_Vector(big*(*npredobs));


  // Create the inverse of the K penalty matrix
  double ld;
  double *Kinv = R_Vector(nb*nb);
  for(b = 0; b < nb; b++){for(bb = 0; bb < nb; bb++){Kinv[b*nb+bb] = K[b*nb+bb];}}
  cholesky(Kinv, nb, &ld);
  inverse_from_cholesky(Kinv, scr1, scr2, nb);

  // stuff I need for tau2h;
  double *thtmp = R_Vector(nb);

  // stuff I need for thetah
  double *sumbeta = R_Vector(nb);
  double *Mstar = R_Vector(nb);
  double *Sstar = R_Vector(nb*nb);
  double *outrmvnorm = R_Vector(nb);

  //stuff I need for lamh
  double olam, nlam, lln, llo, ldo, ldn, llr;
  double *btmp = R_Vector(nb);
  double *nV = R_Vector(nb*nb);
  double *oV = R_Vector(nb*nb);


  // stuff I need for mu
  double *sumtheta = R_Vector(nb);
  double sumtau2;

  // stuff that I need to perform the predictions
  // I am not including prediction at this time. (0.3.1)
  // double *thetatmp = R_VectorInit(nb, 0.0);
  // double *bpred = R_VectorInit(nb, 0.0);
  // double *ppredtmp = R_VectorInit(100,0.0);
  // double lgcon0, lgcat0=0.0;

  // stuff that I need to compute the lpml;
  double llikeval, lpml_iter = 0.0, elppdWAIC;;
  double *CPO = R_VectorInit((*nsubject), 0.0);
  double *llike_iter = R_VectorInit(*nsubject,0.0);
  double *mnlike = R_VectorInit((*nsubject), 0.0);
  double *mnllike = R_VectorInit((*nsubject), 0.0);
  int like0;

  //	Rprintf("npredobs = %d\n", *npredobs);


  if(*balanced==1){
    for(t=0; t < nobs[0]; t++){
      for(kk=0; kk<nb; kk++){
        H[t*(nb) + kk] = Hmat[t*(nb) + kk];
      }
    }  
    for(t=0; t < nobs[0] + (*npredobs); t++){
      for(kk=0; kk<nb; kk++){
        H_pred[t*(nb) + kk] = Hmat_pred[t*(nb) + kk];
      }

    }
    
    mat_transpose(H, tH, nobs[0], nb);
    matrix_product(tH, H, HtH, nb, nb, nobs[0]);

  }


  double *mnmle = R_VectorInit(*ncon,0);
  double *s2mle = R_VectorInit(*ncon,0);
  double sum, sum2;
  if(!(*PPM)){
    for(p = 0; p < *ncon; p++){
      sum = 0.0, sum2=0.0;
      for(j = 0; j < *nsubject; j ++){
        sum = sum + Xcon[j*(*ncon) + p];
        sum2 = sum2 + Xcon[j*(*ncon) + p]*Xcon[j*(*ncon) + p];
      }

      mnmle[p] = sum/((double) *nsubject);
      s2mle[p] = sum2/((double) *nsubject) - mnmle[p]*mnmle[p];
    }
  }


  // ===================================================================================
  //
  // Prior parameter values
  //
  // ===================================================================================

  // prior values for sig2
  //	double asig=100, bsig=1;
  //	double asig=10000.0, bsig=1.0/100.0;
  double Asig = modelPriors[0]; // this is for Metropolis step;


  // priors for mu
  double s2mu = modelPriors[1];

  //	for(j = 0; j<nb; j++) m_mu[j]=bhat[j];

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
  double A = *Aparm;
  Rprintf("A = %f\n", A);

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


//  RprintVecAsMat("sumx = ", sumx, nclus_iter, *ncon);


  // Candidate density sd's for M-H stuff
  double csigLAM = mh[0], csigSIG=mh[1];

  ii = 0;

  GetRNGstate();



  // ===================================================================================
  //
  // start of the mcmc algorithm;
  //
  // ===================================================================================

  for(i = 0; i < *draws; i++){


    if((i+1) % 1000 == 0){
      time_t now;
      time(&now);

      Rprintf("mcmc iter = %d =========================================== \n", i+1);
      Rprintf("number of clusters = %d\n", nclus_iter);
//      Rprintf("%s", ctime(&now));
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
        
        
//        Rprintf("Si_iter[j] = %d\n", Si_iter[j]);
        
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
//          Rprintf("Si_iter[j] = %d\n", Si_iter[j]);

          // The following steps swaps order of cluster specific parameters
          // so that the newly labeled subjects from previous step retain
          // their correct cluster specific parameters
          auxreal = tau2h[iaux-1];
          tau2h[iaux-1] = tau2h[nclus_iter-1];
          tau2h[nclus_iter-1] = auxreal;

          auxreal = lamh[iaux-1];
          lamh[iaux-1] = lamh[nclus_iter-1];
          lamh[nclus_iter-1] = auxreal;

          for(b = 0; b < nb; b++){

            auxreal = thetah[b*(*nsubject) + iaux-1];
            thetah[b*(*nsubject) + iaux-1] = thetah[b*(*nsubject) + nclus_iter-1];
            thetah[b*(*nsubject) + nclus_iter-1] = auxreal;
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
          
//        Rprintf("nclus_iter = %d\n", nclus_iter);
//        Rprintf("Si_iter = %d\n", Si_iter[j]);
//        Rprintf("sumx1 = %f\n", sumx[(nclus_iter-1)*(*ncon) + p]);
//        Rprintf("sumx2 = %f\n", sumx[(Si_iter[j]-1)*(*ncon) + p]);
//        Rprintf("Xcon[j*(*ncon)+p] = %f\n", Xcon[j*(*ncon)+p]);

        // need to reduce the sumx sumx2
        if(!(*PPM)){
          for(p = 0; p < *ncon; p++){
            sumx[(nclus_iter-1)*(*ncon) + p] = sumx[(nclus_iter-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
            sumx2[(nclus_iter-1)*(*ncon) + p] = sumx2[(nclus_iter-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          }

          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){

             nhc[(((nclus_iter-1))*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
                nhc[(((nclus_iter-1))*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;

           }
        }
        
//        Rprintf("sumx1 = %f\n", sumx[(nclus_iter-1)*(*ncon) + p]);
        
        
        // Finally, reduce the number of clusters                
        nclus_iter = nclus_iter - 1;

//        RprintVecAsMat("sumx = ", sumx, nclus_iter, *ncon);

      }


      // The atoms have been relabeled if necessary and now we need to
      // update Si.

      for(b = 0; b < nb; b++){

        btmp[b] = beta_iter[b*(*nsubject) + j];

      }

      // RprintVecAsMat("btmp", btmp, 1, nb);
//      Rprintf("nclus_iter = %d\n", nclus_iter);


      ///////////////////////////////////
      //
      // Begin the cluster probabilities
      //
      //////////////////////////////////
      for(k = 0 ; k < nclus_iter; k++){
//        Rprintf("k = %d\n", k);

        for(b = 0; b < nb; b++){
          thtmp[b] = thetah[b*(*nsubject) + k];
          for(bb = 0; bb < nb; bb++){
            oV[b*nb+bb] = 0.0;
            if(b == bb){
              oV[b*nb+bb] = 1/(lamh[k]*lamh[k]);
            }
          }
        }


        ldo = 2.0*nb*log(lamh[k]);

        lgconY = 0.0;
        lgconN = 0.0;
        lgcatY=0.0;
        lgcatN=0.0;
        
        
       
        if(!(*PPM)){
        
          for(p=0; p<(*ncon); p++){
  
            sumxtmp = sumx[k*(*ncon) + p];
            sumx2tmp = sumx2[k*(*ncon) + p];

            if(*similarity_function==1){ // Auxilliary
              if(*consim==1){ //NN
                lgcont = gsimconNN(m0, v2, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 0, 0, 1);
                lgconN = lgconN + lgcont;
              }
              if(*consim==2){//NNIG
                lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k], 0, 0, 1);
                lgconN = lgconN + lgcont;
              }
  
            }
  
            if(*similarity_function==2){ //Double Dipper
              if(*consim==1){ //NN
                lgcont = gsimconNN(m0, v2, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 1, 0, 1);
                lgconN = lgconN + lgcont;
              }
              if(*consim==2){ //NNIG
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

//            Rprintf("sumxtmp = %f\n", sumxtmp);
//            Rprintf("sumx2tmp = %f\n", sumx2tmp);
  
 
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
              lgcont = gsimconEV(sumxtmp, sumx2tmp, nh[k]+1, alpha,1);
              lgconY = lgconY + lgcont;
            }
  
          }
  
 
      	  // Now calculate similarity for the categorical covariates
          for(p=0; p<(*ncat); p++){
//              Rprintf("p = %d\n", p);
//            RprintVecAsMat("dirweights", dirweights, 1, max_C);

            for(c=0;c<Cvec[p];c++){
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
//            Rprintf("lgcatN = %f\n", lgcatN);

            // include the categorical covariate in the kth cluster
      	  	nhctmp[Xcat[j*(*ncat)+p]] = nhctmp[Xcat[j*(*ncat)+p]] + 1;
  
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
//        Rprintf("lgconY = %f\n", lgconY);
//        Rprintf("lgconN = %f\n", lgconN);
//        Rprintf("lgcatY = %f\n", lgcatY);
//        Rprintf("lgcatN = %f\n", lgcatN);
//        Rprintf("nh[k] = %d\n", nh[k]);
//        Rprintf("dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) = %f\n", dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1));

          ph[k] = dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) +
          log((double) nh[k]) + // Scale parameter from DP
          lgcatY - lgcatN +  // Categorical part only nonzero if PPMx=TRUE
          lgconY - lgconN;   // Continuous part only nonzero if PPMx=TRUE

      }


      // Need to consider allocating subject to new cluster
      lamdraw = runif(0,A);
      tau2draw = 1/rgamma(at,bt); // shape and scale  E(tau2) = atbt


      for(b = 0; b < nb; b++){
        for(bb = 0; bb < nb; bb++){
          oV[b*nb+bb] = 0.0;
          if(b == bb){
            oV[b*nb+bb] = 1/(lamdraw*lamdraw);
          }
          nV[b*nb+bb] = tau2draw*Kinv[b*nb+bb];
        }
      }

      cholesky(nV, nb , &ld);

      ran_mvnorm(mu_iter, nV, nb, scr1, thetadraw);

      ldo = 2.0*nb*log(lamdraw);
      
//      RprintVecAsMat("thetadraw", thetadraw, 1, nb);

      lgcondraw = 0.0;
      lgcatdraw = 0.0;
      if(!(*PPM)){

        // similarity for continuous covariate
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


      ph[nclus_iter] = dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1) +
                           log(Mdp) +
                           lgcondraw +
                           lgcatdraw;



//      RprintVecAsMat("ph", ph, 1, nclus_iter+1);



      /////////////////////////////////////////////////////////////////////////////
      // This is the calibration used when the similarity is standardized by
      // the sum of all cluster similarity values.
      /////////////////////////////////////////////////////////////////////////////
      if(*calibrate==1){

        for(b = 0; b < nb; b++){
          thtmp[b] = thetah[b*(*nsubject) + k];
          for(bb = 0; bb < nb; bb++){
            oV[b*nb+bb] = 0.0;
            if(b == bb){
              oV[b*nb+bb] = 1/(lamh[k]*lamh[k]);
            }
          }
        }

        ldo = 2.0*nb*log(lamh[k]);

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

          ph[k] = dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) +
                     log((double) nh[k]) +
                     lgtilY - lgtilN; //This takes into account both cont and cat vars

        }

      	// calibration for a singleton
        ph[nclus_iter] = dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1) +
                           log(Mdp) +
                           log(gtilN[nclus_iter]) - log(sgN);

      }


      /////////////////////////////////////////////////////////////////////////////
      // End of calibration used when the similarity is normalized
      /////////////////////////////////////////////////////////////////////////////

      maxph = ph[0];
      for(k = 1; k < nclus_iter+1; k++){
      	if(maxph < ph[k]) maxph = ph[k];
      }

      // Rprintf("maxph = %f\n", maxph);

      denph = 0.0;
      for(k = 0; k < nclus_iter+1; k++){
        ph[k] = exp(ph[k] - maxph);
        denph = denph + ph[k];
      }

      for(k = 0; k < nclus_iter+1; k++){
        probh[k] = ph[k]/denph;
      }

//      RprintVecAsMat("probh = ", probh, 1, nclus_iter+1);

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

        tau2h[Si_iter[j]-1] = tau2draw;
        lamh[Si_iter[j]-1] = lamdraw;

        for(b = 0; b < nb; b++){
          thetah[b*(*nsubject) + Si_iter[j]-1] = thetadraw[b];
        }
      }

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
      
//      RprintIVecAsMat("Si", Si_iter, 1, *nsubject);
//      Rprintf("nclus = %d\n", nclus_iter);
//      RprintVecAsMat("sumx", sumx, *ncon, nclus_iter);


    }




    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // update player specific betas, sigma2, and beta0;								//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    csobs=0;
    csobs_pred=0;
    csHobs=0;
    csHobs_pred=0;
    for(j = 0; j < *nsubject; j++){
      //Rprintf("j = %d ====================\n", j);

      for(t = 0; t < nobs[j]; t++){

        y_tmp[t] = y[csobs];
        y_b0[t] = y_tmp[t]-beta0_iter[j];
        csobs = csobs+1;

      }


      if(*balanced!=1){
        for(t=0; t < nobs[j]*(nb); t++){
		  H[t] = Hmat[csHobs + t];
		}
		csHobs = csHobs + nobs[j]*(nb);

        mat_transpose(H, tH, nobs[j], nb);
        matrix_product(tH, H, HtH, nb, nb, nobs[j]);
        
        for(t=0; t < (nobs[j]+(*npredobs))*(nb); t++){
		  H_pred[t] = Hmat_pred[csHobs_pred + t];
		}
		csHobs_pred = csHobs_pred + (nobs[j]+(*npredobs))*(nb);
        
	  }

//      Rprintf("csHobs = %d\n", csHobs);
//      RprintVecAsMat("H", H, nobs[j], nb);

      matrix_product(tH, y_tmp, Hty, nb, 1, nobs[j]);

//      RprintVecAsMat("HtH", HtH, nb, nb);
//      RprintVecAsMat("H'y", Hty, 1, nb);
//      Rprintf("Si = %d\n", Si_iter[j]);
//      Rprintf("lamh[Si] = %f\n", lamh[Si_iter[j]-1]);

      for(b = 0; b < nb; b++){
        for(bb = 0; bb < nb; bb++){

          Sstar[b*nb+bb] = (1/sig2_iter[j])*HtH[b*nb+bb];

          if(b == bb){
            Sstar[b*nb+bb] = (1/sig2_iter[j])*HtH[b*nb+bb] +
              (1/(lamh[Si_iter[j]-1]*lamh[Si_iter[j]-1]));
          }
        }

      }

      cholesky(Sstar, nb, &ld);
      inverse_from_cholesky(Sstar, scr1, scr2, nb);

//      RprintVecAsMat("Sstar", Sstar, nb, nb);
//      RprintVecAsMat("y_b0", y_b0, 1, nobs[j]);
      matrix_product(tH, y_b0, scr2,  nb, 1, nobs[j]);


//      RprintVecAsMat("H'(y-b0)", scr2, 1, nb);

      for(b = 0; b < nb; b++){
        scr3[b] = (1/sig2_iter[j])*scr2[b] +
          (1/(lamh[Si_iter[j]-1]*lamh[Si_iter[j]-1]))*
          thetah[b*(*nsubject) + Si_iter[j]-1];

      }

      matrix_product(Sstar, scr3, Mstar, nb, 1, nb);

//      RprintVecAsMat("Mstar", Mstar, 1, nb);
//      RprintVecAsMat("Sstar", Sstar, nb, nb);

      cholesky(Sstar, nb , &ld);

//      RprintVecAsMat("Mstar", Mstar, 1, nb);

      ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);

//      RprintVecAsMat("betai", outrmvnorm, 1, nb);

       //cholesky(HtH, nb, &ld);
       //inverse_from_cholesky(HtH, scr1, scr2, nb);
       //matrix_product(HtH, Hty, bhat, nb, 1, nb);
//      RprintVecAsMat("bhat", bhat, 1, nb);


//      Rprintf("nb = %d\n", nb);
//      Rprintf("nsubject = %d\n", *nsubject);
      for(b = 0; b < nb; b++){
        btmp[b] = outrmvnorm[b];
        beta_iter[b*(*nsubject) + j] = btmp[b];
      }
   
      matrix_product(H, btmp, Hb, nobs[j], 1, nb);
//      RprintVecAsMat("Hb", Hb, 1, nobs[j]);
      

      sumy_Hb = 0.0;
      for(jj = 0; jj < nobs[j]; jj++){
        sumy_Hb = sumy_Hb + (y_tmp[jj] - Hb[jj]);
      }


      // write the predlines to file straight away so that I don't have to store them 
      // in a matrix
      // out-of-sample predictions for each subject
      matrix_product(H_pred, btmp, Hb_pred, nobs[j]+(*npredobs), 1, nb);
//      RprintVecAsMat("Hb_pred", Hb_pred, 1, nobs[j]+(*npredobs));
      

      if((i >= (*burn)) & (i % (*thin) == 0)){    
       
        for(jj = 0; jj < nobs[j]+(*npredobs); jj++){
          predlines[ii*(N_pred) + csobs_pred] = rnorm(Hb_pred[jj], sqrt(sig2_iter[j]));
          csobs_pred = csobs_pred + 1;
        }
      }

      /////////////////////////////////////////
      //									   //
      // udate sigma2 within the same loop.  //
      //									   //
      /////////////////////////////////////////
      // for(jj = 0; jj < nobs[j]; jj++){
      // 	scr3[jj] = y_b0[jj] - Hb[jj];
      // 	sumy_Hb = sumy_Hb + (y_tmp[jj] - Hb[jj]);
      // }

      //			sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

      //			astar = 0.5*nobs[j] + asig;
      //			bstar = 0.5*sumsq + 1/bsig;


      // these are for the hierarchical variance structure
      //			astar = 0.5*nobs[j] + 0.5*nuh[Si_iter[j]];
      //			bstar = 0.5*sumsq + 0.5;


      //			Rprintf("astar = %f\n", astar);
      //			Rprintf("bstar = %f\n", bstar);
      //bstar is rate and rgamma requires scale hence inverse
      //			sig2_iter[j] = 1/rgamma(astar, 1/bstar);
      //			sig2_iter[j] = 0.000000000000001;

      //			Rprintf("sig2 = %f\n", sig2_iter[j]);


      os = sqrt(sig2_iter[j]);
      ns = rnorm(os, csigSIG);

      if(ns > 0.0){

        for(jj = 0; jj < nobs[j]; jj++){
          scr3[jj] = y_b0[jj] - Hb[jj];
        }

        sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

        llo = -0.5*nobs[j]*log(os*os) - 0.5*(1/(os*os))*sumsq;
        lln = -0.5*nobs[j]*log(ns*ns) - 0.5*(1/(ns*ns))*sumsq;

        llo = llo + dunif(os, 0.0, Asig, 1);
        lln = lln + dunif(ns, 0.0, Asig, 1);

        llr = lln - llo;

        uu = runif(0.0,1.0);
        if(llr > log(uu)) sig2_iter[j] = ns*ns;

      }

      ////////////////////////////////////////////
      //									    //
      // update beta0 within in the same loop;  //
      //									    //
      ////////////////////////////////////////////

      s2star = 1.0/((nobs[j]/sig2_iter[j]) + (1/sig2b0_iter));

      mstar = s2star*((1/sig2_iter[j])*sumy_Hb + (1/sig2b0_iter)*mub0_iter);

      beta0_iter[j] = rnorm(mstar, sqrt(s2star));


      if((i > (*burn-1)) & (i % (*thin) == 0)){
        like0=0;
        llikeval = 0.0;
        for(t = 0; t < nobs[j]; t++){
          llikeval = llikeval + dnorm(y_tmp[t], beta0_iter[j] + Hb[t], sqrt(sig2_iter[j]), 1);
        }

        llike_iter[j] = llikeval;

        // These are needed for WAIC
        mnlike[j] = mnlike[j] + exp(llike_iter[j])/(double) nout;
        mnllike[j] = mnllike[j] + (llike_iter[j])/(double) nout;

        if(like0==0){
          CPO[j] = CPO[j] + (1/(double) nout)*(1/exp(llike_iter[j]));
        }
      }
    }


    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update lam (using MH-step) this is the standard deviation not the variance
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      //			Rprintf("k = %d =================== \n", k);

      olam = lamh[k];
      nlam = rnorm(olam,csigLAM);


      if((nlam > 0) & (nlam < A)){

        for(b = 0; b < nb; b++){

          thtmp[b] = thetah[b*(*nsubject) + k];

          for(bb = 0; bb < nb; bb++){

            oV[b*nb+bb] = 0.0;
            nV[b*nb+bb] = 0.0;

            if(b == bb){

              oV[b*nb+bb] = 1/(olam*olam);
              nV[b*nb+bb] = 1/(nlam*nlam);
            }
          }
        }

        ldo = 2.0*nb*log(olam);
        ldn = 2.0*nb*log(nlam);

        lln = 0.0;
        llo = 0.0;
        for(j = 0; j < *nsubject; j++){

          if(Si_iter[j] == k+1){

            for(b = 0; b < nb; b++){

              btmp[b] = beta_iter[b*(*nsubject) + j];

            }

            llo = llo + dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1);
            lln = lln + dmvnorm(btmp, thtmp, nV, nb, ldn, scr1, 1);

          }

        }


        llo = llo + dunif(olam, 0.0, A, 1);
        lln = lln + dunif(nlam, 0.0, A, 1);

        llr = lln - llo;
        uu = runif(0.0,1.0);

        if(log(uu) < llr) lamh[k] = nlam;

      }

    }

//    RprintVecAsMat("lam", lamh, 1, nclus_iter);


    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate mub0 a global intercept of mean of beta0i     						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    sumb0 = 0.0;
    for(j = 0; j < *nsubject; j++){

      sumb0 = sumb0 + beta0_iter[j];

    }

    s2star = 1.0/((*nsubject)/sig2b0_iter + (1/s2b0));

    mstar = s2star*((1/sig2b0_iter)*(sumb0) + (1/s2b0)*mb0);

    mub0_iter = rnorm(mstar, sqrt(s2star));


    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate sig2b0 a global intercept of mean of beta0i     						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    sumsq = 0.0;
    for(j = 0; j < *nsubject; j++){

      sumsq = sumsq + (beta0_iter[j] - mub0_iter)*(beta0_iter[j] - mub0_iter);

    }

    astar = 0.5*(*nsubject) + ab0;
    bstar = 0.5*sumsq + 1/bb0;

    sig2b0_iter = 1/rgamma(astar, 1/bstar);


    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate thetah each of the cluster specific coefficients;						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      for(b = 0; b < nb; b++){
        for(bb = 0; bb < nb; bb++){

          Sstar[b*nb+bb] = (1/tau2h[k])*K[b*nb+bb];

          if(b == bb) Sstar[b*nb+bb] = ((double) nh[k]/(lamh[k]*lamh[k])) + (1/tau2h[k])*K[b*nb+bb];

        }

        sumbeta[b] = 0.0;

      }


      cholesky(Sstar, nb, &ld);
      inverse_from_cholesky(Sstar, scr1, scr2, nb);

      matrix_product(K, mu_iter, scr1, nb, 1, nb);

      for(j = 0; j < *nsubject; j++){
        if(Si_iter[j] == k+1){
          for(b = 0; b < nb; b++){
            sumbeta[b] = sumbeta[b] + (1/(lamh[k]*lamh[k]))*
                                       beta_iter[b*(*nsubject) + j];
          }
        }
      }

      for(b=0; b<nb; b++) sumbeta[b] = sumbeta[b] + (1/tau2h[k])*scr1[b];


      matrix_product(Sstar, sumbeta, Mstar, nb, 1, nb);


      cholesky(Sstar, nb , &ld);

      ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);
//      RprintVecAsMat("theta", outrmvnorm, 1, nb);

      for(b = 0; b < nb; b++){

        thetah[b*(*nsubject) + k] = outrmvnorm[b];

      }
    }

    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update tau2 for each of the clusters (P-spline smoothing parameter)
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      for(b = 0; b < nb; b++){

        thtmp[b] = thetah[b*(*nsubject) + k] - mu_iter[b];

      }

      sumsq = quform(thtmp,K,nb);

      astar = 0.5*(nb) + at;
      bstar = 1/bt + 0.5*sumsq;

      tau2h[k] = 1/rgamma(astar, 1/bstar);// E(tau2) = astarbstar for gamma.  bstar is scale

    }

    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate mu center of cluster specific thetah;									//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////

    for(b = 0; b < nb; b++) sumtheta[b] = 0.0;
    sumtau2 = 0.0;

    for(k = 0; k < nclus_iter; k++){
      sumtau2 = sumtau2 + (1/tau2h[k]);
      for(b = 0; b < nb; b++){
        sumtheta[b] = sumtheta[b] + (1/tau2h[k])*thetah[b*(*nsubject) + k];
      }
    }

    for(b = 0; b < nb; b++){
      for(bb = 0; bb < nb; bb++){

        Sstar[b*nb+bb] = sumtau2*K[b*nb+bb];

        if(b == bb) Sstar[b*nb+bb] = sumtau2*K[b*nb+bb]+(1/s2mu) ;

      }

    }

    matrix_product(K, sumtheta, scr3, nb, 1, nb);

    cholesky(Sstar, nb, &ld);
    inverse_from_cholesky(Sstar, scr1, scr2, nb);

    matrix_product(Sstar, scr3, Mstar, nb, 1, nb);

    cholesky(Sstar, nb , &ld);

    ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);

//    RprintVecAsMat("mu", outrmvnorm, 1, nb);
    for(b = 0; b < nb; b++){
       mu_iter[b] = outrmvnorm[b];
    }




/*

    //////////////////////////////////////////////////////////////////////////////////
    //
    // Posterior predictives for npred new subjects
    //
    //////////////////////////////////////////////////////////////////////////////////

    if(i > (*burn-1) & i % (*thin) == 0){

      for(pp = 0; pp < *npred; pp++){

        //				Rprintf("pp = %d ================================================ \n", pp+1);



        //				Rprintf("nclus_iter = %d\n", nclus_iter);
        for(k = 0; k < nclus_iter; k++){
          //					Rprintf("k = %d  ========== \n", k);


          lgconN=0.0, lgconY=0.0;
          for(p=0; p<(*ncon); p++){
            //						Rprintf("p = %d ====== \n", p) ;
            nhtmp = 0;
            for(j = 0; j < nobs[j]; j++){
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
            if(*gcontype==1){
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
            if(*gcontype==2){
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
            if(*gcontype==3){
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
            if(*gcontype==1){ // Auxilliary
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
            if(*gcontype==2){ // Double Dipper
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
            if(*gcontype==3){ // Variance
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
            for(j = 0; j < nobs[j]; j++){
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

            if(*gcattype==1){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
              lgcatN = lgcatN + lgcatt;
            }
            if(*gcattype==2){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
              lgcatN = lgcatN + lgcatt;
            }
            if(*gcattype==3){
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

            if(*gcattype==1){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
              lgcatY = lgcatY + lgcatt;
            }
            if(*gcattype==2){
              lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
              lgcatY = lgcatY + lgcatt;
            }
            if(*gcattype==3){// Use entropy
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
          if(*gcontype==1){
            if(*consim==1){
              lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0,1);
              //							lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,1,1);
            }
            if(*consim==2){
              lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,0,1);
              //							lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,1,1);
            }
          }
          if(*gcontype==2){
            if(*consim==1){
              lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,0,1);
              //							lgcon0 = lgcon0 + gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,1,1);
            }
            if(*consim==2){
              lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 0,1);
              //							lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 1,1);
            }
          }
          if(*gcontype==3){
            lgcon0 = lgcon0 + gsimconEV(xcontmp[0],xcontmp[0]*xcontmp[0],1,alpha,1);
          }
        }



        lgcat0 = 0.0;

        for(p=0;p<(*ncat);p++){
          for(c=0;c<Cvec[p];c++) nhc[c] = 0;

          nhc[Xcatp[pp*(*ncat)+p]] = 1;
          //					RprintIVecAsMat("nhc =", nhc, 1, Cvec[p]);

          if(*gcattype==1){
            lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
          }
          if(*gcattype==2){
            lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
          }
          if(*gcattype==3){
            //						lgcat0 = lgcat0 + gsimconEV(xcattmp[0], xcattmp[0]*xcattmp[0],  1,alpha, 1);
            lgcat0 = lgcat0 + -(alpha)*0;

          }


        }


        gtilY[nclus_iter] = lgcat0 + lgcon0;
        gtilN[nclus_iter] = lgcat0 + lgcon0;

        ph[nclus_iter] = log((double) Mdp) + lgcon0 + lgcat0;

        if(*gcontype==4) ph[nclus_iter] = log((double) Mdp) + log(1);


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

            thetatmp[b] = thetah[b*(*nsubject) + iaux-1];

          }


          for(b = 0; b < nb; b++){
            for(bb = 0; bb < nb; bb++){

              Sstar[b*nb+bb] = 0.0;

              // THis is the Cholesky Decomposition as needed in ran_mvnorm
              if(b == bb) Sstar[b*nb+bb] = sqrt(lamh[iaux-1]*lamh[iaux-1]);

            }
          }

          //					RprintVecAsMat("thetatmp = ", thetatmp, 1, nb);
        }else{

          tau2draw = 1/rgamma(at,bt);
          //					Rprintf("tau2draw = %f\n", tau2draw);
          lamdraw = runif(0,A);
          //					Rprintf("lamdraw = %f\n", lamdraw);

          for(b = 0; b < nb; b++){
            for(bb = 0; bb < nb; bb++){

              nV[b*nb+bb] = tau2draw*Kinv[b*nb+bb];

              Sstar[b*nb+bb] = 0.0;
              // THis is the Cholesky Decomposition as needed in ran_mvnorm
              if(b == bb) Sstar[b*nb+bb] = sqrt(lamdraw*lamdraw);

            }


          }


          //					RprintVecAsMat("nV = ", nV, nb, nb);
          cholesky(nV, nb , &ld);

          //					RprintVecAsMat("mu_iter =", mu_iter, 1, nb);
          ran_mvnorm(mu_iter, nV, nb, scr1, thetatmp);




        }


        //				RprintVecAsMat("thetatmp = ", thetatmp, 1, nb);
        //				RprintVecAsMat("Sstar = ", Sstar, nb, nb);


        ran_mvnorm(thetatmp, Sstar, nb, scr1, bpred);


        //				RprintVecAsMat("bpred = ", bpred, 1, nb);

        matrix_product(Hpred, bpred, ppredtmp, *npredobs, 1, nb);

        //				RprintVecAsMat("ppredtmp", ppredtmp, 1, *npredobs);

      }
    }

*/


    //////////////////////////////////////////////////////////////////////////////////
    //																				                                      //
    // Save MCMC iterates															                              //
    //																			  	                                    //
    //////////////////////////////////////////////////////////////////////////////////

    if((i >= (*burn)) & (i % (*thin) == 0)){
//      Rprintf("i = %d\n", i+1);
//      Rprintf("ii = %d\n", ii);


      nclus[ii] = nclus_iter;

      mub0[ii] = mub0_iter;

      sig2b0[ii] = sig2b0_iter;


      for(b = 0; b < nb; b++){
        mu[ii*(nb) + b] = mu_iter[b];
      }



      for(j = 0; j < *nsubject; j++){

        Si[ii*(*nsubject) + j] = Si_iter[j];


        sig2[ii*(*nsubject) + j] = sig2_iter[j];

        beta0[ii*(*nsubject) + j] = beta0_iter[j];

        llike[ii*(*nsubject) + j] = llike_iter[j];

        lam[ii*(*nsubject) + j] = lamh[Si_iter[j]-1];

        tau2[ii*(*nsubject) + j] = tau2h[Si_iter[j]-1];


        for(b = 0; b < nb; b++){
//          Rprintf("beta_iter = %f\n", beta_iter[b*(*nsubject) + j]);
          beta[(ii*(nb) + b)*(*nsubject) + j] = beta_iter[b*(*nsubject) + j];

          theta[(ii*(nb) + b)*(*nsubject) + j] = thetah[b*(*nsubject) + Si_iter[j]-1];

        }

      }

      ii = ii+1;


    }

  }


  lpml_iter=0.0;
  for(j = 0; j < *nsubject; j++){
    lpml_iter = lpml_iter + log(1/CPO[j]);
  }

  lpml[0] = lpml_iter;

  elppdWAIC = 0.0;
  for(j = 0; j < *nsubject; j++){
    elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
  }
  WAIC[0] = -2*elppdWAIC;

  PutRNGstate();

}
