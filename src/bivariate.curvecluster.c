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

* q = the degree of the b-spline or truncated polynomial spline
* nknots = number of values in the preceding vector
* knots = vector containing knot values (typically a uniform array of length 20-30)
*			this vector is the same for all subjects
* K = smoothing matrix for penalized B-splines.

* ncon = number of continuous covariates
* ncat = number of categorical covariates
* Cvec = ncat x 1 vector indicating the number of categories for each categorical covariate.

* Xcon = nsubject x ncon contiguous double vector that contains continuous covariate values
* Xcat = nsubject x ncat contiguous int vector that contains categorical covariate values


* PPM = logical indicating if PPM or PPMx should be used.
* M = double indicating value of M associated with cohesion (scale parameter of DP).

* gcontype = int indicating similarity function for continuous variable
* gcattype = int indicating similarity function for categorical variable
	similarity - an integer indicating which similarity function to use
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


void mcmc_bivariate_curvecluster(int *draws, int *burn, int *thin,
                                 int *nsubject, int *nobs,
                                 double *y1, double *y2, double *z,
                                 int *q, int *nknots, double *knots, double *K,
                                 int *ncon, int *ncat,int *Cvec,
                                 double *Xcon, int *Xcat,
                                 int *PPM, double *M,
                                 int *gcontype, int *gcattype, int *consim,
                                 int *npred, int *npredobs, double *Xconp, int *Xcatp,
                                 double *simParms, double *Aparm,
                                 int *calibrate, double *modelPriors, double *Hmat,
                                 int *balanced, double *mh,
                                 double *beta1, double *beta2, double *beta01, double *beta02,
                                 double *sig21, double *sig22, double *mub01, double *mub02,
                                 double *sig2b01,double *sig2b02,double *lam1, double *lam2,
                                 double *tau21, double *tau22, double *theta1, double *theta2,
                                 double *mu1,double *mu2, int *Si, int *nclus,
                                 double *ppred1,  double *ppred2, int *predclass,
                                 double *llike, double *lpml,double *WAIC){




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

  int c, i, ii, j, t, p, jj, k, kk, pp, b, bb;
  int csobs = 0;
  int nb, N ;

  //	nb = *nknots + *q + 1;  //Number of B-spline basis after filling in all players;
  if(*balanced==0) nb = *nknots + *q + 1;  //Number of B-spline basis my c function;
  if(*balanced==1) nb = *nknots + *q;  //Number of B-spline basis cote;
  int nout = (*draws - *burn)/(*thin);

  Rprintf("nsubject = %d\n", *nsubject);
  Rprintf("nknots = %d\n", *nknots);
  Rprintf("nb = %d\n", nb);
  Rprintf("ncon = %d\n", *ncon);
  Rprintf("ncat = %d\n", *ncat);

//  RprintIVecAsMat("nobs", nobs, 1, *nsubject);
 // RprintIVecAsMat("Cvec", Cvec, 1, *ncat);
  //	RprintVecAsMat("knots", knots, *nknots, *nsubject);


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

  Rprintf("max_C = %d\n", max_C);





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
    Si_iter[j] = 1;
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


  Rprintf("nclus_iter = %d\n", nclus_iter);
  //	RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
  RprintIVecAsMat("nh", nh, 1, nclus_iter);


  // ===================================================================================
  //
  // Initialize arrays of memory needed to update cluster assignments
  //
  // ===================================================================================

  int nhc[(*nobs)*(*ncat)];

  double* xcontmp = R_Vector(*nsubject);



  // ===================================================================================
  //
  // scratch vectors of memory needed to update parameters
  //
  // ===================================================================================

  // These are made particularly big to make sure there is enough memory
  double *scr1 = R_Vector((max_nobs)*(max_nobs));
  double *scr2 = R_Vector((max_nobs)*(max_nobs));
  double *scr3 = R_Vector((max_nobs)*(max_nobs));
  double *scr4 = R_Vector((max_nobs)*(max_nobs));


  // stuff that I need to update Si (cluster labels);
  int iaux, nhtmp;

  double auxt1, auxl1, auxth1, uu;
  double auxt2, auxl2, auxth2;
  double sumsq2, sumsq1, maxph, denph, cprobh;
  double lamdraw1, tau2draw1, lamdraw2, tau2draw2;
  double sumx, sumx2, tmp;

  double *thetadraw1 = R_VectorInit(nb, 0.0);
  double *thetadraw2 = R_VectorInit(nb, 0.0);

  double *ph = R_VectorInit(*nsubject, 0.0);
  double *phtmp = R_VectorInit(*nsubject, 0.0);
  double *probh = R_VectorInit(*nsubject, 0.0);

  double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw;
  double lgcont,lgcatt;
  //	double mgconN, mgconY, mgcatN, mgcatY, sgconN, sgconY, sgcatY, sgcatN;

  double *gtilN = R_VectorInit((*nobs+1),0.0);
  double *gtilY = R_VectorInit((*nobs+1),0.0);
//  double sgY, sgN, lgtilN, lgtilY, maxgtilY, maxgtilN;


  // stuff I need to update sig2 (subject specific), mub0, sig2b0;
  double astar, bstar1, bstar2, os, ns, sumsq;
  double sumb01,sumb02;

  // stuff that I need for player specific beta's;
  //	double *H = R_VectorInit((*nsubject)*(nb)*(max_nobs),0.0);
  double *H = R_VectorInit((*nsubject)*(nb)*(max_nobs),0.0);
  double *tH = R_VectorInit((*nsubject)*(nb)*(max_nobs),0.0);
  double *HtH = R_Vector(max_nobs*max_nobs);
  double *Hty1 = R_Vector(max_nobs*max_nobs);
  double *Hty2 = R_Vector(max_nobs*max_nobs);

  double *z_tmp = R_Vector(max_nobs);
  double *y1_tmp = R_Vector(max_nobs);
  double *y2_tmp = R_Vector(max_nobs);

  double sumy1_Hb, s2star1, mstar1;
  double sumy2_Hb, s2star2, mstar2;
  double *y1_b0 = R_Vector(max_nobs);
  double *y2_b0 = R_Vector(max_nobs);
  double *Hb1 = R_Vector(max_nobs);
  double *Hb2 = R_Vector(max_nobs);


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
  double *thetatmp1 = R_VectorInit(nb, 0.0);
  double *thetatmp2 = R_VectorInit(nb, 0.0);
  double *bpred1 = R_VectorInit(nb, 0.0);
  double *bpred2 = R_VectorInit(nb, 0.0);
  double *ppredtmp1 = R_VectorInit(100,0.0);
  double *ppredtmp2 = R_VectorInit(100,0.0);
  double lgcon0, lgcat0=0.0;


  // stuff that I need to compute the lpml;
  double llikeval, lpml_iter = 0.0, elppdWAIC;;
  double *CPO = R_VectorInit((*nsubject), 0.0);
  double *llike_iter = R_VectorInit(*nsubject,0.0);
  double *mnlike = R_VectorInit((*nsubject), 0.0);
  double *mnllike = R_VectorInit((*nsubject), 0.0);
  int like0,nout_0=0;;



  // 	intial function basis matrix
  // 	B-spline basis
  // 	Also initial marginal likelihood, beta vector, and alpha_star;
  // 	And produce prior mean for mu; (parabola centered a 0.5 with max equal to 0.3)
  //	RprintIVecAsMat("catvar", catvar, 1, *nsubject);
  double *Boundaryknots = R_Vector(2*(*nsubject));
  double *bktmp = R_Vector(2);
  double *knotstmp = R_Vector(*nknots);
  csobs=0;
  for(j = 0; j < *nsubject; j++){

    //		Boundaryknots[0*(*nsubject) + j] = knots[0]-1.0;
    //		Boundaryknots[1*(*nsubject) + j] = knots[(*nknots)*(*nsubject)-1]+1.0;
    Boundaryknots[0*(*nsubject) + j] = -0.05;
    Boundaryknots[1*(*nsubject) + j] = 1.05;

  }

  //	RprintVecAsMat("Boundaryknots", Boundaryknots, 2, *nsubject);
  // Create the H matrix for predictions

  double *Hpred = R_VectorInit((*nsubject)*(nb)*(*npredobs), 0.0);


  // Use same knots used for all players

  for(k = 0; k < *nknots; k++){

    knotstmp[k] = knots[k*(*nsubject)+0];

  }

  bktmp[0] = Boundaryknots[0];
  bktmp[1] = Boundaryknots[2*(*nsubject)-1];

  for(k = 0; k < *npredobs; k++){

    z_tmp[k] = zpred[k]/((double) *npredobs);

  }

  //  bsb(z_tmp, knotstmp, bktmp, Hpred, *npredobs, *nknots, *q);

  if(*balanced==1){
    for(t=0; t < nobs[0]; t++){
      for(kk=0; kk<nb; kk++){
        H[t*(nb) + kk] = Hmat[t*(nb) + kk];
      }
    }
  }

  //	RprintVecAsMat("Hpred", Hpred, *npredobs, nb);
  //	RprintVecAsMat("H", H, *nobs, nb);
  //	RprintVecAsMat("Hmat", Hmat, *nobs, nb);

  //	RprintVecAsMat("z_tmp", z_tmp, 1, *npredobs);
  //	RprintVecAsMat("knotstmp", knotstmp, 1, *nknots);
  //	RprintVecAsMat("bktmp", bktmp, 1, 2);
  //	Rprintf("npredobs = %d\n", *npredobs);
  //	Rprintf("nknots = %d\n", *nknots);
  //	Rprintf("q = %d\n", *q);

  //	bsb(z_tmp, knotstmp, bktmp, Hpred, *npredobs, *nknots, *q);

  //	RprintVecAsMat("Hpred", Hpred, *npredobs, nb);

  double *mnmle = R_Vector(*ncon);
  double *s2mle = R_Vector(*ncon);
  for(p = 0; p < *ncon; p++){
    sumx = 0.0, sumx2=0.0;
    for(j = 0; j < *nobs; j ++){
      sumx = sumx + Xcon[j*(*ncon) + p];
      sumx2 = sumx2 + Xcon[j*(*ncon) + p]*Xcon[j*(*ncon) + p];
    }

    mnmle[p] = sumx/((double) *nobs);
    s2mle[p] = sumx2/((double) *nobs) - mnmle[p]*mnmle[p];
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


      //			Rprintf("j = %d =================== \n", j);

      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
      //			RprintIVecAsMat("nh", nh, 1, *nsubject);
      //			Rprintf("nclus_iter = %d\n", nclus_iter);

      //			RprintVecAsMat("tau2h", tau2h, 1, *nsubject);
      //			RprintVecAsMat("lamh", lamh, 1, *nsubject);
      //			RprintVecAsMat("thetah", thetah, nb, *nsubject);
      //			RprintIVecAsMat("nh", nh, 1, *nsubject);
      //			Rprintf("Si_iter[j] = %d\n", Si_iter[j]);
      //			Rprintf("nh[Si_iter[j]-1] = %d\n", nh[Si_iter[j]-1]);

      if(nh[Si_iter[j]-1] > 1){

        // Observation belongs to a non-singleton ...
        nh[Si_iter[j]-1] = nh[Si_iter[j]-1] - 1;

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
          auxt1 = tau2h1[iaux-1];
          tau2h1[iaux-1] = tau2h1[nclus_iter-1];
          tau2h1[nclus_iter-1] = auxt1;

          auxt2 = tau2h2[iaux-1];
          tau2h2[iaux-1] = tau2h2[nclus_iter-1];
          tau2h2[nclus_iter-1] = auxt2;

          auxl1 = lamh1[iaux-1];
          lamh1[iaux-1] = lamh1[nclus_iter-1];
          lamh1[nclus_iter-1] = auxl1;

          auxl2 = lamh2[iaux-1];
          lamh2[iaux-1] = lamh2[nclus_iter-1];
          lamh2[nclus_iter-1] = auxl2;

          for(b = 0; b < nb; b++){

            auxth1 = thetah1[b*(nclus_iter) + iaux-1];
            thetah1[b*(nclus_iter) + iaux-1] = thetah1[b*(nclus_iter) + nclus_iter-1];
            thetah1[b*(nclus_iter) + nclus_iter-1] = auxth1;

            auxth2 = thetah2[b*(nclus_iter) + iaux-1];
            thetah2[b*(nclus_iter) + iaux-1] = thetah2[b*(nclus_iter) + nclus_iter-1];
            thetah2[b*(nclus_iter) + nclus_iter-1] = auxth2;
          }

          //					RprintVecAsMat("thetah", thetah, nb, *nsubject);

          // the number of members in cluster is also swapped with the last
          nh[iaux-1] = nh[nclus_iter-1];
          nh[nclus_iter-1] = 1;

        }


        // Now remove the ith obs and last cluster;
        nh[nclus_iter-1] = nh[nclus_iter-1] - 1;
        nclus_iter = nclus_iter - 1;


      }

      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);

      //			RprintVecAsMat("tau2h", tau2h, 1, *nsubject);
      //			RprintVecAsMat("lamh", lamh, 1, *nsubject);
      //			RprintVecAsMat("thetah", thetah, nb, *nsubject);
      //			RprintIVecAsMat("nh", nh, 1, *nsubject);

      //			Rprintf("nclus_iter = %d\n", nclus_iter);


      //			RprintIVecAsMat("x", x, 1, *nsubject);
      //			RprintIVecAsMat("x2", x2, 1, *nsubject);


      // The atoms have been relabeled if necessary and now we need to
      // update Si.

      for(b = 0; b < nb; b++){

        btmp1[b] = beta1_iter[b*(*nsubject) + j];
        btmp2[b] = beta2_iter[b*(*nsubject) + j];

      }

      //			RprintVecAsMat("btmp", btmp, 1, nb);


      ///////////////////////////////////
      //
      // Begin the cluster probabilities
      //
      //////////////////////////////////
      for(k = 0 ; k < nclus_iter; k++){

        //				Rprintf("k = %d ==================== \n", k);

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

        //				Rprintf("tau2h = %f\n", tau2h[k]);
        //				Rprintf("lamh = %f\n", lamh[k]);
        //				RprintVecAsMat("thtmp", thtmp, 1, nb);
        //				RprintVecAsMat("btmp", btmp, 1, nb);
        //				RprintVecAsMat("oV", oV, nb, nb);

        ldo1 = 2.0*nb*log(lamh1[k]);
        ldo2 = 2.0*nb*log(lamh2[k]);

        //				Rprintf("nh[k] = %d\n", nh[k]);


        lgconY = 0.0;
        lgconN = 0.0;
        //				RprintIVecAsMat("Si_iter",Si_iter, 1, *nsubject);
        for(p=0; p<(*ncon); p++){

          //					Rprintf("p = %d ====== \n", p) ;
          nhtmp = 0;
          sumx = 0.0;
          sumx2 = 0.0;
          for(jj = 0; jj < *nsubject; jj++){
            if(jj != j){
              if(Si_iter[jj] == k+1){
                tmp = Xcon[jj*(*ncon)+p];

                sumx = sumx + tmp;
                sumx2 = sumx2 + tmp*tmp;

                nhtmp = nhtmp+1;
                //								Rprintf("nhtmp = %d\n", nhtmp);
              }
            }
          }

          //					Rprintf("nh[k] = %d\n", nh[k]);
          //					RprintVecAsMat("xcontmp", xcontmp, 1, nhtmp);

          //					sumx = 0.0;
          //					sumx2 = 0.0;
          //					for(t = 0; t < nhtmp; t++){

          //						sumx = sumx + xcontmp[t];
          //						sumx2 = sumx2 + xcontmp[t]*xcontmp[t];

          //					}

          //					Rprintf("sumx = %f\n", sumx);
          //					Rprintf("sumx2 = %f\n", sumx2);
          //					Rprintf("nhtmp = %d\n", nhtmp);

          if(*gcontype==1){ // Auxilliary
            if(*consim==1){
              lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
              //							lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 1, 1);

              lgconN = lgconN + lgcont;
            }
            if(*consim==2){
              lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
              //							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 1, 1);

              lgconN = lgconN + lgcont;
            }

          }

          if(*gcontype==2){ //Double Dipper
            if(*consim==1){
              lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
              //							lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 1, 1);

              lgconN = lgconN + lgcont;
            }
            if(*consim==2){
              lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
              //							lgonct = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 1, 1);

              lgconN = lgconN + lgcont;
            }
          }

          if(*gcontype==3){ // Cluster Variance
            //						Rprintf("gsimconEV = %f\n", gsimconEV(sumx, sumx2, nhtmp, alpha, 0));
            lgcont = gsimconEV(sumx, sumx2, nhtmp, alpha,1);

            lgconN = lgconN + lgcont;
          }

          //					Rprintf("lgcontN = %f\n", lgcont);
          //					Rprintf("lgconN = %f\n", lgconN);

          // now add jth individual back;
          //					RprintVecAsMat("xconpred", xconpred, 1, *npred);

          //					Rprintf("sumx = %f\n", sumx);
          //					Rprintf("sumx2 = %f\n", sumx2);
          //					Rprintf("Xcon[j*(*ncon)+p] = %f\n", Xcon[j*(*ncon)+p]);

          xcontmp[nhtmp] = Xcon[j*(*ncon)+p];
          sumx = sumx + Xcon[j*(*ncon)+p];
          sumx2 = sumx2 + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          nhtmp = nhtmp+1;



          //					Rprintf("nhtmp = %d\n", nhtmp);
          //					Rprintf("sumx = %f\n", sumx);
          //					Rprintf("sumx2 = %f\n", sumx2);
          if(*gcontype==1){ // Auxilliary
            if(*consim==1){
              lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
              //							lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 1, 1);
              lgconY = lgconY + lgcont;
            }
            if(*consim==2){
              lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
              //							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 1, 1);
              lgconY = lgconY + lgcont;
            }
          }
          if(*gcontype==2){ //Double Dipper
            if(*consim==1){
              lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
              //							lgcont = gsimconNN(m0, v2, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 1, 1);
              lgconY = lgconY + lgcont;
            }
            if(*consim==2){
              lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
              //							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 1, 1);
              lgconY = lgconY + lgcont;
            }
          }
          if(*gcontype==3){ // variance
            //						Rprintf("gsimconEV = %f\n", gsimconEV(sumx, sumx2, nhtmp, alpha, 0));
            lgcont = gsimconEV(sumx, sumx2, nhtmp, alpha,1);
            lgconY = lgconY + lgcont;
          }

          //					Rprintf("lgcontY = %f\n", lgcont);
          //					Rprintf("lgconY = %f\n", lgconY);


        }


        //				RprintIVecAsMat("Si_iter", Si_iter, 1, *nobs);

        lgcatY=0.0;
        lgcatN=0.0;


        for(p=0; p<(*ncat); p++){


          //					Rprintf("p = %d ====== \n", p) ;
          for(c=0;c<Cvec[p];c++){nhc[c]=0;}

          nhtmp = 0;
          for(jj = 0; jj < *nsubject; jj++){
            //						Rprintf("jj = %d\n", jj);
            if(jj != j){
              //							Rprintf("Xcatstd[jj*(*ncat)+p] = %d\n", Xcatstd[jj*(*ncat)+p]);
              //							Rprintf("Si_iter[jj] = %d\n", Si_iter[jj]);

              if(Si_iter[jj]==k+1){
                nhc[Xcat[jj*(*ncat)+p]] = nhc[Xcat[jj*(*ncat)+p]] + 1; // this needs to be a vector
                //                              RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
                nhtmp = nhtmp+1;

              }
            }
            //						RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
          }

          //					RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
          //					Rprintf("nhtmp = %d\n", nhtmp);


          if(*gcattype==1){
            lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
            lgcatN = lgcatN + lgcatt;
          }
          if(*gcattype==2){
            lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
            lgcatN = lgcatN + lgcatt;
          }
          if(*gcattype==3){// Using Entropy instead of variance here
            //						Rprintf("gsimconEV = %f\n", gsimconEV(sumx, sumx2, nhtmp, alpha,  0));
            //						lgcatt = gsimconEV(sumx, sumx2, nhtmp, alpha,1);
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
          //					Rprintf("lgcattN = %f\n", lgcatt);
          //					Rprintf("lgcatN = %f\n", lgcatN);

          //					Rprintf("Xcatstd[j*(*ncat)+p] = %f\n", Xcatstd[j*(*ncat)+p]);

          nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1;
          nhtmp = nhtmp + 1;

          //					RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
          //					Rprintf("nhtmp = %d\n", nhtmp);

          if(*gcattype==1){
            lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
            lgcatY = lgcatY + lgcatt;
          }
          if(*gcattype==2){
            lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
            lgcatY = lgcatY + lgcatt;
          }
          if(*gcattype==3){// Using Entropy instead of variance here
            //						Rprintf("gsimconEV = %f\n", gsimconEV(sumx, sumx2, nhtmp, alpha, 0));
            //						lgcatt = gsimconEV(sumx, sumx2, nhtmp, alpha, 1);
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
          //					Rprintf("lgcattY = %f\n", lgcatt);
          //					Rprintf("lgcatY = %f\n", lgcatY);


        }

        gtilY[k] = lgconY + lgcatY;
        gtilN[k] = lgconN + lgcatN;


        //				RprintVecAsMat("btmp", btmp, 1, nb);
        //				RprintVecAsMat("thtmp", thtmp, 1, nb);

        //				Rprintf("sumsq1 = %f\n", sumsq1);
        //				Rprintf("log(tau2h[k] = %f\n", log(tau2h[k]));
        //				Rprintf("dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) = %f\n",
        //							dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1));
        //				Rprintf("-0.5*(nb*log(tau2h[k]) + (1/tau2h[k])*sumsq) = %f\n",
        //								-0.5*(nb*log(tau2h[k]) + (1/tau2h[k])*sumsq));
        //				Rprintf("log(nh[k]) = %f\n", log(nh[k]));
        //				Rprintf("lgcatY = %f\n", lgcatY);
        //				Rprintf("lgcatN = %f\n", lgcatN);
        //				Rprintf("lgcatY - lgcatN = %f\n", lgcatY - lgcatN);
        //				Rprintf("lgconY = %f\n", lgconY);
        //				Rprintf("lgconN = %f\n", lgconN);
        //				Rprintf("lgconY - lgconN = %f\n", lgconY - lgconN);




        ph[k] = dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1) +
          dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1) +
          log((double) nh[k]) + // Scale parameter from DP
          lgcatY - lgcatN +  // Ratio of categorical variables
          lgconY - lgconN;   // Ratio of continuous variabels


        if(*PPM){
          ph[k] = dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1) +
            dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1) +
            log((double) nh[k]);  // DP part of cohesion function
        }

        //				Rprintf("ph[k] = %f\n", ph[k]);
      }






      //			RprintVecAsMat("ph", ph, 1, nclus_iter);


      lamdraw1 = runif(0,A1);
      lamdraw2 = runif(0,A2);
      tau2draw1 = 1/rgamma(at,bt); // shape and scale  E(tau2) = atbt
      tau2draw2 = 1/rgamma(at,bt); // shape and scale  E(tau2) = atbt

      //			lamdraw = 0.1;
      //			Rprintf("lamdraw = %f\n", lamdraw);
      //			Rprintf("tau2draw = %f\n", tau2draw);

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

      //			RprintVecAsMat("mu_iter =", mu_iter, 1, nb);
      ran_mvnorm(mu1_iter, nV1, nb, scr1, thetadraw1);
      ran_mvnorm(mu2_iter, nV2, nb, scr1, thetadraw2);

      //			RprintVecAsMat("thetadraw", thetadraw, 1, nb);
      //			RprintVecAsMat("btmp", btmp, 1, nb);

      //			Rprintf("K[0]-1 = %f\n", K[0]-1);
      //			Rprintf("1/(K[0]-1) = %f\n", 1/(K[0]-1));

      //			thetadraw[0] = rnorm(mu_iter[0], sqrt(tau2draw*(1/(K[0]-1))));
      //			for(b=0; b<nb-1; b++){
      //				thetadraw[b+1] = rnorm(thetadraw[b] + mu_iter[b], sqrt(tau2draw));
      //			}

      //			RprintVecAsMat("thetadraw", thetadraw, 1, nb);


      ldo1 = 2.0*nb*log(lamdraw1);
      ldo2 = 2.0*nb*log(lamdraw2);


      lgcondraw = 0.0;
      for(p=0;p<(*ncon);p++){
        xcontmp[0] = Xcon[j*(*ncon)+p];
        //				Rprintf("Xcon[j*(*ncon)+p]=%f\n", Xcon[j*(*ncon)+p]);
        if(*gcontype==1){ // Auxilliary
          if(*consim==1){
            lgcont = gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0, 1);
            //						lgcont = gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,1, 1);
            lgcondraw = lgcondraw + lgcont;
          }
          if(*consim==2){
            lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 0,0, 1);
            //						lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 0,1, 1);
            lgcondraw = lgcondraw + lgcont;
          }
        }
        if(*gcontype==2){ // Double Dipper
          if(*consim==1){
            lgcont = gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p], 1, 1, 0, 1);
            //						lgcont = gsimconNN(m0,v2,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p], 1, 1, 1, 1);
            lgcondraw = lgcondraw + lgcont;
          }
          if(*consim==2){
            lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 1, 0, 1);
            //						lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 1, 1, 1);
            lgcondraw = lgcondraw + lgcont;
          }
        }
        if(*gcontype==3){ // Variance
          lgcont = gsimconEV(xcontmp[0], xcontmp[0]*xcontmp[0], 1,alpha,1);
          lgcondraw = lgcondraw + lgcont;
        }

        //				Rprintf("lgcondraw = %f\n", lgcondraw);
      }

      lgcatdraw = 0.0;
      for(p=0;p<(*ncat);p++){
        for(c=0;c<Cvec[p];c++){nhc[c] = 0;}

        nhc[Xcat[j*(*ncat)+p]] = 1;
        //				RprintIVecAsMat("nhc =", nhc, 1, Cvec[p]);

        //				Rprintf("Xcat[j*(*ncon)+p] = %d\n", Xcat[j*(*ncon)+p]);



        if(*gcattype==1){
          lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
          lgcatdraw = lgcatdraw + lgcatt;
        }
        if(*gcattype==2){
          lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
          lgcatdraw = lgcatdraw + lgcatt;
        }
        if(*gcattype==3){
          //					lgcatt = gsimconEV(xcattmp[0], xcattmp[0]*xcattmp[0], 1,alpha, 1);
          //					lgcatdraw = lgcatdraw + lgcatt;
          lgcatdraw = lgcatdraw + -(alpha)*0;
        }

        //				Rprintf("lgcatdraw = %f\n", lgcatdraw);
      }

      gtilY[nclus_iter] = lgcondraw + lgcatdraw;
      gtilN[nclus_iter] = lgcondraw + lgcatdraw;

      //			Rprintf("dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1) = %f\n",
      //						dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1));
      //			Rprintf("-0.5*(nb*log(tau2draw) + (1/tau2draw)*sumsq) = %f\n",
      //						-0.5*(nb*log(tau2draw) + (1/tau2draw)*sumsq));
      //			Rprintf("log(M) = %f\n", log(Mdp));
      //			Rprintf("lgcondraw = %f\n", lgcondraw);
      //			Rprintf("lgcatdraw = %f\n", lgcatdraw);

      ph[nclus_iter] = dmvnorm(btmp1,thetadraw1,oV1,nb,ldo1,scr1,1) +
        dmvnorm(btmp2,thetadraw2,oV2,nb,ldo2,scr2,1) +
        log(Mdp) +
        lgcondraw +
        lgcatdraw;



      if(*PPM){
        ph[nclus_iter] = dmvnorm(btmp1,thetadraw1,oV1,nb,ldo1,scr1,1) +
          dmvnorm(btmp2,thetadraw2,oV2,nb,ldo2,scr2,1) +
          log(Mdp); //DP part
      }

      //			RprintVecAsMat("ph", ph, 1, nclus_iter+1);

      /*
      /////////////////////////////////////////////////////////////////////////////
      // This is the calibration used when the similarity is standardized by
      // the sum of all cluster similarity values.
      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////
      if(*calibrate==1){

      //				Rprintf("k ====== %d\n", k);

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

      //				Rprintf("tau2h = %f\n", tau2h[k]);
      //				Rprintf("lamh = %f\n", lamh[k]);
      //				RprintVecAsMat("thtmp", thtmp, 1, nb);
      //				RprintVecAsMat("btmp", btmp, 1, nb);
      //				RprintVecAsMat("oV", oV, nb, nb);

      ldo1 = 2.0*nb*log(lamh1[k]);
      ldo2 = 2.0*nb*log(lamh2[k]);


      //				Rprintf("k = %d ==================== \n", k);
      //				Rprintf("nh[k] = %d\n", nh[k]);
      maxgtilY = gtilY[0];
      maxgtilN = gtilN[0];
      for(k=1; k < nclus_iter+1; k++){
      if(maxgtilN < gtilN[k]) maxgtilN = gtilY[k];
      if(maxgtilY < gtilY[k]) maxgtilY = gtilY[k];
      }
      //				Rprintf("maxgtilY = %f\n", maxgtilY);
      //				Rprintf("maxgtilN = %f\n", maxgtilN);
      sgY=0.0;
      sgN=0.0;
      for(k=0; k<nclus_iter+1; k++){
      gtilN[k] = exp(gtilN[k] - maxgtilY);
      sgN = sgN + gtilN[k];

      gtilY[k] = exp(gtilY[k] - maxgtilY);
      sgY = sgY + gtilY[k];
      }

      //				Rprintf("sgY = %f\n", sgY);
      //				Rprintf("sgN = %f\n", sgN);
      for(k=0; k<nclus_iter; k++){
      lgtilY = log(gtilY[k]) - log(sgY);
      lgtilN = log(gtilN[k]) - log(sgN);

      ph[k] = dmvnorm(btmp1, thtmp1, oV1, nb, ldo1, scr1, 1) +
      dmvnorm(btmp2, thtmp2, oV2, nb, ldo2, scr2, 1) +
      log((double) nh[k]) +
      lgtilY - lgtilN; //This takes into account both cont and cat vars

      //					Rprintf("lgtilY = %f\n", lgtilY);
      //					Rprintf("lgtilN = %f\n", lgtilN);
      }
      //				Rprintf("mudraw = %f\n", mudraw);
      //				Rprintf("sdraw = %f\n", sdraw);

      // It is not completely clear to me what I need to do with the calibration
      // of g(.) for a new cluster.  Do I use it's weight when standardizing?
      ph[nclus_iter] = dmvnorm(btmp1,thetadraw1,oV1,nb,ldo1,scr1,1) +
      dmvnorm(btmp2,thetadraw2,oV2,nb,ldo2,scr2,1) +
      log(Mdp) +
      log(gtilN[nclus_iter]) - log(sgN);

      //				RprintVecAsMat("ph", ph, 1, nclus_iter+1);

      }
      */


      for(k = 0; k < nclus_iter+1; k++) phtmp[k] = ph[k];


      R_rsort(phtmp,  nclus_iter+1) ;

      //			RprintVecAsMat("phtmp ", phtmp, 1, nclus_iter+1);


      maxph = phtmp[nclus_iter];

      //			Rprintf("maxph = %f\n", maxph);

      denph = 0.0;
      for(k = 0; k < nclus_iter+1; k++){

        ph[k] = exp(ph[k] - maxph);
        //				ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
        denph = denph + ph[k];

      }

      //			RprintVecAsMat("ph", ph, 1, nclus_iter+1);

      for(k = 0; k < nclus_iter+1; k++){

        probh[k] = ph[k]/denph;

      }
      //			Rprintf("denph = %f\n", denph);

      //			RprintVecAsMat("probh", probh, 1, nclus_iter+1);

      uu = runif(0.0,1.0);
      //			Rprintf("uu = %f\n", uu);



      cprobh= 0.0;;
      for(k = 0; k < nclus_iter+1; k++){

        cprobh = cprobh + probh[k];

        if (uu < cprobh){

          iaux = k+1;
          break;
        }
      }


      // 			Rprintf("iaux = %d\n \n \n", iaux);

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

        //				RprintVecAsMat("thetah after", thetah, nb, *nsubject);
      }


      // 			Rprintf("Si_iter = %d\n \n \n", Si_iter[j]);
      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);

      //			RprintVecAsMat("tau2h", tau2h, 1, *nsubject);
      //			RprintVecAsMat("lamh", lamh, 1, *nsubject);
      //			RprintVecAsMat("thetah", thetah, nb, *nsubject);


    }


    //		RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
    //		Rprintf("nclus_iter = %d\n", nclus_iter);







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
    for(j = 0; j < *nsubject; j++){

      for(t = 0; t < nobs[j]; t++){

        y1_tmp[t] = y1[csobs];
        y2_tmp[t] = y2[csobs];

        z_tmp[t] = z[csobs]/((double) nobs[j]);

        y1_b0[t] = y1_tmp[t]-beta01_iter[j];
        y2_b0[t] = y2_tmp[t]-beta02_iter[j];

        csobs = csobs+1;

      }

      for(kk = 0; kk < *nknots; kk++){

        knotstmp[kk] = knots[kk*(*nsubject) + j];

      }

      bktmp[0] = Boundaryknots[0*(*nsubject) + j];
      bktmp[1] = Boundaryknots[1*(*nsubject) + j];


//      if(*balanced==0) bsb(z_tmp, knotstmp, bktmp, H, nobs[j], *nknots, *q);

      mat_transpose(H, tH, nobs[j], nb);

      matrix_product(tH, H, HtH, nb, nb, nobs[j]);

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
