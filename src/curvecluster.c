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
* iCPO =
*
*****************************************************************************************/


void mcmc_curvecluster(	int *draws, int *burn, int *thin,
                        int *nsubject, int *nobs,
                        double *y, double *z,
                        int *q, int *nknots, double *knots, double *K,
                        int *ncon, int *ncat,int *Cvec,
                        double *Xcon, int *Xcat,
                        int *PPM, double *M,
                        int *gcontype, int *gcattype, int *consim,
                        int *npred, int *npredobs, double *Xconp, int *Xcatp,
                        double *simParms, double *Aparm,
                        int *calibrate, double *modelPriors,
                        double *Hmat, int *balanced, double *mh,
                        double *beta, double *beta0, double *sig2, double *mub0,
                        double *sig2b0,double *lam, double *tau2, double *theta,
                        double *mu, int *Si, int *nclus,
                        double *ppred,   int *predclass, double *llike,
                        double *lpml, double *WAIC){

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

  int c, i, ii, j, t, p, jj, k, kk, pp, b, bb;
  int csobs = 0;
  int nb, N ;

  if(*balanced==0) nb = *nknots + *q + 1;  //Number of B-spline basis after filling in all players;
  if(*balanced==1) nb = *nknots + *q;  //Number of B-spline basis after filling in all players;
  int nout = (*draws - *burn)/(*thin);

  Rprintf("nsubject = %d\n", *nsubject);
  Rprintf("nknots = %d\n", *nknots);
  Rprintf("nb = %d\n", nb);
  Rprintf("ncon = %d\n", *ncon);
  Rprintf("ncat = %d\n", *ncat);

  //	RprintIVecAsMat("nobs", nobs, 1, *nsubject);
  RprintIVecAsMat("Cvec", Cvec, 1, *ncat);
  //	RprintVecAsMat("knots", knots, *nknots, *nsubject);


  double *zpred = R_Vector(*npredobs);
  for(pp = 0; pp < *npredobs; pp++) zpred[pp] = (double) (pp+1);
  //	RprintVecAsMat("zpred", zpred, 1, *npredobs);
  //	nstrtkntsindx = ((int) floor((double) *nstartknots));

  //	Rprintf("nstrtkntsindx = %d\n", nstrtkntsindx);
  int  max_nobs;

  double *max_z = R_Vector(*nsubject);
  double *min_z = R_Vector(*nsubject);

  double *ypy = R_VectorInit(*nsubject, 0.0);
  double *sumy = R_VectorInit(*nsubject,0.0);

  max_nobs = nobs[0];
  csobs = 0;
  for(j = 0; j < (*nsubject); j++){
    //		Rprintf("j = %d\n", j);
    max_z[j] = z[csobs];
    min_z[j] = z[csobs];

    for(t = 0; t < nobs[j]; t++){

      //			Rprintf("t = %d\n", t);

      ypy[j] = ypy[j] + y[csobs]*y[csobs];
      sumy[j] = sumy[j] + y[csobs];

      if(max_z[j] < z[csobs]) max_z[j] = z[csobs];
      if(min_z[j] > z[csobs]) min_z[j] = z[csobs];

      csobs = csobs + 1;

    }
    //		Rprintf("min_z = %f\n",min_z[j]);
    //		Rprintf("max_z = %f\n",max_z[j]);
    if(max_nobs < nobs[j]) max_nobs = nobs[j];


  }


  N = csobs;


  //	RprintVecAsMat("lnobs", lnobs, 1, *nsubject);

  //	RprintVecAsMat("y", y, 1, N);
  //	RprintVecAsMat("z", z, 1, N);

  //	RprintVecAsMat("Xcon", Xcon, *nsubject, *ncon);
  //	RprintIVecAsMat("Xcat", Xcat, *nsubject, *ncat);

//  RprintVecAsMat("ypy = ", ypy, 1, *nsubject);
//  RprintVecAsMat("sumy", sumy, 1, *nsubject);

  //	RprintVecAsMat("max_z = ", max_z, 1, *nsubject);
  //	RprintVecAsMat("min_z = ", min_z, 1, *nsubject);
  //	Rprintf("max_nobs = %d\n", max_nobs);
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
  // Calculate intelligent starting values
  //
  // ===================================================================================
  // Calculate the average career length (in years from first game played) and
  // average number of games played using all data;
  double mnNi=0.0;
  for(j = 0; j < *nsubject; j++){
    mnNi = mnNi + nobs[j]*(1/(double) *nsubject);
  }




  // ===================================================================================
  //
  // Memory vectors to hold MCMC iterates for non cluster specific parameters
  //
  // ===================================================================================

  double *sig2_iter = R_VectorInit(*nsubject, ((modelPriors[0]-0.0)/2.0)*((modelPriors[0]-0.0)/2.0));
  double *beta0_iter = R_VectorInit(*nsubject, 0.0);

  double *beta_iter = R_VectorInit(nb*(*nsubject),0.0);

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
    //		Si_iter[j] = j+1;
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
  RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
  RprintIVecAsMat("nh", nh, 1, *nsubject);


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
  double *scr3 = R_Vector((max_nobs));


  // stuff that I need to update Si (cluster labels);
  int iaux, nhtmp;

  double auxt, auxl, auxth, uu, tmp;
  double sumsq, maxph, denph, cprobh;
  double lamdraw, tau2draw;
  double sumx, sumx2;

  double *thetadraw = R_VectorInit(nb, 0.0);

  double *ph = R_VectorInit(*nsubject, 0.0);
  double *phtmp = R_VectorInit(*nsubject, 0.0);
  double *probh = R_VectorInit(*nsubject, 0.0);

  double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw;
  double lgcont,lgcatt;
  //	double mgconN, mgconY, mgcatN, mgcatY, sgconN, sgconY, sgcatY, sgcatN;

  double *gtilN = R_VectorInit((*nobs+1),0.0);
  double *gtilY = R_VectorInit((*nobs+1),0.0);
  double sgY, sgN,  lgtilN, lgtilY, maxgtilY, maxgtilN;


  // stuff I need to update sig2 (player specific), mub0, sig2b0;
  double astar, bstar, os, ns, sumb0;

  // stuff that I need for player specific beta's;
  double *H = R_VectorInit((*nsubject)*(nb)*(max_nobs),0.0);
  double *tH = R_VectorInit((*nsubject)*(nb)*(max_nobs),0.0);
  double *HtH = R_Vector(max_nobs*max_nobs);

  double *z_tmp = R_Vector(max_nobs);
  double *y_tmp = R_Vector(max_nobs);

  double sumy_Hb, s2star, mstar;
  double *y_b0 = R_Vector(max_nobs);
  double *Hb = R_Vector(max_nobs);
  double *sumH = R_VectorInit(max_nobs,0.0);


  // Create the inverse of the K penalty matrix
  double ld;
  double *Kinv = R_Vector(nb*nb);
  for(b = 0; b < nb; b++){for(bb = 0; bb < nb; bb++){Kinv[b*nb+bb] = K[b*nb+bb];}}
  cholesky(Kinv, nb, &ld);
  inverse_from_cholesky(Kinv, scr1, scr2, nb);

  //	RprintVecAsMat("Kinv", Kinv, nb,nb);
  //	RprintVecAsMat("K", K, nb,nb);
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
  int like0,nout_0=0;;



  // 	intial function basis matrix
  // 	B-spline basis
  // 	Also initial marginal likelihood, beta vector, and alpha_star;
  // 	And produce prior mean for mu; (parabola centered a 0.5 with max equal to 0.3)
  //	RprintIVecAsMat("catvar", catvar, 1, *nsubject);
  double *Boundaryknots = R_Vector(2*(*nsubject));
  double *bktmp = R_Vector(2);
  double *knotstmp = R_Vector(*nknots);
  //	double *zmu_tmp = R_Vector(101);
  //	double *ymu_tmp = R_Vector(101);
  csobs=0;
  for(j = 0; j < *nsubject; j++){

    //		Boundaryknots[0*(*nsubject) + j] = knots[0]-1.0;
    //		Boundaryknots[1*(*nsubject) + j] = knots[(*nknots)*(*nsubject)-1]+1.0;
    Boundaryknots[0*(*nsubject) + j] = -0.05;
    Boundaryknots[1*(*nsubject) + j] = 1.05;

  }

  //	RprintVecAsMat("Boundaryknots", Boundaryknots, 2, *nsubject);
  // Create the H matrix for predictions

  // double *Hpred = R_VectorInit((*nsubject)*(nb)*(*npredobs), 0.0);



  // Use same knots used for all players

  for(k = 0; k < *nknots; k++){

    knotstmp[k] = knots[k*(*nsubject)+0];

  }

  bktmp[0] = Boundaryknots[0];
  bktmp[1] = Boundaryknots[2*(*nsubject)-1];

  for(k = 0; k < *npredobs; k++){

    z_tmp[k] = zpred[k]/((double) *npredobs);

  }

  //	RprintVecAsMat("z_tmp", z_tmp, 1, *npredobs);
  //	RprintVecAsMat("knotstmp", knotstmp, 1, *nknots);

  //	RprintVecAsMat("bktmp", bktmp, 1, 2);
  //	Rprintf("npredobs = %d\n", *npredobs);
  //	Rprintf("nknots = %d\n", *nknots);
  //	Rprintf("q = %d\n", *q);

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
          auxt = tau2h[iaux-1];
          tau2h[iaux-1] = tau2h[nclus_iter-1];
          tau2h[nclus_iter-1] = auxt;

          auxl = lamh[iaux-1];
          lamh[iaux-1] = lamh[nclus_iter-1];
          lamh[nclus_iter-1] = auxl;

          for(b = 0; b < nb; b++){

            auxth = thetah[b*(*nsubject) + iaux-1];
            thetah[b*(*nsubject) + iaux-1] = thetah[b*(*nsubject) + nclus_iter-1];
            thetah[b*(*nsubject) + nclus_iter-1] = auxth;
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

        btmp[b] = beta_iter[b*(*nsubject) + j];

      }

      // RprintVecAsMat("btmp", btmp, 1, nb);


      ///////////////////////////////////
      //
      // Begin the cluster probabilities
      //
      //////////////////////////////////
      for(k = 0 ; k < nclus_iter; k++){

        //				Rprintf("k = %d ==================== \n", k);

        for(b = 0; b < nb; b++){

          thtmp[b] = thetah[b*(*nsubject) + k];


          for(bb = 0; bb < nb; bb++){

            oV[b*nb+bb] = 0.0;

            if(b == bb){

              oV[b*nb+bb] = 1/(lamh[k]*lamh[k]);
              //							oV[b*nb+bb] = 1/(lamh[k]);
            }
          }
        }

        // Rprintf("tau2h = %f\n", tau2h[k]);
//         Rprintf("lamh = %f\n", lamh[k]);
//         RprintVecAsMat("thtmp", thtmp, 1, nb);
//         RprintVecAsMat("btmp", btmp, 1, nb);
        // RprintVecAsMat("oV", oV, nb, nb);

        ldo = 2.0*nb*log(lamh[k]);

        //Rprintf("nh[k] = %d\n", nh[k]);


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

        //				Rprintf("lamh[k] = %f\n", lamh[k]);
        //				Rprintf("log(tau2h[k] = %f\n", log(tau2h[k]));
        //				Rprintf("dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) = %f\n",
        //							dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1));
        //				Rprintf("log(nh[k]) = %f\n", log(nh[k]));
        //				Rprintf("lgcatY = %f\n", lgcatY);
        //				Rprintf("lgcatN = %f\n", lgcatN);
        //				Rprintf("lgcatY - lgcatN = %f\n", lgcatY - lgcatN);
        //				Rprintf("lgconY = %f\n", lgconY);
        //				Rprintf("lgconN = %f\n", lgconN);
        //				Rprintf("lgconY - lgconN = %f\n", lgconY - lgconN);




        ph[k] = dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) +
          log((double) nh[k]) + // Scale parameter from DP
          lgcatY - lgcatN +  // Ratio of categorical variables
          lgconY - lgconN;   // Ratio of continuous variabels

        if(*PPM){
          ph[k] = dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) +
            log((double) nh[k]);  // DP part of cohesion function
        }

//        Rprintf("ph[k] = %f\n", ph[k]);
      }






      //			RprintVecAsMat("ph", ph, 1, nclus_iter);


      lamdraw = runif(0,A);
      tau2draw = 1/rgamma(at,bt); // shape and scale  E(tau2) = atbt

//       Rprintf("lamdraw = %f\n", lamdraw);
      // Rprintf("tau2draw = %f\n", tau2draw);

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

//      RprintVecAsMat("mu_iter =", mu_iter, 1, nb);
      ran_mvnorm(mu_iter, nV, nb, scr1, thetadraw);

//      RprintVecAsMat("thetadraw", thetadraw, 1, nb);
//      RprintVecAsMat("btmp", btmp, 1, nb);

      //			Rprintf("K[0]-1 = %f\n", K[0]-1);
      //			Rprintf("1/(K[0]-1) = %f\n", 1/(K[0]-1));




      ldo = 2.0*nb*log(lamdraw);
      //			ldo = nb*log(lamdraw);


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
      //			Rprintf("log(M) = %f\n", log(Mdp));
      //			Rprintf("lgcondraw = %f\n", lgcondraw);
      //			Rprintf("lgcatdraw = %f\n", lgcatdraw);

      ph[nclus_iter] = dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1) +
        log(Mdp) +
        lgcondraw +
        lgcatdraw;



      if(*PPM){
        ph[nclus_iter] = dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1) +
          log(Mdp); //DP part
      }

//      RprintVecAsMat("ph", ph, 1, nclus_iter+1);


      /////////////////////////////////////////////////////////////////////////////
      // This is the calibration used when the similarity is standardized by
      // the sum of all cluster similarity values.
      /////////////////////////////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////
      if(*calibrate==1){

        //				Rprintf("k ====== %d\n", k);

        for(b = 0; b < nb; b++){

          thtmp[b] = thetah[b*(*nsubject) + k];


          for(bb = 0; bb < nb; bb++){

            oV[b*nb+bb] = 0.0;

            if(b == bb){

              oV[b*nb+bb] = 1/(lamh[k]*lamh[k]);
              //							oV[b*nb+bb] = 1/(lamh[k]);
            }
          }
        }

        //				Rprintf("tau2h = %f\n", tau2h[k]);
        //				Rprintf("lamh = %f\n", lamh[k]);
        //				RprintVecAsMat("thtmp", thtmp, 1, nb);
        //				RprintVecAsMat("btmp", btmp, 1, nb);
        //				RprintVecAsMat("oV", oV, nb, nb);

        ldo = 2.0*nb*log(lamh[k]);


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

          ph[k] = dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) +
            log((double) nh[k]) +
            lgtilY - lgtilN; //This takes into account both cont and cat vars

          //					Rprintf("lgtilY = %f\n", lgtilY);
          //					Rprintf("lgtilN = %f\n", lgtilN);
        }
        //				Rprintf("mudraw = %f\n", mudraw);
        //				Rprintf("sdraw = %f\n", sdraw);

        // It is not completely clear to me what I need to do with the calibration
        // of g(.) for a new cluster.  Do I use it's weight when standardizing?
        ph[nclus_iter] = dmvnorm(btmp,thetadraw,oV,nb,ldo,scr1,1) +
          log(Mdp) +
          log(gtilN[nclus_iter]) - log(sgN);

        //				RprintVecAsMat("ph", ph, 1, nclus_iter+1);

      }


      //			RprintVecAsMat("ph ", ph, 1, nclus_iter+1);

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

        tau2h[Si_iter[j]-1] = tau2draw;
        lamh[Si_iter[j]-1] = lamdraw;

        for(b = 0; b < nb; b++){
          thetah[b*(*nsubject) + Si_iter[j]-1] = thetadraw[b];
        }

        //				RprintVecAsMat("thetah after", thetah, nb, *nsubject);
      }


      // 			Rprintf("Si_iter = %d\n \n \n", Si_iter[j]);
      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);

      //			RprintVecAsMat("tau2h", tau2h, 1, *nsubject);
      //			RprintVecAsMat("lamh", lamh, 1, *nsubject);
      //			RprintVecAsMat("thetah", thetah, nb, *nsubject);


    }


    //		RprintVecAsMat("lamh", lamh, 1, *nsubject);
    //		RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);
    //		Rprintf("nclus_iter = %d\n", nclus_iter);





    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // update player specific betas, sigma2, and beta0;								//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    csobs=0;
    for(j = 0; j < *nsubject; j++){
      //			Rprintf("j = %d ====================\n", j);


      //			Rprintf("nobs = %d\n", nobs[j]);
      //			Rprintf("ni_iter = %d\n", ni_iter[j]);
      for(t = 0; t < nobs[j]; t++){

        y_tmp[t] = y[csobs];
        z_tmp[t] = z[csobs]/((double) nobs[j]);
        y_b0[t] = y_tmp[t]-beta0_iter[j];
        csobs = csobs+1;

      }


      for(kk = 0; kk < *nknots; kk++){

        knotstmp[kk] = knots[kk*(*nsubject) + j];

      }

      bktmp[0] = Boundaryknots[0*(*nsubject) + j];
      bktmp[1] = Boundaryknots[1*(*nsubject) + j];

      //			RprintVecAsMat("y_tmp", y_tmp, 1, nobs[j]);
      //			RprintVecAsMat("z_tmp", z_tmp, 1, nobs[j]);
      //			RprintVecAsMat("knots", knotstmp, 1, *nknots);
      //			RprintVecAsMat("Boundaryknots", bktmp, 1, 2);
      //			Rprintf("nobs = %d\n", nobs[j]);
      //			Rprintf("nknots = %d\n", *nknots);
      //			Rprintf("q = %d\n", *q);
      //			Rprintf("nb = %d\n", nb);

//      if(*balanced==0) bsb(z_tmp, knotstmp, bktmp, H, nobs[j], *nknots, *q);

      for(t = 0; t < nobs[j]; t++){
        sumH[t] = 0.0;
      }

      for(t = 0; t < nobs[j]; t++){
        for(b = 0; b < nb; b++){
          sumH[t] = sumH[t] + H[t*nb + b];
        }
      }

      //			RprintVecAsMat("H", H, nobs[j], nb);
      //			RprintVecAsMat("sumH", sumH, 1, nobs[j]);


      mat_transpose(H, tH, nobs[j], nb);

      matrix_product(tH, H, HtH, nb, nb, nobs[j]);

      //			matrix_product(tH, y_tmp, Hty, nb, 1, nobs[j]);

      //			RprintVecAsMat("HtH", HtH, nb, nb);
      //			RprintVecAsMat("H'y", Hty, 1, nb);
      //			Rprintf("Si = %d\n", Si_iter[j]);
      //			Rprintf("lamh[Si] = %f\n", lamh[Si_iter[j]-1]);

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

      // RprintVecAsMat("Sstar", Sstar, nb, nb);
      // RprintVecAsMat("y_b0", y_b0, 1, nobs[j]);
      matrix_product(tH, y_b0, scr2,  nb, 1, nobs[j]);

      // RprintVecAsMat("H'(y-b0)", scr2, 1, nb);

      for(b = 0; b < nb; b++){
	  // Rprintf("thetah[b*(*nsubject) + Si_iter[j]-1] = %f\n", thetah[b*(*nsubject) + Si_iter[j]-1]);
        scr3[b] = (1/sig2_iter[j])*scr2[b] +
          (1/(lamh[Si_iter[j]-1]*lamh[Si_iter[j]-1]))*
          thetah[b*(*nsubject) + Si_iter[j]-1];

      }

      matrix_product(Sstar, scr3, Mstar, nb, 1, nb);

	  // RprintVecAsMat("Mstar", Mstar, 1, nb);
      // RprintVecAsMat("Sstar", Sstar, nb, nb);

      cholesky(Sstar, nb , &ld);


      // RprintVecAsMat("Mstar", Mstar, 1, nb);


      ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);

      // RprintVecAsMat("betai", outrmvnorm, 1, nb);
      // RprintVecAsMat("bhat", bhat, 1, nb);

      for(b = 0; b < nb; b++){

        beta_iter[b*(*nsubject) + j] = outrmvnorm[b];
        btmp[b] = beta_iter[b*(*nsubject) + j];

        // beta_iter[b*(*nsubject) + j] = bhat[b];
        // btmp[b] = bhat[b];

      }

      // RprintVecAsMat("beta_iter", beta_iter, nb, *nsubject);

      // cholesky(HtH, nb, &ld);
      // inverse_from_cholesky(HtH, scr1, scr2, nb);
      // matrix_product(HtH, Hty, bhat, nb, 1, nb);


      matrix_product(H, btmp, Hb, nobs[j], 1, nb);


      sumy_Hb = 0.0;
      for(jj = 0; jj < nobs[j]; jj++){

        sumy_Hb = sumy_Hb + (y_tmp[jj] - Hb[jj]);
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
      //			Rprintf("os = %f\n", os);
      //			Rprintf("ns = %f\n", ns);

      if(ns > 0.0){

        for(jj = 0; jj < nobs[j]; jj++){
          scr3[jj] = y_b0[jj] - Hb[jj];
        }

        sumsq = inner_product(scr3, 1, scr3, 1, nobs[j]);

        llo = -0.5*nobs[j]*log(os*os) - 0.5*(1/(os*os))*sumsq;
        lln = -0.5*nobs[j]*log(ns*ns) - 0.5*(1/(ns*ns))*sumsq;

        //				Rprintf("llo = %f\n", llo);
        //				Rprintf("lln = %f\n", lln);

        llo = llo + dunif(os, 0.0, Asig, 1);
        lln = lln + dunif(ns, 0.0, Asig, 1);

        //				Rprintf("llo = %f\n", llo);
        //				Rprintf("lln = %f\n", lln);


        llr = lln - llo;

        //				Rprintf("llr = %f\n", llr);

        uu = runif(0.0,1.0);
        if(llr > log(uu)) sig2_iter[j] = ns*ns;

      }


      //			sig2_iter[j] = 0.000000000000001;

      //			Rprintf("sig2 = %f\n", sig2_iter[j]);

      ////////////////////////////////////////////
      //									      //
      // update beta0 within in the same loop;  //
      //									      //
      ////////////////////////////////////////////

      s2star = 1.0/((nobs[j]/sig2_iter[j]) + (1/sig2b0_iter));

      mstar = s2star*((1/sig2_iter[j])*sumy_Hb + (1/sig2b0_iter)*mub0_iter);

      //			Rprintf("mstar = %f\n", mstar);
      //			Rprintf("s2star = %f\n", s2star);

      beta0_iter[j] = rnorm(mstar, sqrt(s2star));

      //			Rprintf("beta0 = %f\n", beta0_iter[j]);

      //			beta0_iter[j] = 0.0;


      if((i > (*burn-1)) & (i % (*thin) == 0)){
        like0=0;
        llikeval = 0.0;
        for(t = 0; t < nobs[j]; t++){
          //					Rprintf("t = %d\n", t);
          //					Rprintf("y_tmp[t] = %f\n", y_tmp[t]);
          //					Rprintf("beta0_iter[j] = %f\n", beta0_iter[j]);
          //					Rprintf("Hb[t] = %f\n", Hb[t]);
          //					Rprintf("sqrt(sig2_iter[j]) = %f\n", sqrt(sig2_iter[j]));
          //					Rprintf("dnorm(y_tmp[t], beta0_iter[j] + Hb[t], sqrt(sig2_iter[j]), 0) = %f\n", dnorm(y_tmp[t], beta0_iter[j] + Hb[t], sqrt(sig2_iter[j]), 0));
          llikeval = llikeval + dnorm(y_tmp[t], beta0_iter[j] + Hb[t], sqrt(sig2_iter[j]), 1);

          //					Rprintf("like = %f\n", likeval);
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





    // RprintVecAsMat("sig2_iter", sig2_iter, 1.0, *nsubject);
    // RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);


    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update lam (using MH-step) this is the standard deviation not the variance
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      //			Rprintf("k = %d =================== \n", k);

      olam = lamh[k];
      nlam = rnorm(olam,csigLAM);

      //			Rprintf("olam = %f\n", olam);
      //			Rprintf("nlam = %f\n", nlam);
      //			Rprintf("A = %f\n", A);

      if((nlam > 0) & (nlam < A)){

        for(b = 0; b < nb; b++){

          thtmp[b] = thetah[b*(*nsubject) + k];


          for(bb = 0; bb < nb; bb++){

            oV[b*nb+bb] = 0.0;
            nV[b*nb+bb] = 0.0;

            if(b == bb){

              oV[b*nb+bb] = 1/(olam*olam);
              nV[b*nb+bb] = 1/(nlam*nlam);
              //							oV[b*nb+bb] = 1/(olam);
              //							nV[b*nb+bb] = 1/(nlam);


            }
          }
        }

        ldo = 2.0*nb*log(olam);
        ldn = 2.0*nb*log(nlam);

        //				ldo = nb*log(olam);
        //				ldn = nb*log(nlam);

        //				RprintVecAsMat("thtmp", thtmp, 1, nb);

        //				RprintVecAsMat("oV", oV, nb, nb);
        //				RprintVecAsMat("nV", nV, nb, nb);

        //				Rprintf("ldo = %f\n", ldo);
        //				Rprintf("ldn = %f\n", ldn);

        lln = 0.0;
        llo = 0.0;
        for(j = 0; j < *nsubject; j++){

          if(Si_iter[j] == k+1){

            for(b = 0; b < nb; b++){

              btmp[b] = beta_iter[b*(*nsubject) + j];

            }

            //						RprintVecAsMat("btmp", btmp, 1, nb);
            //						Rprintf("dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1) = %f\n", dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1));

            llo = llo + dmvnorm(btmp, thtmp, oV, nb, ldo, scr1, 1);
            lln = lln + dmvnorm(btmp, thtmp, nV, nb, ldn, scr1, 1);

          }

        }


        llo = llo + dunif(olam, 0.0, A, 1);
        lln = lln + dunif(nlam, 0.0, A, 1);


        //				Rprintf("llo = %f\n", llo);
        //				Rprintf("lln = %f\n", lln);

        llr = lln - llo;
        uu = runif(0.0,1.0);

        //				Rprintf("llr = %f\n", llr);
        //				Rprintf("log(uu) = %f\n", log(uu));

        if(log(uu) < llr) lamh[k] = nlam;
        //				lamh[k] = 1000;

      }

    }




    //		RprintVecAsMat("lamh", lamh, 1, *nsubject);
    //		RprintVecAsMat("lamh", lamh, 1, nclus_iter);




    //		RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);

    //		Rprintf("nclus_iter = %d\n", nclus_iter);
    //		RprintVecAsMat("tau2h", tau2h,1,nclus_iter);
    //		RprintVecAsMat("lamh", lamh,1,nclus_iter);




    //		RprintVecAsMat("beta", beta_iter, nb, *nsubject);
    //		RprintVecAsMat("sig2", sig2_iter, 1, *nsubject);
    //		RprintVecAsMat("beta0", beta0_iter, 1, *nsubject);


    //		mnb0 = 0.0;
    //		for(j=0; j < *nsubject; j++){
    //			mnb0 = mnb0 + (beta0_iter[j]/(*nsubject));
    //		}
    //		for(j=0; j < *nsubject; j++) beta0_iter[j] = beta0_iter[j] - mnb0;
    //		RprintVecAsMat("beta0", beta0_iter, 1, *nsubject);

    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate mub0 a global intercept of mean of beta0i     						//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////
    sumb0 = 0.0;
    for(j = 0; j < *nsubject; j++){

      sumb0 = sumb0 + beta0_iter[j];

    }

    //		Rprintf("sumb0 = %f\n", sumb0);

    s2star = 1.0/((*nsubject)/sig2b0_iter + (1/s2b0));

    mstar = s2star*((1/sig2b0_iter)*(sumb0) + (1/s2b0)*mb0);

    //		Rprintf("mstar = %f\n", mstar);
    //		Rprintf("s2star = %f\n", s2star);

    mub0_iter = rnorm(mstar, sqrt(s2star));


    //		Rprintf("mub0 = %f\n", mub0);
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


    //		Rprintf("sig2b0 = %f\n", sig2b0);

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

          Sstar[b*nb+bb] = (1/tau2h[k])*K[b*nb+bb];

          if(b == bb) Sstar[b*nb+bb] = ((double) nh[k]/(lamh[k]*lamh[k])) +
            (1/tau2h[k])*K[b*nb+bb];
          //				 if(b == bb) Sstar[b*nb+bb] = (nh[k]/lamh[k]) + (1/tau2h[k])*K[b*nb+bb];

        }

        sumbeta[b] = 0.0;

      }

      //			RprintVecAsMat("Sstar",Sstar, nb, nb);

      cholesky(Sstar, nb, &ld);
      inverse_from_cholesky(Sstar, scr1, scr2, nb);

      //			RprintVecAsMat("Sstar",Sstar, nb, nb);
      //			RprintVecAsMat("mu_iter = ", mu_iter, 1, nb);

      matrix_product(K, mu_iter, scr1, nb, 1, nb);

      //			RprintVecAsMat("Kmu",scr1, 1, nb);
      //			RprintIVecAsMat("Si_iter", Si_iter, 1, *nsubject);

      for(j = 0; j < *nsubject; j++){

        if(Si_iter[j] == k+1){

          for(b = 0; b < nb; b++){

            sumbeta[b] = sumbeta[b] + (1/(lamh[k]*lamh[k]))*
              beta_iter[b*(*nsubject) + j];
          }
        }
      }

      for(b=0; b<nb; b++) sumbeta[b] = sumbeta[b] + (1/tau2h[k])*scr1[b];


      //			RprintVecAsMat("sumbeta", sumbeta, 1, nb);


      matrix_product(Sstar, sumbeta, Mstar, nb, 1, nb);

      //			RprintVecAsMat("Mstar", Mstar, 1, nb);
      //			RprintVecAsMat("Sstar", Sstar, nb, nb);


      cholesky(Sstar, nb , &ld);



      ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);


//      RprintVecAsMat("thetah", outrmvnorm, 1, nb);

      for(b = 0; b < nb; b++){

        thetah[b*(*nsubject) + k] = outrmvnorm[b];


      }


    }


    //		RprintVecAsMat("thetah", thetah, nb, *nsubject);
    //		RprintVecAsMat("thetah", thetah, nb, nclus_iter);



    //		RprintVecAsMat("thetah", thetah, nb, *nsubject);
    //		RprintVecAsMat("thetah", thetah, nb, nclus_iter);
    //		RprintVecAsMat("thetah", thetah, 1, nb*(*nsubject));
    //////////////////////////////////////////////////////////////////////////////////
    //
    // Update tau2 for each of the clusters (P-spline smoothing parameter)
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < nclus_iter; k++){

      //			Rprintf("k = %d =================== \n", k);

      for(b = 0; b < nb; b++){

        thtmp[b] = thetah[b*(*nsubject) + k] - mu_iter[b];

      }

      //			RprintVecAsMat("thtmp", thtmp, 1, nb);
      //			RprintVecAsMat("mu_iter", mu_iter, 1, nb);

      sumsq = quform(thtmp,K,nb);

      //			Rprintf("sumsq = %f\n", sumsq);

      astar = 0.5*(nb) + at;
      bstar = 1/bt + 0.5*sumsq;

      //			Rprintf("astar = %f\n", astar);
      //			Rprintf("bstar = %f\n", bstar);

      tau2h[k] = 1/rgamma(astar, 1/bstar);// E(tau2) = astarbstar for gamma.  bstar is scale

      //			Rprintf("tau2h = %f\n", tau2h[k]);

      //			tau2h[k] = 0.01;

    }


    //		RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);



    //		RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);
    //		Rprintf("nclus_iter = %d\n", nclus_iter);







    //////////////////////////////////////////////////////////////////////////////////
    //																				//
    // udpate mu center of cluster specific thetah;									//
    //																				//
    //////////////////////////////////////////////////////////////////////////////////

    //		Rprintf("nb = %d\n", nb);
    //		RprintIVecAsMat("nc", nc, 1, nclus_iter);
    //		RprintVecAsMat("tau2h", tau2h, 1, nclus_iter);
    for(b = 0; b < nb; b++) sumtheta[b] = 0.0;
    sumtau2 = 0.0;

    for(k = 0; k < nclus_iter; k++){
      sumtau2 = sumtau2 + (1/tau2h[k]);
      for(b = 0; b < nb; b++){
        //				Rprintf("thetah[b*(*nsubject) + k] = %f\n", thetah[b*(*nsubject) + k]);
        sumtheta[b] = sumtheta[b] + (1/tau2h[k])*thetah[b*(*nsubject) + k];
      }
    }

    //		Rprintf("sumtau2 = %f\n", sumtau2);
    //		RprintVecAsMat("sumtheta",sumtheta, 1, nb);


    for(b = 0; b < nb; b++){
      for(bb = 0; bb < nb; bb++){

        Sstar[b*nb+bb] = sumtau2*K[b*nb+bb];

        if(b == bb) Sstar[b*nb+bb] = sumtau2*K[b*nb+bb]+(1/s2mu) ;

      }

    }

    //		RprintVecAsMat("sumtheta",sumtheta, 1, nb);
    //		RprintVecAsMat("Sstar",Sstar, nb, nb);

    matrix_product(K, sumtheta, scr3, nb, 1, nb);

    //		RprintVecAsMat("sumthetaK",scr3, 1, nb);

    cholesky(Sstar, nb, &ld);
    inverse_from_cholesky(Sstar, scr1, scr2, nb);

    //		RprintVecAsMat("Sstar",Sstar, nb, nb);


    matrix_product(Sstar, scr3, Mstar, nb, 1, nb);

    //		RprintVecAsMat("Mstar", Mstar, 1, nb);
    //		RprintVecAsMat("Sstar", Sstar, nb, nb);

    cholesky(Sstar, nb , &ld);

    ran_mvnorm(Mstar, Sstar, nb, scr1, outrmvnorm);


    //		RprintVecAsMat("mu", outrmvnorm, 1, nb);

    for(b = 0; b < nb; b++){
       mu_iter[b] = outrmvnorm[b];
//      mu_iter[b] = 0;
    }

    //		RprintVecAsMat("thetah", thetah, nb, *nsubject);



    //		RprintVecAsMat("mu_iter", mu_iter, 1, nb);

/*

    //////////////////////////////////////////////////////////////////////////////////
    //
    // Posterior predictives.
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

    if((i > (*burn-1)) & (i % (*thin) == 0)){


      nclus[ii] = nclus_iter;


      mub0[ii] = mub0_iter;

      sig2b0[ii] = sig2b0_iter;


      for(b = 0; b < nb; b++){
        mu[b*(nb) + ii] = mu_iter[b];
      }



      for(j = 0; j < *nsubject; j++){

        Si[ii*(*nsubject) + j] = Si_iter[j];


        sig2[ii*(*nsubject) + j] = sig2_iter[j];

        beta0[ii*(*nsubject) + j] = beta0_iter[j];

        llike[ii*(*nsubject) + j] = llike_iter[j];

        lam[ii*(*nsubject) + j] = lamh[Si_iter[j]-1];

        tau2[ii*(*nsubject) + j] = tau2h[Si_iter[j]-1];


        for(b = 0; b < nb; b++){
          beta[(ii*(nb) + b)*(*nsubject) + j] = beta_iter[b*(*nsubject) + j];

          theta[(ii*(nb) + b)*(*nsubject) + j] = thetah[b*(*nsubject) + Si_iter[j]-1];

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
