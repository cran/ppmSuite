/****************************************************************************************
 * Copyright (c) 2018 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit Ordinal data model with PPMx prior.
 *
 * letting rho denote a partition of spatial locations we
 * have the following model
 *
 *	Pr(rho) \propto \prod c(S_j)g(x*_j) with x_i p-dimsional

 *	Y_i | mu*_j, c_i, sig2*_j ~ N(mu*_{c_i}, sig2*_{c_i})
 *
 *	Priors:
 *
 *  mu*_j ~ N(mu0, sig20) with mu0 ~ N(m, s^2) and sig0 ~ UN(0,a1)
 *
 *  sig*_j ~ UN(0, m)
 *
 ****************************************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

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
* Cvec = ncat x 1 vector indicating the number of categories for each categorical covariate.
* PPM = logical indicating if PPM or PPMx should be used.
* cohesion = integer indicating which cohesion function to use
*   1 - DP stype
*   2 - uniform
* similarity_function = int indicating similarity function to be used
		1 - auxiliary model
		2 - double dipper
		3 - a*exp(-variance)
		4 - exp(-sum_ijdissimilarity)
* consim = 1 or 2.  1 implies sim for con var is N-N.  2 implies sim is N-NIG
* M = double indicating value of M associated with cohesion.
* y = nobs x 1 vector that contains ordinal response values
* Xcon = nobs x ncon contiguous double vector that contains continuous covariate values
* Xcat = nobs x ncat contiguous int vector that contains categorical covariate values
* npred = integer indicating number of out of sample predictions
* Xconp = npred x ncon matrix of continuous covariates to make predictions
* Xcatp = npred x ncat matrix of categorical covariates to make predictions
* simparms = vector containing similarity functions that are ordered in the following way;
*           c(m0, s20, v, k0, nu0, a0, alpha)
* dissimtn = nobs x nobs pairwise dissimilarity matrix using gower's distance.
* dissimtt = nobs x nobs dissimilarity between predictors and model covariate
* calibrate = integer that determines if similarity functions are calibrated.
* totgowdist = integer that determines if use the average (0) or total (1) pairwise gower distance .
* co = vector containing fixed cut-off points for associated with latent variable
* modelPriors = vector containing model parameter priors
* 		  m,  s2,    smin, smax, s0min, s0max,
*       c(0, 10*10,  0,     0.1,    0,     10)
*
* Output:
*
*****************************************************************************************/

void ordinal_ppmx(int *draws, int *burn, int *thin, int *nobs, int *ncon, int *ncat,
                    int *Cvec, int *PPM, int *cohesion, int *similarity_function, int *consim,
                    double *M,
                    int *y, double *Xcon, int *Xcat,
                    int *npred, int *nordcat, double *Xconp, int *Xcatp,
                    double *simParms, double *dissimtn, double *dissimtt,
                    int *calibrate, double *co, double *modelPriors, int *verbose, double *mh,
                    double *mu, double *sig2, double *mu0, double *sig20,
                    int *Si, int *nclus, double *zi, double *like, double *WAIC, double *lpml,
                    double *ispred, int *isordpred, double *ppred, int *predclass, int *ordppred,
                    double *rbpred, int *rbordpred, double *predclass_prob){


	// i - MCMC index
	// j - individual index
	// jj - second individual index (for double for loops)
	// jjj - third individual index
	// c - categorical variable index
	// p - number of covariates index
	// pp - second covariate index
	// k - cluster index
	// t - subset of covariates index

	int i, j, jj, jjj, c, p, pp, k, t, ii;

	double max_C, nout, sumx, sumx2;
	nout = (*draws-*burn)/(*thin);
	int ncov = *ncon + *ncat;

	max_C = Cvec[0];

	for(p = 0; p < (*ncat); p++){
		if(max_C < Cvec[p]) max_C = Cvec[p];
	}
	if(*ncat == 0) max_C = 1.0;

  // MLEs needed for a specific similarity function
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

	// =============================================================================================
	//
	// Memory vectors to hold a single MCMC iterate
	//
	// =============================================================================================

	double mu0_iter = 0.0;
	double sig20_iter =1.0;

	int nclus_iter = 0;

	int iaux;

	int Si_iter[*nobs];
	int nh[*nobs];
	int nhc[(*nobs)*(*ncat)];

	double* zi_iter = R_VectorInit(*nobs,0.0);

	double* xcontmp = R_Vector(*nobs);


  // To begin assign everybody to the same cluster
	for(j = 0; j < *nobs; j++){
		Si_iter[j] = 1;
		nh[j] = 0;
	}
  // Create vector of cluster sizes
	for(j = 0; j < *nobs; j++){
		for(jj = 0; jj < *nobs; jj++){
			if(Si_iter[j] == jj+1) nh[jj] = nh[jj] + 1;
		}
	}
  // Count the number of clusters
	for(j = 0; j < *nobs; j++){
		if(nh[j] > 0) nclus_iter = nclus_iter + 1;
	}


	double *sig2h = R_VectorInit(*nobs, 0.05);
	double *muh = R_VectorInit(*nobs, 0.0);

	// Stuff to compute the posterior predictive
    double *ispred_iter = R_VectorInit(*nobs, 0.0);
    int isordpred_iter[*nobs];

    double *ppred_iter = R_Vector((*npred)*(*nobs));
    double *rbpred_iter = R_Vector((*npred)*(*nobs));
    int ordppred_iter[*npred];
    int rbordpred_iter[*npred];

    int predclass_iter[*nobs];

    double *predclass_prob_iter = R_Vector((*npred)*(*nobs));

	// ===================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// ===================================================================================
	// stuff that I need to update Si (cluster labels);
	int nhtmp;
	double auxm, auxs2, tmp, sumdis,npdN,npdY,npd;
	double mudraw, sdraw, maxph, denph, cprobh, uu;
	double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw;
	double lgcont,lgcatt;

	double *ph = R_VectorInit(*nobs, 0.0);
	double *probh = R_VectorInit(*nobs, 0.0);

	double *gtilN = R_VectorInit((*nobs+1),0.0);
	double *gtilY = R_VectorInit((*nobs+1),0.0);
	double *lgtilN = R_VectorInit((*nobs+1),0.0);
	double *lgtilY = R_VectorInit((*nobs+1),0.0);
	double sgY, sgN,  lgtilNk, lgtilYk, maxgtilY, maxgtilN;


	// stuff I need to update z (latent parameter)
	double mntmp, s2tmp;

	//stuff I need for muh
	double mstar,s2star, sumz;

	//suff I need to upate mu0
	double summu;

	//stuff I need to update sig20
	double llo, lln, llr,  os0,ns0;

	//stuff I need to update sig2
	double osig,nsig;

	// Stuff for out of sample predictions
	double lgcon0, lgcat0, mupred, sig2pred;

	// Stuff to compute lpml, likelihood, WAIC, and Rao-Blackwellized density values
	double lpml_iter, elppdWAIC, sdens;
	double *CPOinv = R_VectorInit(*nobs, 0.0);
	double *like_iter = R_VectorInit(*nobs, 0.0);
	double *mnlike = R_VectorInit(*nobs, 0.0);
	double *mnllike = R_VectorInit(*nobs, 0.0);

    // priors for mu0
    double m = modelPriors[0]; double s2 = modelPriors[1];

    // prior values for sig2h and sig20;
    double smin=0, smax=modelPriors[2];
    double s0min=0, s0max=modelPriors[3];


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

    // M-H tuning parameters
    double csigSIG0=mh[0], csigSIG=mh[1];

  if(*verbose){
	  Rprintf("nobs = %d\n", *nobs);
	  Rprintf("ncon = %d\n", *ncon);
	  Rprintf("ncat = %d\n", *ncat);
	  Rprintf("ncov = %d\n", ncov);
	  Rprintf("npred = %d\n", *npred);
	  RprintVecAsMat("co", co, 1, *nordcat);
	  if(*ncat > 0 ) RprintIVecAsMat("Cvec",Cvec, 1, *ncat);
	  Rprintf("Prior values: m = %.2f, s2 = %.2f smin = %.2f, smax = %.2f, s0min = %.2f, s0max = %.2f\n", m, s2, smin, smax, s0min, s0max);
	  Rprintf("M = %f\n", Mdp);
	  Rprintf("Simlarity values: m0 = %.2f, s20 = %.2f v = %.2f, k0 = %.2f, nu0 = %.2f, a0 = %.2f, alpha = %.2f\n", m0, s20, v, k0, nu0, simParms[5], alpha);
	  if(*nobs > 50) RprintIVecAsMat("First fifty response values",y, 1, 50);
  }



	GetRNGstate();

	// =============================================================================================
	//
	// start of the mcmc algorithm;
	//
	// =============================================================================================
	ii = 0;

	for(i = 0; i < *draws; i++){

		if(*verbose & ((i+1) % 100 == 0)){
			time_t now;
			time(&now);

			Rprintf("mcmc iter = %d ===================================================== \n", i+1);
			Rprintf("%s", ctime(&now));
		}

		//////////////////////////////////////////////////////////////////////////////////
		//
		// update the the latent variables z_i using a truncated normal
		//
		//////////////////////////////////////////////////////////////////////////////////

		for(j = 0; j < *nobs; j++){

			mntmp = muh[Si_iter[j]-1];
			s2tmp = sig2h[Si_iter[j]-1];

			zi_iter[j] = rtnorm(mntmp, sqrt(s2tmp), co[y[j]], co[y[j]+1]);

		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// update the cluster labels using the polya urn scheme of
		// algorithm 8 found in  Radford Neal's
		//	"Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
		//	paper.
		//
		//////////////////////////////////////////////////////////////////////////////////

		for(j = 0; j < *nobs; j++){

			if(nh[Si_iter[j]-1] > 1){

				// Observation belongs to a non-singleton ...
				nh[Si_iter[j]-1] = nh[Si_iter[j]-1] - 1;

			}else{

				// Observation is a member of a singleton cluster ...

				iaux = Si_iter[j];

				if(iaux < nclus_iter){

					// Need to relabel clusters.  I will do this by swapping cluster labels
					// Si_iter[j] and nclus_iter along with cluster specific parameters;


					// All members of last cluster will be assigned subject i's cluster label
					for(jj = 0; jj < *nobs; jj++){

						if(Si_iter[jj] == nclus_iter){

							Si_iter[jj] = iaux;

						}

					}


					Si_iter[j] = nclus_iter;

					// The following steps swaps order of cluster specific parameters
					// so that the newly labeled subjects from previous step retain
					// their correct cluster specific parameters

					auxm = muh[iaux-1];
					muh[iaux-1] = muh[nclus_iter-1];
					muh[nclus_iter-1] = auxm;

					auxs2 = sig2h[iaux-1];
					sig2h[iaux-1] = sig2h[nclus_iter-1];
					sig2h[nclus_iter-1] = auxs2;


					// the number of members in cluster is also swapped with the last
					nh[iaux-1] = nh[nclus_iter-1];
					nh[nclus_iter-1] = 1;

				}


				// Now remove the ith obs and last cluster;
				nh[nclus_iter-1] = nh[nclus_iter-1] - 1;
				nclus_iter = nclus_iter - 1;


			}

			// The atoms have been relabeled if necessary and now we need to update Si.

			// Begin the cluster probabilities

			for(k=0; k<nclus_iter; k++){

				lgconY = 0.0;
				lgconN = 0.0;

				for(p=0; p<(*ncon); p++){
					nhtmp = 0;
					sumx = 0.0;
					sumx2 = 0.0;
					sumdis = 0.0;
					for(jj = 0; jj < *nobs; jj++){
						if(jj != j){
							if(Si_iter[jj] == k+1){
								tmp = Xcon[jj*(*ncon)+p];

								sumx = sumx + tmp;
								sumx2 = sumx2 + tmp*tmp;

								nhtmp = nhtmp+1;
							}
						}
					}

					if(*similarity_function==1){ // Auxilliary
						if(*consim==1){ // NN
							lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
							lgconN = lgconN + lgcont;
						}
						if(*consim==2){// NNIG
							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
							lgconN = lgconN + lgcont;
						}

					}
					if(*similarity_function==2){ //Double Dipper
						if(*consim==1){
							lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
							lgconN = lgconN + lgcont;
						}
						if(*consim==2){
							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
							lgconN = lgconN + lgcont;
						}
					}
					if(*similarity_function==3){ // variance
						lgcont = gsimconEV(sumx, sumx2, nhtmp, alpha,1);
						lgconN = lgconN + lgcont;
					}


					// now add jth individual back;

					xcontmp[nhtmp] = Xcon[j*(*ncon)+p];
					sumx = sumx + Xcon[j*(*ncon)+p];
					sumx2 = sumx2 + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
					nhtmp = nhtmp+1;

					if(*similarity_function==1){ // Auxilliary
						if(*consim==1){
							lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
							lgconY = lgconY + lgcont;
						}
						if(*consim==2){
							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
							lgconY = lgconY + lgcont;
						}
					}
					if(*similarity_function==2){ //Double Dipper
						if(*consim==1){
							lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
							lgconY = lgconY + lgcont;
						}
						if(*consim==2){
							lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
							lgconY = lgconY + lgcont;
						}
					}
					if(*similarity_function==3){ // variance
						lgcont = gsimconEV(sumx, sumx2, nhtmp, alpha,1);
						lgconY = lgconY + lgcont;
					}

				}


        // Now calculate similarity for the categorical covariates
				lgcatY=0.0;
				lgcatN=0.0;
				for(p=0; p<(*ncat); p++){
					for(c=0;c<Cvec[p];c++){nhc[c]=0;}

					nhtmp = 0;
					for(jj = 0; jj < *nobs; jj++){
						if(jj != j){

							if(Si_iter[jj]==k+1){
								nhc[Xcat[jj*(*ncat)+p]] = nhc[Xcat[jj*(*ncat)+p]] + 1; // this needs to be a vector
								nhtmp = nhtmp+1;

							}
						}
					}


					if(*similarity_function==1){ // Auxiliary
						lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
						lgcatN = lgcatN + lgcatt;
					}
					if(*similarity_function==2){// Double dipper
						lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
						lgcatN = lgcatN + lgcatt;
					}
					if(*similarity_function==3){// Using Entropy instead of variance here
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


					nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1;
					nhtmp = nhtmp + 1;

					if(*similarity_function==1){
						lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
						lgcatY = lgcatY + lgcatt;
					}
					if(*similarity_function==2){
						lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
						lgcatY = lgcatY + lgcatt;
					}
					if(*similarity_function==3){// Using Entropy instead of variance here
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

						if((Si_iter[jj] == k+1) & (jj != j)){
							lgconY = lgconY + dissimtn[jj*(*nobs) + j];
							for(jjj = 0; jjj < jj; jjj++){
								if((Si_iter[jjj] == k+1) & (jjj != j)){
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
					// lgconN = -(alpha)*lgconN/(npdN);
					// lgconY = -(alpha)*lgconY/(npdY);

					// Just use cluster-total gower dissimilarity
					lgconN = -(alpha)*lgconN;
					lgconY = -(alpha)*lgconY;

				}
				//////////////////////////////////////////////////////////
				// End of Gower similarity
				//////////////////////////////////////////////////////////


        // Compute the unnormalized cluster probabilities
				ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
				        	log((double) nh[k]) +
				        	lgcatY - lgcatN +
							    lgconY - lgconN;


				if(*PPM){
					ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
						    	log((double) nh[k]);  // DP part of cohesion function

				}

				if(*calibrate == 2){
					ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
				            log((double) nh[k]) +
				            (1/((double)*ncon + (double)*ncat))*(lgcatY + lgconY - lgcatN - lgconN);
				}
				// Uniform cohesion
				if(*cohesion==2){
					ph[k] = ph[k] - log((double) nh[k]);
				}


			}

      // Need to consider allocating subject to new cluster
			mudraw = rnorm(mu0_iter, sqrt(sig20_iter));
			sdraw = runif(smin, smax);

			lgcondraw = 0.0;
			for(p=0;p<(*ncon);p++){
				xcontmp[0] = Xcon[j*(*ncon)+p];
				if(*similarity_function==1){ // Auxilliary
					if(*consim==1){
						lgcont = gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0, 1);
						lgcondraw = lgcondraw + lgcont;
					}
					if(*consim==2){
						lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 0,0, 1);
						lgcondraw = lgcondraw + lgcont;
					}
				}
				if(*similarity_function==2){ // Double Dipper
					if(*consim==1){
						lgcont = gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p], 1, 1, 0, 1);
						lgcondraw = lgcondraw + lgcont;
					}
					if(*consim==2){
						lgcont = gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 1, 0, 1);
						lgcondraw = lgcondraw + lgcont;
					}
				}
				if(*similarity_function==3){ // Variance
					lgcont = gsimconEV(xcontmp[0], xcontmp[0]*xcontmp[0], 1,alpha,1);
					lgcondraw = lgcondraw + lgcont;
				}
			}
			if(*similarity_function==4){ // Dissimilarity
				lgcondraw = -(alpha)*0;
			}

			lgcatdraw = 0.0;
			for(p=0;p<(*ncat);p++){
				for(c=0;c<Cvec[p];c++){nhc[c] = 0;}

				nhc[Xcat[j*(*ncat)+p]] = 1;


				if(*similarity_function==1){
					lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
					lgcatdraw = lgcatdraw + lgcatt;
				}
				if(*similarity_function==2){
					lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
					lgcatdraw = lgcatdraw + lgcatt;
				}
				if(*similarity_function==3){
					lgcatdraw = lgcatdraw + -(alpha)*0;
				}

			}

			if(*similarity_function==4){
				lgcatdraw = -(alpha)*0;
			}

			gtilY[nclus_iter] = lgcondraw + lgcatdraw;
			gtilN[nclus_iter] = lgcondraw + lgcatdraw;

			ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
			                 	log(Mdp) +
			                 	lgcondraw +
			                 	lgcatdraw;


			if(*PPM){
				ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
				                 log(Mdp); //DP part
			}
			if(*cohesion==2){
				ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
			                 	lgcondraw +
			                 	lgcatdraw;
			}


			if(*calibrate==2){
				ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
				                 log(Mdp) +
				                 (1/((double)*ncon + (double)*ncat))*(lgcondraw + lgcatdraw);

			  if(*cohesion==2){
				  ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
			                 	(1/((double)*ncon + (double)*ncat))*(lgcondraw + lgcatdraw);
			  }

			}


			/////////////////////////////////////////////////////////////////////////////
			// This is the calibration used when the similarity is normalized
			/////////////////////////////////////////////////////////////////////////////
			if((*calibrate==1) & (*PPM != 1)){
				maxgtilN = gtilN[0];
				maxgtilY = gtilY[0];
				for(k=1; k < nclus_iter+1; k++){

					if(maxgtilN < gtilN[k]) maxgtilN = gtilN[k];

					if(k < nclus_iter){
						if(maxgtilY < gtilY[k]) maxgtilY = gtilY[k];
					}
				}

				sgY=0.0;
				sgN=0.0;
				for(k=0; k<nclus_iter+1; k++){

					lgtilN[k] = gtilN[k] - maxgtilN;
					sgN = sgN + exp(lgtilN[k]);

					if(k < nclus_iter){// If x is included in an existing cluster in cannot be a singleton
						lgtilY[k] = gtilY[k] - maxgtilY;
						sgY = sgY + exp(lgtilY[k]);
					}
				}

				// Calibrate the unnormalized cluster probabilities
				for(k=0; k<nclus_iter; k++){
					lgtilNk = lgtilN[k] - log(sgN);
					lgtilYk = lgtilY[k] - log(sgY);

					ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
			        	    log((double) nh[k]) +  // Cohesion part
							      lgtilYk - lgtilNk; //This takes into account both cont and cat vars



				}

				// calibration for a singleton
				ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
			                 	    log(Mdp) +
								            lgtilN[nclus_iter] - log(sgN);

			}
			/////////////////////////////////////////////////////////////////////////////
			// End of calibration used when the similarity is normalized
			/////////////////////////////////////////////////////////////////////////////


			maxph = ph[0];
			for(k = 1; k < nclus_iter+1; k++){
				if(maxph < ph[k]) maxph = ph[k];
			}

			denph = 0.0;
			for(k = 0; k < nclus_iter+1; k++){
				ph[k] = exp(ph[k] - maxph);
				denph = denph + ph[k];
			}

			for(k = 0; k < nclus_iter+1; k++){
				probh[k] = ph[k]/denph;
			}

			uu = runif(0.0,1.0);

			cprobh= 0.0;;
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

				muh[Si_iter[j]-1] = mudraw;
				sig2h[Si_iter[j]-1] = sdraw*sdraw;
			}

			// Compute the CPO and lpml using the mixture
			sdens = 0.0;
			for(k = 0; k < nclus_iter; k++){
				s2tmp = sig2h[k];
				mntmp = muh[k];
				sdens = sdens + probh[k]*(pnorm(co[y[j]+1], mntmp, sqrt(s2tmp), 1, 0) -
							              pnorm(co[y[j]],   mntmp, sqrt(s2tmp), 1, 0));
			}
			sdens = sdens + probh[nclus_iter]*(pnorm(co[y[j]+1], mudraw, sdraw, 1, 0) -
								               pnorm(co[y[j]],   mudraw, sdraw, 1, 0));

			like_iter[j] = sdens;

			if((i > (*burn-1)) & (i % (*thin) == 0)){

				// These are needed for WAIC
				mnlike[j] = mnlike[j] + (like_iter[j])/(double) nout;
				mnllike[j] = mnllike[j] + log(like_iter[j])/(double) nout;

				CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/like_iter[j]);
			}
		}



		//////////////////////////////////////////////////////////////////////////////////
		//
		// update muh's cluster specific means
		//
		//////////////////////////////////////////////////////////////////////////////////

		for(k = 0; k < nclus_iter; k++){
			sumz = 0.0;
			for(j = 0; j < *nobs; j++){
				if(Si_iter[j] == k+1){
					sumz = sumz + zi_iter[j];
				}
			}

			s2star = 1/((double) nh[k]/sig2h[k] + 1/sig20_iter);
			mstar = s2star*( (1/sig2h[k])*sumz + (1/sig20_iter)*mu0_iter);

			muh[k] = rnorm(mstar, sqrt(s2star));
		}


		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update mu0  prior mean of muh
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		summu = 0.0;
		for(k = 0; k < nclus_iter; k++){
			summu = summu + muh[k];
		}

		s2star = 1/(((double) nclus_iter/sig20_iter) + (1/s2));
		mstar = s2star*((1/sig20_iter)*summu + (1/s2)*m);

		mu0_iter = rnorm(mstar, sqrt(s2star));

		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update sig20  prior variance of muh
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		os0 = sqrt(sig20_iter);
		ns0 = rnorm(os0,csigSIG0);

		if(ns0 > 0){

			lln = 0.0;
			llo = 0.0;
			for(k = 0; k < nclus_iter; k++){

				llo = llo + dnorm(muh[k], mu0_iter, os0,1);
				lln = lln + dnorm(muh[k], mu0_iter, ns0,1);
			}

			llo = llo + dunif(os0, s0min, s0max, 1);
			lln = lln + dunif(ns0, s0min, s0max, 1);

			llr = lln - llo;
			uu = runif(0,1);

			if(log(uu) < llr){
				sig20_iter = ns0*ns0;
			}

		}


		//////////////////////////////////////////////////////////////////////////////////
		//
		// Update sig2h  cluster specific variance parameters with metropolis and
		// Uniform prior.
		//
		//////////////////////////////////////////////////////////////////////////////////
		for(k = 0; k < nclus_iter; k++){
			osig = sqrt(sig2h[k]);
			nsig = rnorm(osig,csigSIG);

			if((nsig > 0) & (nsig < smax)){

				lln = 0.0;
				llo = 0.0;
				for(j = 0; j < *nobs; j++){
					if(Si_iter[j] == k+1){
						llo = llo + dnorm(zi_iter[j], muh[k], osig,1);
						lln = lln + dnorm(zi_iter[j], muh[k], nsig,1);
					}
				}

				llo = llo + dunif(osig, smin, smax, 1);
				lln = lln + dunif(nsig, smin, smax, 1);

				llr = lln - llo;
				uu = runif(0,1);

				if(log(uu) < llr){
					sig2h[k] = nsig*nsig;
				}
			}
		}



		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// in sample prediction to assess model fit
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & (i % (*thin) == 0)){

			for(j = 0; j < *nobs; j++){
				ispred_iter[j] = rnorm(muh[Si_iter[j]-1], sqrt(sig2h[Si_iter[j]-1]));

            for(c=0; c < *nordcat-1; c++){
              if((ispred_iter[j] > co[c]) & (ispred_iter[j] < co[c+1])) isordpred_iter[j] = c;
            }
          }
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// out of sample prediction using posterior predictive?
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & (i % (*thin) == 0)){
			for(pp = 0; pp < *npred; pp++){

				for(k = 0; k < nclus_iter; k++){

					lgconN=0.0, lgconY=0.0;
					for(p=0; p<(*ncon); p++){
						nhtmp = 0;
						for(j = 0; j < *nobs; j++){
							if(Si_iter[j] == k+1){
								xcontmp[nhtmp] = Xcon[j*(*ncon)+p]; //create cluster specific x-vector
								nhtmp = nhtmp+1;
							}
						}

						sumx = 0.0;
						sumx2 = 0.0;
						for(t = 0; t < nhtmp; t++){
							sumx = sumx + xcontmp[t];
							sumx2 = sumx2 + xcontmp[t]*xcontmp[t];
						}

						if(*similarity_function==1){
							if(*consim==1){
								lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
								lgconN = lgconN + lgcont;
							}
							if(*consim==2){
								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
								lgconN = lgconN + lgcont;
							}
						}
						if(*similarity_function==2){
							if(*consim==1){
								lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
								lgconN = lgconN + lgcont;
							}
							if(*consim==2){
								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
								lgconN = lgconN + lgcont;
							}
						}
						if(*similarity_function==3){
							lgcont = gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
							lgconN = lgconN + lgcont;
						}

					  // now add ppth prediction to cluster;
						xcontmp[nhtmp] = Xconp[pp*(*ncon)+p];
						sumx = sumx + Xconp[pp*(*ncon)+p];
						sumx2 = sumx2 + Xconp[pp*(*ncon)+p]*Xconp[pp*(*ncon)+p];
						nhtmp = nhtmp + 1;

						if(*similarity_function==1){ // Auxilliary
							if(*consim==1){
								lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
								lgconY = lgconY + lgcont;
							}
							if(*consim==2){
								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
								lgconY = lgconY + lgcont;
							}
						}
						if(*similarity_function==2){ // Double Dipper
							if(*consim==1){
								lgcont = gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
								lgconY = lgconY + lgcont;
							}
							if(*consim==2){
								lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
								lgconY = lgconY + lgcont;
							}
						}
						if(*similarity_function==3){ // Variance
							lgcont = gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
							lgconY = lgconY + lgcont;
						}
					}


					lgcatY=0.0, lgcatN=0.0;
					for(p=0; p<(*ncat); p++){
						for(c=0;c<Cvec[p];c++){nhc[c]=0;}
						nhtmp=0;
						for(j = 0; j < *nobs; j++){

							if(Si_iter[j]==k+1){
								nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1; // this needs to be a vector
								nhtmp = nhtmp+1;
							}
						}

						if(*similarity_function==1){
							lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
							lgcatN = lgcatN + lgcatt;
						}
						if(*similarity_function==2){
							lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
							lgcatN = lgcatN + lgcatt;
						}
						if(*similarity_function==3){
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

						nhc[Xcatp[pp*(*ncat)+p]] = nhc[Xcatp[pp*(*ncat)+p]] + 1;
						nhtmp=nhtmp + 1;

						if(*similarity_function==1){
							lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
							lgcatY = lgcatY + lgcatt;
						}
						if(*similarity_function==2){
							lgcatt = gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
							lgcatY = lgcatY + lgcatt;
						}
						if(*similarity_function==3){// Use entropy
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
					}

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
							if(Si_iter[j] == k+1){

								lgconY = lgconY + dissimtt[pp*(*nobs) + j];

								for(jj = 0; jj < j; jj++){
									if(Si_iter[jj] == k+1){

										lgconN = lgconN + dissimtn[j*(*nobs) + jj];
										lgconY = lgconY + dissimtn[j*(*nobs) + jj];
									}
								}
							}
						}

						npdN = nh[k]*(nh[k]-1)/2;
						if(npdN==0)npdN=1;
						npdY = (nh[k]+1)*(nh[k])/2;

            // This is cluster-mean Gower dissimilarity
						// lgconN = -(alpha)*lgconN/(npdN);
						// lgconY = -(alpha)*lgconY/(npdY);

						  // Use the cluster-total Gower dissimilarity
							lgconN = -(alpha)*lgconN;
							lgconY = -(alpha)*lgconY;
					}
			    //////////////////////////////////////////////////////////
				  // End of Gower similarity values for gower dissimilarity
				  //////////////////////////////////////////////////////////



					ph[k] = log((double) nh[k]) +
				         	lgcatY - lgcatN +
							    lgconY - lgconN;

					if(*PPM) ph[k] = log((double) nh[k]);


					if(*calibrate == 2){
						ph[k] = log((double) nh[k]) +
						           (1/((double)*ncon + (double)*ncat))*
						            (lgcatY + lgconY - lgcatN - lgconN);
					}

					if(*cohesion==2) ph[k] = ph[k] - log((double) nh[k]);
				}


				lgcon0=0.0;
				for(p=0;p<*ncon;p++){ // 0 means that data are not missing
					xcontmp[0] = Xconp[pp*(*ncon)+p];
					if(*similarity_function==1){
						if(*consim==1){
							lgcon0 = lgcon0 + gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0,1);
						}
						if(*consim==2){
							lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,0,1);
						}
					}
					if(*similarity_function==2){
						if(*consim==1){
							lgcon0 = lgcon0 + gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,0,1);
						}
						if(*consim==2){
							lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 0,1);
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

					if(*similarity_function==1){
						lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
					}
					if(*similarity_function==2){
						lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
					}
					if(*similarity_function==3){
						lgcat0 = lgcat0 + -(alpha)*0;

					}


				}


				gtilY[nclus_iter] = lgcat0 + lgcon0;
				gtilN[nclus_iter] = lgcat0 + lgcon0;


        // Unormalized predictive probabilities that now include singleton
				ph[nclus_iter] = log((double) Mdp) + lgcon0 + lgcat0;

				if(*similarity_function==4) ph[nclus_iter] = log((double) Mdp) + log(1);

				if(*PPM) ph[nclus_iter] = log((double) Mdp);

				if(*cohesion==2) ph[nclus_iter] = lgcon0 + lgcat0;

				if(*calibrate==2){
					ph[nclus_iter] = log(Mdp) +
									(1/((double)*ncon + (double)*ncat))*(lgcon0 + lgcat0);

				  if(*cohesion==2) ph[nclus_iter] = (1/((double)*ncon + (double)*ncat))*(lgcon0 + lgcat0);
				}

				//////////////////////////////////////////////////////////////////////////
				// This is the calibration used when the similarity is standardized by
				//////////////////////////////////////////////////////////////////////////

				if(*calibrate==1){

					maxgtilN = gtilN[0];
					maxgtilY = gtilY[0];
					for(k=1; k < nclus_iter+1; k++){
						if(maxgtilN < gtilN[k]) maxgtilN = gtilN[k];

						if(k < nclus_iter){
							if(maxgtilY < gtilY[k]) maxgtilY = gtilY[k];
						}
					}


					sgY=0.0;
					sgN=0.0;

					for(k=0; k<nclus_iter+1; k++){

						lgtilN[k] = gtilN[k] - maxgtilN;
						sgN = sgN + exp(lgtilN[k]);

						if(k < nclus_iter){
							lgtilY[k] = gtilY[k] - maxgtilY;
							sgY = sgY + exp(lgtilY[k]);
						}
					}


					for(k=0; k<nclus_iter; k++){
						lgtilNk = lgtilN[k] - log(sgN);
						lgtilYk = lgtilY[k] - log(sgY);

						ph[k] = log((double) nh[k]) + lgtilYk - lgtilNk; //This takes into account both cont and cat vars

						if(*cohesion==2) ph[k] = ph[k] - log((double) nh[k]);

					}

									// calibration for a singleton
				  ph[nclus_iter] =  log(Mdp) +
								            lgtilN[nclus_iter] - log(sgN);

				  if(*cohesion==2){// Note with a uniform cohesion, for a new cluster
				                 // the value of log(c({nclus_iter}}) = log(1) = 0;
				    ph[nclus_iter] = ph[nclus_iter] - log(Mdp);
			    }

				}
				//////////////////////////////////////////////////////////////////////////
				// End of calibration used when the similarity is standardized by
				//////////////////////////////////////////////////////////////////////////


                // Normalize cluster predictive probabilities
				maxph = ph[0];
				for(k = 1; k < nclus_iter+1; k++){
					if(ph[k] > maxph) maxph=ph[k];
				}

				denph = 0.0;
				for(k = 0; k < nclus_iter+1; k++){

					ph[k] = exp(ph[k] - maxph);
					denph = denph + ph[k];

				}

				for(k = 0; k < nclus_iter+1; k++){
					probh[k] = ph[k]/denph;
				}

				uu = runif(0.0,1.0);
				cprobh= 0.0;

				for(k = 0; k < nclus_iter+1; k++){
					cprobh = cprobh + probh[k];
					if (uu < cprobh){
						iaux = k+1;
						break;
					}
				}

				if(iaux <= nclus_iter){
					mupred = muh[(iaux-1)];
					sig2pred = sig2h[(iaux-1)];
				}else{
					mupred = rnorm(mu0_iter,sqrt(sig20_iter));
					sig2pred = runif(smin, smax);
					sig2pred = sig2pred*sig2pred;
				}

				ppred_iter[pp] = rnorm(mupred, sqrt(sig2pred));
				predclass_iter[pp] = iaux;

                // Carry out Rao-Blackwellized prediction
                rbpred_iter[pp] = 0.0;
                for(k=0;k<nclus_iter;k++){
                  rbpred_iter[pp] = rbpred_iter[pp] + probh[k]*muh[k];
                  predclass_prob_iter[pp*(*nobs) + k] = probh[k];
                }
                rbpred_iter[pp] = rbpred_iter[pp] + probh[nclus_iter]*rnorm(mu0_iter,sqrt(sig20_iter));
                predclass_prob_iter[pp*(*nobs) + nclus_iter] = probh[nclus_iter];

                // convert the auxiliary variable to the ordinal scale
                for(c=0; c < *nordcat-1; c++){
                  if((ppred_iter[pp] > co[c]) & (ppred_iter[pp] < co[c+1])) ordppred_iter[pp] = c;
                  if((rbpred_iter[pp] > co[c]) & (rbpred_iter[pp] < co[c+1])) rbordpred_iter[pp] = c;
                }
			}
		}



		////////////////////////////////////////////
		//																				//
		// Save MCMC iterates											//
		//																				//
		////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin ==0)){

			mu0[ii] = mu0_iter;
			sig20[ii] = sig20_iter;

			nclus[ii] = nclus_iter;

			for(j = 0; j < *nobs; j ++){

				mu[ii*(*nobs) + j] = muh[Si_iter[j]-1];
				sig2[ii*(*nobs) + j] = sig2h[Si_iter[j]-1];
				Si[ii*(*nobs) + j] = Si_iter[j];
				zi[ii*(*nobs) + j] = zi_iter[j];

				like[ii*(*nobs) + j] = like_iter[j];
				ispred[ii*(*nobs) + j] = ispred_iter[j];
				isordpred[ii*(*nobs) + j] = isordpred_iter[j];

			}

			for(pp = 0; pp < *npred; pp++){

				ppred[ii*(*npred) + pp] = ppred_iter[pp];
				predclass[ii*(*npred) + pp] = predclass_iter[pp];
				ordppred[ii*(*npred) + pp] = ordppred_iter[pp];

                rbpred[ii*(*npred) + pp] = rbpred_iter[pp];
                rbordpred[ii*(*npred) + pp] = rbordpred_iter[pp];

            }

            for(pp = 0; pp < (*nobs)*(*npred); pp++){
                predclass_prob[ii*((*nobs)*(*npred)) + pp] = predclass_prob_iter[pp];
            }

			ii = ii+1;

		}
	}


	PutRNGstate();


	//////////////////////////////////////////////////////////////////////////////////
	// calculate LPML
	//////////////////////////////////////////////////////////////////////////////////
	lpml_iter=0.0;

	for(jj = 0; jj < *nobs; jj++){
		lpml_iter = lpml_iter + log(1/CPOinv[jj]);
	}
	lpml[0] = lpml_iter;

	////////////////////////////////////////////////////////////////////////////////////////////
	// Computing WAIC  (see Gelman article in lit review folder)
	////////////////////////////////////////////////////////////////////////////////////////////
	elppdWAIC = 0.0;
	for(j = 0; j < *nobs; j++){
	  elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
	}

	WAIC[0] = -2*elppdWAIC;

}
