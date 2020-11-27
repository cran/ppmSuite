/****************************************************************************************
 * Copyright (c) 2014 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit PPMx models when covariate values are missing and
 * with an ordinal response
 *
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
* The inputs of the function that are read from R are the following
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
* y = nobs x 1 vector that contains ordinal response values
* co = vector containing fixed cut-off points for associated with latent variable
* Xcon = nobs x ncon contiguous double vector that contains continuous covariate values
* Xcat = nobs x ncat contiguous int vector that contains categorical covariate values
* npred = integer indicating number of out of sample predictions
* nordcat = integer indication number of orgdinal categories
* Xconp = npred x ncon matrix of continuous covariates to make predictions
* Xcatp = npred x ncat matrix of categorical covariates to make predictions
* Mconp = npred x ncon matrix of logicals indicating if continuous covariate is missing
* Mcatp = npred x ncat matrix of logicals indicating if categorical covariate is missing
* calibrate = integer indicating how similarity function should be calibrated
    0 - no calibration
    1 - "calibration" similarity function
    2 - "coarsened" similarity function
* simparms = vector containing similarity functions that are ordered in the following way;
* modelPriors =
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
*
******************************************************************************************/

void ordinal_missing_ppmx(int *draws, int *burn, int *thin, int *nobs, int *ncon, int *ncat,
              int *Cvec, int *PPM, int *cohesion, int *similarity_function, int *consim,
              double *M,
              int *y, double *co, double *Xcon, int *Xcat, int *Mcon, int *Mcat,
			        int *npred, int *nordcat, double *Xconp, int *Xcatp, int *Mconp, int *Mcatp,
		          double *simParms, double *dissimtn, double *dissimtt,
              int *calibrate, double *modelPriors, int *verbose, double *mh,
			        double *mu, double *sig2, double *mu0, double *sig20,
			        int *Si, int *nclus, double *zi, double *like, double *WAIC, double *lpml,
			        double *ispred, int *isordpred, double *ppred, int *predclass, int *ordppred,
			        double *rbpred, int *rbordpred, double *predclass_prob){



	// i - MCMC index
	// j - response index
	// jj - second response index (for double for loops)
	// c - categorical variable index - and ordinal category
	// p - number of covariates index
	// pp - prediction index
	// k - cluster index
	// t - obs per cluster indx

	int i, j, jj, jjj, c, p, pp, k, t, ii;
	int ncov = *ncon + *ncat;

	double nout, sumx, sumx2;
	nout = (*draws-*burn)/(*thin);

	int max_C = Cvec[0];

	for(p = 0; p < (*ncat); p++){
		if(max_C < Cvec[p]) max_C = Cvec[p];
	}
	if(*ncat == 0) max_C = 1.0;

    // Compute mles to be used as if gcontype=3 possible similarity function
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
	double* xcattmp = R_Vector(*nobs);


	for(j = 0; j < *nobs; j++){
		Si_iter[j] = 1;
		nh[j] = 0;
	}

	for(j = 0; j < *nobs; j++){

		for(k = 0; k < *nobs; k++){

			if(Si_iter[j] == k+1) nh[k] = nh[k] + 1;
		}
	}

	for(j = 0; j < *nobs; j++){

		if(nh[j] > 0) nclus_iter = nclus_iter + 1;

	}



	double *sig2h = R_VectorInit(*nobs, 0.1);
	double *muh = R_VectorInit(*nobs, 0.0);

	// Things to produce predictions compute the posterior predictive
	double *ispred_iter = R_VectorInit(*nobs, 0.0);
	int isordpred_iter[*nobs];

	double *ppred_iter = R_Vector((*npred));
	double *rbpred_iter = R_Vector((*npred));
	int ordppred_iter[*npred];
	int rbordpred_iter[*npred];

	int predclass_iter[*nobs];

	double *predclass_prob_iter = R_Vector((*npred)*(*nobs));




	// =============================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// =============================================================================================


	// stuff that I need to update Si (cluster labels);
	int nhtmp, pc;
	double auxm, auxs2, npdN,npdY,npd;
	double mudraw, sdraw=1.0, maxph, denph, cprobh, uu;
	double lgconN,lgconY,lgcatN,lgcatY,lgcondraw,lgcatdraw,lgcatt;

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
	double sumz,mstar,s2star;

	//suff I need to upate mu0
	double summu;

	//stuff I need to update sig20
	double llo, lln, llr,  os0,ns0;

	//stuff I need to update sig2
	double osig,nsig;


	// Stuff for out of sample predictions
	double lgcon0, lgcat0, mupred, sig2pred;

	// Stuff to compute lpml, likelihood, and WAIC
	double lpml_iter, elppdWAIC, sdens;
	double *CPOinv = R_VectorInit(*nobs, 0.0);
	double *like_iter = R_VectorInit(*nobs, 0.0);
	double *mnlike = R_VectorInit(*nobs, 0.0);
	double *mnllike = R_VectorInit(*nobs, 0.0);


	// =============================================================================================
	//
	// Prior parameter values
	//
	// =============================================================================================

 	// priors for mu0
	double m = modelPriors[0]; double s2 = modelPriors[1];

    // prior values for sig2h and sig20;
	double smin=0, smax=modelPriors[2];
	double s0min=0, s0max=modelPriors[3];

    Rprintf("m = %f, s2 = %f\n", m, s2);
    Rprintf("smin = %f, smax = %f, s0min = %f, s0max = %f\n", smin, smax, s0min, s0max);

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


    if(*verbose){
	  Rprintf("nobs = %d\n", *nobs);
	  Rprintf("ncon = %d\n", *ncon);
	  Rprintf("ncat = %d\n", *ncat);
	  Rprintf("ncov = %d\n", ncov);
	  Rprintf("npred = %d\n", *npred);
	  RprintVecAsMat("co", co, 1, 3);
	  if(*ncat > 0 ) RprintIVecAsMat("Cvec",Cvec, 1, *ncat);
	  Rprintf("Prior values being used are: \n m = %.1f, s2 = %.1f \n smin = %.1f, smax = %.1f,\n s0min = %.1f, s0max = %.1f\n", m, s2, smin, smax, s0min, s0max);
	  Rprintf("M = %f\n", Mdp);
	  Rprintf("Simlarity values: m0 = %.2f, s20 = %.2f v = %.2f, k0 = %.2f, nu0 = %.2f, a0 = %.2f, alpha = %.2f\n", m0, s20, v, k0, nu0, simParms[5], alpha);
	  if(*nobs > 50) RprintIVecAsMat("First fifty response values",y, 1, 50);
    }

	// M-H tuning parameters
	double csigSIG0=mh[0], csigSIG=mh[1];

	GetRNGstate();

	ii = 0;


	// =============================================================================================
	//
	// start of the mcmc algorithm;
	//
	// =============================================================================================

	double calc_time = 0.0;
	clock_t  begin = clock();

	for(i = 0; i < *draws; i++){

	    if(*verbose & ((i+1) % 1000 == 0)){
//	      time_t now;
//	      time(&now);

//	      Rprintf("mcmc iter = %d ===================================================== \n", i+1);
//        Rprintf("%s", ctime(&now));

	    }

	    if(1){
          clock_t ith_iterate = clock();
	      calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

          Rprintf("Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*draws))*100.0, calc_time);
//          fflush(stdout);
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
		// To accomodate missing here, if a covariate is missing for individual i, then
		// I am simply not evaluating the similarity function for that individual.
		// This is done by dragging around a logical vector that indicates missing or not.
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


			// The atoms have been relabeled if necessary and now we need to
			// update Si.

			// Begin the cluster probabilities

			for(k=0; k<nclus_iter; k++){

				lgconY = 0.0;
				lgconN = 0.0;

				// Evaluate the similarity function for each continuous covariate
				for(p=0; p<(*ncon); p++){
					nhtmp = 0;
					for(jj = 0; jj < *nobs; jj++){
						if(jj != j){
							if((Si_iter[jj] == k+1) & (Mcon[jj*(*ncon)+p] == 0)){
								xcontmp[nhtmp] = Xcon[jj*(*ncon)+p];
								nhtmp = nhtmp+1;
							}
						}
					}

					sumx = 0.0;
					sumx2 = 0.0;
					for(t = 0; t < nhtmp; t++){

						sumx = sumx + xcontmp[t];
						sumx2 = sumx2 + xcontmp[t]*xcontmp[t];

					}

                    // Evaluate similarity function for continuous covariate on cluster *excluding* the ith individual
					if(nhtmp > 0){
						if(*similarity_function==1){ // Auxilliary similarity function
							if(*consim==1){ // NN
							  lgconN = lgconN + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
							}
							if(*consim==2){ // NNIG
								lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
							}

						}
						if(*similarity_function==2){ // Double Dipper similarity function
							if(*consim==1){// NN
								lgconN = lgconN + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
							}
							if(*consim==2){// NNIG
								lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
							}
						}
						if(*similarity_function==3){ // Cluster Variance similarity function
							lgconN = lgconN + gsimconEV(sumx, sumx2, nhtmp, alpha,1);
						}

					}

					// now add jth individual back;
					if(Mcon[j*(*ncon)+p] == 0){ // 0 indicates it is not missing
						xcontmp[nhtmp] = Xcon[j*(*ncon)+p];
						sumx = sumx + Xcon[j*(*ncon)+p];
						sumx2 = sumx2 + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
						nhtmp = nhtmp+1;
					}

					// Evaluate similiarty function for continuous covariate on cluster *including* ith individual
					if(nhtmp > 0){
						if(*similarity_function==1){ // Auxilliary Similarity function
							if(*consim==1){ // NN
								lgconY = lgconY + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
							}
							if(*consim==2){ // NNIG
								lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
							}
						}
						if(*similarity_function==2){ //Double Dipper Similarity function
							if(*consim==1){// NN
								lgconY = lgconY + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
							}
							if(*consim==2){// NNIG
								lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
							}
						}
						if(*similarity_function==3){ // Cluster Variance similarity function
							lgconY = lgconY + gsimconEV(sumx, sumx2, nhtmp, alpha,1);
						}
					}

				}


                // Now add the categorical covaraite similarity function
				lgcatY=0.0;
				lgcatN=0.0;
				for(p=0; p<(*ncat); p++){
					for(c=0;c<Cvec[p];c++){nhc[c]=0;}

					nhtmp = 0;
					for(jj = 0; jj < *nobs; jj++){
						if(jj != j){

							if((Si_iter[jj]==k+1) & (Mcat[jj*(*ncat)+p] == 0)){ // 0 indicates not missing
								nhc[Xcat[jj*(*ncat)+p]] = nhc[Xcat[jj*(*ncat)+p]] + 1; // this needs to be a vector
								xcattmp[nhtmp] = Xcat[jj*(*ncat)+p];
								nhtmp = nhtmp+1;

							}
						}
					}

					// Evaluate similiarty function for categorical covariate on cluster *excluding* ith individual
					if(nhtmp >0){
					  if(*similarity_function==1){// Auxilliary Similarity function
							lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
					  }
					  if(*similarity_function==2){//Double Dipper Similarity function
							lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
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
					}


					if(Mcat[j*(*ncat)+p] == 0){
						xcattmp[nhtmp] = Xcat[j*(*ncat)+p];
						nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1;

						nhtmp = nhtmp + 1;
					}

					// Evaluate similiarty function for categorical covariate on cluster *including* ith individual

					if(nhtmp > 0){
						if(*similarity_function==1){// Auxilliary Similarity function
							lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
						}
						if(*similarity_function==2){//Double Dipper Similarity function
							lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
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
				}

                // I need to save similarity values to employ if calibration is used
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
					 lgconN = -(alpha)*lgconN/(npdN);
					 lgconY = -(alpha)*lgconY/(npdY);

					// Just use cluster-total gower dissimilarity
					lgconN = -(alpha)*lgconN;
					lgconY = -(alpha)*lgconY;

				}
			   	//////////////////////////////////////////////////////////
				// End of Gower similarity
				//////////////////////////////////////////////////////////








                // Calculate unormalized cluster probabilities
				ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
				        log((double) nh[k]) +
				        lgcatY - lgcatN +
						lgconY - lgconN;

                // If PPM is used rather than PPMx remove the covariate part
				if(*PPM){
					ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
						    	log((double) nh[k]);  // DP part of cohesion function
				}

				// If calibrate through coarsening the use the following
				if(*calibrate == 2){

					ph[k] = dnorm(zi_iter[j], muh[k], sqrt(sig2h[k]), 1) +
				            log((double) nh[k]) +
				            (1/((double)*ncon + (double)*ncat))*(lgcatY + lgconY - lgcatN - lgconN);
				}


				// If a uniform cohesion is used rather than the DP type cohesion, remove the cluster size
				if(*cohesion==2){
					ph[k] = ph[k] - log((double) nh[k]);
				}


			}

			// Now need to consider ith observation being assigned to own cluster
			mudraw = rnorm(mu0_iter, sqrt(sig20_iter));
			sdraw = runif(smin, smax);

            // Evaluate similarity function for *continuous* covariate when ith observation is in singleton cluster
			lgcondraw = 0.0;
			for(p=0;p<(*ncon);p++){
				if(Mcon[j*(*ncon)+p] == 0){
					xcontmp[0] = Xcon[j*(*ncon)+p];

					if(*similarity_function==1){ // Auxilliary Similarity
						if(*consim==1){ // NN
							lgcondraw = lgcondraw + gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0, 1);
						}
						if(*consim==2){ // NNIG
							lgcondraw = lgcondraw + gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 0,0, 1);
						}
					}
					if(*similarity_function==2){ // Double Dipper Similarity
						if(*consim==1){ // NN
							lgcondraw = lgcondraw + gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p], 1, 1, 0, 1);
						}
						if(*consim==2){ // NNIG
							lgcondraw = lgcondraw + gsimconNNIG(m0, k0, nu0, s20, xcontmp[0], xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p], 1, 1, 0, 1);
						}
					}
					if(*similarity_function==3){ // Variance Similarity function
						lgcondraw = lgcondraw + gsimconEV(xcontmp[0], xcontmp[0]*xcontmp[0], 1,alpha,1);
					}
				}
			}


            // Evaluate similarity function for *categorical* covariate when ith observation is in singleton cluster
			lgcatdraw = 0.0;
			for(p=0;p<(*ncat);p++){
				for(c=0;c<Cvec[p];c++){nhc[c] = 0;}
	            if(Mcat[j*(*ncat)+p] == 0){

				  nhc[Xcat[j*(*ncat)+p]] = 1;

			      if(*similarity_function==1){ // Auxilliary Similarity
					  lgcatdraw = lgcatdraw + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
				  }
				  if(*similarity_function==2){// Double Dipper Similarity
					  lgcatdraw = lgcatdraw + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
				  }
				  if(*similarity_function==3){// Entropy similarity
				    lgcatdraw = lgcatdraw + -(alpha)*0;
				  }
			   }
			}

            // Need to add this to vector to calibrate later
			gtilY[nclus_iter] = lgcondraw + lgcatdraw;
			gtilN[nclus_iter] = lgcondraw + lgcatdraw;


			// unormalized singleton cluster probability
			ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
			                 log(Mdp) +
			                 lgcondraw +
			                 lgcatdraw;

            // If similarity is calibrated via coarsening
			if(*calibrate==2){
				ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
				                 log(Mdp) +
				                 (1/((double)*ncon + (double)*ncat))*(lgcondraw + lgcatdraw);
			}

			// If PPM rather than PPMx remove the covariate part
			if(*PPM){
				ph[nclus_iter] = dnorm(zi_iter[j],mudraw,sdraw,1) +
				                 log(Mdp); //DP part
			}


            // If a uniform cohesion is used instead
			if(*cohesion==2){
				ph[nclus_iter] = ph[nclus_iter] - log(Mdp);
			}

			/////////////////////////////////////////////////////////////////////////////
			// This is the calibration used when the similarity is normalized
			/////////////////////////////////////////////////////////////////////////////
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

			// Now to normalize the cluster probabilities
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

			// Select cluster label based on normalized probabilities
			cprobh= 0.0;
			pc=nclus_iter+1;
			for(k = 0; k < nclus_iter+1; k++){
				cprobh = cprobh + probh[k];
				if (uu < cprobh){
					pc = k+1;
					break;
				}
			}

			// If assigned to an existing cluster add one to cluster size
			if(pc <= nclus_iter){

				Si_iter[j] = pc;
				nh[Si_iter[j]-1] = nh[Si_iter[j]-1] + 1;

			}else{

			  // If assigned to singleton cluster, increase number of clusters by one
			  // and create the new cluster's mean and variance.
				nclus_iter = nclus_iter + 1;
				Si_iter[j] = nclus_iter;
				nh[Si_iter[j]-1] = 1;

				muh[Si_iter[j]-1] = mudraw;
				sig2h[Si_iter[j]-1] = sdraw*sdraw;

			}

			// Compute the CPO and lpml using the mixture and Gaussian cdf
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
		// update muh's cluster specific intercepts
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
		if((i > (*burn-1)) & ((i+1) % *thin ==0)){

			for(j = 0; j < *nobs; j++){

				ispred_iter[j] = rnorm(muh[Si_iter[j]-1], sqrt(sig2h[Si_iter[j]-1]));

                for(c=0; c < *nordcat-1; c++){
                  if((ispred_iter[j] > co[c]) & (ispred_iter[j] < co[c+1])) isordpred_iter[j] = c;
                }
			}
		}



		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// out of sample prediction using posterior predictive
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		if((i > (*burn-1)) & ((i+1) % *thin ==0)){
			for(pp = 0; pp < *npred; pp++){

				for(k = 0; k < nclus_iter; k++){

                    // evaluate similarity of cluster if "new" subject's covariate value included
					lgconN=0.0, lgconY=0.0;
					for(p=0; p<(*ncon); p++){
						nhtmp = 0;
						for(j = 0; j < *nobs; j++){
							if((Si_iter[j] == k+1) & (Mcon[j*(*ncon)+p] == 0)){
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

                        // First evalute similarity of cluster *excluding* the new individual
						if(nhtmp > 0){
							if(*similarity_function==1){// Auxilliary Similarity Function
								if(*consim==1){// NN
									lgconN = lgconN + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconN = lgconN + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
								}
							}
							if(*similarity_function==2){// Double Dipper Similarity
								if(*consim==1){// NN
									lgconN = lgconN + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconN = lgconN +  gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
								}
							}
							if(*similarity_function==3){// Cluster variance similarity
								lgconN = lgconN + gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
							}
						}

                        // Next evalute the similarity of cluster *including* the new individual
						if(Mconp[pp*(*ncon)+p] == 0){
							xcontmp[nhtmp] = Xconp[pp*(*ncon)+p];
							sumx = sumx + Xconp[pp*(*ncon)+p];
							sumx2 = sumx2 + Xconp[pp*(*ncon)+p]*Xconp[pp*(*ncon)+p];
							nhtmp = nhtmp + 1;
						}

						if(nhtmp > 0){
							if(*similarity_function==1){ // Auxilliary Similarity function
								if(*consim==1){// NN
									lgconY = lgconY + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 0, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
								}
							}
							if(*similarity_function==2){ // Double Dipper Similarity function
								if(*consim==1){// NN
									lgconY = lgconY + gsimconNN(m0, v, s20, sumx, sumx2, mnmle[p], nhtmp, 1, 0, 1);
								}
								if(*consim==2){// NNIG
									lgconY = lgconY + gsimconNNIG(m0, k0, nu0, s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
								}
							}
							if(*similarity_function==3){ // Cluster variance similarity function
								lgconY = lgconY + gsimconEV(sumx, sumx2, nhtmp,alpha, 1);
							}
						}

					}

					// Now need to add the categorical covariate similarity value
					lgcatY=0.0, lgcatN=0.0;
					for(p=0; p<(*ncat); p++){

						for(c=0;c<Cvec[p];c++){nhc[c]=0;}

						nhtmp=0;
						for(j = 0; j < *nobs; j++){

							if((Si_iter[j]==k+1) & (Mcat[j*(*ncat)+p] == 0)){
								nhc[Xcat[j*(*ncat)+p]] = nhc[Xcat[j*(*ncat)+p]] + 1; // this needs to be a vectore
								xcattmp[nhtmp] = Xcat[j*(*ncat)+p];
								nhtmp = nhtmp+1;
							}
						}

                        // First evaluate similarity function for categorical variable *excluding* the "new" observation
						if(nhtmp > 0){
							if(*similarity_function==1){// Auxilliary Similarity
								lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
							}
							if(*similarity_function==2){// Double Dipper Similarity
								lgcatN = lgcatN + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
							}
							if(*similarity_function==3){// Entropy similarity
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
						}



                        // Now evaluate similarity function for categorical variable *including* the "new" observation
						if(Mcatp[pp*(*ncat)+p] == 0){
								nhc[Xcatp[pp*(*ncat)+p]] = nhc[Xcatp[pp*(*ncat)+p]] + 1;
								nhtmp=nhtmp + 1;
						}

						if(nhtmp > 0){
							if(*similarity_function==1){// Auxilliary simlarity
								lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
							}
							if(*similarity_function==2){// Double Dipper similarity
								lgcatY = lgcatY + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
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


						  // Use the cluster-total Gower dissimilarity
						  lgconN = -(alpha)*lgconN;
						  lgconY = -(alpha)*lgconY;
					}
			        //////////////////////////////////////////////////////////
				    // End of Gower similarity values for gower dissimilarity
				    //////////////////////////////////////////////////////////



					// calculate the unormalized predictive probability
					ph[k] = log((double) nh[k]) +
				         	  lgcatY - lgcatN +
							      lgconY - lgconN;

                    // If calibration through coarsening
					if(*calibrate == 2){
						ph[k] = log((double) nh[k]) +
						           (1/((double)*ncon + (double)*ncat))*
						            (lgcatY + lgconY - lgcatN - lgconN);

					}

                    // If PPM rather than PPMx then only cluster size influentials predictive distribution
					if(*PPM) ph[k] = log((double) nh[k]);

                    // If uniform cohesion, then remove cluster size
					if(*cohesion==2) ph[k] = lgcatY - lgcatN + lgconY - lgconN;

			  }

              // Now need to consider allocating "new" individual to own cluster
              lgcon0=0.0;
              for(p=0;p<*ncon;p++){ // 0 means that data are not missing

                if(Mconp[pp*(*ncon)+p] == 0){
        	      xcontmp[0] = Xconp[pp*(*ncon)+p];

        		  if(*similarity_function==1){ // Auxilliary
        		    if(*consim==1){ // NN
        				lgcon0 = lgcon0 + gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,0,0,1);
        			}
        			if(*consim==2){// NNIG
        				lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],s2mle[p],1,0,0,1);
        			}
        		  }
        		  if(*similarity_function==2){ // Double Dipper
        			if(*consim==1){ // NN
        				lgcon0 = lgcon0 + gsimconNN(m0,v,s20,xcontmp[0],xcontmp[0]*xcontmp[0], mnmle[p],1,1,0,1);
        			}
        			if(*consim==2){// NNIG
        				lgcon0 = lgcon0 + gsimconNNIG(m0, k0, nu0, s20,xcontmp[0],xcontmp[0]*xcontmp[0],mnmle[p],s2mle[p],1, 1, 0,1);
        			}
        		  }
        		  if(*similarity_function==3){ // Variance
        			lgcon0 = lgcon0 + gsimconEV(xcontmp[0], xcontmp[0]*xcontmp[0], 1,alpha,1);
        		  }

        	    }
              }

             // Now evaluate categorical similarity if "new" subject is allocated to singleton cluster

             lgcat0 = 0.0;
             for(p=0;p<(*ncat);p++){
               for(c=0;c<Cvec[p];c++) nhc[c] = 0;

                 if(Mcatp[pp*(*ncat)+p] == 0){

                   nhc[Xcatp[pp*(*ncat)+p]] = 1;

        	       if(*similarity_function==1){
        	         lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 0, 1);
        		   }
        		   if(*similarity_function==2){
        			 lgcat0 = lgcat0 + gsimcatDM(nhc, dirweights, Cvec[p], 1, 1);
        		   }
        		   if(*similarity_function==3){
        			 lgcat0 = lgcat0 +  -(alpha)*0;
        		   }

        	    }

              }

			  gtilY[nclus_iter] = lgcat0 + lgcon0;
			  gtilN[nclus_iter] = lgcat0 + lgcon0;


              // Unormalized predictive probabilities that now include singleton
              ph[nclus_iter] = log((double) Mdp) + lgcon0 + lgcat0;

			  // Gower similarity
			  if(*similarity_function==4) ph[nclus_iter] = log((double) Mdp) + log(1);


              // Calibrate through coarsening
              if(*calibrate==2){
        	     ph[nclus_iter] = log(Mdp) +
        					(1/((double)*ncon + (double)*ncat))*(lgcon0 + lgcat0);

              }

              // Here if PPM is selected, remove covariate information.  Singleton weight based on M parameter
              if(*PPM) ph[nclus_iter] = log((double) Mdp);

              // If uniform cohesion selected, remove M parameter from DP cohesion
              if(*cohesion==2) ph[nclus_iter] = ph[nclus_iter] - log((double) Mdp);


			//////////////////////////////////////////////////////////////////////////
			// This is the calibration used when the similarity is standardized by
			//////////////////////////////////////////////////////////////////////////

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

					if(k < nclus_iter){
						lgtilY[k] = gtilY[k] - maxgtilY;
						sgY = sgY + exp(lgtilY[k]);
					}
				}


				for(k=0; k<nclus_iter; k++){
					lgtilNk = lgtilN[k] - log(sgN);
					lgtilYk = lgtilY[k] - log(sgY);

					ph[k] = log((double) nh[k]) + lgtilYk - lgtilNk; //This takes into account both cont and cat vars

				  	if(*cohesion==2){
			        	ph[k] = ph[k] - log((double) nh[k]);
			    	  	}
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

        // Need to normalize the predictive probabilities
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

        // Probabilistically allocate new individual to a cluster
        cprobh= 0.0;
        pc = nclus_iter+1;
        for(k = 0; k < nclus_iter+1; k++){
          cprobh = cprobh + probh[k];
          if (uu < cprobh){
        	  pc = k+1;
        	  break;
        	}

        }

        // If allocated to an existing cluster, then predicted mean is cluster mean
        if(pc <= nclus_iter){
        	mupred = muh[(pc-1)];
        	sig2pred = sig2h[(pc-1)];
        }else{
        	mupred = rnorm(mu0_iter,sqrt(sig20_iter));
        	sig2pred = runif(smin, smax);
        	sig2pred = sig2pred*sig2pred;
        }

        // Take a draw from the posterior predictive of the auxiliary variable
        ppred_iter[pp] = rnorm(mupred, sqrt(sig2pred));
		predclass_iter[pp] = pc;

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
          if((rbpred_iter[pp] > co[c]) & (rbpred_iter[pp] < co[c+1]) ) rbordpred_iter[pp] = c;
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
	// (see page 81 Christensen Hansen and Johnson)
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
