/*************************************************************
 * Copyright (c) 2014 Garritt Leland Page
 *
 * This file contains C code for an MCMC algorithm constructed
 * to fit the spatial PPM model that contains spatial stucture
 * only in partition model.
 *
 * letting rho denote a partition of spatial locations we
 * have the following model
 *
 *	Pr(rho) \propto \prod c(S_j)g(x*_j)

 *	Y(s_i) | mu^*, sig2^*, rho ~ N(m^*_{c_i}, sig2^*_{c_i})
 *
 *	Priors:
 *
 *  mu_j ~ N(m, S)
 *  sig2_j ~ Unif(0, ms)
 *
 * The definition of c(S_j) is key to producing spatially
 * "pleasing" clusters
 *************************************************************/

#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>



/***************************************************************************************************
* The following are the inputs of the function that are read from R
*
* draws = total number of MCMC draws
* burn = number of MCMC draws discarded as burn-in
* thin = indicates how much MCMC chain should be thinned
* nobs = number of spatial locations
* npred = number of locations at which to make predictions
* Cohesion = number indicating what spatial cohesion is employed.
* cParm = cohesion parameters

* modelPriors = priors associated with data model
* PPM = logical indicating wether to use PPM or sPPM default is sPPM
* M = value of DP precisions parameter
* y = vector containing nobs response values at each spatial location
* alpha = weight parameter associated with Cohesion 1
* s1 = vector containing nobs vertical location value
* s2 = vector containing nobs horizontal location value
* s1p = vector containing npred vertical location value
* s2p = vector containing npred horizontal location value
*
* output - vectors
* 	mu -
* 	sig2 -
* 	Si -
* 	nclus -
* 	mu0 -
* 	sig20 -
* 	like -
* 	lpml -
* 	WAIC -
* 	ispred -
* 	ppred -
***************************************************************************************************/



void mcmc_sppm(int *draws, int *burn, int *thin, int *nobs, int *npred, int *cohesion,
			  int *PPM, double *M, double *y,
			  double *s1, double *s2,double *s1p, double *s2p,
			  double *modelPriors, double *cParms, double *mh,
			  double *mu, double *sig2, double *mu0, double *sig20,
			  int *Si, int *nclus, double *like, double *WAIC, double *lpml,
			  double *ispred, double *ppred){



	// i - MCMC iterate
	// ii - MCMC save iterate
	// j - spatial location iterate
	// jj - second spatial location iterate
	// k - cluster iterate
	// pp - prediction iterate

	int i, ii, j, jj, k, pp;


	int nout = (*draws - *burn)/(*thin);
	Rprintf("nout = %d\n", nout);


	double uu;

	Rprintf("nobs = %d\n", *nobs);

	Rprintf("npred = %d\n", *npred);
	Rprintf("cohesion = %d\n", *cohesion);

	double mny;
	mny = 0.0;
	for(j=0; j<*nobs; j++){
		mny = mny + y[j]/(*nobs);
	}

	Rprintf("mny = %f\n", mny);


	// =============================================================================================
	//
	// Memory vectors to hold a single MCMC iterate
	//
	// =============================================================================================

	double sig20_iter=modelPriors[3]*0.5, mu0_iter=0.0;
	int nclus_iter = 0;

	int iaux=1;


	int Si_iter[*nobs];
	int nh[*nobs];

	double *ppred_iter = R_VectorInit(*npred, 0.0);

	// Cluster specific parameters
	double *muh = R_VectorInit((*nobs), 0.0);
	double *sig2h = R_VectorInit((*nobs), modelPriors[2]*0.5);

	// Start each observation in the same cluster
	for(j=0; j < *nobs; j++){
		Si_iter[j] = 1;
		nh[j] = 0;
	}

	for(j = 0; j < *nobs; j++){

		for(jj = 0; jj < *nobs; jj++){

			if(Si_iter[j] == jj+1) nh[jj] = nh[jj] + 1;
		}
	}

	for(j = 0; j < *nobs; j++){

		if(nh[j] > 0) nclus_iter = nclus_iter + 1;

	}



	// =============================================================================================
	//
	// scratch vectors of memory needed to update parameters
	//
	// =============================================================================================


	// stuff I need for sig2 and s2b0; (tried gibbs and MH
	double osig, nsig, llo, lln, llr, os0, ns0;

	// stuff I need for mu;
	double sumy, s2star, mstar, summu;


	// stuff that I need to update Si (cluster labels);
	int indx;
	double auxm, auxs2, lCo=0, lCn=0, lCn_1=0;
	double mudraw, sigdraw, maxph, denph, cprobh;

	double *ph = R_VectorInit(*nobs, 0.0);
	double *phtmp = R_VectorInit(*nobs, 0.0);
	double *probh = R_VectorInit(*nobs, 0.0);

	double *s1o = R_Vector(*nobs);
	double *s2o = R_Vector(*nobs);
	double *s1n = R_Vector(*nobs);
	double *s2n = R_Vector(*nobs);


	// Stuff to do predictions
	double lC1=log(1), mupred=0, sig2pred=1.0;

	// Stuff to compute lpml, likelihood, and WAIC
	double lpml_iter, elppdWAIC;
	double *CPO = R_VectorInit(*nobs, 0.0);
	double *like_iter = R_VectorInit(*nobs, 0.0);
	double *mnlike = R_VectorInit(*nobs, 0.0);
	double *mnllike = R_VectorInit(*nobs, 0.0);




	// prior for mu0 and s20
	double m0 = modelPriors[0], s20 = modelPriors[1];

	double ms0=modelPriors[3]; // max value of sig20 for MH
	Rprintf("ms0 = %f\n", ms0);

	// prior values for sig2
	double ms=modelPriors[2];  //max value of sig2 when MH is used
	Rprintf("ms = %f\n", ms);




	// Cohesion parameters
	// DP Cohesion function parameters (use DP analog)
	double Mdp=*M;
	Rprintf("M = %f\n", Mdp);

	// Epsilon value for Cohesion 1
	double eC1 = cParms[0];
	Rprintf("eC1 = %f\n", eC1);
	// distance cut-off parameter of Cohesion 2
	double aC2 = cParms[1];
	Rprintf("aC2 = %f\n", aC2);

	// Cohesion auxiliary model paramaters for Cohesions 3 and 4
	double k0=cParms[3], v0=cParms[4];
	double *mu0c = R_VectorInit(2,cParms[2]);
	double *L0 = R_VectorInit(2*2,0);
	for(j = 0; j < 2; j++){
		for(jj = 0; jj < 2; jj++){
			if(jj==j) L0[jj*2+j] = cParms[5];
		}
	}

	Rprintf("v0 = %f\n", v0);
	Rprintf("k0 = %f\n", k0);
	RprintVecAsMat("mu0c", mu0c, 1, 2);
	RprintVecAsMat("L0", L0, 2, 2);


	// stuff for M-H
	double csigSIG=mh[0], csigSIG0=mh[1];


	ii = 0;

	GetRNGstate();


	// =============================================================================================
	//
	// start of the mcmc algorithm;
	//
	// =============================================================================================

	for(i = 0; i < *draws; i++){


		if((i+1) % 10000 == 0){
			time_t now;
			time(&now);

			Rprintf("mcmc iter = %d ===================================================== \n", i+1);
			Rprintf("%s", ctime(&now));
		}



		//////////////////////////////////////////////////////////////////////////////////
		//
		// update the cluster labels using the polya urn scheme of
		// algorithm 8 found in  Radford Neal's
		// "Markov Chain Sampling Methods for Dirichlet Process Mixture Models"
		// paper.
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


					// All members of last cluster will be assigned subject i's cluster lable
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
			for(k = 0; k < nclus_iter; k++){


				indx = 0;
				for(jj = 0; jj < *nobs; jj++){

					if((Si_iter[jj] == k+1) & (j != jj)){

						s1o[indx] = s1[jj];
						s2o[indx] = s2[jj];

						s1n[indx] = s1[jj];
						s2n[indx] = s2[jj];

						indx = indx+1;
					}
					if(j == jj){
						s1n[nh[k]] = s1[jj];
						s2n[nh[k]] = s2[jj];

					}


				}


				if(*cohesion==1){
					lCo = Cohesion1(s1o, s2o, eC1, nh[k],1);
					lCn = Cohesion1(s1n, s2n, eC1, nh[k]+1,1);
				}

				if(*cohesion==2){
					lCo = Cohesion2(s1o, s2o, aC2, nh[k],1);
					lCn = Cohesion2(s1n, s2n, aC2, nh[k]+1,1);
				}

				if(*cohesion==3){
					lCo = Cohesion3_4(s1o, s2o, mu0c, k0, v0, L0, nh[k], 3, 1);
					lCn = Cohesion3_4(s1n, s2n, mu0c, k0, v0, L0, nh[k]+1,3, 1);
				}

				if(*cohesion==4){

					lCo = Cohesion3_4(s1o, s2o, mu0c, k0, v0, L0, nh[k], 4, 1);
					lCn = Cohesion3_4(s1n, s2n, mu0c, k0, v0, L0, nh[k]+1,4, 1);

				}



				ph[k] = dnorm(y[j], muh[k], sqrt(sig2h[k]), 1) +
					        lCn - lCo + // Spatial part of cohesion function
					    	log((double) nh[k]);  // DP part of cohesion function

				if(*PPM){
					ph[k] = dnorm(y[j], muh[k], sqrt(sig2h[k]), 1) +
					    		log((double) nh[k]);  // DP part of cohesion function
				}



			}



			// Now to evaluate probability of new cluster.

			mudraw = rnorm(mu0_iter, sqrt(sig20_iter));
			sigdraw = runif(0,ms);
//			Rprintf("b0draw = %f\n", b0draw);

			s1o[0] = s1[j];
			s2o[0] = s2[j];
			if(*cohesion==1) lCn_1 = log(1); // This is because set the cohesion to M with singletons
			if(*cohesion==2) lCn_1 = log(1);
			if(*cohesion==3) lCn_1 = Cohesion3_4(s1o, s2o, mu0c, k0, v0, L0, 1, 3, 1);
			if(*cohesion==4) lCn_1 = Cohesion3_4(s1o, s2o, mu0c, k0, v0, L0, 1, 4, 1);

			ph[nclus_iter] = dnorm(y[j], mudraw, sigdraw, 1) +
			                 lCn_1 +  //Spatial part
			                 log(Mdp); // DP part

			if(*PPM){

				ph[nclus_iter] = dnorm(y[j], mudraw, sigdraw, 1) +
					                 log(Mdp); // DP part

			}



			for(k = 0; k < nclus_iter+1; k++) phtmp[k] = ph[k];


			R_rsort(phtmp,  nclus_iter+1) ;

			maxph = phtmp[nclus_iter];

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
			iaux=nclus_iter+1;
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
				sig2h[Si_iter[j]-1] = sigdraw*sigdraw;


			}



		}



		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update mu (this is done with a Gibbs Sampler)
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		for(k = 0; k < nclus_iter; k++){

			sumy = 0.0;
			for(j = 0; j < *nobs; j++){
				if(Si_iter[j] == k+1){
					sumy = sumy + y[j];
				}

			}


			s2star = 1/((double) nh[k]/sig2h[k] + 1/sig20_iter);
			mstar = s2star*( (1/sig2h[k])*sumy + (1/sig20_iter)*mu0_iter);

			muh[k] = rnorm(mstar, sqrt(s2star));

		}


		//////////////////////////////////////////////////////////////////////////////////
		//																				//
		// Update sig2  cluster specific Use metropolis and Uniform prior.				//
		//																				//
		//////////////////////////////////////////////////////////////////////////////////
		for(k = 0; k < nclus_iter; k++){

			osig = sqrt(sig2h[k]);
			nsig = rnorm(osig,csigSIG);


			if(nsig > 0){
				lln = 0.0;llo = 0.0;
				for(j = 0; j < *nobs; j++){
					if(Si_iter[j] == k+1){

						llo = llo + dnorm(y[j], muh[k], osig,1);
						lln = lln + dnorm(y[j], muh[k], nsig,1);
					}
				}

				llo = llo + dunif(osig, 0, ms, 1);
				lln = lln + dunif(nsig, 0, ms, 1);

				llr = lln - llo;
				uu = runif(0,1);

				if(log(uu) < llr){
					sig2h[k] = nsig*nsig;
				}
			}
		}


		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update mu0  prior mean of mu
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		summu= 0.0;
		for(k = 0; k < nclus_iter; k++){

			summu = summu + muh[k];

		}


		s2star = 1/(((double) nclus_iter/sig20_iter) + (1/s20));
		mstar = s2star*((1/sig20_iter)*summu + m0*(1/s20));

		mu0_iter = rnorm(mstar, sqrt(s2star));




		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Update sig20  prior variance of mu
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


			llo = llo + dunif(os0, 0, ms0, 1);
			lln = lln + dunif(ns0, 0, ms0, 1);


			llr = lln - llo;
			uu = runif(0,1);

			if(log(uu) < llr){
				sig20_iter = ns0*ns0;
			}
		}



		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// evaluating likelihood that will be used to calculate LPML?
		// (see page 81 Christensen Hansen and Johnson)
		//
		////////////////////////////////////////////////////////////////////////////////////////////

		if((i > (*burn-1)) & (i % (*thin) == 0)){
			for(j = 0; j < *nobs; j++){

				like_iter[j] = dnorm(y[j], muh[Si_iter[j]-1], sqrt(sig2h[Si_iter[j]-1]), 0);

				// These are needed for WAIC
				mnlike[j] = mnlike[j] + (like_iter[j])/(double) nout;
				mnllike[j] = mnllike[j] + log(like_iter[j])/(double) nout;

				CPO[j] = CPO[j] + (1/(double) nout)*
				                  (1/like_iter[j]);

			}

		}


		lpml_iter=0.0;
		if(i == (*draws-1)){

			for(jj = 0; jj < *nobs; jj++){

				lpml_iter = lpml_iter + log(1/CPO[jj]);

			}

			lpml[0] = lpml_iter;

		}

		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// Computing WAIC  (see Gelman article in lit review folder)
		// An issue with WAIC is that considering spatially structure data
		//
		////////////////////////////////////////////////////////////////////////////////////////////
		elppdWAIC = 0.0;
		if(i == (*draws - 1)){

			for(j = 0; j < *nobs; j++){

				elppdWAIC = elppdWAIC + (2*mnllike[j] - log(mnlike[j]));
			}

			WAIC[0] = -2*elppdWAIC;
		}

		////////////////////////////////////////////////////////////////////////////////////////////
		//
		// out of sample prediction using posterior predictive?
		//
		////////////////////////////////////////////////////////////////////////////////////////////

		if((i > (*burn-1)) & (i % (*thin) == 0)){
			for(pp = 0; pp < *npred; pp++){

				for(k = 0; k < nclus_iter; k++){

					indx = 0;
					for(jj = 0; jj < *nobs; jj++){

						if(Si_iter[jj] == k+1){

							s1o[indx] = s1[jj];
							s2o[indx] = s2[jj];

							s1n[indx] = s1[jj];
							s2n[indx] = s2[jj];

							indx = indx+1;

						}

					}
					s1n[nh[k]] = s1p[pp];
					s2n[nh[k]] = s2p[pp];

					if(*cohesion==1){
						lCo = Cohesion1(s1o, s2o, eC1, nh[k], 1);
						lCn = Cohesion1(s1n, s2n, eC1, nh[k]+1, 1);
						lC1 = log(1); // this is because we set this cohesion to M for singleton clusters;
					}

					if(*cohesion==2){
						lCo = Cohesion2(s1o, s2o, aC2, nh[k], 1);
						lCn = Cohesion2(s1n, s2n, aC2, nh[k]+1, 1);
						lC1 = Cohesion2(&s1p[pp],&s2p[pp], aC2, 1, 1);
					}

					if(*cohesion==3){
						lCo = Cohesion3_4(s1o, s2o, mu0c, k0, v0, L0, nh[k], 3, 1);
						lCn = Cohesion3_4(s1n, s2n, mu0c, k0, v0, L0, nh[k]+1,3, 1);
						lC1 = Cohesion3_4(&s1p[pp],&s2p[pp], mu0, k0, v0, L0, 1, 3, 1);
					}

					if(*cohesion==4){
						lCo = Cohesion3_4(s1o, s2o, mu0c, k0, v0, L0, nh[k], 4, 1);
						lCn = Cohesion3_4(s1n, s2n, mu0c, k0, v0, L0, nh[k]+1,4, 1);
						lC1 = Cohesion3_4(&s1p[pp],&s2p[pp], mu0c, k0, v0, L0, 1, 4, 1);
					}




					ph[k] = ((double) nh[k])*exp(lCn-lCo);

					if(*PPM) ph[k] = ((double) nh[k]);

				}

				ph[nclus_iter] = (Mdp)*exp(lC1);

				if(*PPM) ph[nclus_iter] = (Mdp);

				denph = 0.0;
				for(k = 0; k < nclus_iter+1; k++){

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
					sig2pred = runif(0.0,ms);
				}



				ppred_iter[pp] = rnorm(mupred, sqrt(sig2pred));



			}

		}






		/////////////////////////////////////////////
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

				like[ii*(*nobs) + j] = like_iter[j];
				ispred[ii*(*nobs) + j] = mu[Si_iter[j]-1];

			}

			for(pp = 0; pp < *npred; pp++){

				ppred[ii*(*npred) + pp] = ppred_iter[pp];

			}


			ii = ii+1;

		}

/**/

	}




	PutRNGstate();

}


