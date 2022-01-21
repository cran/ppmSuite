/****************************************************************************************
*
*
* My attempt to produce and MCMC algorith for a finite gaussian mixture using 
* the .Call function
*
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

#include <math.h>
#include <stdio.h>
#include <time.h>



// inputs
// y - data vector
// n - length of y
// N - number of mixture components
// m - prior mean of the mean from mixture components
// v - prior variance of the mean from mixture components
// a - prior shape of the variance from mixture components
// b - prior rate of the variance from mixture components
// alpha - prior shape parameters for mixture component weights (dirichlet)
// niter - number of MCMC iterates
// nburn - number of MCMC iterates to be discarded
// nthin - amount of thinning applied to the MCMC chain
//
// outputs
// mu - matrix containing MCMC draws for component means
// sigma2 - matrix containing MCMC draws for component variances
// pi - matrix containing MCMC draws for component weights
// z - matrix containing MCMC draws for component labels


static void ordinal_ppmx(
                  int *y, double *co, int *nobs, int *nordcat,
                  double *Xcon, int *ncon,
                  int *Xcat,  int *ncat, int *Cvec,
                  int *npred,
                  double *Xconp,
                  int *Xcatp,
                  double *M, int *meanModel, double *modelPriors, double *simParms,
                  int *PPM, int *cohesion, int *similarity_function, int *consim,
                  double *dissimtn, double *dissimtt, int *calibrate, double * mh,
                  int *verbose, int *draws, int *burn, int *thin,
                  int *Si, int *nclus, double *mu, double *sig2, double *beta, double *zi,
                  double *mu0, double *sig20,
                  double *like, double *waic, double *lpml,
                  double *ispred, int *isordpred, double *ppred, int *ordppred, int *predclass,
                  double *rbpred, int *rbordpred,
                  double *predclass_prob){ 



  // i - MCMC index
  // j - individual index
  // jj - second individual index (for double for loops)
  // jjj - third individual index
  // c - categorical variable index
  // p - number of covariates index
  // pp - second covariate index
  // k - cluster index
  // ii - save MCMC iterates index
  // b - xb index
  // bb - second xb index


  int i, j, jj, jjj, c, p, pp, k, ii, b, bb;
  int nout = (*draws - *burn)/(*thin);
  int ncov = (*ncon) + (*ncat);
 
  
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
  // zero mean and unit standard deviation
  double *mnmle = R_VectorInit(*ncon, 0.0);
  double *s2mle = R_VectorInit(*ncon, 1.0);
  double sum, sum2;
  if(!(*PPM)){
    for(p = 0; p < *ncon; p++){
      sum = 0.0, sum2=0.0;
      for(j = 0; j < *nobs; j ++){
        sum = sum + Xcon[j*(*ncon) + p];
        sum2 = sum2 + Xcon[j*(*ncon) + p]*Xcon[j*(*ncon) + p];
      }
      for(pp = 0; pp < *npred; pp++){
        sum = sum + Xconp[pp*(*ncon) + p];
        sum2 = sum2 + Xconp[pp*(*ncon) + p]*Xconp[pp*(*ncon) + p];
      }
      mnmle[p] = sum/((double) (*nobs)+(*npred));
      s2mle[p] = sum2/((double) (*nobs)+(*npred)) - mnmle[p]*mnmle[p];
    }
  }

  // Create Xmat to update the covariates from the
  // global regression that is included in the data model
  // Only used if *meanModel == 2;
  double *fullXmat = R_VectorInit((*nobs)*(ncov), 0.0);
  double *fullXmatp = R_VectorInit((*npred)*(ncov), 0.0);

  if(!(*PPM)){
    for(j = 0; j < *nobs; j++){
	  for(b = 0; b < ncov; b++){
	    if(b < *ncon){
          fullXmat[j*ncov+b] = Xcon[j*(*ncon)+b];
//          Xcon[j*(*ncon) + b] = (Xcon[j*(*ncon)+b] - mnmle[b])/sqrt(s2mle[b]);
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
//          Xconp[pp*(*ncon) + b] = (Xconp[pp*(*ncon)+b] - mnmle[b])/sqrt(s2mle[b]);
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
  double *_beta = R_VectorInit(ncov,0.0);
  double *_zi = R_VectorInit(*nobs,0.0);
  double *_muh = R_VectorInit(*nobs,0.0);
  double *_sig2h = R_VectorInit(*nobs,0.5*modelPriors[2]);
  double *_like = R_VectorInit(*nobs,0.0);

  int _nclus=0;
  int _Si[*nobs], nh[*nobs];
    
  for(j = 0; j < *nobs; j++){
    _zi[j] = rtnorm(0, 1, co[y[j]], co[y[j]+1]);
    _Si[j] = rbinom(2, 0.25) + 1;
    nh[j] = 0; // initialize the number of units in each cluster
  }
       
//  RprintIVecAsMat("_Si", _Si, 1, *nobs);       
  
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
  int _isordpred[*nobs];
  
  double *_ppred = R_Vector((*npred));
  double *_rbpred = R_Vector((*npred));
  int _ordppred[*npred];
  int _rbordpred[*npred];
  int _predclass[*npred];
  double *_predclass_prob = R_VectorInit((*npred)*(*nobs+1), 0.0);


  
  // ===================================================================================
  //
  // scratch vectors of memory needed to update parameters
  //
  // ===================================================================================

  // stuff that I need to update Si (cluster labels);
  int iaux, auxint;
  int nhctmp[max_C];
  double auxreal, sumxtmp, sumx2tmp, npdN,npdY,npd, mn, xcontmp, uu, xb;
  double mudraw, sdraw, maxph, denph, cprobh, sdens, s2tmp;
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
  int  nhc[(*nobs)*(*ncat)*max_C];

  // stuff I need to update muh
  double *sumz = R_VectorInit((*nobs),0.0);
  double *sumz2 = R_VectorInit((*nobs),0.0);
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
        sumx[(_Si[j]-1)*(*ncon) + p] = sumx[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p];
        sumx2[(_Si[j]-1)*(*ncon) + p] = sumx2[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
      }
      for(p=0; p<*ncat; p++){
        nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
          nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] + 1;
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
    RprintVecAsMat("co", co, 1, *nordcat);
    if(*ncat > 0 ) RprintIVecAsMat("Cvec",Cvec, 1, *ncat);
    Rprintf("Prior values: m = %.2f, s2 = %.2f smin = %.2f, smax = %.2f, s0min = %.2f, s0max = %.2f\n", m, s2, smin, smax, s0min, s0max);
    Rprintf("M = %f\n", Mdp);
    Rprintf("Simlarity values: m0 = %.2f, s20 = %.2f v = %.2f, k0 = %.2f, nu0 = %.2f, a0 = %.2f, alpha = %.2f\n", m0, s20, v, k0, nu0, simParms[5], alpha);
    if(*nobs > 50) RprintIVecAsMat("First fifty response values",y, 1, 50);
  }
  
  
  

  ii = 0;
  
  // ===================================================================================
  //
  // Beginning of the MCMC loop
  //
  // ===================================================================================
  double calc_time = 0.0;
  clock_t  begin = clock();

  for(i=0; i<*draws; i++){
  
//    Rprintf("i = %d\n", i);

    if(*verbose){
      clock_t ith_iterate = clock();
      calc_time = (ith_iterate - begin)/CLOCKS_PER_SEC;

      Rprintf("Progress:%.1f%%, Time:%.1f seconds\r", ((double) (i+1) / (double) (*draws))*100.0, calc_time);
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


      if(nh[_Si[j]-1] > 1){
        
        // Observation belongs to a non-singleton ...
        nh[_Si[j]-1] = nh[_Si[j]-1] - 1;
        
        if(!(*PPM)){
          // need to reduce the sumx sumx2 to
          for(p = 0; p < *ncon; p++){
            sumx[(_Si[j]-1)*(*ncon) + p] = sumx[(_Si[j]-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
            sumx2[(_Si[j]-1)*(*ncon) + p] = sumx2[(_Si[j]-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          }
        
          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){
        
            nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
               nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;
        
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
            for(p = 0; p < *ncon; p++){
              auxreal = sumx[(iaux-1)*(*ncon) + p];
              sumx[(iaux-1)*(*ncon) + p] = sumx[(_nclus-1)*(*ncon) + p];
              sumx[(_nclus-1)*(*ncon) + p] = auxreal;
      
              auxreal = sumx2[(iaux-1)*(*ncon) + p];
              sumx2[(iaux-1)*(*ncon) + p] = sumx2[(_nclus-1)*(*ncon) + p];
              sumx2[(_nclus-1)*(*ncon) + p] = auxreal;
             
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
      
  
        // need to reduce the sumx sumx2
        if(!(*PPM)){
          for(p = 0; p < *ncon; p++){
            sumx[(_nclus-1)*(*ncon) + p] = sumx[(_nclus-1)*(*ncon) + p] - Xcon[j*(*ncon)+p];
            sumx2[(_nclus-1)*(*ncon) + p] = sumx2[(_nclus-1)*(*ncon) + p] - Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
          }
         
          // need to reduce the nhc
          for(p = 0; p < *ncat; p++){
      
             nhc[((_nclus-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] =
                nhc[((_nclus-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] - 1;
         
           }
        }
  
  		// finally reduce the number of clusters
      	_nclus = _nclus - 1;
      
      }
      
      if(*meanModel==2){
        xb = 0.0;
        for(b = 0; b < ncov; b++){
          xb = xb + fullXmat[j*ncov+b]*_beta[b];
        }
      }
      
      
      // The atoms have been relabeled if necessary and now we need to update Si.
      
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

      	    if(*similarity_function==1){ // Auxilliary
      	      if(*consim==1){ // NN
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 0, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
      	      if(*consim==2){// NNIG
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k], 0, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
          
      	    }
      	    if(*similarity_function==2){ //Double Dipper
      	      if(*consim==1){
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 1, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
      	      if(*consim==2){
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k], 1, 0, 1);
      	        lgconN = lgconN + lgcont;
      	      }
      	    }
      	    if(*similarity_function==3){ // variance
      	      lgcont = gsimconEV(sumxtmp, sumx2tmp, nh[k], alpha,1);
      	      lgconN = lgconN + lgcont;
      	    }
          
          
      	    // now add jth individual back;
      	    sumxtmp = sumxtmp + Xcon[j*(*ncon)+p];
      	    sumx2tmp = sumx2tmp + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];

          
      	    if(*similarity_function==1){ // Auxilliary
      	      if(*consim==1){
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k]+1, 0, 0, 1);
      	        lgconY = lgconY + lgcont;
      	      }
      	      if(*consim==2){
      	        lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k]+1, 0, 0, 1);
      	        lgconY = lgconY + lgcont;
      	      }
      	    }
      	    if(*similarity_function==2){ //Double Dipper
      	      if(*consim==1){
      	        lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k]+1, 1, 0, 1);
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
        if(*meanModel==2) mn = _muh[k] + xb;
      
      
      	// Compute the unnormalized cluster probabilities
      	// Note that if PPMx = FALSE then 
      	// lgcatY = lgcatN = lgconY = lgconN = 0;
      	ph[k] = dnorm(_zi[j], mn, sqrt(_sig2h[k]), 1) +
      		    	log((double) nh[k]) + // cohesion part
      		        lgcatY - lgcatN + // Categorical part only nonzero if PPMx=TRUE
      				lgconY - lgconN;  // Continuous part only nonzero if PPMx=TRUE
      
      
      
      	if(*calibrate == 2){
      		ph[k] = dnorm(_zi[j], mn, sqrt(_sig2h[k]), 1) +
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
        
        if(*similarity_function==4){
          lgcatdraw = -(alpha)*0;
        }
        
        gtilY[_nclus] = lgcondraw + lgcatdraw;
        gtilN[_nclus] = lgcondraw + lgcatdraw;
      } // THIS ENDS THE PPMX PART.
      
      mn = mudraw;
      if(*meanModel == 2) mn = mudraw + xb;
      
      // Note that if PPMx = FALSE, then 
      // lgcondraw = lgcondraw =  0;
      
      
      ph[_nclus] = dnorm(_zi[j],mn,sdraw,1) +
                       	log(Mdp) +
                       	lgcondraw +
                       	lgcatdraw;
      
      if(*calibrate==2){
      	ph[_nclus] = dnorm(_zi[j],mn,sdraw,1) +
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
      		if(*meanModel == 2) mn = _muh[k] + xb;
      
      		ph[k] = dnorm(_zi[j], mn, sqrt(_sig2h[k]), 1) +
              	    	log((double) nh[k]) +  // Cohesion part
      				    lgtilYk - lgtilNk; //This takes into account both cont and cat vars
      
      		if(*cohesion==2){
      	    	ph[k] = ph[k] - log((double) nh[k]);
          	}
      
      	}
      
      	// calibration for a singleton
      	mn = mudraw;
      	if(*meanModel == 2) mn = mudraw + xb;
      
      	ph[_nclus] = dnorm(_zi[j],mn,sdraw,1) +
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
          sumx[(_Si[j]-1)*(*ncon) + p] = sumx[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p];
          sumx2[(_Si[j]-1)*(*ncon) + p] = sumx2[(_Si[j]-1)*(*ncon) + p] + Xcon[j*(*ncon)+p]*Xcon[j*(*ncon)+p];
        }
      
        // need to now add the xcat to the cluster to which it was assigned; 
        for(p = 0; p < *ncat; p++){
      
          nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] = 
             nhc[((_Si[j]-1)*(*ncat) + p)*(max_C) + Xcat[j*(*ncat)+p]] + 1;
      
        }
      } 
      
      
      if((i >= (*burn)) & ((i) % *thin ==0)){
        /////////////////////////////////////////////
        //
        // Compute the CPO and lpml using the mixture
        //
        /////////////////////////////////////////////
        sdens = 0.0;
        for(k = 0; k < _nclus; k++){
      	  s2tmp = _sig2h[k];
      
      	  mn = _muh[k];
      	  if(*meanModel==2) mn = _muh[k] + xb;
      
      	  sdens = sdens + probh[k]*(pnorm(co[y[j]+1], mn, sqrt(s2tmp), 1, 0) -
      				                pnorm(co[y[j]],   mn, sqrt(s2tmp), 1, 0));
        }
      
      
        mn = mudraw;
        if(*meanModel==2) mn = mudraw + xb;
      
        sdens = sdens + probh[_nclus]*(pnorm(co[y[j]+1], mn, sdraw, 1, 0) -
      					               pnorm(co[y[j]],   mn, sdraw, 1, 0));
      
        _like[j] = sdens;
      
      
      	// These are needed for WAIC
      	mnlike[j] = mnlike[j] + (_like[j])/(double) nout;
      	mnllike[j] = mnllike[j] + log(_like[j])/(double) nout;
      
      	CPOinv[j] = CPOinv[j] + (1/(double) nout)*(1/_like[j]);
      } 
      
      
    }
    

//    RprintIVecAsMat("Si", _Si, 1, *nobs);
//    RprintIVecAsMat("nh", nh, 1, _nclus);
//    Rprintf("nclus = %d\n", _nclus);
    
    //////////////////////////////////////////////////////////////////////////////////
    //
    // update the the latent variables z_i using a truncated normal
    //
    //////////////////////////////////////////////////////////////////////////////////

    for(k = 0; k < _nclus; k++){
      sumz[k] = 0.0;
      sumz2[k] = 0.0;
    }


    for(j = 0; j < *nobs; j++){
    

      if(*meanModel==1){
        mn = _muh[_Si[j]-1];
      }
      if(*meanModel==2){
        xb = 0.0;
        for(b = 0; b < ncov; b++){
          xb = xb + fullXmat[j*ncov+b]*_beta[b];
        }
        mn =  _muh[_Si[j]-1] + xb;
      }

      _zi[j] = rtnorm(mn, sqrt(_sig2h[_Si[j]-1]), co[y[j]], co[y[j]+1]);
//      _zi[j] = r_truncnorm(co[y[j]], co[y[j]+1], _muh[_Si[j]-1], sqrt(_sig2h[_Si[j]-1]));

      // sum with-in this for loop to save some computing time;
      // This are used to update muh which is next
      if(*meanModel==1){
        sumz[_Si[j]-1] = sumz[_Si[j]-1] + _zi[j];
        sumz2[_Si[j]-1] = sumz2[_Si[j]-1] + _zi[j]*_zi[j];
      }
      if(*meanModel==2){
        sumz[_Si[j]-1] = sumz[_Si[j]-1] + (_zi[j] - xb);
        sumz2[_Si[j]-1] = sumz2[_Si[j]-1] + (_zi[j] - xb)*(_zi[j] - xb);
      
      }
    }

//    RprintVecAsMat("sumz", sumz, 1, _nclus);

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
        
        llo = -(nh[k])*log(osig) - 1/(2*osig*osig)*(sumz2[k] - 2*_muh[k]*sumz[k] + nh[k]*_muh[k]*_muh[k]);
        lln = -(nh[k])*log(nsig) - 1/(2*nsig*nsig)*(sumz2[k] - 2*_muh[k]*sumz[k] + nh[k]*_muh[k]*_muh[k]);
        
        llo = llo + dunif(osig, smin, smax, 1);
        lln = lln + dunif(nsig, smin, smax, 1);
        
        llr = lln - llo;
        uu = runif(0,1);
        
        if(log(uu) < llr){
        	_sig2h[k] = nsig*nsig;
        }
      }
    }


//    RprintVecAsMat("sig2h", _sig2h, 1, _nclus);

    //////////////////////////////////////////////////////////////////////////////////
    //
    // update muh's cluster specific means
    //
    //////////////////////////////////////////////////////////////////////////////////
    summu = 0.0, summu2=0.0;
    for(k = 0; k < _nclus; k++){
    
      s2star = 1/((double) nh[k]/_sig2h[k] + 1/_sig20);
      mstar = s2star*( (1/_sig2h[k])*sumz[k] + (1/_sig20)*_mu0);

      _muh[k] = rnorm(mstar, sqrt(s2star));
      
      summu = summu + _muh[k];  // This is used in the updating of mu0 and sig20 
      summu2 = summu2 + _muh[k]*_muh[k];  // This is used in the updating of sig20
    }

//    RprintVecAsMat("muh", _muh, 1, _nclus);


    //////////////////////////////////////////////////////////////////////////////////////
    //
    // Update mu0  prior mean of muh
    //
    //////////////////////////////////////////////////////////////////////////////////////
    
    s2star = 1/(((double) _nclus/_sig20) + (1/s2));
    mstar = s2star*((1/_sig20)*summu + (1/s2)*m);
    
    _mu0 = rnorm(mstar, sqrt(s2star));
    

//    Rprintf("mu0 = %f\n", _mu0);

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

//    Rprintf("sig20 = %f\n", _sig20);


  
//    Rprintf("meanModel = %d\n", *meanModel);

//    RprintVecAsMat("Sstar", Sstar, ncov, ncov);

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
      
      	  sumyxbt[b] = sumyxbt[b] +  (1/_sig2h[_Si[j]-1])*fullXmat[j*ncov+b]*
      					                                  (_zi[j] - _muh[_Si[j]-1]);
      
      	  for(bb = 0; bb < ncov; bb++){
      
            sumXXp[b*ncov+bb] = sumXXp[b*ncov+bb] + (1/_sig2h[_Si[j]-1])*
      						                         fullXmat[j*ncov+b]*fullXmat[j*ncov+bb];
      
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
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    //
    // in sample prediction to assess model fit
    //
    ////////////////////////////////////////////////////////////////////////////////////////////
    if((i >= (*burn)) & ((i) % *thin ==0)){
      
      for(j = 0; j < *nobs; j++){
        
        mn = _muh[_Si[j]-1];
        if(*meanModel == 2){
          xb = 0.0;
          for(b = 0; b < ncov; b++){
            xb = xb + fullXmat[j*ncov+b]*_beta[b];
          }
          
          mn = _muh[_Si[j]-1] + xb;
        }
  
  
        _ispred[j] = rnorm(mn, sqrt(_sig2h[_Si[j]-1]));
  
        for(c=0; c < *nordcat-1; c++){
          if((_ispred[j] > co[c]) & (_ispred[j] < co[c+1])) _isordpred[j] = c;
        }
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
      	  lgcatY=0.0, lgcatN=0.0;
        
      	  if(!(*PPM)){
      	  	
      	    for(p=0; p<(*ncon); p++){
              sumxtmp = sumx[k*(*ncon) + p];
              sumx2tmp = sumx2[k*(*ncon) + p];
          
          
      	  	  if(*similarity_function==1){
      	  	    if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k], 0, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumx2tmp, sumx2[k], mnmle[p], s2mle[p], nh[k], 0, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==2){
      	  	  	if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumx2tmp, sumx2tmp, mnmle[p], nh[k], 1, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k], 1, 0, 1);
      	  	      lgconN = lgconN + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==3){
      	  	    lgcont = gsimconEV(sumxtmp, sumx2tmp, nh[k], alpha, 1);
      	  	  	lgconN = lgconN + lgcont;
      	  	  }
          
      	  	  // now add ppth prediction to cluster;
      	  	  sumxtmp = sumxtmp + Xconp[pp*(*ncon)+p];
      	  	  sumx2tmp = sumx2tmp + Xconp[pp*(*ncon)+p]*Xconp[pp*(*ncon)+p];
          
      	  	  if(*similarity_function==1){ // Auxilliary
      	  	  	if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k]+1, 0, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k]+1, 0, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==2){ // Double Dipper
      	  	  	if(*consim==1){
      	  	      lgcont = gsimconNN(m0, v, s20, sumxtmp, sumx2tmp, mnmle[p], nh[k]+1, 1, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  	if(*consim==2){
      	  	      lgcont = gsimconNNIG(m0, k0, nu0, s20, sumxtmp, sumx2tmp, mnmle[p], s2mle[p], nh[k]+1, 1, 0, 1);
      	  	      lgconY = lgconY + lgcont;
      	  	  	}
      	  	  }
      	  	  if(*similarity_function==3){ // Variance
      	  	    lgcont = gsimconEV(sumxtmp, sumx2tmp, nh[k]+1,alpha, 1);
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
      	  	  nhctmp[Xcatp[pp*(*ncat)+p]] = nhctmp[Xcatp[pp*(*ncat)+p]] + 1;
          
          
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
        
      	  }
      	  	
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
          
         
          for(p=0;p<(*ncat);p++){
            for(c=0;c<Cvec[p];c++){nhctmp[c] = 0;}
          
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
        iaux = _nclus+1;
        for(k = 0; k < _nclus+1; k++){
          cprobh = cprobh + probh[k];
          if (uu < cprobh){
            iaux = k+1;
            break;
          }
        }
        
       
        if(iaux <= _nclus){
          mudraw = _muh[(iaux-1)];
          sdraw = sqrt(_sig2h[(iaux-1)]);
        }else{
          mudraw = rnorm(_mu0,sqrt(_sig20));
          sdraw = runif(smin, smax);
        }
        
        mn = mudraw;
        if(*meanModel==2){
          xb = 0.0;
          for(b = 0; b < ncov; b++){
            xb = xb + fullXmatp[pp*ncov+b]*_beta[b];
          }
          mn = mudraw + xb;
        }
        
        _ppred[pp] = rnorm(mn, sdraw);
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
    
        // convert the auxiliary variable to the ordinal scale
        for(c=0; c < *nordcat-1; c++){
          if((_ppred[pp] > co[c]) & (_ppred[pp] < co[c+1])) _ordppred[pp] = c;
          if((_rbpred[pp] > co[c]) & (_rbpred[pp] < co[c+1])) _rbordpred[pp] = c;
        }
        
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
        mu[ii + nout*j] = _muh[_Si[j]-1];
        sig2[ii + nout*j] = _sig2h[_Si[j]-1];
        Si[ii + nout*j] = _Si[j];
        zi[ii + nout*j] = _zi[j];


        like[ii + nout*j] = _like[j];
        ispred[ii + nout*j] = _ispred[j];
        isordpred[ii + nout*j] = _isordpred[j];
        
      }

      for(b = 0; b < ncov; b++){
        beta[ii + nout*b] = _beta[b];
      }      

      for(pp = 0; pp < *npred; pp++){
        ppred[ii + nout*pp] = _ppred[pp];
        predclass[ii + nout*pp] = _predclass[pp];
		ordppred[ii + nout*pp] = _ordppred[pp];
        
        rbpred[ii + nout*pp] = _rbpred[pp];
        rbordpred[ii + nout*pp] = _rbordpred[pp];
      }
      for(pp = 0; pp < (*nobs)*(*npred); pp++){
        predclass_prob[ii + nout*pp] = _predclass_prob[pp];
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




SEXP ORDINAL_PPMX(SEXP y, SEXP co, SEXP nobs, SEXP nordcat,
                  SEXP Xcon, SEXP ncon,
                  SEXP Xcat, SEXP ncat, SEXP Cvec,
                  SEXP npred,
                  SEXP Xconp,
                  SEXP Xcatp,
                  SEXP M, SEXP meanModel, SEXP modelPriors, SEXP simParms,
                  SEXP PPM, SEXP cohesion, SEXP similarity_function, SEXP consim,
                  SEXP dissimtn, SEXP dissimtt, SEXP calibrate, SEXP mh,
                  SEXP verbose, SEXP draws, SEXP burn, SEXP thin) {


  int nprot = 0;
  int _nobs = asInteger(nobs);
  int _ncon = asInteger(ncon);
  int _ncat = asInteger(ncat);
  int _npred = asInteger(npred);
  int _nordcat = asInteger(nordcat);

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
  
  y = PROTECT(coerceVector(y, INTSXP)); nprot++;
  co = PROTECT(coerceVector(co, REALSXP)); nprot++;
  Xcon =  PROTECT(coerceVector(Xcon, REALSXP)); nprot++;
  Xcat =  PROTECT(coerceVector(Xcat, INTSXP)); nprot++;
  Xconp =  PROTECT(coerceVector(Xconp, REALSXP)); nprot++;
  Xcatp =  PROTECT(coerceVector(Xcatp, INTSXP)); nprot++;
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
  SEXP Zi = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++; 
  SEXP NCLUS = PROTECT(allocMatrix(INTSXP, nout, 1)); nprot++;
  SEXP MU0 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;
  SEXP SIG20 = PROTECT(allocMatrix(REALSXP, nout, 1)); nprot++;

  SEXP ISPRED = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++;
  SEXP ISORDPRED = PROTECT(allocMatrix(INTSXP, nout, _nobs)); nprot++;
  SEXP PPRED = PROTECT(allocMatrix(REALSXP, nout, _npred)); nprot++;
  SEXP ORDPPRED = PROTECT(allocMatrix(INTSXP, nout, _npred)); nprot++;
  SEXP RBPRED = PROTECT(allocMatrix(REALSXP, nout, _npred)); nprot++;
  SEXP ORDRBPRED = PROTECT(allocMatrix(INTSXP, nout, _npred)); nprot++;
  SEXP PREDCLASS = PROTECT(allocMatrix(INTSXP, nout, _npred)); nprot++;
  SEXP PREDCLASS_PROB = PROTECT(allocMatrix(REALSXP, nout, (_nobs)*(_npred))); nprot++;


  SEXP LIKE = PROTECT(allocMatrix(REALSXP, nout, _nobs)); nprot++;
  SEXP WAIC = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;
  SEXP LPML = PROTECT(Rf_allocVector(REALSXP, 1)); nprot++;

  double *MUout, *SIG2out, *Ziout, *MU0out, *SIG20out, *BETAout, *LIKEout, *WAICout, *LPMLout;
  double *ISPREDout, *PPREDout, *RBPREDout, *PREDCLASS_PROBout;
  int *Siout, *NCLUSout, *PREDCLASSout, *ISORDPREDout, *ORDPPREDout, *ORDRBPREDout;  

  MUout = REAL(MU);
  SIG2out = REAL(SIG2);
  MU0out = REAL(MU0);
  SIG20out = REAL(SIG20);
  Ziout = REAL(Zi);
  BETAout = REAL(BETA);
  Siout = INTEGER(Si);
  NCLUSout = INTEGER(NCLUS);
  
  LIKEout = REAL(LIKE);
  WAICout = REAL(WAIC);
  LPMLout = REAL(LPML);

  ISPREDout = REAL(ISPRED);
  PPREDout = REAL(PPRED);
  RBPREDout = REAL(RBPRED);
  ISORDPREDout = INTEGER(ISORDPRED);
  ORDPPREDout = INTEGER(ORDPPRED);
  ORDRBPREDout = INTEGER(ORDRBPRED);
  PREDCLASS_PROBout = REAL(PREDCLASS_PROB);
  PREDCLASSout = INTEGER(PREDCLASS);



  GetRNGstate();

  ordinal_ppmx(INTEGER(y), REAL(co), &_nobs, &_nordcat,
               REAL(Xcon), &_ncon,
               INTEGER(Xcat),  &_ncat, INTEGER(Cvec),
               &_npred,
               REAL(Xconp),
               INTEGER(Xcatp),
               &_M, &_meanModel, REAL(modelPriors),  REAL(simParms),
               &_PPM, &_cohesion, &_similarity_function, &_consim,
               REAL(dissimtn), REAL(dissimtt), &_calibrate, REAL(mh),
               &_verbose, &_niter, &_nburn, &_nthin,
               Siout, NCLUSout, MUout, SIG2out, BETAout, Ziout,
               MU0out, SIG20out,
               LIKEout, WAICout, LPMLout,
               ISPREDout, ISORDPREDout, PPREDout, ORDPPREDout, PREDCLASSout,
               RBPREDout, ORDRBPREDout,
               PREDCLASS_PROBout);

  PutRNGstate();



  SEXP ans = PROTECT(allocVector(VECSXP, 19)); nprot++;
  SET_VECTOR_ELT(ans, 0, MU);
  SET_VECTOR_ELT(ans, 1, SIG2);
  SET_VECTOR_ELT(ans, 2, BETA);
  SET_VECTOR_ELT(ans, 3, Si);
  SET_VECTOR_ELT(ans, 4, Zi);
  SET_VECTOR_ELT(ans, 5, MU0);
  SET_VECTOR_ELT(ans, 6, SIG20);
  SET_VECTOR_ELT(ans, 7, NCLUS);
  SET_VECTOR_ELT(ans, 8, LIKE);
  SET_VECTOR_ELT(ans, 9, WAIC);
  SET_VECTOR_ELT(ans, 10, LPML);
  SET_VECTOR_ELT(ans, 11, ISPRED);
  SET_VECTOR_ELT(ans, 12, PPRED);
  SET_VECTOR_ELT(ans, 13, PREDCLASS);
  SET_VECTOR_ELT(ans, 14, RBPRED);
  SET_VECTOR_ELT(ans, 15, PREDCLASS_PROB);
  SET_VECTOR_ELT(ans, 16, ISORDPRED);
  SET_VECTOR_ELT(ans, 17, ORDPPRED);
  SET_VECTOR_ELT(ans, 18, ORDRBPRED);


  SEXP nm = allocVector(STRSXP, 19);
  setAttrib(ans, R_NamesSymbol, nm);
  SET_STRING_ELT(nm, 0, mkChar("mu"));
  SET_STRING_ELT(nm, 1, mkChar("sig2"));
  SET_STRING_ELT(nm, 2, mkChar("beta"));
  SET_STRING_ELT(nm, 3, mkChar("Si"));
  SET_STRING_ELT(nm, 4, mkChar("zi"));
  SET_STRING_ELT(nm, 5, mkChar("mu0"));
  SET_STRING_ELT(nm, 6, mkChar("sig20"));
  SET_STRING_ELT(nm, 7, mkChar("nclus"));
  SET_STRING_ELT(nm, 8, mkChar("like"));
  SET_STRING_ELT(nm, 9, mkChar("WAIC"));
  SET_STRING_ELT(nm, 10, mkChar("lpml"));
  SET_STRING_ELT(nm, 11, mkChar("fitted.values"));
  SET_STRING_ELT(nm, 12, mkChar("ppred"));
  SET_STRING_ELT(nm, 13, mkChar("predclass"));
  SET_STRING_ELT(nm, 14, mkChar("rbpred"));
  SET_STRING_ELT(nm, 15, mkChar("predclass_prob"));
  SET_STRING_ELT(nm, 16, mkChar("ord.fitted.values"));
  SET_STRING_ELT(nm, 17, mkChar("ord.ppred"));
  SET_STRING_ELT(nm, 18, mkChar("ord.rbpred"));

  UNPROTECT(nprot);
  return(ans);
}



