#include "matrix.h"
#include "Rutil.h"

#include <R_ext/Lapack.h>
#include <R.h>
#include <Rmath.h>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>


void rppmx(int *N, int *similarity_function, int *simparm, double *M, int *ncon, int *ncat,
           double *xcon, int *xcat, int *Cvec, int *ppm,
           double *m0, double *k0, double *v0, double *s20, double *v, double *dirweights,
           double *alpha, int *Si, int *nk, int *nh){
/**************************************************************************************************
 * Function that generates draws from a product partition model prior
 *
 * Inputs:
 *
 * similarity - an integer indicating which similarity function to use
	1 - auxiliary model
	2 - double dipper
	3 - alpha*exp(-variance)
 * simparm indicates which parametrization of similarity to use
  1 - NN model
  2 - NNIG model
 * M - DP scale parameter
 * N - integer that indicates number of observations
 * ncon - integer indicating number of continuous covariates
 * ncat - integer indicating number of categorical covariates
 * xcon - m x ncon matrix containing continuous covariate values
 * xcat - m x ncat integer matrix containing categorical covariate values
 * Cvec = 1 x ncat integer indicating the number of categories for each categorical variable
 * ppm = logical indicating if data comes from PPM (ppm=1) or PPMx (ppm=0)
 * m0 -  double holding mean of NN NIG aux and DD
 * k0 - double holding variance scale parameter (number of obs apriori)
 * v0 - double holding degrees of freedom (number of obs apriori)
 * s20 - double holding scale parameter
 * v - double for "likelihood" variance when NN is used and variance known
 * dirweights = max(cvec) x 1 vector of 0.1 as prior for Dir-Mult similarity
 * alpha - tuning parameter associated with alpha*exp(-variance) similarity
 *
 * Outputs:
 * Si - nx1 scratch array of contiguous memory that holds partition of n objects
 * nk - an integer indicating number of clusters
 * nh - an nkx1 scratch vector that holds number of subjects per cluster.
 *
 *************************************************************************************************/


	int i, ii, k, p, c, iaux, nhtmp;
	double maxph, cprobh, denph, uu, xi;
	double sumx0, sumx20, sumx, sumx2;
	double lgconY, lgconN, lgcatY, lgcatN, lgcon1=0.0, lgcat1=0.0;
	double *ph = R_Vector(*N);
	double *probh = R_Vector(*N);
	int nhc0[*N], nhc[*N], nhc1[*N];
	for(i=0;i<*N;i++){
		Si[i]=0;
		nh[i]=0;
	}
	Si[0] = 1;
	nh[0] = 1;
	nk[0] = 1; // identifies the number of clusters
//	RprintVecAsMat("s1", s1, 1, m);
//	RprintVecAsMat("s2", s2, 1, m);

	double *mnmle = R_Vector(*ncon);
	double *s2mle = R_Vector(*ncon);
	for(p = 0; p < *ncon; p++){
		sumx = 0.0, sumx2=0.0;
		for(ii = 0; ii < *N; ii++){
			sumx = sumx + xcon[ii*(*ncon) + p];
			sumx2 = sumx2 + xcon[ii*(*ncon) + p]*xcon[ii*(*ncon) + p];
		}

		mnmle[p] = sumx/((double) *N);
		s2mle[p] = sumx2/((double) *N) - mnmle[p]*mnmle[p];
	}



//	RprintVecAsMat("xcatstd", xcatstd, *N, *ncat);

	for(i = 1; i < *N; i++){
//		Rprintf("i = %d ====================\n", i);

		for(k=0; k < nk[0]; k++){
//			Rprintf("k = %d =========\n", k);


			lgconY = 0.0;
			lgconN = 0.0;
			lgcon1 = 0.0;


			for(p=0; p<(*ncon); p++){
//				Rprintf("p = %d ====== \n", p) ;
				nhtmp = 0;
				sumx0 = 0.0;
				sumx20 = 0.0;
//				Rprintf("nhtmp = %d\n", nhtmp);
				for(ii = 0; ii < i; ii++){
//					Rprintf("ii = %d ===\n", ii);
//					Rprintf("Si = %d\n", Si[ii]);
					if(Si[ii] == k+1){
						sumx0 = sumx0 + xcon[ii*(*ncon)+p];
						sumx20 = sumx20 + xcon[ii*(*ncon)+p]*xcon[ii*(*ncon)+p];
						nhtmp = nhtmp+1;
//						Rprintf("nhtmp = %d\n", nhtmp);
					}
				}

				xi = xcon[i*(*ncon)+p];

				sumx = sumx0 + xi;
				sumx2 = sumx20 + xi*xi;

				if(nhtmp > 0){
					if(*similarity_function==1){ // Auxilliary
					  if(*simparm==1){// NN model
  						lgconN = lgconN + gsimconNN(*m0, *v, *s20, sumx0, sumx20, mnmle[p], nhtmp, 0, 0, 1);
	  					lgconY = lgconY + gsimconNN(*m0, *v, *s20, sumx, sumx2, mnmle[p], nhtmp+1, 0, 0, 1);
		  				lgcon1 = lgcon1 + gsimconNN(*m0, *v, *s20, xi, xi*xi, mnmle[p], 1, 0, 0, 1);
					  }
					  if(*simparm==2){// NNIG
  						lgconN = lgconN + gsimconNNIG(*m0, *k0, *v0, *s20, sumx0, sumx20, mnmle[p], s2mle[p], nhtmp, 0, 0, 1);
	  					lgconY = lgconY + gsimconNNIG(*m0, *k0, *v0, *s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp+1, 0, 0, 1);
		  				lgcon1 = lgcon1 + gsimconNNIG(*m0, *k0, *v0, *s20, xi, xi*xi, mnmle[p],s2mle[p], 1, 0, 0, 1);
					  }

					}
					if(*similarity_function==2){ //Double Dipper
					  if(*simparm==1){// NN model
  						lgconN = lgconN + gsimconNN(*m0, *v, *s20, sumx0, sumx20, mnmle[p], nhtmp, 1, 0, 1);
	  					lgconY = lgconY + gsimconNN(*m0, *v, *s20, sumx, sumx2, mnmle[p], nhtmp+1, 1, 0, 1);
		  				lgcon1 = lgcon1 + gsimconNN(*m0, *v, *s20, xi, xi*xi, mnmle[p], 1, 1, 0, 1);
					  }

					  if(*simparm==2){// NNIG
  						lgconN = lgconN + gsimconNNIG(*m0, *k0, *v0, *s20, sumx0, sumx20, mnmle[p], s2mle[p], nhtmp, 1, 0, 1);
	  					lgconY = lgconY + gsimconNNIG(*m0, *k0, *v0, *s20, sumx, sumx2, mnmle[p], s2mle[p], nhtmp+1, 1, 0, 1);
		  				lgcon1 = lgcon1 + gsimconNNIG(*m0, *k0, *v0, *s20, xi, xi*xi, mnmle[p],s2mle[p], 1, 1, 0, 1);
					  }

					}
					if(*similarity_function==3){ // variance
						lgconN = lgconN + gsimconEV(sumx0, sumx20, nhtmp, *alpha,1);
						lgconY = lgconY + gsimconEV(sumx, sumx2, nhtmp+1, *alpha,1);
						lgcon1 = lgcon1 + gsimconEV(xi, xi*xi, 1, *alpha,1);
					}

//					Rprintf("lgconN = %f\n", lgconN);
//					Rprintf("lgconY = %f\n", lgconY);
				}

			}


			lgcatY=0.0;
			lgcatN=0.0;
			lgcat1=0.0;

			for(p=0; p<(*ncat); p++){
//				Rprintf("p = %d ====== \n", p) ;
				for(c=0;c<Cvec[p];c++){
					nhc0[c]=0;
					nhc1[c]=0;
					nhc[c]=0;
				}

				nhtmp = 0;
				for(ii = 0; ii < i; ii++){
//					Rprintf("jj = %d\n", jj);
					if(Si[ii] == k+1){
						nhc0[xcat[ii*(*ncat)+p]] = nhc0[xcat[ii*(*ncat)+p]] + 1; // this needs to be a vectore
						nhc[xcat[ii*(*ncat)+p]] = nhc[xcat[ii*(*ncat)+p]] + 1; // this needs to be a vectore
						nhtmp = nhtmp+1;

					}
				}

				nhc[xcat[i*(*ncat)+p]] =  nhc[xcat[i*(*ncat)+p]] + 1;
				nhc1[xcat[i*(*ncat)+p]] =  nhc1[xcat[i*(*ncat)+p]] + 1;


//				Rprintf("sumx = %f\n", sumx);
//				Rprintf("sumx2 = %f\n", sumx2);
//				RprintIVecAsMat("nhc", nhc, 1, Cvec[p]);
				if(nhtmp >0){// Auxilliary
					if(*similarity_function==1){
						lgcatN = lgcatN + gsimcatDM(nhc0, dirweights, Cvec[p], 0, 1);
						lgcatY = lgcatY + gsimcatDM(nhc,  dirweights, Cvec[p], 0, 1);
						lgcat1 = lgcat1 + gsimcatDM(nhc1,  dirweights, Cvec[p], 0, 1);
					}
					if(*similarity_function==2){// Double Dipper
						lgcatN = lgcatN + gsimcatDM(nhc0, dirweights, Cvec[p], 1, 1);
						lgcatY = lgcatY + gsimcatDM(nhc,  dirweights, Cvec[p], 1, 1);
						lgcat1 = lgcat1 + gsimcatDM(nhc1,  dirweights, Cvec[p], 1, 1);
					}
					if(*similarity_function==3){// Use entropy
					  lgcatN = lgcatN + 0.0;
					  lgcatY = lgcatY + 0.0;

					  for(c=0;c<Cvec[p];c++){
					    lgcat1 = lgcat1 + -(*alpha)*0;
              if(nhc0[c]==0){
                 lgcatN = lgcatN + 0.0;
              }else{
                 lgcatN = lgcatN +  -((double) nhc0[c]/(double) nhtmp)*(
                      	            log((double) nhc0[c]/(double) nhtmp)/log(2));
              }

              if(nhc[c]==0){
                 lgcatN = lgcatY + 0.0;
              }else{
                 lgcatN = lgcatY +  -((double) nhc[c]/(double) (nhtmp+1))*(
                      	            log((double) nhc[c]/(double) (nhtmp+1))/log(2));
              }

					  }
				  }

			  }
			}

			ph[k] = log((double) nh[k]) +
		        	lgcatY - lgcatN +
					    lgconY - lgconN;

			if(*ppm) ph[k] = log((double) nh[k]);

//			RprintVecAsMat("ph", ph, 1, k+1);

		}

		ph[nk[0]] = log(*M) + lgcon1 + lgcat1;

		if(*ppm) ph[nk[0]] = log(*M);

//		RprintVecAsMat("ph", ph, 1, nk[0]+1);

		maxph = ph[0];
		for(k = 1; k < nk[0]+1; k++){
			if(maxph < ph[k]) maxph = ph[k];
		}
//		Rprintf("maxph = %f\n", maxph);


//		RprintVecAsMat("ph", ph, 1, nk[0]+1);

		denph = 0.0;
		for(k = 0; k < nk[0]+1; k++){

			ph[k] = exp(ph[k] - maxph);
//			ph[k] = pow(exp(ph[k] - maxph), (1 - exp(-0.0001*(i+1))));
			denph = denph + ph[k];

		}

//		RprintVecAsMat("ph", ph, 1, nk[0]+1);

		for(k = 0; k < nk[0]+1; k++){

			probh[k] = ph[k]/denph;

		}
//		Rprintf("denph = %f\n", denph);

//		RprintVecAsMat("probh", probh, 1, nk[0]+1);

		uu = runif(0.0,1.0);
//		Rprintf("uu = %f\n", uu);

		cprobh= 0.0;
    iaux=nk[0]+1;
		for(k = 0; k < nk[0]+1; k++){

			cprobh = cprobh + probh[k];

			if (uu < cprobh){

				iaux = k+1;
				break;
			}
		}

		if(iaux <= nk[0]){

			Si[i] = iaux;
			nh[Si[i]-1] = nh[Si[i]-1] + 1;

		}else{

			nk[0] = nk[0] + 1;
			Si[i] = nk[0];
			nh[Si[i]-1] = 1;

		}
//		RprintIVecAsMat("Si", Si, 1, *N);
//		Rprintf("nk[0] = %d\n", nk[0]);
//		RprintIVecAsMat("nh", nh, 1, nk[0]);
	}



}
