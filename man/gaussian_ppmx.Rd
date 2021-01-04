\name{gaussian_ppmx}
\alias{gaussian_ppmx}
\title{Function that fits Gaussian PPMx model}
\description{
\code{gaussian_ppmx} is the main function used to fit Gaussian PPMx model.
}
\usage{
gaussian_ppmx(y, X=NULL, Xpred=NULL,
                  meanModel=1,
                  cohesion=1,
                  M=1,
                  similarity_function=1,
                  consim=1,
                  calibrate=0,
                  simParms=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1),
                  modelPriors=c(0, 100^2, 1, 1),
                  mh=c(0.5, 0.5),
                  draws=1100,burn=100,thin=1,
                  verbose=FALSE)

}

\arguments{
\item{y}{numeric vector for the response variable}
\item{X}{a data frame whose columns consist of covariates that will be incorporated in the partition model. Those with class of "character" or "factor" will be treated as categorical covariates.  All others will be treated as continuous covariates.  If NULL, then a PPM model is fit.}

\item{Xpred}{a data frame containing covariate values for which out-of-sample predictions are desired.  The format of and order of Xpred must be the same as that found in X.}

\item{meanModel}{Type of mean model included in the likelihood that is to be used.  Options are 1 or 2 with
  \itemize{

    \item  1 - cluster-specific means with no covariates in likelihood.

    \item 2 - cluster-specific intercepts and a global regression of the type Xbeta is included in the likelihood.
  }
}

\item{cohesion}{Type of cohesion function to use in the PPMx prior.  Options are 1 or 2 with
  \itemize{

    \item  1 - Dirichlet process style of cohesion c(S) = M x (|S| - 1)!

    \item  2 - Uniform cohesion c(S) = 1
  }
}
\item{M}{Precision parameter.  Default is 1.}

\item{similarity_function}{Type of similarity function that is employed for the PPMx prior on partitions.  Options are 1-4 with

  \itemize{
    \item  1 - Auxilliary similarity

    \item  2 - Double dipper similarity

    \item  3 - Cluster variance or entropy for categorical covariates

    \item 4 - Mean Gower dissimilarity (Gower dissimilarity is not available if missing values are present in X)
  }
}
\item{consim}{If similarity_function is either set to 1 or 2, then this specifies the type of marginal likelihood used as the similarity function. Options are 1 or 2 with
  \itemize{
     \item 1 - N-N(m0, s20, v) (v variance of ``likelihood'', m0 and s20 ``prior'' parameters),

    \item 2 - N-NIG(m0,s20, k0, nu0) (m0 and s20 center and scale of Gaussian, k0 and nu0 shape and scale of Inverse Gamma)
  }
}
\item{calibrate}{Indicates if the similarity should be calibrated.  Options are 0-2 with
\itemize{
    \item 0 - no calibration

    \item 1 - standardize similarity value for each covariate

    \item 2 - coarsening is applied so that each similarity is raised to the 1/p power}
}
\item{simParms}{Vector of parameter values employed in the similarity function of the PPMx. Entries of the vector correspond to
\itemize{
    \item m0 - center continuous similarity with default 0,

    \item s20 - spread of 'prior' continuous similarity with default 1,

    \item v2 - spread of 'likelihood' for conitnuous similarity (smaller values place more weight partitions with clusters that contain homogeneous  covariate values)

    \item k0 - degrees of freedom upper for v (only used for N-NIG similarity model)

    \item nu0 - scale for v (only used for N-NIG similarity model)

    \item a0 - dirichlet weight for categorical similarity with default of 0.1 (smaller more weight placed on this variable)

    \item alpha - weight associated with cluster-variance and Gower disimilarity
 }
}
\item{modelPriors}{Vector of prior parameter values used in the PPMx prior.
  \itemize{
    \item m - prior mean for mu0 with default equal to 0,

    \item s2 - prior variance mu0 with default equal to 100^2,

    \item A - upper bound on sigma2*_j with default equal to 10

   \item B - upper bound on sig20 with default equal to 10
}
}
\item{mh}{two dimensional vector containing values for tunning parameter associated with MH update for sigma2 and sigma20}

\item{draws}{number of MCMC iterates to be collected. default is 1100}

\item{burn}{number of MCMC iterates discared as burn-in. default is 100}

\item{thin}{number by which the MCMC chain is thinne. default is 1.  Thin must be selected so that it is a multilple of (draws - thin)}

\item{verbose}{Logical indicating if information regarding data and MCMC iterate should be printed to screen}

}

\details{
This function is able to fit a Gaussian PPM or PPMx model as detailed in (Mueller, Quintana, and Rosner, 2011).  The data model is a Gaussian distribution with cluster-specific means and variances.  If meanModel = 2, then a ``global'' regression component is added to the mean.  Conjugate priors are used for cluster-specific means while uniform priors are used for variance components.  A variety of options associated with the similarity function of the PPMx are available.  See Page, Quintana 2018; Mueller, Quintana, Rosner 2011 for more details.

If covariate matrix contains missing values, then the approach described in  Page, Quintana, Mueller (2020) is automatically employed.

All continuous X's are standardized to have mean 0 and unit standard deviation before being passed to the PPMx partition model.  However, for meanModel = 2 the regression coefficients are estimated on the original scale and are ordered such that the continuous covariates appear first and the categorical covariates come after.

The computational implementation of the model is based on algorithm 8 found in Neal 2000.


}

\value{
The function returns a list containing arrays filled with MCMC iterates corresponding to model parameters and model fit metrics. In order to provide more detail, in what follows let

"T" - be the number of MCMC iterates collected,

"N" - be the number of observations,

"P" - be the number of predictions.

"C" - be the total number of covariates

The output list contains the following
\itemize{
\item mu - a matrix  of dimension (T, N) containing MCMC iterates associated with each subjects mean parameter (mu*_{c_i}).

\item sig2 - a matrix of dimension (T, N) containing MCMC iterates associated with each subjects variance parameter (sigma2*_{c_i})

\item beta - if meanModel = 2, then this is a matrix of dimension (T,C) containing MCMC iterates associated coefficients in the global regression

\item Si - a matrix of dimension (T, N) containing MCMC iterates assocated with each subjects cluster label.

\item mu0 - vector of length T containing MCMC iterates for mu0 parameter

\item sig20 - vector of length T containing MCMC iterates for sig20

\item nclus - vector of length T containing number of clusters at each MCMC iterate

\item like - a matrix of dimension (T, N) containing likelihood values at each MCMC iterate.

\item WAIC - scalar containing the WAIC value

\item lpml - scalar containing lpml value

\item fitted.values - a matrix of dimension (T, N) containing fitted (or in sample predictions) for each subject at each MCMC iterate

\item ppred - a matrix of dimension (T, P) containing out of sample preditions for each "new" subject at each MCMC iterate of the posterior predictive distribution

\item predclass - a matrix of dimension (T, P) containing MCMC iterates of cluster two which "new" subject is allocated

\item rbpred - a matrix of dimension (T, P) containing out of sample preditions for each "new" subject at each MCMC iterate based on the rao-blackwellized prediction

\item predclass_prob - a matrix of dimension (T, P*N) that contains the cluster allocation probabilities.  They are organized so that each row corresponds to an MCMC iterate.  Letting nclus represent the number of components at the t-th MCMC sample, the the first nclus columns of row t correspond to probabilities for new subject 1.  Then columns (N+1):(N+nclus) correspond to probabilities associated with new subject 2, etc.
}
}
\examples{

\donttest{
data(bear)

# plot length, sex, and weight of bears
ck <- c(4,3,2)
pairs(bear[,ck])


# response is length
Y <- bear$weight

# Continuous Covariate is chest
# Categorical covariate is sex
X <- bear[,c("length", "sex")]
X$sex <- as.factor(X$sex)

# Randomly partition data into 44 training and 10 testing
set.seed(1)
trainObs <- sample(1:length(Y),44, replace=FALSE)

Ytrain <- Y[trainObs]
Ytest <- Y[-trainObs]

Xtrain <- X[trainObs,,drop=FALSE]
Xtest <- X[-trainObs,,drop=FALSE]

simParms <- c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1)
modelPriors <- c(0, 100^2, 0.5*sd(Y), 100)
M <- 1.0

niter <- 100000
nburn <- 50000
nthin <- 50

nout <- (niter - nburn)/nthin

mh <- c(1,10)

# Run MCMC algorithm for Gaussian PPMx model
out1 <- gaussian_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, M=M,
              meanModel = 1,
		          similarity_function=1,
		          consim=1,
		          calibrate=0,
		          simParms=simParms,
		          modelPriors = modelPriors,
		          draws=niter, burn=nburn, thin=nthin,
		          mh=mh)



# plot a select few posterior distributions
plot(density(out1$mu[,1])) # first observation's mean
plot(density(out1$sig2[,1])) # first observation's variance
plot(table(out1$nc)/nout,type='h') # distribution
plot(density(out1$mu0), type='l')
plot(density(out1$sig20))


# The first partition iterate is used for plotting
# purposes only. We recommended using the salso
# R-package to estimate the partition based on Si
pairs(bear[trainObs,ck],col=out1$Si[1,], pch=out1$Si[1,])


# Compare fit and predictions when covariates are not included
# in the partition model.  That is, refit data with PPM rather than PPMx
out2 <- gaussian_ppmx(y=Ytrain, X=NULL, Xpred=Xtest, M=M,
		          similarity_function=1,
		          consim=1,
		          calibrate=0,
		          simParms=simParms,
		          modelPriors = modelPriors,
		          draws=niter, burn=nburn, thin=nthin,
		          mh=mh)



oldpar <- par(no.readonly = TRUE)

par(mfrow=c(1,2))
plot(Xtrain[,1], Ytrain, ylab="weight", xlab="length", pch=20)
points(Xtrain[,1], apply(out2$fitted,2,mean), col='red',pch=2, cex=1)
points(Xtrain[,1], apply(out1$fitted,2,mean), col='blue',pch=3, cex=1)
legend(x="topleft",legend=c("Observed","PPM","PPMx"),
          col=c("black","red","blue", "green"),pch=c(20,2,3,4))


plot(Xtest[,1], Ytest, ylab="weight", xlab="length",pch=20)
points(Xtest[,1], apply(out2$ppred,2,mean), col='red',pch=2, cex=1)
points(Xtest[,1], apply(out1$ppred,2,mean), col='blue',pch=3, cex=1)
legend(x="topleft",legend=c("Observed","PPM","PPMx"),
          col=c("black","red","blue","green"),pch=c(20,2,3,4))



par(oldpar)


}



}
