\name{gaussian_ppmx}
\alias{gaussian_ppmx}
\title{Function that fits Gaussian PPMx model}
\usage{
gaussian_ppmx(y, X=NULL, Xpred=NULL,
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
\item{X}{a data frame whose columns consist of covariates that will be incorporated in the partition model. Those with class of "character" or "factor" will be treated as categorical covaraites.  All others will be treated as continuous covariates.  If NULL, then a PPM model  is fit.}

\item{Xpred}{a data frame containing covariate values for which out of sample predictions are desired.The format of Xpred must be the same as for X.}

\item{cohesion}{Type of cohesion function to use in the PPMx prior.

    1 - Dirichlet process style of cohesion c(S) = M x (|S| - 1)!

    2 - Uniform cohesion c(S) = 1}

\item{M}{Precision parameter.  Default is 1.}

\item{similarity_function}{Type of similarity function that is employed for the PPMx prior on partitions.  Options are

    1 - Auxilliary similarity

    2 - Double dipper similarity

    3 - Cluster variance or entropy for categorical covariates

    4 - Mean Gower disimilarity (this one not available if missing values are present in X)}

\item{consim}{If similarity_function is either set to 1 or 2, then this specifies the type of marginal likelihood used as the similarity function. Options are

    1 - N-N(m0, s20, v) (v variance of ''likelihood'', m0 and s20 ''prior'' parameters),

    2 - N-NIG(m0,s20, k0, nu0) (m0 and s20 center and scale of Gaussian, k0 and nu0 )}

\item{calibrate}{Whether the similarity should be calibrated.  Options are

    0 - no calibration

    1 - standardize similarity value for each covariate

    2 - coarsening is applied so that each similarity is raised to the 1/p power}

\item{simParms}{Vector of parameter values employed in the similarity function of the PPMx. Entries of the vector correspond to

    m0 - center continuous similarity with default 0,

    s20 - spread of 'prior' continuous similarity with default 1,

    v2 - spread of 'likelihood' for conitnuous similarity (smaller values place more weight partitions with clusters that contain homogeneous  covariate values)

    k0 - degrees of freedom upper for v (only used for N-NIG similarity model)

    nu0 - scale for v (only used for N-NIG similarity model)

    a0 - dirichlet weight for categorical similarity with default of 0.1 (smaller more weight placed on this variable)

    alpha - weight associated with cluster-variance and Gower disimilarity
}

\item{modelPriors}{Vector of prior parameter values used in the PPMx prior.

    m - prior mean for mu0 with default equal to 0,

    s2 - prior variance mu0 with default equal to 100^2,

    A - upper bound on sigma2*_j with default equal to 10

    B - upper bound on sig20 with default equal to 10
}

\item{mh}{two dimensional vector containing values for tunning parameter associated with MH update for sigma2 and sigma20}

\item{draws}{number of MCMC iterates to be collected. default is 1100}

\item{burn}{number of MCMC iterates discared as burn-in. default is 100}

\item{thin}{number by which the MCMC chain is thinne. default is 1.  Thin must be selected so that it is a multilple of (draws - thin)}

\item{verbose}{Logical indicating if information regarding data and MCMC iterate should be printed to screen}

}

\details{
This generic function fits a Gaussian PPMx model (Muller, Quintana, and Gosner, 2011):
  \deqn{y_i | \mu^*, \sigma^{2*}, c_i \sim N(\mu_{c_i}^*, \sigma^{2*}_{c_i}), i=1,\ldots,n}{y_i ~ N(mu*_{c_i}, sigma2*_{c_i}), i=1,\ldots,n}
  \deqn{\mu^*_j | \mu_0, \sigma^2_0 \sim N(\mu_0,\sigma^2_0)}{mu*_j | mu0, sig20 ~ N(mu0,sig20)}
  \deqn{\sigma^*_j | A \sim UN(0,  A)}{sigma*_j | A  ~ UN(0, A)}
  \deqn{\rho|M,\xi  \sim PPMx(M, \xi)}{rho | M, xi ~ PPMx(M, xi)}
To complete the model specification, the folloing hyperpriors are assumed,
  \deqn{\mu_0 | m, s^2 \sim N(m,s^2)}{mu0 | m, s2 ~ N(m,s2)}
  \deqn{\sigma_0 | B  \sim UN(0,B)}{sigma0 | B ~ UN(0,B)}

Note that we employ uniform prior distributions on variance components as suggest in Gelman's 2006 Bayesian paper.  The PPMx(M, xi) denotes the following
\eqn{Pr(\rho | x, M, \xi) \propto \prod_{j=1}^k M(|S_j| - 1)! g(x^*_j|\xi)}

The computational implementation of the model is based algorithm 8 found in Neal's 2000 JCGS paper.
}

\value{
The function returns a list containing arrays filled with MCMC iterates corresponding to model parameters and model fit metrics. In order to provide more detail, in what follows let

"T" - be the number of MCMC iterates collected,

"N" - be the number of observations,

"P" - be the number of predictions.


The output list contains the following

mu - a matrix  of dimension (T, N) containing MCMC iterates associated with each subjects mean parameter (mu*_{c_i}).

sig2 - a matrix of dimension (T, N) containing MCMC iterates associated with each sujbects variance parameter (sigma2*_{c_i})

Si - a matrix of dimension (T, N) containing MCMC iterates assocated with each subjects cluster label.

like - a matrix of dimension (T, N) containing likelihood values at each MCMC iterate.

fitted - a matrix of dimension (T, N) containing fitted (or in sample predictions) for each subject at each MCMC iterate

ppred - a matrix of dimension (T, P) containing out of sample preditions for each "new" subject at each MCMC iterate

mu0 - vector of length T containing MCMC iterates for mu0 parameter

sig20 - vector of length T containing MCMC iterates for sig20

nclus - vector of length T containing number of clusters at each MCMC iterate

WAIC - scalar containing the WAIC value

lpml - scalar containing lpml value
}
\description{
\code{ppmx} is the main function used to fit Gaussian PPMx model.
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
out <- gaussian_ppmx(y=Ytrain, X=Xtrain, Xpred=Xtest, M=M,
		          similarity_function=1,
		          consim=1,
		          calibrate=0,
		          simParms=simParms,
		          modelPriors = modelPriors,
		          draws=niter, burn=nburn, thin=nthin,
		          mh=mh)



# plot MCMC iterats
plot(density(out$mu[,1:10]),type='l')
plot(density(out$sig2[,1:10]),type='l')
plot(density(out$nc),type='l')
plot(density(out$mu0), type='l')
plot(density(out$sig20), type='l')


# The first partition iterate is used for plotting
# purposes only. We recommended using the salso
# R-package to estimate the partition based on Si
pairs(bear[trainObs,ck],col=out$Si[1,], pch=out$Si[1,])


# To compare fit and predictions when covariates not included
# in the partition model, refit data with PPM rather than PPMx
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
points(Xtrain[,1], apply(out$fitted,2,mean), col='blue',pch="+", cex=1.5)
points(Xtrain[,1], apply(out2$fitted,2,mean), col='red',pch=2, cex=1)
legend(x="topleft",legend=c("Observed","PPM","PPMx"), col=c("black","red","blue"),pch=c(20,3,2))

plot(Xtest[,1], Ytest, ylab="weight", xlab="length",pch=20)
points(Xtest[,1], apply(out$ppred,2,mean), col='blue',pch="+", cex=1.5)
points(Xtest[,1], apply(out2$ppred,2,mean), col='red',pch=2, cex=1)
legend(x="topleft",legend=c("Observed","PPM","PPMx"), col=c("black","red","blue"),pch=c(20,3,2))

par(oldpar)


}



}
