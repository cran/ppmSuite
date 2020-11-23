\name{ordinal_ppmx}

\alias{ordinal_ppmx}

\title{Function that fits Ordinal probit model with a PPMx as a prior on partitions}

\description{
\code{ordinal_ppmx} is the main function used to fit ordinal probit model with a PPMx as a prior on partitions.
}

\usage{
ordinal_ppmx(y, co, X=NULL,Xpred=NULL,
              cohesion=1,
              M=1,
              similarity_function=1,
              consim=1,
              calibrate=0,
              simParms=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1),
              modelPriors=c(0, 10, 1, 1),
              mh=c(0.5, 0.5),
              draws=1100,burn=100,thin=1,
              verbose=FALSE)}


\arguments{
\item{y}{Response vector containing ordinal categories that have been mapped to natural numbers beginning with 0}

\item{co}{Vector specifying the boundaries associated with auxiliary variables of the probit model.  If the number of ordinal categories is c, then the dimension of this vector must be c+1.}

\item{X}{a data frame whose columns consist of covariates that will be incorporated in the partition model.
         Those with class of "character" or "factor" will be treated as categorical covaraites.  All others
         will be treated as continuous covariates.  If NULL, then a PPM is fit.  All continuous covariates are
         standardized to have mean 0 and variance 1 before employing the PPMx model.}

\item{Xpred}{a data frame containing covariate values for which out of sample predictions are desired.
             The format of Xpred must be the same as for X.}

\item{cohesion}{Type of cohesion function to use in the PPMx prior.

    1 - Dirichlet process style of cohesion c(S) = M x (|S| - 1)!

    2 - Uniform cohesion c(S) = 1}

\item{M}{Precision parameter of the PPMx if a DP style cohesion is used.  See above.  Default is 1.}

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
An ordinal probit model is fit.  If covariates containin missing values, then approach employed by page et al is automatically employed.
The computational implementation of the model is based algorithm 8 found in Neal's 2000 JCGS paper
}

\value{
The function returns A list containing arrays filled with MCMC iterates corresponding to model parameters and model fit metrics. In order to provide more detail, in what follows let

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
\examples{

\donttest{



# Continuous Covariate
X1 <- runif(100, 0,1)

# Binary Covariate
X2 <- rbinom(100, 1, 0.5)

pi <- exp(2*X1 + -2*X2)/(exp(2*X1 + -2*X2) + 1)

# Binary response
Y <- rbinom(100, 1, pi)


keep <- sample(1:100, 75, replace=FALSE)

X <- cbind(X1, X2)

Xtn <- X[keep,]
ytn <- Y[keep]
Xtt <- X[-keep,]
ytt <- Y[-keep]


# Since Multinomial need to select boundaries of "latent states".
# For co below the latent states are < 0 and > 0
co <- c(-100000, 0, 100000)
# Their selection is arbitrary and doesn't impact things
# See Thanasis work on this.

#             m0   s20  v    k0   n0   a0
simParms <- c(0.0, 1.0, 2.0, 1.0, 2.0, 0.1)
#                m  s2  s  s  s0 s0
modelPriors <- c(0, 10, 0, 1, 0, 1)


draws <- 50000
burn <- 25000
thin <- 25
nout <- (draws - burn)/thin


# Takes about 15 seconds to run
fit <- ordinal_ppmx(y = ytn, co=co, X=Xtn, Xpred=Xtt,
                     similarity_function=1, consim=1,
                     calibrate=0,
                     simParms=c(0.0, 1.0, 0.5, 1.0, 2.0, 0.1, 1),
                     modelPriors=c(0, 1, 0, 0.5, 0, 10),
                     draws=draws, burn=burn, thin=thin, verbose=FALSE)

# The first partition iterate is used for plotting
# purposes only. We recommended using the salso
# R-package to estimate the partition based on Si
pairs(cbind(Y, X), col=fit$Si[1,])

}

}
