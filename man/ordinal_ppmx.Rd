\name{ordinal_ppmx}

\alias{ordinal_ppmx}

\title{Function that fits ordinal probit model with a PPMx as a prior on partitions}

\description{
\code{ordinal_ppmx} is the main function used to fit ordinal probit model with a PPMx as a prior on partitions.
}

\usage{
ordinal_ppmx(y, co, X=NULL,Xpred=NULL,
              meanModel=1,
              cohesion=1,
              M=1,
              PPM = FALSE,
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
         Those with class of "character" or "factor" will be treated as categorical covariates.  All others
         will be treated as continuous covariates. }

\item{Xpred}{a data frame containing covariate values for which out of sample predictions are desired.
             The format of Xpred must be the same as for X.}

\item{meanModel}{Type of mean model included in the likelihood that is to be used
\itemize{
    \item 1 - cluster-specific means with no covariates in likelihood.

    \item 2 - cluster-specific intercepts and a global regression of the type Xbeta is included in the likelihood.}
}

\item{cohesion}{Type of cohesion function to use in the PPMx prior.
\itemize{
    \item 1 - Dirichlet process style of cohesion c(S) = M x (|S| - 1)!

    \item 2 - Uniform cohesion c(S) = 1}
}

\item{M}{Precision parameter of the PPMx if a DP style cohesion is used.  See above.  Default is 1.}

\item{PPM}{Logical argument that indicates if the PPM or PPMx partition model should be employed.  If PPM = FALSE, then an X matrix must be supplied.}


\item{similarity_function}{Type of similarity function that is employed for the PPMx prior on partitions.  Options are
\itemize{
    \item 1 - Auxilliary similarity

    \item 2 - Double dipper similarity

    \item 3 - Cluster variance or entropy for categorical covariates

    \item 4 - Mean Gower disimilarity (this one not available if missing values are present in X)}
}
\item{consim}{If similarity_function is set to either 1 or 2, then consim specifies the type of marginal likelihood used as the similarity function. Options are (see simparms argument for more details)
\itemize{
    \item 1 - N-N(m0, s20, v) (v variance of ''likelihood'', m0 and s20 ''prior'' parameters),

    \item 2 - N-NIG(m0, k0, nu0, s20) (m0 and k0 center and inverse scalar of a Gaussian, and nu0 and s20 are the number of prior observations and prior variance guess of a Inverse-Chi-Square distribution.)}
}
\item{calibrate}{This argument determines if the similarity should be calibrated.  Options are
\itemize{
    \item 0 - no calibration

    \item 1 - standardize similarity value for each covariate

    \item 2 - coarsening is applied so that each similarity is raised to the 1/p power}
}
\item{simParms}{Vector of parameter values employed in the similarity function of the PPMx. Entries of the vector correspond to
\itemize{
    \item m0 - center continuous similarity with default 0,

    \item s20 - spread of continuous similarity with default 1 if consim=1.  For consim=2 guess of x's variance,

    \item v2 - spread of 'likelihood' for conitnuous similarity (smaller values place more weight on partitions with clusters that contain homogeneous  covariate values)

   \item  k0 - inverse scale for v (only used for N-NIG similarity model)

    \item nu0 - prior number of x "observations" (only used for N-NIG similarity model)

    \item a0 - dirichlet weight for categorical similarity with default of 0.1 (smaller values place more weight on partitions with individuals that are in the same category.)

   \item  alpha - weight associated with cluster-variance and Gower disimilarity
}
}
\item{modelPriors}{Vector of prior parameter values for priors assigned to parameters of the Gaussian latent model.
\itemize{
    \item m - prior mean for mu0 with default equal to 0,

    \item s2 - prior variance mu0 with default equal to 100^2,

    \item A - upper bound on sigma2*_j with default equal to 10

    \item A0 - upper bound on sig20 with default equal to 10
}
}
\item{mh}{two dimensional vector containing values for tunning parameter associated with MH update for sigma2 and sigma20}

\item{draws}{number of MCMC iterates to be collected. default is 1100}

\item{burn}{number of MCMC iterates discared as burn-in. default is 100}

\item{thin}{number by which the MCMC chain is thinne. default is 1.  Thin must be selected so that it is a multilple of (draws - thin)}

\item{verbose}{Logical indicating if information regarding data and MCMC iterate should be printed to screen}

}

\details{
This function fits an ordinal probit model with either a PPM or PPMx prior on partitions.  For details on the ordinal probit model see Kottas et al (2005) and Page, Quintana, Rosner (2020).  Cutpoints listed in the ``co'' argument can be arbitrarily selected, but values that are too far from zero will result in an algorithm that will require more burn-in.   Based on these cutpoints latent variables are introduced.  The latent variables are assumed to follow a Gaussian distribution with cluster-specific means and variances. If meanModel = 2, then a ``global'' regression component is added to the mean resulting in a model with cluster-specific parallel regression lines.  Commonly used conjugate priors are then employed in the regression component.

If covariates contain missing values, then the approach developed in Page, Quintana, Mueller (2022) is automatically employed.  Missing values must be denoted using "NA".  Currently, NAs cannot be accommodated if a ``global'' regression is desired (i.e., meanMode = 2).

We recommend standardizing covariates so thay they have mean zero and standard deviation one.  This makes the default values provided for the similarity function reasonable in most cases.  If covariates are standardized and meanModel = 2 the regression coefficients are estimated on the original scale and are ordered such that the continuous covariates appear first and the categorical covariates come after.

The MCMC algorithm used to sample from the joint posterior distribution is based on algorithm 8 found in Neal 2000.
}

\value{
The function returns a list containing arrays filled with MCMC iterates corresponding to model parameters and also returns a couple of model fit metrics. In order to provide more detail, in what follows let

"T" - be the number of MCMC iterates collected,

"N" - be the number of observations,

"P" - be the number of predictions.

"C" - be the total number of covariates

The output list contains the following
\itemize{
\item mu - a matrix  of dimension (T, N) containing MCMC iterates associated with each subjects mean parameter (mu*_{c_i}).

\item sig2 - a matrix of dimension (T, N) containing MCMC iterates associated with each sujbects variance parameter (sigma2*_{c_i})

\item beta - available only if meanModel = 2, then this is a matrix of dimension (T,C) containing MCMC iterates associated coefficients in the global regression

\item Si - a matrix of dimension (T, N) containing MCMC iterates assocated with each subjects cluster label.

\item zi - a matrix of dimension (T, N) containing MCMC iterates assocated with each subjects latent variable.

\item mu0 - vector of length T containing MCMC iterates for mu0 parameter

\item sig20 - vector of length T containing MCMC iterates for sig20

\item nclus - vector of length T containing number of clusters at each MCMC iterate

\item like - a matrix of dimension (T, N) containing likelihood values at each MCMC iterate.

\item WAIC - scalar containing the WAIC value

\item lpml - scalar containing lpml value

\item fitted.values - a matrix of dimension (T, N) containing fitted values at the latent variable level for each subject at each MCMC iterate

\item ppred - a matrix of dimension (T, P) containing out of sample preditions at the latent variable level for each ``new'' subject at each MCMC iterate

\item predclass - a matrix of dimension (T, P) containing MCMC iterates of cluster two which "new" subject is allocated

\item rbpred - a matrix of dimension (T, P) containing out of sample preditions at the latent variable level for each "new" subject at each MCMC iterate based on the rao-blackwellized prediction

\item ord.fitted.values - a matrix of dimension (T, N) containing fitted values on the ordinal variable scale for each subject at each MCMC iterate.

\item ord.ppred - a matrix of dimension (T, P) containing out of sample preditions on the ordinal variable scale for each ``new'' subject at each MCMC iterate from the posterior predictive distribution.

\item ord.rbpred - a matrix of dimension (T, P) containing out of sample preditions on the ordinal variable scale  for each "new" subject at each MCMC iterate based on the rao-blackwellized prediction.

}
}
\examples{
\donttest{
n <- 100
# Continuous Covariate
X1 <- runif(n, -1,1)

# Binary Covariate
X2 <- rbinom(n, 1, 0.5)

pi <- exp(2*X1 + -2*X2)/(exp(2*X1 + -2*X2) + 1)

# Binary response
Y <- rbinom(n, 1, pi)


keep <- 1:(n-25)

# standardize X1 to have mean zero and sd 1.
X <- data.frame(X1=scale(X1), X2=as.factor(X2))

Xtn <- X[keep,]
ytn <- Y[keep]
Xtt <- X[-keep,]
ytt <- Y[-keep]


# Since we have a binary response there are two "latent states".
# The boundaries of the latent states can be selected arbitrarily.
# Below I essentially use (-Inf, 0, Inf) to define the two latent spaces.
co <- c(-1e5, 0, 1e5)


#             m0   s20  v  k0   n0   a0   alpha
simParms <- c(0, 1.0, 0.5, 1.0, 2.0, 0.1, 1)
#                m  s2  s  s0
modelPriors <- c(0, 10, 1, 1)


draws <- 50000
burn <- 25000
thin <- 25
nout <- (draws - burn)/thin


# Takes about 20 seconds to run
fit <- ordinal_ppmx(y = ytn, co=co, X=Xtn, Xpred=Xtt,
                     meanModel=1,
                     similarity_function=1, consim=1,
                     calibrate=0,
                     simParms=simParms,
                     modelPriors=modelPriors,
                     draws=draws, burn=burn, thin=thin, verbose=FALSE)

# The first partition iterate is used for plotting
# purposes only. We recommended using the salso
# R-package to estimate the partition based on Si
pairs(cbind(Y, X), col=fit$Si[1,])

# in-sample confusion matrix
table(ytn, apply(fit$ord.fitted.values, 2, median))

# out-of-sample confusion matrix based on posterior predictive samples
table(ytt, apply(fit$ord.ppred, 2, median))
}

}
