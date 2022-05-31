\name{ccp_ppm}

\alias{ccp_ppm}

\title{Function that fits a multivariate correlated product partition change
point model}

\description{
\code{ccp_ppm} is a function that fits a Bayesian product partition change
point model, where the set of change point indicators between time series are correlated.
}

\usage{
ccp_ppm(ydata,
         nu0, mu0, sigma0,
         mltypes, thetas,
         devs,
         nburn, nskip, nsave,
         verbose = FALSE)
}

\arguments{
\item{ydata}{An \eqn{L \times n} data matrix, where \eqn{L} is the number of
time series and \eqn{n}, the number of time points.}

\item{nu0}{Degrees of freedom of the multivariate Student's
t-distribution (see section Details).}

\item{mu0}{Location vector of dimension \eqn{L} (see section Details).}

\item{sigma0}{Positive definite scale matrix of order \eqn{L \times L} (see
section Details).}

\item{mltypes}{Type of marginal likelihood. Currently only available is:
\itemize{

    \item{\code{mltypes = 1}.} Observations within a block are conditionally
    independent \eqn{Normal(\mu, \sigma^2)} variates with mean \eqn{\mu} and
    variance \eqn{\sigma^2}. The desired marginal likelihood is obtained after
    integrating \eqn{(\mu, \sigma^2)} with respect to a
    \eqn{Normal-Inverse-Gamma(\mu_{0}, \kappa_{0}, \alpha_{0}, \beta_{0})}
    prior.
    }
}

\item{thetas}{An \eqn{L \times q} matrix containing hyperparameters associated
with the marginal likelihood. The number of rows \eqn{(L)} corresponds to the
number of series. The number of columns \eqn{(q)} depend on the marginal
likelihood:
  \itemize{

    \item  If \code{mltypes = 1}, then \eqn{q = 4} and \code{thetas} equals
    the hyperparameter \eqn{(\mu_{0}, \kappa_{0}, \alpha_{0}, \beta_{0})} of
    the Normal-Inverse-Gamma prior.
  }
}

\item{devs}{An \eqn{L \times (n - 1)} matrix containing the standard deviations
of the candidate density associated with the random walk Metropolis-Hastings
steps for updating change point probabilities.}

\item{nburn}{The number of initial MCMC iterates to be discarded as burn-in.}

\item{nskip}{The amount to thinning that should be applied to the MCMC chain.}

\item{nsave}{Then number of MCMC iterates to be stored.}

\item{verbose}{Logical indicating whether to print to screen the MCMC progression. The default value is \code{verbose = FALSE}.}

}

\details{
As described in Quinlan et al. (add cite), for each time series
\eqn{\boldsymbol{y}_{i} = (y_{i,1}, \ldots , y_{i,n})'}:

\deqn{\boldsymbol{y}_{i} \mid \rho_{i} \sim
\prod_{j = 1}^{b_{i}}\mathcal{F}(\boldsymbol{y}_{i,j} \mid
\boldsymbol{\theta}_{i})}

\deqn{\rho_{i} \mid (p_{i,1}, \ldots , p_{i,n-1})' \sim
\prod_{t \in T_{i}}p_{i,t} \prod_{t \notin T_{i}}(1 - p_{i,t}) :
T_{i} = \{\tau_{i,1}, \ldots, \tau_{i,b_{i} - 1}\}}

\deqn{(p_{1,t}, \ldots , p_{L,t})' \sim logit-t(\nu_{0}, \boldsymbol{\mu}_{0},
\boldsymbol{\Sigma}_{0}).}

Here, \eqn{\rho_{i} = \{S_{i,1}, \ldots , S_{i,b_{i}}\}} is a partition of
the set \eqn{\{1, \ldots , n\}} into \eqn{b_{i}} contiguous blocks, and
\eqn{\boldsymbol{y}_{i,j} = (y_{i,t} : t \in S_{i,j})'}. Also,
\eqn{\tau_{i,j} = \max(S_{i,j})} and
\eqn{\mathcal{F}( \cdot \mid \boldsymbol{\theta}_{i})} is a marginal
likelihood function which depends on the nature of \eqn{\boldsymbol{y}_{i}},
indexed by a hyperparameter \eqn{\boldsymbol{\theta}_{i}}. In addition,
\eqn{logit-t(\nu_{0}, \boldsymbol{\mu}_{0}, \boldsymbol{\Sigma}_{0})}
is the logit of a multivariate Student's t-distribution with degrees of
freedom \eqn{\nu_{0}}, location vector \eqn{\boldsymbol{\mu}_{0}} and scale
matrix \eqn{\boldsymbol{\Sigma}_{0}}.}

\value{
The function returns a list containing arrays filled with MCMC iterates
corresponding to model parameters. In order to provide more detail, in what
follows let \eqn{M} be the number of MCMC iterates collected. The output list
contains the following:
\itemize{

  \item{C.} An \eqn{M \times \{L(n - 1)\}} matrix containing MCMC iterates
  associated with each series indicators of a change point. The \eqn{m}th
  row in \code{C} is divided into \eqn{L} blocks; the first \eqn{(n - 1)}
  change point indicators for time series 1, the next \eqn{(n - 1)} change
  point indicators for time series 2, and so on.

  \item{P.} An \eqn{M \times \{L(n - 1)\}} matrix containing MCMC iterates
  associated with each series probability of a change point. The \eqn{m}th
  row in \code{P} is divided into \eqn{L} blocks; the first \eqn{(n - 1)}
  change point probabilities for time series 1, the next \eqn{(n - 1)} change
  point probabilities for time series 2, and so on.
  }
}

\examples{
\donttest{
# Generate data that has two series, each with 100 observations
y1 <- replicate(25, rnorm(4, c(-1, 0, 1, 2), c(0.1, 0.25, 0.5, 0.75)))
y2 <- replicate(25, rnorm(4, c(2, 1, 0, -2), c(0.1, 0.25, 0.5, 0.75)))
y <- rbind(c(t(y1)), c(t(y2)))

# Marginal likelihood parameters
thetas <- matrix(1, nrow = 2,ncol = 4)
thetas[1,] <- c(0, 1, 2, 1)
thetas[2,] <- c(0, 1, 2, 1)

# M-H candidate density standard deviations
devs = matrix(0.1, nrow = 2, ncol = (dim(y)[2] - 1))

# Prior parameters for logit-t distribution
L <- nrow(y)
pivar <- 10
picorr <- 0.9
pimu <- rep(-6, L) # mean associated with logit of p_i
piSigma <- pivar*picorr*(rep(1, L) \%*\% t(rep(1, L))) +
           pivar*(1 - picorr)*diag(L)
nu0 = 3
mu0 = pimu
sigma0 = piSigma

# Fit the bayesian ppm change point model
fit <- ccp_ppm(nburn = 1000, nskip = 1, nsave = 1000, ydata = y, nu0 = nu0,
               mu0 = mu0, sigma0 = sigma0, mltypes = c(1, 1), thetas = thetas,
               devs = devs)
}
}
