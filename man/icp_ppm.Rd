\name{icp_ppm}

\alias{icp_ppm}

\title{Function that fits the multivariate independent product partition
change point model}

\description{
\code{icp_ppm} is a function that fits a Bayesian product partition change
point model. Each series is treated independently.
}

\usage{
icp_ppm(ydata,
         a0, b0,
         mltypes,
         thetas,
         nburn, nskip, nsave,
         verbose = FALSE)
}
\arguments{
\item{ydata}{An \eqn{ascii}{L \times n} data matrix, where \eqn{L} is the number of
time series and \eqn{n}, the number of time points.}

\item{a0}{Vector of dimension \eqn{L} with shape 1 Beta parameters (see section
Details).}

\item{b0}{Vector of dimension \eqn{L} with shape 2 Beta parameters (see section
Details).}

\item{mltypes}{Type of marginal likelihood. Currently only available is:
\itemize{

    \item{\code{mltypes = 1}.} Observations within a block are conditionally
    independent \eqn{Normal(\mu, \sigma^2)} variates with mean \eqn{mu} and
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

\item{nburn}{The number of initial MCMC iterates to be discarded as burn-in.}

\item{nskip}{The amount to thinning that should be applied to the MCMC chain.}

\item{nsave}{Then number of MCMC iterates to be stored.}

\item{verbose}{Logical indicating whether to print to screen the MCMC
progression. The default value is \code{verbose = FALSE}.}

}

\details{
As described in Loschi (add cite), for each time series
\eqn{\boldsymbol{y}_{i} = (y_{i,1}, \ldots , y_{i,n})'}:

\deqn{\boldsymbol{y}_{i} \mid \rho_{i} \sim
\prod_{j = 1}^{b_{i}}\mathcal{F}(\boldsymbol{y}_{i,j} \mid
\boldsymbol{\theta}_{i})}

\deqn{\rho_{i} \mid p_{i} \sim p_{i}^{b_{i} - 1}(1 - p_{i})^{n - b_{i}}}

\deqn{p_{i} \sim Beta(a_{i,0}, b_{i,0}).}

Here, \eqn{\rho_{i} = \{S_{i,1}, \ldots , S_{i,b_{i}}\}} is a partition of
the set \eqn{\{1, \ldots , n\}} into \eqn{b_{i}} contiguous blocks, and
\eqn{\boldsymbol{y}_{i,j} = (y_{i,t} : t \in S_{i,j})'}. Also,
\eqn{\mathcal{F}( \cdot \mid \boldsymbol{\theta}_{i})} is a marginal
likelihood function which depends on the nature of \eqn{\boldsymbol{y}_{i}},
indexed by a hyperparameter \eqn{\boldsymbol{\theta}_{i}}. Notice that
\eqn{p_{i}} is the probability of observing a change point in series \eqn{i},
at each time \eqn{t \in \{2, \ldots , n\}}.
}

\value{
The function returns a list containing arrays filled with MCMC iterates
corresponding to model parameters. In order to provide more detail, in what
follows let \eqn{M} be the number of MCMC iterates collected. The output list
contains the following:
\itemize{

  \item{C.} An \eqn{M \times \{L * (n - 1)\}} matrix containing MCMC iterates
  associated with each series indicators of a change point. The \eqn{m}th
  row in \code{C} is divided into \eqn{L} blocks; the first \eqn{(n - 1)}
  change point indicators for time series 1, the next \eqn{(n - 1)} change
  point indicators for time series 2, and so on.

  \item{P.} An \eqn{M \times \{L * (n - 1)\}} matrix containing MCMC iterates
  associated with each series probability of a change point. The \eqn{m}th
  row in \code{P} is divided into \eqn{L} blocks; the first \eqn{(n - 1)}
  change point probabilities for time series 1, the next \eqn{ascii}{(n - 1)} change
  point probabilities for time series 2, and so on.
  }
}

\examples{
\donttest{
# Generate data that has two series, each with 100 observations
y1 <- replicate(25, rnorm(4, c(-1, 0, 1, 2), c(0.1, 0.25, 0.5, 0.75)))
y2 <- replicate(25, rnorm(4, c(2, 1, 0, -2), c(0.1, 0.25, 0.5, 0.75)))
y <- rbind(c(t(y1)), c(t(y2)))
n <- ncol(y)
# Marginal likelihood parameters
thetas <- matrix(1, nrow = 2, ncol = 4)
thetas[1,] <- c(0, 1, 2, 1)
thetas[2,] <- c(0, 1, 2, 1)

# Fit the Bayesian ppm change point model
fit <- icp_ppm(ydata = y,
               a0 = c(1, 1),
               b0 = c(1, 1),
               mltypes = c(1, 1),
               thetas = thetas,
               nburn = 1000, nskip = 1, nsave = 1000)

cpprobsL <- matrix(apply(fit$C,2,mean), nrow=n-1, byrow=FALSE)


}
}
