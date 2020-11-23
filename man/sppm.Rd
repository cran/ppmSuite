\name{sppm}
\alias{sppm}
\title{Function that fits spatial product partition model with Gaussian likelihood}
\usage{
sppm(y,s,
    s.pred=NULL,
    cohesion,
    M=1,
    modelPriors=c(0, 100^2, 10, 10),
    cParms=c(1, 1.5, 0, 1, 2, 2),
    mh=c(0.5, 0.5),
    draws=1100,burn=100,thin=1)
}

\description{
\code{sppm} is the main function used to fit model with Guassian likelihood and spatial PPM as prior on partitions.

}

\arguments{
\item{y}{numeric vector containing response variable}

\item{s}{Two-column matrix containing spatial locations (i.e., longitude and lattitude).}

\item{s.pred}{Two-column matrix containing spatial locations at which out-of-sample predictions will be collected.}

\item{cohesion}{Scalar that indicates which cohesion to use.
\itemize{
    \item 1 - distance from centroids \cr
	  \item 2 - upper bound \cr
	  \item 3 - auxiliary similarity \cr
	  \item 4 - double dipper similarity \cr
  }
}

\item{M}{Parameter related to Dirichlet process scale or dispersion parameter.}

\item{modelPriors}{Vector containing model prior values (see below for more details)}

\item{cParms}{Vector containing partition model prior values (see below for more details)}

\item{mh}{Tuning standard deviations for metropolis updates for sigma2 and sigma20}

\item{draws}{Number of MCMC samples to collect}

\item{burn}{Number of the MCMC samples discarded in the burn-in phase of the sampler}

\item{thin}{The amount of thinning desired for the chain}

}
\details{

The vector \verb{modelPriors = c(m0, s20, ms, ms0)} where each prior parameter is listed in the model description below.

The cParm vector contains values associated with the cohesion function.\cr
\verb{
cParm = c(
  epsilon value - cohesion 1 only,
  distance bound - cohesion 2 only,
  mu0 - center of NNIG for cohesion 3 and 4
  k0 - scale parm of gaussian part of NNIG for cohesion 3 and 4
  v0 - degrees of freedom IG part of NNIG for cohesion 3 and 4
  L0 - scale parm (scalar of identity matrix) IG part of NNIG for cohesion 3 and 4).}

The model this function fits is Gaussian likelihood model using the sPPM prior on partitions (Page and Quintana, 2016).  Specific model details are
  \deqn{y_i | \mu^*, \sigma^{2*}, c_i \sim N(\mu_{c_i}^*, \sigma^{2*}_{c_i}), i=1,\ldots,n}{y_i ~ N(mu*_{c_i}, sigma2*_{c_i}), i=1,\ldots,n}
  \deqn{\mu^*_j | \mu_0, \sigma^2_0 \sim N(\mu_0,\sigma^2_0)}{mu*_j | mu0, sig20 ~ N(mu0,sig20)}
  \deqn{\sigma^*_j | A \sim UN(0,  ms)}{sigma*_j | A  ~ UN(0, ms)}
  \deqn{\rho|M,\xi  \sim sPPM}{rho ~ sPPM}

To complete the model specification, the folloing hyperpriors are assumed,
  \deqn{\mu_0 | m, s^2 \sim N(m0,s0^2)}{mu0 | m0, s02 ~ N(m0,s02)}
  \deqn{\sigma_0 | B  \sim UN(0,ms0)}{sigma0 | B ~ UN(0,ms0)}

Note that we employ uniform prior distributions on variance components as suggest in Gelman's 2006 Bayesian paper.  "sPPM" in the model specificaiton denotes the the spatial product partition model.  The computational implementation of the model is based algorithm 8 found in Neal's 2000 JCGS paper.
}

\value{
This function returns in a list all MCMC interates for each model parameter, posterior predictive, and fitted values.  In addition the LPML model fit metric is provided.
}
\examples{


\donttest{


data(scallops)

Y<-log(scallops[,5]+1)
s_coords <- scallops[,3:4] #lat and long
m <- dim(s_coords)[1]



# standardize spatial coordinates
smn <- apply(s_coords,2,mean)
ssd <- apply(s_coords,2,sd)
s_std <- t((t(s_coords) - smn)/ssd)

# Create a grid of prediction locations
np <- 10

sp <- expand.grid(seq(min(s_coords[,1]), max(s_coords[,1]),length=np),
                   seq(min(s_coords[,2]), max(s_coords[,2]), length=np))

sp_std <- t((t(sp) - smn)/ssd) # standardized prediction spatial coordinates


niter <- 20000
nburn <- 10000
nthin <- 10
nout <- (niter - nburn)/nthin


out <- sppm(y=Y,s=s_std,s.pred=sp_std,cohesion=4, M=1, draws=niter, burn=nburn, thin=nthin)

# fitted values
fitted.values <- out$fitted
fv.mn <- apply(fitted.values, 2,mean)
mean((Y - fv.mn)^2) # MSE
out$lpml #lpml value

ppred <- out$ppred
predmn <- apply(ppred,2,mean)

# The first partition iterate is used for plotting
# purposes only. We recommended using the salso
# R-package to estimate the partition based on Si
Si <- out$Si
plot(s_coords[,1], s_coords[,2], col=Si[1,])

}



}
