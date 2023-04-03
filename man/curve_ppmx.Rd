\name{curve_ppmx}
\alias{curve_ppmx}
\title{Gaussian PPMx Model for Functional Realizations}
\description{
\code{curve_ppmx} is the main function used to fit Functional Gaussian PPMx model.
}
\usage{
curve_ppmx(y, z, subject,
            Xcon=NULL,Xcat=NULL,
            Xconp=NULL,Xcatp=NULL,
            PPM, M,
            q=3, rw_order=1, balanced=1,
            nknots,knots,npredobs,
            Aparm, modelPriors,
            similarity_function=1,
            consim, calibrate,
            simParms,
            mh=c(1,1),
            draws=1100,burn=100,thin=1)

}

\arguments{
\item{y}{numeric vector or a matrix with two columns that contains measured functional response in long format}
\item{z}{numeric vector contains time points at which functional response is measured in long format}
\item{subject}{vector of the same length as z that identies the subject to which each measurement in y corresponds.  }

\item{Xcon}{a data-frame with number of rows being equal to the number of subjects and whose columns consist of continuous covariates.  These covariates are included in the PPMx model and therefore influence clusters and thus are only used if the PPM argument is FALSE.  This argument is set to NULL by default}
\item{Xcat}{a data-frame  with nsubject number of rows and whose columns consist of categorical covariates.  These covariates are included in the PPMx model and therefore influence clusters and thus are only used if the PPM argument is FALSE.  This argument is set to NULL by default}

\item{Xconp}{a data-frame with the number of rows corresponding to the number of out-of-sample predictions that are desired and columns consist of continuous covariates that are contained in Xcon.}
\item{Xcatp}{a data-frame with the number of rows corresponding to the number of out-of-sample predictions that are desired and columns consist of categorical covariates that are contained in Xcat.}

\item{PPM}{Logical argument that indicates if the PPM or PPMx partition model should be employed.  If PPM = FALSE, then at least one of Xcon and Xcat must be supplied.}

\item{M}{Scale parameter connected to the dispersion parameter of a Dirichlet process.  Default is 1.}

\item{q}{Degree of B-spline employed to fit curves}

\item{rw_order}{Order of the random walk.  This specifies the type of penalty matrix employed in the penalized B-splines.}

\item{balanced}{scalar with 1 - indicating the design was balanced in the sense that all subjects measurements occur at the same time and 0 - indicating design was not balanced.}

\item{nknots}{scalar indicating the number of evenly spaced knots to be used.  This is overidden if the vector of knots is provided (see next argument)}

\item{knots}{vector containing the knots to be used when formulating the spline basis matrix.}

\item{npredobs}{number of time predictions to be made when predicting curve of a new subject.}

\item{Aparm}{Upper bound parameter for lambda wich regulates the similarity of curves with in a cluster.  Larger values result in clusters with curves that can be more dissimilar.}

\item{modelPriors}{Vector of prior parameter values for priors assigned to parameters of the Gaussian Functional data model.
  \itemize{
    \item A - upper bound on sigma*_j.
    \item s2mu - prior variance for mu the mean vector of theta which are cluster specific spline coefficients,

    \item mb0 - prior mean for beta0_i (subject specific intercept)
    \item s2b0 - prior variance for beta0_i (subject specific intercept)

    \item as2b0 - prior shape associated with IG prior on variance of beta0_i (subject specific intercept)
    \item bs2b0 - prior scale associated with IG prior on variacne of beta0_i (subject specific intercept)

    \item at - prior shape associated with IG prior on tau (smoothing parameter for theta)
    \item bt - prior scale associated with IG prior on tau (smoothing parameter for theta)
}
}


\item{similarity_function}{Type of similarity function that is employed for the PPMx prior on partitions.  Options are 1-4 with

  \itemize{
    \item  1 - Auxilliary similarity

    \item  2 - Double dipper similarity

    \item  3 - Cluster variance or entropy for categorical covariates

    \item 4 - Mean Gower dissimilarity (Gower dissimilarity is not available if missing values are present in X)
  }
}
\item{consim}{If similarity_function is set to either 1 or 2, then consim specifies the type of marginal likelihood used as the similarity function. Options are 1 or 2 with
  \itemize{
     \item 1 - N-N(m0, s20, v) (v variance of ``likelihood'', m0 and s20 ``prior'' parameters),

    \item 2 - N-NIG(m0, k0, nu0, s20) (m0 and k0 center and inverse scalar of a Gaussian, and nu0 and s20 are the number of prior observations and prior variance guess of a Inverse-Chi-Square distribution.)
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

    \item s20 - spread of continuous similarity with default 1 if consim=1.  For consim=2 guess of x's variance,

    \item v2 - spread of 'likelihood' for conitnuous similarity (smaller values place more weight on partitions with clusters that contain homogeneous  covariate values)

   \item  k0 - inverse scale for v (only used for N-NIG similarity model)

    \item nu0 - prior number of x "observations" (only used for N-NIG similarity model)

    \item a0 - dirichlet weight for categorical similarity with default of 0.1 (smaller values place more weight on partitions with individuals that are in the same category.)

   \item  alpha - weight associated with cluster-variance and Gower disimilarity
}
}
\item{mh}{two dimensional vector containing values for tunning parameter associated with MH update for sigma2 and sigma20}

\item{draws}{number of MCMC iterates to be collected. default is 1100}

\item{burn}{number of MCMC iterates discared as burn-in. default is 100}

\item{thin}{number by which the MCMC chain is thinne. default is 1.  Thin must be selected so that it is a multilple of (draws - thin)}

}

\details{
This function is able to fit a Gaussian PPM or PPMx model as detailed in (Mueller, Quintana, and Rosner, 2011) to functional data.
}

\value{
The function returns a list containing arrays filled with MCMC iterates corresponding to model parameters and model fit metrics. In order to provide more detail, in what follows let

"T" - be the number of MCMC iterates collected,

"N" - be the number of subjects/units,

"P" - be the number of knots + degree of spline.


The output list contains the following
\itemize{
\item Si - a matrix of dimension (T, N) containing MCMC iterates assocated with each subjects cluster label.

\item nclus - a matrix  of dimension (T, 1) containing MCMC iterates associated with the number of clusters

\item beta - an array of dimension (N, P, T) containing the MCMC iterates assocated with each subjects P-dimensional B-spline coefficients

\item theta - an array of dimension (N, P, T) containing the MCMC iterates assocated with the cluster specific P-dimensional B-spline coefficients.  Each subjects theta value is reported.

\item sig2 - a matrix of dimension (T, N) containing MCMC iterates associated with each subjects variance parameter (sigma2*_{c_i})

\item tau2 - a matrix of dimension (T, N) containing MCMC iterates associated with each the cluster-specific smoothing parameter for theta

\item mu - a matrix of dimension (T, P) containing MCMC iterates for the the P-dimensional B-spline coefficients associated with the global mean.

\item lam - a matrix of dimension (T, N) containing MCMC iterates for the cluster-specific lambda parameter that dictates the similarity of curves within a cluster

\item beta0 - a matrix of dimension (T, N) containing MCMC iterates for the subject-specific intercepts


\item mub0 - vector of length T containing MCMC iterates for mean of beta0

\item sig2b0 - vector of length T containing MCMC iterates for variance of beta0

\item like - a matrix of dimension (T, N) containing likelihood values at each MCMC iterate.

\item WAIC - scalar containing the WAIC value

\item lpml - scalar containing lpml value

\item Hmat - a spline design matrix of dimension (N, P)
}
}

\examples{

\donttest{

# generate data for two clusters with 10 subjects each

nobs <- 100
nsubject <- 2*10

set.seed(101)
xx <- seq(0,2*pi, length=nobs)
y <- cbind(replicate(n=10, sin(xx) + rnorm(nobs,0,0.5)),
           replicate(n=10, cos(xx) + rnorm(nobs,0,0.5)))

dat <- data.frame(y=c(y),
                  z=rep(1:nobs, times=nsubject),
                  Name=rep(1:nsubject, each=nobs))

subject_obs_vec <- dat$Name

nknots <- 15



# Small number of iterates for illustrative purposes only
niter <- 5000
nburn <- 2000
nthin <- 3
nout <- (niter-nburn)/nthin

z <- dat$z

## the order here is c(mu0, s20, v, k0, nu0, a0, alpha)
## If simularity is N-NIG then k0 and nu0 are used but v is not
## If simularity is N-N then v is used but no k0 and nu0
simparms <- c(0.0, 1.0, 0.1, 1.0, 1.0, 0.1, 1)

fits <- list()

# fit vgrf only
y <- dat$y

modelPriors <- c(0.5,      # Asig
                 1000^2,   # s2_mu
                 0,        # mb0
                 1000^2,   # s2b0
                 1,        # as2b0
                 1,        # bs2b0
                 1,        # at
                 1.0/0.05) # bt


fit <- curve_ppmx(y=cbind(y), z=z,
             subject=subject_obs_vec,
             Xcon = NULL, Xcat = NULL,
             Xconp=NULL, Xcatp=NULL,
             PPM=TRUE, M=1,
             q=3, rw_order=1, balanced=1,
             nknots=nknots,  npredobs=1, knots=1,
             Aparm=100,
             modelPriors=modelPriors,
             similarity_function=1,
             consim=1, calibrate=0,
             simParms=simparms,
             mh=c(0.1, 1e-4),
             draws=niter,
             burn=nburn,
             thin=nthin)


Hmat <- fit$Hmat

# For a point estimate of partition, take first MCMC interate
# This is done only for illustrative purposes.  Recommend using
# the salso R package.

p.est <- fit$Si[1,]

nc <- length(unique(p.est))


oldpar <- par(no.readonly = TRUE)


# Plot individual subject fits.

tmp <- c(1,6,11,16)
par(mfrow=c(2,2))
for(j in tmp){
  bmn <- apply(fit$beta[j,,],1,mean)
  b0mn <- mean(fit$beta0[,j])

  ytmp <- y[dat$Name==j]

  b0vec <- rep(b0mn, nobs)

  plot(1:nobs,c(ytmp),
		  type='n',ylab="Response",
		  xlab="Time")

  points(1:nobs,ytmp)
  lines(1:nobs,  b0vec+Hmat\%*\%bmn, col=p.est[j],lwd=2)

}

# plot all curves in one plot
par(mfrow=c(1,1))

plot(dat$z, dat$y, type="n",ylab="",xlab="Time")

for(j in 1:nsubject){

  bmn <- apply(fit$beta[j,,],1,mean)
  b0mn <- mean(fit$beta0[,j])

  b0vec <- rep(b0mn, nobs)

  lines((1:nobs),  b0vec+Hmat\%*\%bmn, col=p.est[j],lwd=0.5)

}


par(oldpar)


}



}
