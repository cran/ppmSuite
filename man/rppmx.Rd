\name{rppmx}
\alias{rppmx}
\title{Function generates random realizations from a PPM or PPMx}
\description{
\code{rppmx} Employes the ploya urn sampling scheme to randomly generate a partition from the PPM or PPMx.
}

\usage{
rppmx(m, X=NULL,
      similarity,
      simparm,
      M=1,
      m0=0,s20=1,v=2,k0=10,v0=1,alpha=1)
}



\arguments{
\item{m}{Number of unites that are allocated to partitions}
\item{X}{a data frame whose columns consist of covariates that will be incorporated in the partition model.
         Those with class of "character" or "factor" will be treated as categorical covaraites.  All others
         will be treated as continuous covariates.  If NULL, then a PPM partition is produced.}

\item{similarity}{Type of similarity function that is employed for covariates.  Options are

    1 - Auxilliary similarity,

    2 - Double dipper similarity

    3 - variance similarity}

\item{simparm}{Type of similarty model employed for continuous covariates. Options are

    1 - N-N(m0, s20, v) (v variance of ''likelihood'', m0 and s20 ''prior'' parameters),

    2 - N-NIG(m0,k0, k0, v0, s20) (m0 and k0 center and scale of Gaussian, n0 and s20 shape and scale of IG )}


\item{M}{Precision parameter.  Default is 1.}

\item{m0}{Continuous similarity function value (see above)}
\item{s20}{Continuous similarity function value (see above)}
\item{v}{Continuous similarity function value (see above)}
\item{k0}{Continuous similarity function value (see above)}
\item{v0}{Continuous similarity function value (see above)}
\item{alpha}{Penalty value when using the variance similarity}
}

\details{
Use polya urn scheme to sample from the PPM or the PPMx
}

\value{
The function returns randomly generated partition
}
\examples{

X <- cbind(rnorm(100), rbinom(100,1,0.5))
p <- rppmx(m=100, X=X, similarity=1, simparm=1, M=1)
p


}




