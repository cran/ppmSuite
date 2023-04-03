# Function that computes ``robust'' lpmls and waic
lpml.robust <- function(llike){

	omega <- 1/exp(llike)

	lpml <- sum(log(1/(apply(omega,2,mean))))

	omegabar <- apply(omega,2,mean)

	omegatil <- matrix(NA, nrow=nrow(omega), ncol=ncol(omega))
	for(i in 1:nrow(omega)){
		for(j in 1:ncol(omega)){

			omegatil[i,j] <- min(omega[i,j], sqrt(nrow(omega)*omegabar[j]))
		}
	}

	lpml.robust <- sum(log(1/(apply(omegatil,2,mean))))

	# Compute lpml by droping iteration that produces inf
	indx <- which(llike < -709, arr.ind=TRUE)[,1]
	CPO <- apply(llike[-indx,], 2, function(x) mean(1/exp(x)))
	if(length(indx) == 0){
		CPO <- apply(llike, 2, function(x) mean(1/exp(x)))
	}
	lpml.drop.inf <- sum(log(1/CPO))

	# Compute lpml by removing all llike values less than
  # the 5th quantile
  mn995 <- function(x) mean(x[x < quantile(x, 0.995)])

	lpml.quant <- sum(log(1/(apply(omega,2,mn995))))


	# Comput lpml by summing off min (see paper for details)
	nobs <- ncol(llike)
	niter <- nrow(llike)
	minll <- apply(llike,2,min)
	tmp <- numeric()
	for(i in 1:nobs){

		tmp[i] <- sum(exp(minll[i] - llike[,i]))

	}

	lpml.log.scale <- sum(log(niter) + minll - log(tmp))

	# compute the waic
  mnllike <- apply(llike, 2, mean)
  mnlike <- apply(exp(llike), 2, mean)
  elppWAIC <- sum(2*mnllike - log(mnlike))
  waic <- -2*elppWAIC

	out <- c(lpml, lpml.robust, lpml.log.scale, lpml.drop.inf, lpml.quant, waic)
	names(out) <- c("lpml",
                    "lpml.stable.gelman",
                    "lpml.stable.log",
                    "lpml.drop.inf",
                    "lpml.drop.quant995",
                    "waic")
	out
}
