
## wrapper for sppm model with no covariates
sppm <- function(y,s,s.pred=NULL,cohesion,M=1,
					modelPriors=c(0, 100^2, 10, 10), cParms=c(1, 1.5, 0, 1, 2, 2),
					mh=c(0.5, 0.5),draws=1100,burn=100,thin=1){
	# Cohesion  1 - distance from centroids
	#			      2 - upper bound
	#			      3 - auxiliary similarity
	#			      4 - double dipper similarity

	# modelPriors = c(m0, s20, ms, ms0)

	# cParm = c(epsilon value - cohesion 1
	#        	distance bound - cohesion 2,
	#			    mu0 - center of NNIG for cohesion 3 and 4
	#			    k0 - scale parm of gaussian part of NNIG for cohesion 3 and 4
	#			    v0 - degrees of freedom IG part of NNIG for cohesion 3 and 4
	#			    L0 - scale parm (scalar of identity matrix) IG part of NNIG for cohesion 3 and 4


  out <- NULL

	nobs <- length(y)
	cat("nobs = ", nobs, "\n")

  s1 <- s[,1]; s2 <- s[,2]
  s1p <- s.pred[,1]; s2p <- s.pred[,2]

	nout <- (draws - burn)/thin
	npred <- ifelse(is.null(s1p), 0, length(s1p))
	cat("npred = ", npred, "\n")

  mu <- sig2 <- Si <- like <- ispred <- matrix(1,nrow=nout,ncol=nobs)
  mu0 <- sig20 <- nclus <- rep(1,nout)
  ppred <- matrix(1, nrow=nout, ncol=npred)
  WAIC <- lpml <- rep(1,1)


  PPM <- FALSE

	C.out <- .C("mcmc_sppm",
              	as.integer(draws), as.integer(burn), as.integer(thin),
              	as.integer(nobs), as.integer(npred), as.integer(cohesion),
              	as.integer(PPM),as.double(M), as.double(y),
              	as.double(s1), as.double(s2),as.double(s1p), as.double(s2p),
              	as.double(modelPriors), as.double(cParms), as.double(mh),
              	mu.out=as.double(mu), sig2.out=as.double(sig2),
              	mu0.out=as.double(mu0), sig20.out=as.double(sig20),
              	Si.out=as.integer(Si), nclus.out=as.integer(nclus), like.out=as.double(like),
              	WAIC.out=as.double(WAIC), lpml.out=as.double(lpml), ispred.out=as.double(ispred),
              	ppred.out=as.double(ppred))

    out$mu <- matrix(C.out$mu.out, nrow=nout, byrow=TRUE)
  	out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
  	out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
  	out$like <- matrix(C.out$like.out, nrow=nout, byrow=TRUE)
  	out$fitted <- matrix(C.out$ispred.out, nrow=nout, byrow=TRUE)
  	out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
  	out$mu0 <- C.out$mu0.out
  	out$sig20 <- C.out$sig20.out
  	out$nclus <- C.out$nclus.out
  	out$WAIC <- C.out$WAIC.out
  	out$lpml <- C.out$lpml.out
  	out


}


