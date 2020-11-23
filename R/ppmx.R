# Wrapper the executes the C code with the .C function
# I am standardizing the x's.

gaussian_ppmx <- function(y,X=NULL,Xpred=NULL,
                           cohesion=1, M=1,
                           similarity_function=1, consim=1,
                           calibrate=0,
                           simParms=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1.0),
                           modelPriors = c(0, 100^2, 1, 1),
                           mh=c(0.5, 0.5),
                           draws=1100,burn=100,thin=1,
                           verbose=FALSE){

  # X - data.frame whose columns are
  # gcontype - similarity function (1 - Auxilliary, 2 - double dipper)
  # consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
  # consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)
  # modelPriors = (mu0, s^2, a1, m)
  # simParms = (m0, s20, v2, k0, nu0, dir, alpha)
  # mh - tuning parameters for MCMC updates of sig2 and sig20
  out <- NULL

  PPM <- ifelse(is.null(X), TRUE, FALSE)

  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)
  if(!is.data.frame(Xpred) & !is.null(Xpred)){
    Xpred <- data.frame(Xpred);
    colnames(Xpred) <- colnames(X)
  }

  nout <- (draws-burn)/thin

  nobs <- length(y)
  nxobs <- ifelse(is.null(X), 0, nrow(X))
  npred <- ifelse(is.null(Xpred), 0, nrow(Xpred))

  Xall <- rbind(X, Xpred)
  print(Xall)
  nmissing <- sum(is.na(Xall))

  if(nmissing > 0){
    Mall <- 1*is.na(Xall) # create the missing matrix with 1 = missing
  }

  # Function that relabels categorical variables to begin with 0
  relab <- function(x) as.numeric(as.factor(as.character(x))) - 1

  classes <- sapply(Xall, class)
  catvars <- classes %in% c("factor","character")

  # standardize continuous covariates
  if(nxobs > 0){
    if(sum(!catvars) > 0){
  	  Xconstd <- apply(Xall[,!catvars, drop=FALSE], 2, scale)
  	  Xcon <- Xconstd[1:nobs,,drop=FALSE];
  	  ncon <- ncol(Xcon)
  	  if(nmissing > 0) Mcon <- Mall[1:nobs, !catvars, drop=FALSE]
    }else{
      Xcon <- cbind(rep(0,nobs));
      ncon <- 0
      if(nmissing > 0) Mcon <- cbind(rep(0,nobs))
    }


    if(sum(catvars) > 0){
  	  # Change the factors or characters into integers with category starting at 0.
  	  Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
  	  if(nmissing > 0) Mcat <- Mall[1:nobs, catvars, drop=FALSE]

  	  Xcat <- Xcatall[1:nobs,,drop=FALSE];
  	  Cvec <- apply(Xcat,2,function(x)length(unique(x)))
  	  ncat <- ncol(Xcat)

    }else{
  	  Xcat <- cbind(rep(0,nobs));
      if(nmissing > 0) Mcat <- cbind(rep(0,nobs))
  	  Cvec <- 0
  	  ncat <- 0
    }
  } else {
	  Xcon <- cbind(rep(0,1));
	  Xcat <- cbind(rep(0,1));
  	if(nmissing > 0) Mcon <- cbind(rep(0,1));
  	if(nmissing > 0) Mcat <- cbind(rep(0,1));
  }


  # Create the matrices of continuous and categorical variables
  # that will be used in th posterior predictive distributions
  if(npred > 0){
  	if(sum(!catvars) > 0){
  	  Xconstd <- apply(Xall[,!catvars, drop=FALSE], 2, scale)
  	  Xconp <- Xconstd[(nrow(Xall)-npred+1):nrow(Xall),,drop=FALSE];
  	  if(nmissing > 0) Mconp <- Mall[(nrow(Xall)-npred+1):nrow(Xall), !catvars, drop=FALSE]
  	  ncon <- ncol(Xconp)
  	} else {
  	  Xconp <- cbind(rep(0,npred));
  	  if(nmissing > 0) Mconp <- cbind(rep(0,npred));
      ncon <- 0
	  }

	  if(sum(catvars) > 0){
  	  Xcatall <- apply(Xall[, catvars,drop=FALSE], 2, relab)
	    Xcatp <- Xcatall[(nrow(Xall)-npred+1):nrow(Xall),,drop=FALSE];
  	  if(nmissing > 0) Mcatp <- Mall[(nrow(Xall)-npred+1):nrow(Xall), catvars, drop=FALSE]
	    ncat <- ncol(Xcatp)
  	  Cvec <- apply(Xcatall,2,function(x)length(unique(x)))
	  } else {
	    Xcatp <- cbind(rep(0,npred));
	    ncat <- 0
	    Cvec <- 0
  	  if(nmissing > 0) Mcatp <- cbind(rep(0,npred));
	  }

  } else {
	  Xconp <- cbind(rep(0,1));
	  Xcatp <- cbind(rep(0,1));
  	if(nmissing > 0) Mconp <- cbind(rep(0,1));
  	if(nmissing > 0) Mcatp <- cbind(rep(0,1));
  }

  # Need to compute Gower distance
  if(nmissing == 0 & similarity_function == 4){
	  dissim <- as.matrix(cluster::daisy(Xall, metric="gower"))
    dissimtn <- dissim[1:nobs, 1:nobs]
    dissimtt <- dissim[-c(1:nobs), 1:nobs]
  }else{
    dissimtn <- 0
    dissimtt <- 0
  }


  # Create empty vectors that will hold MCMC iterates
  mu <- sig2 <- Si <- like <- ispred <- zi <- isordpred <- matrix(1,nrow=nout,ncol=nobs)
  mu0 <- sig20 <- nclus <- rep(1,nout)
  ppred <- predclass <-  rbpred <-  matrix(1, nrow=nout, ncol=npred)
  predclass_prob <- matrix(1, nrow=nout, ncol=npred*nobs)
  WAIC <- lpml <- rep(1,1)

  if(nmissing > 0){

    Xcon[is.na(Xcon)] <- 999;Xconp[is.na(Xconp)] <- 999;
    Xcat[is.na(Xcat)] <- 999;Xcatp[is.na(Xcatp)] <- 999;

    message("There are", nmissing, "missing covariate values.
        The missing covariate values will be
        accommodated using extentions to the
        ppmx model detailed in Page et. al (2020)")

		C.out <- .C("mcmc_missing",
        as.integer(draws), as.integer(burn), as.integer(thin),
        as.integer(nobs),as.integer(ncon), as.integer(ncat),
        as.integer(Cvec), as.integer(PPM), as.integer(cohesion),
        as.integer(similarity_function), as.integer(consim),
        as.double(M),
        as.double(y),
        as.double(t(Xcon)),as.integer(t(Xcat)),
        as.integer(t(Mcon)), as.integer(t(Mcat)),
        as.integer(npred),
        as.double(t(Xconp)),as.integer(t(Xcatp)),
        as.integer(t(Mconp)), as.integer(t(Mcatp)),
        as.double(simParms),
        as.double(dissimtn),as.double(dissimtt),
        as.integer(calibrate),as.double(modelPriors),
        as.integer(verbose), as.double(mh),
        mu.out=as.double(mu), sig2.out=as.double(sig2),
        mu0.out=as.double(mu0), sig20.out=as.double(sig20),
        Si.out=as.integer(Si), nclus.out=as.integer(nclus),
        like.out=as.double(like), WAIC.out=as.double(WAIC),
        lpml.out=as.double(lpml),ispred.out=as.double(ispred),
        ppred.out=as.double(ppred),predclass.out=as.integer(predclass),
        rbpred.out=as.double(rbpred),
        predclass_prob.out=as.double(predclass_prob))


  } else {
   C.out <- .C("mcmcppmx",
        as.integer(draws), as.integer(burn), as.integer(thin),
        as.integer(nobs),as.integer(ncon), as.integer(ncat),
        as.integer(Cvec), as.integer(PPM), as.integer(cohesion),
        as.integer(similarity_function),
        as.integer(consim), as.double(M), as.double(y),
        as.double(t(Xcon)),as.integer(t(Xcat)),
        as.integer(npred),
        as.double(t(Xconp)),as.integer(t(Xcatp)),
        as.double(simParms),
        as.double(t(dissimtn)),as.double(t(dissimtt)),
        as.integer(calibrate),
        as.double(modelPriors),as.integer(verbose),
        as.double(mh),
        mu.out=as.double(mu), sig2.out=as.double(sig2),
        mu0.out=as.double(mu0), sig20.out=as.double(sig20),
        Si.out=as.integer(Si), nclus.out=as.integer(nclus),
        like.out=as.double(like), WAIC.out=as.double(WAIC),
        lpml.out=as.double(lpml),ispred.out=as.double(ispred),
        ppred.out=as.double(ppred),predclass.out=as.integer(predclass),
        rbpred.out=as.double(rbpred),
        predclass_prob.out=as.double(predclass_prob))

  }

  out$mu <- matrix(C.out$mu.out, nrow=nout, byrow=TRUE)
  out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
  out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
  out$like <- matrix(C.out$like.out, nrow=nout, byrow=TRUE)
  out$fitted <- matrix(C.out$ispred.out, nrow=nout, byrow=TRUE)
  out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
  out$predclass <- matrix(C.out$predclass.out, nrow=nout, byrow=TRUE)
  out$predclass_prob <- matrix(C.out$predclass_prob.out, nrow=nout, byrow=TRUE)
  out$rbpred <- matrix(C.out$rbpred.out, nrow=nout, byrow=TRUE)
  out$mu0 <- C.out$mu0.out
  out$sig20 <- C.out$sig20.out
  out$nclus <- C.out$nclus.out
  out$WAIC <- C.out$WAIC.out
  out$lpml <- C.out$lpml.out
  out
}

#-----------------------------------------------------------------------------------------------------
