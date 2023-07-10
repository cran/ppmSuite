# Wrapper the executes the C code with the .C function
# I am standardizing the x's.

gaussian_ppmx <- function(y,X=NULL,Xpred=NULL,
                           meanModel=1,
                           cohesion=1, M=1,
                           PPM = FALSE,
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
  # meanModel = 1 => no covariates in the likelihood
  # meanModel = 2 => global regression included in the likelihood
  # If there is missing, then only meanModel 1 is available.
  out <- NULL

  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)
  if(!is.data.frame(Xpred) & !is.null(Xpred)){
    Xpred <- data.frame(Xpred);
    colnames(Xpred) <- colnames(X)
  }

  cnames <- colnames(X)


  nout <- (draws-burn)/thin
  nobs <- length(y)

  # If X and Xpred are both NULL I need to
  # initialize all vectors to zero.
  Xcon <- cbind(rep(0,1));
  Xconp <- cbind(rep(0,1));
  ncon <- 0
  Mcon <- cbind(rep(0,1))
  Mconp <- cbind(rep(0,1));

  Xcat <- cbind(rep(0,1));
  Xcatp <- cbind(rep(0,1));
  Cvec <- 0
  ncat <- 0
  Mcat <- cbind(rep(0,1))
  Mcatp <- cbind(rep(0,1))

  nmissing <- 0
  npred <- 0

  # If at least one of X and or Xpred are not NULL
  if(!(is.null(X) & is.null(Xpred))){
    nxobs <- ifelse(is.null(X), 0, nrow(X))
    npred <- ifelse(is.null(Xpred), 0, nrow(Xpred))

    Xall <- rbind(X, Xpred)
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
        mn <- apply(Xall[,!catvars, drop=FALSE],2,mean, na.rm=TRUE)
        sd <- apply(Xall[,!catvars, drop=FALSE],2,sd, na.rm=TRUE)
        if(sum(round(mn + sd,10) == 1) != length(mn)){
          message("Note: It appears that continuous covariates are not standardized using both training and testing observations.  The default similarity values may not be appropriate.")
        }
  	    Xconstd <- Xall[,!catvars, drop=FALSE]
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
    }

    # Now consider the case when number of covariates for prediction are greater than zero
    if(npred > 0){
  	  if(sum(!catvars) > 0){
  	    Xconstd <- Xall[,!catvars, drop=FALSE]
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

    }

  }



  # if Gower's dissimilarity is selected
  if(nmissing > 0 & similarity_function == 4){
    stop("Gower's dissimilarity similarity function cannot be used with missing covariates")
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



  if(nmissing > 0){

    Xcon[is.na(Xcon)] <- 999;Xconp[is.na(Xconp)] <- 999;
    Xcat[is.na(Xcat)] <- 999;Xcatp[is.na(Xcatp)] <- 999;

    message("Note: There are a total of ", nmissing, " missing covariate values. They will be accommodated using extensions to the ppmx model detailed in Page et. al (2022)")

    run <- .Call("GAUSSIAN_PPMX_MISSING",
              as.double(y), as.integer(nobs),
              as.double(t(Xcon)), as.integer(t(Mcon)), as.integer(ncon),
              as.integer(t(Xcat)), as.integer(t(Mcat)),as.integer(ncat), as.integer(Cvec),
              as.integer(npred),
              as.double(t(Xconp)), as.integer(t(Mconp)),
              as.integer(t(Xcatp)), as.integer(t(Mcatp)),
              as.double(M), as.integer(meanModel), as.double(modelPriors), as.double(simParms),
              as.integer(PPM), as.integer(cohesion), as.integer(similarity_function),as.integer(consim),
              as.double(dissimtn), as.double(dissimtt), as.integer(calibrate), as.double(mh),
              as.integer(verbose), as.integer(draws), as.integer(burn), as.integer(thin))


  } else {
    run <- .Call("GAUSSIAN_PPMX",
                  as.double(y), as.integer(nobs),
                  as.double(t(Xcon)), as.integer(t(Xcat)), as.integer(ncon),
                  as.integer(ncat), as.integer(Cvec),
                  as.double(t(Xconp)), as.integer(t(Xcatp)), as.integer(npred),
                  as.integer(meanModel), as.double(modelPriors), as.double(mh),
                  as.integer(PPM), as.integer(cohesion), as.integer(similarity_function),
                  as.integer(consim), as.double(M), as.double(simParms),
                  as.double(dissimtn), as.double(dissimtt), as.integer(calibrate),
                  as.integer(verbose),
                  as.integer(draws), as.integer(burn), as.integer(thin))


  }
  if(meanModel == 2) colnames(run$beta) <- c(cnames[!catvars], cnames[catvars])
  if(meanModel == 1) run <- run[-3]
  out <- run

  if(nmissing > 0) out$Missmat <- Mall

  out
}

#-----------------------------------------------------------------------------------------------------
