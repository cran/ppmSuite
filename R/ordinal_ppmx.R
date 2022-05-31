## I am standardizing the x's.  For
ordinal_ppmx  <- function(y, co, X=NULL,Xpred=NULL,
                      meanModel=1,
                  		cohesion=1, M=1,
                      PPM = FALSE,
                  		similarity_function=1, consim=1,
                  		calibrate=0,
                  		simParms=c(0.0, 1.0, 0.1, 1.0, 2.0, 0.1, 1),
                  		modelPriors=c(0, 10, 1, 1),
                      mh=c(0.5, 0.5),
                  		draws=1100,burn=100,thin=1, verbose=FALSE){

  # X - is a data.frame whose columns that are factors or characters are
  #     considered categorical
  #
  # consim argument
  # 1 - similarity function of continuous covariate is N-N model (v_j is fixed)
  # 2 - similarity function of continuous covariate is N-NIG model (v_j is unknown)
  #
  # cohesion argument:
  # 1 - DP style
  # 2 - uniform (i.e., c(S_j) = 1)
  #
  # similarity_function argument:
  # 1 - Auxiliary
  # 2 - Double Dipper
  # 3 - Variance
  # 4 - Gower's dissimilarity.  This one only available if no missing
  #
  # calibrate:
  # 0 - no calibration
  # 1 - calibration using "calibrating" similarity
  # 2 - calibration using "coarsening" similarity
  #
  # meanModel = 1 => no covariates in the likelihood
  # meanModel = 2 => global regression included in the likelihood
  #
  # If there is missing, then only meanModel 1 is available.
  #
  # simparms: parameter values for the similarity functions
  #           all are used for marginal likelihood similarities
  #           save the last entry which is alpha and is use
  #           the variance or gower dissimilarity similarity functions
  # c(m0, s20, v, k0, nu0, a0, alpha)
  #
  out <- NULL

  if(!is.data.frame(X) & !is.null(X)) X <- data.frame(X)
  if(!is.data.frame(Xpred) & !is.null(Xpred)){
    Xpred <- data.frame(Xpred);
    colnames(Xpred) <- colnames(X)
  }

  cnames <- colnames(X)

  if(is.null(X) & meanModel != 1){
    stop("No training covariates are provided and a regression is included in the mean model")
  }

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
          message("Note: It appears that continuous covariates are not standardized using both training and testing observations.")
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


  nordcat <- length(co)

  if(nmissing > 0){
    Xcon[is.na(Xcon)] <- 999;Xconp[is.na(Xconp)] <- 999;
    Xcat[is.na(Xcat)] <- 999;Xcatp[is.na(Xcatp)] <- 999;

    message(" There are a total of ", nmissing, " missing covariate values. \n They will be accommodated using extentions to the ppmx model detailed in Page et. al (2020).")

    run <- .Call("ORDINAL_PPMX_MISSING",
                  as.integer(y), as.double(co), as.integer(nobs), as.integer(nordcat),
                  as.double(t(Xcon)), as.integer(t(Mcon)), as.integer(ncon),
                  as.integer(t(Xcat)), as.integer(t(Mcat)), as.integer(ncat), as.integer(Cvec),
                  as.integer(npred),
                  as.double(t(Xconp)), as.integer(t(Mconp)),
                  as.integer(t(Xcatp)), as.integer(t(Mcatp)),
                  as.double(M), as.integer(meanModel), as.double(modelPriors), as.double(simParms),
                  as.integer(PPM), as.integer(cohesion), as.integer(similarity_function), as.integer(consim),
                  as.double(dissimtn), as.double(dissimtt), as.integer(calibrate), as.double(mh),
                  as.integer(verbose), as.integer(draws), as.integer(burn), as.integer(thin))


  } else {
    run <- .Call("ORDINAL_PPMX",
                  as.integer(y), as.double(co), as.integer(nobs), as.integer(nordcat),
                  as.double(t(Xcon)), as.integer(ncon),
                  as.integer(t(Xcat)), as.integer(ncat), as.integer(Cvec),
                  as.integer(npred),
                  as.double(t(Xconp)),
                  as.integer(t(Xcatp)),
                  as.double(M), as.integer(meanModel), as.double(modelPriors), as.double(simParms),
                  as.integer(PPM), as.integer(cohesion), as.integer(similarity_function), as.integer(consim),
                  as.double(dissimtn), as.double(dissimtt), as.integer(calibrate), as.double(mh),
                  as.integer(verbose), as.integer(draws), as.integer(burn), as.integer(thin))


  }


  if(meanModel == 2) colnames(run$beta) <- c(cnames[!catvars], cnames[catvars])
  if(meanModel == 1) run <- run[-3]
  out <- run


  if(nmissing > 0) out$Missmat <- Mall

  out


}

