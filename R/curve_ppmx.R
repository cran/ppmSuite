
#-----------------------------------------------------------------------------------------------------


## R wrapper that accesses C code to fit curve cluster model

curve_ppmx <- function(y, z, subject,
                       Xcon=NULL,Xcat=NULL,
                       Xconp=NULL,Xcatp=NULL,
                       PPM, M,
                       q=3, rw_order=1, balanced=1,
                       nknots,npredobs,
                       Aparm, modelPriors,
                       similarity_function=1,
                       consim,calibrate,
                       simParms,
                       mh=c(1,1),
                       draws=1100,burn=100,thin=1){



	## consim=1 => that similarity function of continuous covariate is N-N model (v_j is fixed)
	## consim=2 => that similarity functino of continuous covariate is N-NIG model (v_j is unknown)

  out <- NULL
  nobs <- tapply(subject, subject, length)
  nsubject <- length(nobs)

  if(is.null(Xcon)){
	  ncon <- 0; Xcon <- cbind(rep(0,nsubject))
  }else{
    ncon <- ncol(Xcon)
  }
  cat("ncon = ", ncon, "\n")

  if(is.null(Xcat)){
    ncat <- 0; Xcat <- cbind(rep(0,nsubject))
    Cvec <- 0
  }else{
    ncat <- ncol(Xcat)
  	Cvec <- apply(Xcat,2,function(x)length(unique(x)))
  }
  cat("ncat = ", ncat, "\n")


  if(is.null(Xconp) & is.null(Xcatp)){
  	npred <- 0
  } else {
	  npred <- ifelse(is.null(Xconp),nrow(Xcatp), nrow(Xconp))
  }

  cat("npred = ", npred, "\n")

  if(is.null(Xconp)){
    Xconp <- cbind(rep(0,npred))
  }
  if(is.null(Xcatp)){
    Xcatp <- cbind(rep(0,npred))
  }


  ndx <- nknots  # number of segments (inner knots)

  # Best to feed in each subjects basis created from Cote's code rather than
  # creating code in C. Not the Hmat matrix depends on the balanced argument
  #
  # balanced 1 - All subjects have same number of measurements and time at which
  #              they were measure and so have same design matrix
  #		       0 - Subjects do not have same number of measurements and they are
  #		           measured at different time points so each subject needs their
  #              own design matrix

  # Unique time points
  tt <- sort(unique(z))

  # Create the basis
  G <- bbase(tt, ndx = ndx, bdeg = q)

  Hmat <- G

  # For prediction
  tpred <- seq(min(z), max(z)+npredobs, length = length(tt)+npredobs)
  Gp <- predict.bbase(G, tpred)
  Hmat_pred <- Gp

  D <- diff(diag(ncol(G)), differences = rw_order)
  K <- crossprod(D)
  if(rw_order != 1) stop("random walk orders not equal to 1 are not currently available.")

  K[1,1] <- 1+1

  nb <- ncol(K)

  if(!balanced) {
    Hmat <- NULL
    Hmat_pred <- NULL
    cnobs <- 0
    for(j in 1:nsubject){
      ttmp <- z[(cnobs+1):(cnobs + nobs[j])]
      bsb <- predict.bbase(G, ttmp)
      Hmat <- c(Hmat, c(t(bsb)))

      ttmp <- c(z[(cnobs+1):(cnobs + nobs[j])],
                (z[(cnobs + nobs[j])]+1):(z[(cnobs + nobs[j])]+npredobs))
      bsb <- predict.bbase(G, ttmp)
      Hmat_pred <- c(Hmat_pred, c(t(bsb)))

      cnobs <- cnobs + nobs[j]

    }
  }



  nout <- (draws-burn)/thin;

  if(ncol(y) == 1){

    Si <- llike <- ispred <- lam <- matrix(1,nrow=nout,ncol=nsubject)
  	sig2 <- beta0 <-  matrix(1,nrow=nout,ncol=nsubject)
  	iCPO <- tau2 <- matrix(1,nrow=nout,ncol=nsubject)

  	beta <- theta <- matrix(0, nrow = nout, ncol = nb*nsubject)
  	mu <- matrix(0, nrow=nout, ncol=nb)
  	nclus <- mub0 <- sig2b0 <- rep(1, nout)
  	predclass <- matrix(1, nrow=nout, ncol=npred)
  	ppred <- matrix(0, nrow=nout, ncol=npred*npredobs)
  	predlines <- matrix(0, nrow=nout, ncol=sum(nobs+npredobs))
  	lpml <- WAIC <- rep(0,1)


  	C.out <- .C("mcmc_curvecluster",
  	              as.integer(draws), as.integer(burn),as.integer(thin),
  	              as.integer(nsubject), as.integer(nobs),
  	              as.double(y), as.double(t(z)),
  	              as.double(t(K)), as.integer(nb),
  	              as.integer(ncon), as.integer(ncat), as.integer(Cvec),
  	              as.double(t(Xcon)), as.integer(t(Xcat)),
  	              as.integer(PPM), as.double(M),
  	              as.integer(similarity_function), as.integer(consim),
  	              as.integer(npred),as.integer(npredobs),
  	              as.double(t(Xconp)),as.integer(t(Xcatp)),
  	              as.double(simParms), as.double(Aparm),
  	              as.integer(calibrate), as.double(modelPriors),
  	              as.double(t(Hmat)), as.integer(balanced), as.double(mh),
  	              as.double(t(Hmat_pred)),
  	              beta.out= as.double(beta),beta0.out= as.double(beta0),
  	              sig2.out= as.double(sig2), mub0.out= as.double(mub0),
  	              sig2b0.out= as.double(sig2b0), lam.out= as.double(lam),
  	              tau2.out= as.double(tau2), theta.out= as.double(theta),
  	              mu.out= as.double(mu),Si.out= as.integer(Si),
  	              nclus.out=as.integer(nclus),
  	              ppred.out=as.double(ppred),
  	              predclass.out=as.integer(predclass),
  	              predlines.out=as.double(predlines),
  	              llike.out=as.double(llike),
  	              lpml.out=as.double(lpml),WAIC.out=as.double(WAIC))

    out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
  	out$nclus <- matrix(C.out$nclus.out, nrow=nout, byrow=TRUE)
  	out$beta <- array(C.out$beta.out, c(nsubject, nb, nout))
  	out$theta <- array(C.out$theta.out, c(nsubject, nb, nout))
  	out$sig2 <- matrix(C.out$sig2.out, nrow=nout, byrow=TRUE)
  	out$tau2 <- matrix(C.out$tau2.out, nrow=nout, byrow=TRUE)
  	out$mu <- matrix(C.out$mu.out, nrow=nout, byrow=TRUE)
  	out$lam <- matrix(C.out$lam.out, nrow=nout, byrow=TRUE)
  	out$beta0 <- matrix(C.out$beta0.out, nrow=nout, byrow=TRUE)
  	out$sig2b0 <- matrix(C.out$sig2b0.out, nrow=nout, byrow=TRUE)
  	out$mub0 <- matrix(C.out$mub0.out, nrow=nout, byrow=TRUE)
  	out$llike <- matrix(C.out$llike.out, nrow=nout, byrow=TRUE)
  	#	out$fitted <- matrix(C.out$ispred.out, nrow=nout, byrow=TRUE)
  	#	out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
  	#	out$predclass <- matrix(C.out$predclass.out, nrow=nout, byrow=TRUE)
  	out$pred_curve <- matrix(C.out$predlines.out, nrow=nout, byrow=TRUE)
  	out$WAIC <- C.out$WAIC.out
  	out$lpml <- C.out$lpml.out

  }


  if(ncol(y) == 2){
    y1 <- y[,1]
    y2 <- y[,2]

    Si <- llike <- ispred <- lam1 <- lam2 <-  matrix(1,nrow=nout,ncol=nsubject)
    sig21 <- sig22 <- beta01 <- beta02 <-  matrix(1,nrow=nout,ncol=nsubject)
    iCPO <- tau21 <- tau22 <- matrix(1,nrow=nout,ncol=nsubject)

    beta1 <- beta2 <- theta1 <- theta2 <- matrix(0, nrow = nout, ncol = nb*nsubject)
    mu1 <- mu2 <- matrix(0, nrow=nout, ncol=nb)
    nclus <- mub01 <- mub02 <- sig2b01 <- sig2b02 <- rep(1, nout)
    predclass <- matrix(1, nrow=nout, ncol=npred)
    ppred1 <- ppred2 <- matrix(0, nrow=nout, ncol=npred*npredobs)
    lpml <- WAIC <- rep(0,1)

    C.out <- .C("mcmc_bivariate_curvecluster",
                  as.integer(draws), as.integer(burn),as.integer(thin),
                  as.integer(nsubject),as.integer(nobs),
                  as.double(y1), as.double(y2), as.double(t(z)),
                  as.double(t(K)), as.integer(nb),
                  as.integer(ncon), as.integer(ncat),as.integer(Cvec),
                  as.double(t(Xcon)),as.integer(t(Xcat)),
                  as.integer(PPM), as.double(M),
                  as.integer(similarity_function), as.integer(consim),
                  as.integer(npred),as.integer(npredobs),
                  as.double(t(Xconp)),as.integer(t(Xcatp)),
                  as.double(simParms), as.double(Aparm),
                  as.integer(calibrate),as.double(modelPriors),
                  as.double(t(Hmat)), as.integer(balanced),  as.double(mh),
                  beta1.out= as.double(beta1), beta2.out= as.double(beta2),
                  beta01.out= as.double(beta01), beta02.out= as.double(beta02),
                  sig21.out= as.double(sig21), sig22.out= as.double(sig22),
                  mub01.out= as.double(mub01), mub02.out= as.double(mub02),
                  sig2b01.out= as.double(sig2b01), sig2b02.out= as.double(sig2b02),
                  lam1.out= as.double(lam1), lam2.out= as.double(lam2),
                  tau21.out= as.double(tau21), tau22.out= as.double(tau22),
                  theta1.out= as.double(theta1), theta2.out= as.double(theta2),
                  mu1.out= as.double(mu1), mu2.out= as.double(mu2),
                  Si.out= as.integer(Si), nclus.out=as.integer(nclus),
                  ppred1.out=as.double(ppred1),ppred2.out=as.double(ppred2),
                  predclass.out=as.integer(predclass), llike.out=as.double(llike),
                  lpml.out=as.double(lpml),WAIC.out=as.double(WAIC))

    out$Si <- matrix(C.out$Si.out, nrow=nout, byrow=TRUE)
    out$nclus <- matrix(C.out$nclus.out, nrow=nout, byrow=TRUE)
    out$beta1 <- array(C.out$beta1.out, c(nsubject, nb, nout))
    out$beta2 <- array(C.out$beta2.out, c(nsubject, nb, nout))
    out$theta1 <- array(C.out$theta1.out, c(nsubject, nb, nout))
    out$theta2 <- array(C.out$theta2.out, c(nsubject, nb, nout))
    out$sig21 <- matrix(C.out$sig21.out, nrow=nout, byrow=TRUE)
    out$sig22 <- matrix(C.out$sig22.out, nrow=nout, byrow=TRUE)
    out$tau21 <- matrix(C.out$tau21.out, nrow=nout, byrow=TRUE)
    out$tau22 <- matrix(C.out$tau22.out, nrow=nout, byrow=TRUE)
    out$mu1 <- matrix(C.out$mu1.out, nrow=nout, byrow=TRUE)
    out$mu2 <- matrix(C.out$mu2.out, nrow=nout, byrow=TRUE)
    out$lam1 <- matrix(C.out$lam1.out, nrow=nout, byrow=TRUE)
    out$lam2 <- matrix(C.out$lam2.out, nrow=nout, byrow=TRUE)
    out$beta01 <- matrix(C.out$beta01.out, nrow=nout, byrow=TRUE)
    out$beta02 <- matrix(C.out$beta02.out, nrow=nout, byrow=TRUE)
    out$sig2b01 <- matrix(C.out$sig2b01.out, nrow=nout, byrow=TRUE)
    out$sig2b02 <- matrix(C.out$sig2b02.out, nrow=nout, byrow=TRUE)
    out$mub01 <- matrix(C.out$mub01.out, nrow=nout, byrow=TRUE)
    out$mub02 <- matrix(C.out$mub02.out, nrow=nout, byrow=TRUE)
    out$llike <- matrix(C.out$llike.out, nrow=nout, byrow=TRUE)
    #	out$fitted <- matrix(C.out$ispred.out, nrow=nout, byrow=TRUE)
    #	out$ppred <- matrix(C.out$ppred.out, nrow=nout, byrow=TRUE)
    #	out$predclass <- matrix(C.out$predclass.out, nrow=nout, byrow=TRUE)
    out$WAIC <- C.out$WAIC.out
    out$lpml <- C.out$lpml.out


  }

  out$Hmat <- Hmat
  out$Hmat_pred <- Hmat_pred
  out$number.of.basis <- nb
  out
}



