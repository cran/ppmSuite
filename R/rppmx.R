#-----------------------------------------------------------------------------------------------------
# rppm.R
# This produces one draw from a Product Partition Distribution.
# The arguments are
# m - number of objects that will be assigned to clusters
# M - is a parameter that influens co-clustering probabilities

rppmx <- function(m, X=NULL, similarity, simparm,M=1,m0=0,s20=1,v=2,k0=10,v0=1,alpha=1){



  out <- NULL

  ppm <- ifelse(is.null(X), TRUE, FALSE)

  if(!ppm){
    nobs <- nrow(X)
    if(!is.data.frame(X)) X <- data.frame(X)

    classes <- sapply(X, class)
    catvars <- classes %in% c("factor","character")

    # standardize continuous covariates
    if(sum(!catvars) > 0){
  	  xcon <- apply(X[,!catvars, drop=FALSE], 2, scale)
  	  ncon <- ncol(xcon)
    }else{
      xcon <- cbind(rep(0,nobs));
      ncon <- 0
    }


    # Function that relabels categorical variables to begin with 0
    relab <- function(x) as.numeric(as.factor(as.character(x))) - 1

    if(sum(catvars) > 0){
  	  # Change the factors or characters into integers with category starting at 0.
  	  xcat <- apply(X[, catvars,drop=FALSE], 2, relab)
  	  Cvec <- apply(xcat,2,function(x)length(unique(x)))
  	  ncat <- ncol(xcat)
    }else{
  	  xcat <- cbind(rep(0,nobs));
  	  Cvec <- 0
  	  ncat <- 0
    }
  } else { # If draw from PPM don't need the following.
    similarity = simparm = dirweights = Cvec = 0;
    ncon = ncat = xcon = xcat = 0;
    m0 = k0 = v0 = s20 = v = alpha = 0;
  }
	nk <- 1
	nh <- rep(0,m)
	Si <- rep(0,m)

	dirweights <- rep(0.1, length=max(Cvec));
	N <- m

	C.out <- .C("rppmx",
                as.integer(N),as.integer(similarity),
	              as.integer(simparm),as.double(M),
              	as.integer(ncon), as.integer(ncat),
              	as.double(t(xcon)), as.integer(t(xcat)),
              	as.integer(Cvec), as.integer(ppm),
              	as.double(m0),as.double(k0),
              	as.double(v0), as.double(s20),
              	as.double(v),as.double(dirweights),
              	as.double(alpha),
              	Si.out = as.integer(Si),
              	nclus.out = as.integer(nk),
              	nh.out = as.integer(nh))

	Si <- C.out$Si.out
	nclus <- C.out$nclus.out
	nh <- C.out$nh.out
	out$Si <- Si
	out$nclus <- nclus
	out$nh <- nh[1:nclus]
	out

}



