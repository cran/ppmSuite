


K1 <- function(dim){
  K <- Matrix::bandSparse(dim, k=-c(1), diagonals=list(rep(-1,dim)), symmetric=TRUE)
  diag(K) <- c(1,rep(2,dim-2),1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}

K2 <- function(dim){
  tmp <- list(c(-2,rep(-4,dim-3),-2),rep(1,dim))
  K <- Matrix::bandSparse(dim,k=-c(1:2), diagonals=tmp,symmetric=TRUE)
  diag(K) <- c(1,5,rep(6,dim-4),5,1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}

K3 <- function(dim){
  tmp <- list(c(-3,-12,rep(-15,dim-5), -12,-3),
              c(3,rep(6,dim-4),3),
              rep(-1,dim))
  K <- Matrix::bandSparse(dim,k=-c(1:3), diagonals=tmp,symmetric=TRUE)
  diag(K) <- c(1,10,19,rep(20,dim-6),19,10,1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}

K4 <- function(dim){
  tmp <- list(c(-4,-28,-52,rep(-56,dim-6-1),-52,-28,-4),
              c(6,22,rep(28,dim-4-2),22,6),
              c(-4,rep(-8,dim-2-3),-4),
              rep(1,dim-4))
  K <- Matrix::bandSparse(dim,k=-c(1:4), diagonals=tmp,symmetric=TRUE)
  diag(K) <- c(1,17,53,69, rep(70,dim-8),69,53,17,1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}


K5 <- function(dim){
  tmp <- list(c(-5,-55,-155,-205,rep(-210,dim-8-1),-205,-155,-55,-5),
              c(10,60,110,rep(120,dim-6-2),110,60,10),
              c(-10,-35,rep(-45,dim-4-3),-35,-10),
              c(5,rep(10,dim-2-4),5),
              rep(-1,dim-5))
  K <- Matrix::bandSparse(dim,k=-c(1:5), diagonals=tmp,symmetric=TRUE)
  diag(K) <- c(1,26,126,226,251, rep(252,dim-10),251,226,126,26,1)
  K <- matrix(as.numeric(K),nrow=dim, byrow=TRUE)
  K
}


tpower <- function(x, t, p) {
  # Function for truncated p-th power function
  return((x - t) ^ p * (x > t))
}
bbase <- function(x, xl = min(x), xr = max(x), ndx = 10, bdeg = 3, eps = 1e-5) {
  # Function for B-spline basis
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B[B < eps] = 0
  res <- list(B = B, knots = knots)
  res
}
