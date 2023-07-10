tpower <- function(x, t, p) {
  # Function for truncated p-th power function
  return((x - t) ^ p * (x > t))
}
bbase <- function(x, ndx, bdeg = 3, eps = 1e-5) {
  xl = min(x)
  xr = max(x)
  dx <- (xr - xl)/ndx
  knots <- seq(xl - bdeg*dx, xr + bdeg*dx, by=dx)
  P <- outer(x, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1) / (gamma(bdeg + 1) * dx ^ bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B[B < eps] <- 0
  attr(B,"knots") <- knots
  attr(B,"bdeg") <- bdeg
  attr(B,"eps") <- eps
  class(B) <- c("bbase")
  B
}


# Prediction function
predict.bbase <- function(object, newx) {
  knots <- attr(object,"knots")
  bdeg <- attr(object,"bdeg")
  eps <- attr(object,"eps")

  dx <- diff(knots)[1]
  P <- outer(newx, knots, tpower, bdeg)
  n <- dim(P)[2]
  D <- diff(diag(n), diff = bdeg + 1)/(gamma(bdeg + 1)*dx^bdeg)
  B <- (-1) ^ (bdeg + 1) * P %*% t(D)
  B[B < eps] = 0
  B
}

