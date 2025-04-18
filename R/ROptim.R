ROptim <- function(
    S,
    R,
    Rinv,
    lambda,
    tol.inner = 1e-4,
    tol.outer = 1e-4,
    max.inner.iter = 10,
    max.outer.iter = 100) {
  if (is.null(Rinv)) {
    Rinv <- solve(R)
  }

  # TODO: tol.inner, max.inner.iter

  p <- dim(R)[1]
  lambda_matrix <- (matrix(rep(1, p * p), ncol = p) - diag(p)) * lambda
  ans <- ROptim_to_fortran(
    S, rho = lambda_matrix, thr = tol.outer, maxIt = max.outer.iter,
    start = "warm", w.init = Rinv, wi.init = R
  )

  R <- ans$wi
  Rinv <- ans$w
  outer.count <- ans$niter
  loglik <- NA

  return(list(R = R, Rinv = Rinv, outer.count = outer.count, loglik = loglik))
}


# GLASSO algorithm of Friedman et al. 2008 with FORTRAN implementation of Sustik and Calderhead 2012.
# Ported to R by J. Clavel <julien.clavel@hotmail.fr> / <clavel@biologie.ens.fr> - 2017.

#' Fast optimization for the corelation matrix
#'
#' Part of the [pcglassoFast()] function.
#' This implementation is a small modification of the `glassoFast` function from `glassoFast` package.
#'
#' @param S Covariance matrix (a p by p symmetric matrix).
#' @param rho Regularization parameter (a non-negative value or a p by p matrix).
#' @param thr Threshold for convergence. Default is 1e-4.
#' @param maxIt Maximum number of iterations. Default is 10,000.
#' @param start Type of start: \code{"cold"} or \code{"warm"}.
#' @param w.init Optional starting values for the estimated covariance matrix (p x p). Used only for warm starts.
#' @param wi.init Optional starting values for the inverse covariance matrix (p x p). Used only for warm starts.
#' @param trace Logical. If \code{TRUE}, prints iteration info.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{w}{Estimated covariance matrix}
#'   \item{wi}{Estimated inverse covariance matrix}
#'   \item{errflag}{Memory allocation error flag: 0 = no error}
#'   \item{niter}{Number of iterations}
#' }
#'
#' @keywords glasso covariance matrix regularization penalized likelihood
ROptim_to_fortran <- function(S, rho, thr = 1.0e-4, maxIt = 1e4, start = c("cold", "warm"), w.init = NULL, wi.init = NULL, trace = FALSE) {
  n <- nrow(S) # dimension of S
  if (is.matrix(rho)) {
    if (length(rho) != n * n) stop("The input matrix for \"rho\" must be of size ", n, " by ", n)
    L <- rho # matrix of regularization parameters
  } else {
    L <- matrix(rho, n, n) # matrix of regularization parameters
  }

  # cold or warm start
  start.type <- match.arg(start)
  if (start.type == "cold") {
    is <- 0
    W <- X <- matrix(0, nrow = n, ncol = n)
  }
  if (start.type == "warm") {
    is <- 1
    if (is.null(w.init) | is.null(wi.init)) {
      stop("Warm start specified: w.init and wi.init must be non-null")
    }
    W <- w.init
    X <- wi.init
  }

  Wd <- WXj <- numeric(n)

  msg <- 1 * trace
  info <- 0
  mode(n) <- "integer"
  mode(S) <- "double"
  mode(L) <- "double"
  mode(thr) <- "double"
  mode(maxIt) <- "integer"
  mode(msg) <- "integer"
  mode(is) <- "integer"
  mode(X) <- "double"
  mode(W) <- "double"
  mode(info) <- "integer"


  LASSO <- .Fortran(
    "roptim",
    n,
    S,
    L,
    thr,
    maxIt,
    msg,
    is,
    X,
    W,
    Wd,
    WXj,
    info
  )

  results <- list(w = LASSO[[9]], wi = LASSO[[8]], errflag = LASSO[[12]], niter = LASSO[[5]])
  return(results)
}
