ROptim <- function(
    S,
    R,
    Rinv,
    lambda,
    tol = 1e-4,
    max_inner_iter = 10,
    max_outer_iter = 100) {
  if (is.null(Rinv)) {
    Rinv <- solve(R)
  }

  p <- dim(R)[1]
  lambda_matrix <- matrix(lambda, p, p); diag(lambda_matrix) <- 0
  ans <- ROptim_to_fortran(
    S,
    rho = lambda_matrix, thr = tol, maxIt = max_outer_iter,
    maxItLasso = max_inner_iter, start = "warm", Rinv_init = Rinv, R_init = R
  )

  if (ans$niter == max_outer_iter + 1) {
    rlang::warn(paste0("Optimization of R matrix reached the max number of iterations (", max_outer_iter, "). Consider increasing the `max_iter_R_outer` parameter in `pcglassoFast()` function."))
  }

  smallest_eigen_value <- eigen(ans$wi, TRUE, TRUE)$values[p]
  if (smallest_eigen_value < 0) {
    rlang::warn("Optimization of R matrix resulted in a matrix that is not positive definite. We highly encourage to increse the `max_iter_R_outer` parameter in `pcglassoFast()` function.")

    desired_smallest_eigen_value <- 0.01
    x <- (1-desired_smallest_eigen_value) / (1-smallest_eigen_value)
    ans$wi <- x*ans$wi + diag(1-x, p)
    ans$w <- solve(ans$wi)
  }

  list(
    R           = ans$wi,
    Rinv        = ans$w,
    outer.count = ans$niter,
    loglik      = NA
  )
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
#' @param Rinv_init Optional starting values for the estimated covariance matrix (p x p). Used only for warm starts.
#' @param R_init Optional starting values for the inverse covariance matrix (p x p). Used only for warm starts.
#' @param trace Logical. If \code{TRUE}, prints iteration info.
#'
#' @noRd
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
ROptim_to_fortran <- function(S, rho, thr = 1.0e-4, maxIt = 1e4, maxItLasso = 500, start = c("cold", "warm"), Rinv_init = NULL, R_init = NULL, trace = FALSE) {
  n <- nrow(S) # dimension of S
  if (is.matrix(rho)) {
    if (length(rho) != n * n) stop("The input matrix for \"rho\" must be of size ", n, " by ", n)
    L <- rho # matrix of regularization parameters
  } else {
    L <- matrix(rho, n, n) # matrix of regularization parameters
  }

  # cold or warm start
  switch(match.arg(start),
    cold = {
      is_flag <- 0
      W <- X <- matrix(0, n, n)
    },
    warm = {
      is_flag <- 1
      if (is.null(Rinv_init) || is.null(R_init)) {
        stop("Warm start specified: Rinv_init and R_init must be non-null")
      }
      W <- Rinv_init
      X <- R_init
    }
  )

  Wd <- WXj <- numeric(n)

  msg <- 1 * trace
  info <- 0
  mode(n) <- "integer"
  mode(S) <- "double"
  mode(L) <- "double"
  mode(thr) <- "double"
  mode(maxIt) <- "integer"
  mode(maxItLasso) <- "integer"
  mode(msg) <- "integer"
  mode(is_flag) <- "integer"
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
    maxItLasso,
    msg,
    is_flag,
    X,
    W,
    Wd,
    WXj,
    info
  )

  results <- list(w = LASSO[[10]], wi = LASSO[[9]], errflag = LASSO[[13]], niter = LASSO[[5]])
  return(results)
}
