######
#' blockwise optimization for pcglasso
#'
#' @param S (p x p matrix) empirical covariance matrix derived from the data.
#' @param lambda,alpha (double \eqn{0\le\lambda}, double \eqn{0\le\alpha < 1})
#'         Parameters of the method. See Details section below.
#'   * \eqn{\lambda} is a penalty for off-diagonal
#'   * \eqn{\alpha} is a penalty for on-diagonal
#' @param R,R_inv (p x p matrix, unit diagonal) initial estimation of
#'         precision matrix (it is recommended to left it default).
#' @param D (vector of length p) diagonal of initial estimation of diagonal
#'         matrix (it is recommended to left it default).
#' @param max.iter (integer) maximum number of iterations.
#' @param tolerance (double) tolerance for convergence.
#' @param R.tol.inner,R.tol.outer,R.max.inner.iter,R.max.outer.iter
#'         Parameters passed to [ROptim()] function.
#' @param D.tol,D.max.starting.iter,D.max.outer.iter
#'         Parameters passed to [DOptim()] function.
#'
#' @details
#' The function maximizes the
#' \eqn{f(R, D) = log(det(R)) + (1-\alpha)log(det(D^2)) - tr(DSDR) - \lambda ||R||_1}
#' function, where \eqn{||R||_1} is only for off-diagonal elements.
#'
#' The function employs coordinate descent,
#' also known as blockwise optimization,
#' to iteratively optimize the variables `R` and `D`
#' while fixing the other variable.
#' It continues this process until convergence or
#' until the maximum number of iterations is reached.
#'
#' @return list of three elements:
#' * "R" - found correlation matrix
#' * "D" - found diagonal matrix
#' * "n_iters" - number of iterations of the outer loop
#' @md
#'
#' @export
#' @examples
#' p <- 7
#' R.true <- toeplitz(c(1, -0.5, 0, 0, 0, 0, 0))
#' D.true <- sqrt(rchisq(p, 3))
#' S_inv.true <- diag(D.true) %*% R.true %*% diag(D.true)
#'
#' S <- solve(S_inv.true) # data
#'
#' alpha <- 4 / 20 # 4 / n, as in PCGLASSO paper
#'
#' pcglassoFast(S, 0.11, alpha, max.iter = 1000)
pcglassoFast <- function(
    S, lambda, alpha,
    R = diag(dim(S)[1]), R_inv = NULL, D = rep(1, dim(S)[1]),
    max.iter = 100, tolerance = 1e-6,
    R.tol.inner = 1e-2, R.tol.outer = 1e-3,
    R.max.inner.iter = 10, R.max.outer.iter = 100,
    D.tol = 1e-4,
    D.max.starting.iter = 500, D.max.outer.iter = 100) {
  D <- rep(1, dim(S)[1])
  stop_loop <- FALSE
  i <- 1
  loss_vec <- rep(0, max.iter)
  loss_D <- rep(0, max.iter)
  loss_old <- function_to_optimize(R, D, S, lambda, alpha)
  loss_vec[i] <- loss_old
  loss_D[i] <- loss_old
  while (!stop_loop & i < max.iter) {
    resD <- DOptim(
      A = R * S,
      D0 = D,
      tol = D.tol,
      max.starting.iter = D.max.starting.iter,
      max.outer.iter = D.max.outer.iter,
      alpha = alpha
    )
    D <- resD$D
    loss_D[i + 1] <- function_to_optimize(R, D, S, lambda, alpha)

    resR <- ROptim(
      S = sweep(sweep(S, 1, D, "*"), 2, D, "*"),
      R = R,
      lambda = lambda,
      Rinv = R_inv,
      tol.inner = R.tol.inner,
      tol.outer = R.tol.outer,
      max.inner.iter = R.max.inner.iter,
      max.outer.iter = R.max.outer.iter
    )

    R <- resR$R
    R_inv <- resR$Qinv

    loss_new <- function_to_optimize(R, D, S, lambda, alpha)
    stop_loop <- (loss_new - loss_old < tolerance)
    loss_old <- loss_new

    i <- i + 1
    loss_vec[i] <- loss_old
  }

  loss_vec <- loss_vec[1:i]
  loss_D <- loss_D[1:i]
  D <- as.vector(D)
  return(list(
    "Sinv" = diag(D) %*% R %*% diag(D),
    "R" = R,
    "D" = D,
    "R_inv" = R_inv,
    "n_iters" = i,
    "loss" = loss_vec,
    "loss_D" = loss_D
  ))
}


#' Function that the `pcglassoFast()` is maximizing
#'
#' In `ROptim()`, this function is maximized with respect to `R`.
#' In `DOptim()`, this function is maximized with respect to `D`.
#'
#' \eqn{function\_to\_optimize(R, D) = log(det(R)) + (1-\alpha)log(det(D^2)) - tr(DSDR) - \lambda ||R||_1}
#' function, where \eqn{||R||_1} is only for off-diagonal elements.
function_to_optimize <- function(R, d, S, lambda, alpha) {
  my_norm_1 <- function(my_matrix) {
    diag(my_matrix) <- 0 # no diagonal elements
    sum(abs(my_matrix))
  }

  2 * sum(log(diag(chol(R)))) + (1 - alpha) * 2 * sum(log(d)) - trace_DSDR(d, S, R) - lambda * my_norm_1(R)
}

#' compute tr(DSDR) where D are diagonal matrices
trace_DSDR <- function(d, S, R) {
  # return(sum(c(S)*rep(d,times=length(d))*rep(d,each=length(d))*c(R)))
  return(sum(d * (S * R) %*% d))
}
