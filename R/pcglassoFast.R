#' blockwise optimization for pcglasso
#'
#' @param S (p x p matrix) empirical covariance matrix derived from the data.
#' @param lambda,alpha (double \eqn{0\le\lambda}, double \eqn{\alpha \in \mathbb{R}})
#'         Parameters of the method. See Details section below.
#'   * \eqn{\lambda} is a penalty for off-diagonal
#'   * \eqn{\alpha} is a penalty for on-diagonal
#' @param R,R_inv (p x p matrix, unit diagonal) initial estimation of
#'         precision matrix (it is recommended to left it default).
#' @param D (vector of length p) diagonal of initial estimation of diagonal
#'         matrix (it is recommended to left it default).
#' @param max_iter (integer) maximum number of iterations.
#' @param tolerance (double) tolerance for convergence.
#' @param tol_R,max_iter_R_inner,max_iter_R_outer
#'         Parameters passed to [ROptim()] function.
#' @param tol_D,max_iter_D_newton,max_iter_D_ls
#'         Parameters passed to [DOptim()] function.
#' @param verbose (boolian) print of loss.
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
#' @seealso [pcglassoPath()] to compute a full λ‑path of solutions
#'
#' @export
#' @examples
#' p <- 7
#' R.true <- toeplitz(c(c(1, -0.5), rep(0, p - 2)))
#' D.true <- sqrt(rchisq(p, 3))
#' S_inv.true <- diag(D.true) %*% R.true %*% diag(D.true)
#'
#' S <- solve(S_inv.true) # data
#'
#' alpha <- 4 / 20 # 4 / n, as in Carter's paper
#'
#' pcglassoFast(S, 0.11, alpha, max_iter = 15, diagonal_Newton = TRUE, verbose = TRUE)
pcglassoFast <- function(
    S, lambda, alpha,
    R = diag(dim(S)[1]), R_inv = solve(R), D = rep(1, dim(S)[1]),
    max_iter = 100, tolerance = 1e-6,
    tol_R = 1e-4,
    max_iter_R_inner = 10, max_iter_R_outer = 100,
    tol_D = 1e-4,
    max_iter_D_newton = 500, max_iter_D_ls = 100,
    diagonal_Newton = TRUE,
    verbose = FALSE) {
  stopifnot(
    is.matrix(S), nrow(S) == ncol(S),
    is.numeric(lambda), lambda >= 0,
    is.numeric(alpha),
    max_iter >= 1, tolerance > 0,
    is.matrix(R), is.matrix(R_inv)
  )

  stop_loop <- FALSE
  loss_history <- function_to_optimize(R, D, S, lambda, alpha)
  if (verbose) {
    print(paste0(round(loss_history, 4), ", starting loss"))
  }
  while (!stop_loop && (length(loss_history)/2) < max_iter) {
    # D step
    resD <- DOptim(
      A = R * S,
      d0 = D,
      tol = tol_D,
      max_newton_iter = max_iter_D_newton,
      max_ls_steps = max_iter_D_ls,
      alpha = alpha, diagonal_Newton = diagonal_Newton
    )
    proposed_loss <- function_to_optimize(R, resD$D, S, lambda, alpha)
    if (verbose) {
      print(paste0(round(proposed_loss, 4), ", after ", resD$iter, " iters of D optim"))
    }
    if (proposed_loss <= loss_history[length(loss_history)] - (tol_D * 2)) {
      rlang::warn("D optimization decreased the goal. This should not occur. We recommend to decrease the `tol_D` parameter.")
      break
    }
    D <- resD$D
    loss_history <- c(loss_history, proposed_loss)

    # R step
    resR <- ROptim(
      S = S * (D %o% D),
      R = R,
      Rinv = R_inv,
      lambda = lambda,
      tol = tol_R,
      max_inner_iter = max_iter_R_inner,
      max_outer_iter = max_iter_R_outer
    )
    proposed_loss <- function_to_optimize(resR$R, D, S, lambda, alpha)
    if (verbose) {
      print(paste0(round(proposed_loss, 4), ", after ", resR$outer.count, " iters of R optim"))
    }
    if (proposed_loss <= loss_history[length(loss_history)] - (tol_R * 2)) {
      rlang::warn("R optimization decreased the goal. This should not occur. We recommend to decrease the `tol_R` parameter.")
      break
    }
    R <- resR$R
    R_inv <- resR$Rinv
    loss_history <- c(loss_history, proposed_loss)

    # loop
    stop_loop <- (loss_history[length(loss_history)] - loss_history[length(loss_history) - 2] < tolerance)
  }

  D <- as.vector(D)
  list(
    "Sinv" = R * (D %o% D),
    "S" = R_inv * ((1 / D) %o% (1 / D)),
    "R" = R,
    "D" = D,
    "R_inv" = R_inv,
    "n_iters" = floor(length(loss_history)/2),
    "loss" = loss_history
  )
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
  sum(d * ((S * R) %*% d))
}
