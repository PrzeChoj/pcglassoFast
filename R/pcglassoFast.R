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
#' pcglassoFast(S, 0.12, alpha, max_iter = 15, diagonal_Newton = TRUE, verbose = TRUE)
pcglassoFast <- function(
    S, lambda, alpha,
    R = diag(dim(S)[1]), R_inv = solve(R), D = rep(1, dim(S)[1]),
    max_iter = 100, tolerance = 1e-3,
    tol_R = 1e-4,
    max_iter_R_inner = 500, max_iter_R_outer = 100,
    tol_D = 1e-4,
    max_iter_D_newton = 500, max_iter_D_ls = 100,
    diagonal_Newton = TRUE,
    verbose = FALSE) {
  stopifnot(
    is.matrix(S), nrow(S) == ncol(S),
    is.numeric(lambda), lambda >= 0,
    is.numeric(alpha),
    max_iter >= 1, tolerance > 0,
    is.matrix(R), is.matrix(R_inv),
    length(diagonal_Newton) == 1, is.logical(diagonal_Newton), !is.na(diagonal_Newton),
    all(diag(R) == rep(1, nrow(S)))
  )
  C <- cov2cor(S)

  stop_loop <- FALSE
  loss_history <- function_to_optimize(R, D, C, lambda, alpha)
  if (verbose) {
    print(paste0("starting loss: ", round(loss_history, 4)))
  }

  tol_D_curr <- 10
  tol_R_curr <- 1
  while (!stop_loop && (length(loss_history)/2) < max_iter) {
    # D step
    A <- C * R
    resD <- DOptim(
      A = A,
      d0 = D,
      tol = tol_D_curr,
      max_newton_iter = max_iter_D_newton,
      max_ls_steps = max_iter_D_ls,
      alpha = alpha, diagonal_Newton = diagonal_Newton
    )
    proposed_loss <- function_to_optimize(R, resD$D, C, lambda, alpha)
    if (verbose) {
      print(paste0("loss: ", round(proposed_loss, 4), ", after ", resD$iter, " iters of D optim"))
    }
    if (length(loss_history) > 1) {
      improvement_D <- proposed_loss - loss_history[length(loss_history)]
      improvement_R <- loss_history[length(loss_history)] - loss_history[length(loss_history) - 1]

      if (improvement_D < improvement_R * 2) {
        tol_D_curr <- max(tol_D, tol_D_curr / 2)
        message(paste0("decrease tol_D_curr to ", tol_D_curr))
      }
    }

    if (proposed_loss <= loss_history[length(loss_history)] - (tol_D * 2)) {
      break
    }
    D <- resD$D
    loss_history <- c(loss_history, proposed_loss)

    # R step
    R_curr <- R
    R_inv_curr <- R_inv
    iterations_done <- 0
    S_for_Fortran <- C * (D %o% D)
    repeat {
      if (verbose) {
        print("======================== Before R optimization ========================")
        print(paste0("dual norm: ", determinant(R_inv_curr)$modulus - sum(diag(R_inv_curr))))
        print(paste0("lambda = ", lambda, "; biggest err = ", max(abs(R_inv_curr - S_for_Fortran) - diag(diag(abs(R_inv_curr - S_for_Fortran))))))
      }
      resR <- ROptim(
        S = S_for_Fortran,
        R = R_curr,
        Rinv = R_inv_curr,
        lambda = lambda,
        tol = tol_R_curr,
        max_inner_iter = max_iter_R_inner,
        max_outer_iter = max_iter_R_outer - iterations_done
      )
      if (any(is.nan(resR$Rinv)) | any(is.nan(resR$R))) {
        warn("NaNs introduced in Fortran calculations")
        # TODO: Return with code error
        return()
      }
      if (verbose) {
        print("======================== After R optimization ========================")
        print(paste0("dual norm: ", determinant(resR$Rinv)$modulus - sum(diag(resR$Rinv))))
        print(paste0("lambda = ", lambda, "; biggest err = ", max(abs(resR$Rinv - S_for_Fortran) - diag(diag(abs(resR$Rinv - S_for_Fortran))))))
      }
      iterations_done <- iterations_done + resR$outer.count
      proposed_loss <- function_to_optimize(resR$R, D, C, lambda, alpha)
      if (verbose) {
        print(paste0("loss: ", round(proposed_loss, 4), ", after ", resR$outer.count, " iters of R optim"))
      }

      used_all_iterations <- (iterations_done >= max_iter_R_outer)
      tolerance_cannot_be_smaller <- (tol_R_curr <= tol_R)
      loss_is_better <- (proposed_loss > loss_history[length(loss_history)] - (tol_R * 2))
      dual_constraint_satisfied <- {
        error_matrix <- resR$Rinv - S_for_Fortran
        diag(error_matrix) <- 0 # dual constraint involve only off-diagonal
        max(abs(error_matrix)) <= lambda * 1.1 # 0.1 for approximate satisfaction
      }

      if (used_all_iterations | tolerance_cannot_be_smaller) {
        # no improvement can be made
        break
      }

      if (loss_is_better & dual_constraint_satisfied) {
        # no improvement is needed, we've converged
        break
      }

      tol_R_curr <- max(tol_R_curr / 2, tol_R)
      message(paste0("decrease tol_R_curr to ", tol_R_curr))
      R_curr <- resR$R
      R_inv_curr <- resR$Rinv
    }

    if (proposed_loss <= loss_history[length(loss_history)] - (tol_R * 2)) {
      if (verbose) {
        print("ending optimization as loss increased")
      }
      break
    }
    R <- resR$R
    R_inv <- resR$Rinv
    loss_history <- c(loss_history, proposed_loss)

    # loop
    stop_loop <- (loss_history[length(loss_history)] - loss_history[length(loss_history) - 2] < tolerance)
    if (verbose && stop_loop) {
      print(paste0("ending optimization as loss decreased less than a tolerance (", tolerance, ")"))
    }

    if (verbose && (length(loss_history)/2) >= max_iter) {
      print(paste0("ending optimization as number of iterations reached max_iter (", max_iter, ")"))
    }
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
#' \eqn{function\_to\_optimize(R, D) = log(det(R)) + (1-\alpha)log(det(D^2)) - tr(DCDR) - \lambda ||R||_1}
#' function, where \eqn{||R||_1} is only for off-diagonal elements.
function_to_optimize <- function(R, d, C, lambda, alpha) {
  my_norm_1 <- function(my_matrix) {
    diag(my_matrix) <- 0 # no diagonal elements
    sum(abs(my_matrix))
  }

  log_det_R <- determinant(R)[["modulus"]]
  attributes(log_det_R) <- NULL
  log_det_R + (1 - alpha) * 2 * sum(log(d)) - trace_DCDR(d, C, R) - lambda * my_norm_1(R)
}

#' compute tr(DCDR) where D are diagonal matrices
trace_DCDR <- function(d, C, R) {
  sum(d * ((C * R) %*% d))
}
