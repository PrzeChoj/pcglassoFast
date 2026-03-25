R_step_primal <- function(
    C, D, lambda, alpha, R_curr, R_inv_curr,
    tolerance_full_optimization, times_tol_decrease, tol_R, tol_R_curr,
    max_iter_R, max_iter_R_outer, prev_objective, verbose, iteration_number) {
  digits_to_print <- max(0, -floor(log10(tolerance_full_optimization)))

  resR <- ROptimPrimal(
    S = C * (D %o% D),
    R = R_curr, Rinv = R_inv_curr,
    lambda = lambda, tol = tol_R_curr, max_outer_iter = max_iter_R
  )

  proposed_objective <- function_to_optimize(resR$R_symetric, D, C, lambda, alpha)
  iterations_done <- length(resR$loglik)

  if (verbose >= 2) {
    print(paste0("Iteration ", iteration_number, ". Objective: ", round(proposed_objective, digits_to_print), ". Objective diff: ", round(proposed_objective - prev_objective, digits_to_print), ", after ", iterations_done, " iters of R optim"))
  }

  list(
    R = resR$R,
    R_symetric = resR$R_symetric,
    R_inv = resR$Rinv,
    proposed_objective = proposed_objective,
    iterations_done = iterations_done
  )
}


#' compute log likelihood
#' @param S (p x p) empirical covariance matrix
#' @param Q (p x p) precision matrix
#' @return loglik (double)
#' @noRd
loglik <- function(S, Q) {
  R <- chol(Q)
  return((sum(log(diag(R))) - 0.5 * sum(diag(S %*% Q))))
}

#' @useDynLib pcglassoFast, .registration = TRUE
#' @import Rcpp
ROptimPrimal <- function(
    S,
    R,
    Rinv,
    lambda,
    tol = 1e-8,
    max_outer_iter = 100) {
  stopifnot(!is.null(Rinv))
  p <- dim(R)[1]
  Q <- R
  Qinv <- Rinv
  lambda <- lambda/2
  shr <- sum(abs(S - diag(diag(S))))
  if (shr == 0) {
    # S is diagonal
    list(
      R           = diag(p),
      R_symetric  = diag(p),
      Rinv        = diag(p),
      outer.count = 0,
      loglik      = loglik(S, diag(p))
    )
  }
  tol.outer <- tol*shr/(p-1)
  tol.inner <- tol.outer / p
  max.outer.iter <- max_outer_iter
  max.inner.iter <- max.outer.iter * p

  p <- dim(Q)[1]
  loglik_old <- loglik(S, Q) - lambda * (sum(abs(Q)) - p) # diagonal of Q is 1
  loglik_vec <- rep(NA, max.outer.iter)
  loglik_vec[1] <- loglik_old
  loglik <- loglik_old
  outer.count <- 0
  crit.outer <- TRUE
  Q[1,1] <- Q[1,1] + 0
  Qinv[1,1] <- Qinv[1,1] + 0
  while (crit.outer) {
    # loop over vector
    loglik <- updateLoopCpp(S, Q, Qinv, loglik, lambda, tol.inner, max.inner.iter)

    err.outer <- loglik - loglik_old
    crit.outer <- ((err.outer > tol.outer) & (outer.count < max.outer.iter))

    outer.count <- outer.count + 1
    loglik_vec[outer.count+1] <- loglik
    loglik_old <- loglik
  }
  loglik_vec <- loglik_vec[!is.na(loglik_vec)]

  list(
    R           = Q,
    R_symetric  = (Q + t(Q)) / 2,
    Rinv        = Qinv,
    outer.count = outer.count,
    loglik      = loglik_vec
  )
}
