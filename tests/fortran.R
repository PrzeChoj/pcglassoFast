# This is the same code as in Fortran

#####
# Fortran’s sign(b,a): abs(b) with the sign of a
sign_fortran <- function(b, a) {
  if (a >= 0) abs(b) else -abs(b)
}

roptim_fortran_in_R <- function(
    n, S, L, thr, maxIt, maxItLasso, msg, warm, X, W
  ) {
  Wd  <- numeric(n)
  WXj <- numeric(n)

  info      <- 0
  EPS       <- 1.1e-16
  iter      <- 0
  converged <- FALSE

  shr <- sum(abs(S) - diag(diag(abs(S))))
  shr <- thr*shr/(n-1)
  thrLasso <- shr/n
  if (thrLasso < 2*EPS)
    thrLasso <- 2*EPS

  for (i in 1:n) {
    X[, i] <- -X[, i] / X[i, i]
    X[i, i] <- 0
  }

  for (iter in seq_len(maxIt)) {
    if(iter == 2){
      #browser()
      X_iter_2 <- X
      W_iter_2 <- W
    }
    if (msg != 0) cat("iter:", iter, "\n")
    dw <- 0.0

    for (j in seq_len(n)) {
      innerIter <- 0
      WXj[] <- 0
      for (i in seq_len(n)) {
        if (X[i,j] != 0) {
          WXj <- WXj + W[,i] * X[i,j]
        }
      }

      ## --- inner LASSO loop ----
      repeat {
        dlx <- 0
        for (i in seq_len(n)) if (i != j) {
          a     <- S[i,j] - WXj[i] + W[i,i] * X[i,j]
          b     <- abs(a) - L[i,j]
          c_new <- if (b > 0) sign_fortran(b, a) / W[i,i] else 0
          # c_new == sign(a)*max(abs(a)-lambda, 0) / W[i,i]
          delta <- c_new - X[i,j]
          if (delta != 0) {
            X[i,j] <- c_new
            WXj    <- WXj + delta * W[,i]
            dlx    <- max(dlx, abs(delta))
          }
        }
        if (dlx < thrLasso || innerIter >= maxItLasso) break
        innerIter <- innerIter + 1
      }

      WXj[j] <- W[j,j]
      dw      <- max(dw, sum(abs(WXj - W[,j])))

      W[, j]    <- WXj
      W[j, ]    <- WXj
      dw      <- max(dw, abs(W[j, j] - (1 + sum(X[,j] * W[,j]))))
      W[j, j]   <- 1 + sum(X[,j] * W[,j])
    }

    if (dw < shr) {
      converged <- TRUE
      break
    }
  }

  # --- final cleanup ---
  for (i in seq_len(n)) {
    X[1:n, i] <- -X[1:n, i]
    X[i, i]   <- 1
  }

  # Symmetrize X
  X <- (X + t(X))/2

  list(
    X          = X,
    W          = W,
    iterations = iter,
    converged  = converged,
    info       = info
  )
}

ROptim_in_R <- function(
    S,
    R,
    Rinv,
    lambda,
    tol_inner = 1e-4,
    tol_outer = 1e-4,
    max_inner_iter = 10,
    max_outer_iter = 100) {
  if (is.null(Rinv)) {
    Rinv <- solve(R)
  }

  p <- dim(R)[1]
  lambda_matrix <- matrix(lambda, p, p); diag(lambda_matrix) <- 0
  ans <- roptim_fortran_in_R(
    n = p, S = S, L = lambda_matrix,
    thr = tol_outer, maxIt = max_outer_iter, maxItLasso = max_inner_iter,
    msg = 0, warm = 1, X = Rinv, W = R
  )

  if (ans$iterations == max_outer_iter + 1) {
    rlang::warn(paste0("Optimization of R matrix reached the max number of iterations (", max_outer_iter, "). Consider increasing the `max_iter_R_outer` parameter in `pcglassoFast()` function."))
  }

  list(
    R           = ans$X,
    Rinv        = ans$W,
    outer.count = ans$iterations,
    loglik      = NA
  )
}


#####
# Test of the code:
# S <- structure(c(2.146, 2.307, 2.307, 2.590), dim = c(2L, 2L))
# p <- 50; S <- structure(rnorm(p*p), dim = c(p, p)); S <- S %*% t(S)
#
# lambda <- 1
#
# p <- dim(S)[1]
# lambda_matrix <- matrix(lambda, p, p); diag(lambda_matrix) <- 0
# R <- R_inv <- diag(p)
#
#
# R_time <- Sys.time()
# R_fortran_ans <- roptim_fortran_in_R(
#   n = p, S = S, L = lambda_matrix,
#   thr = 1.0e-4, maxIt = 120, maxItLasso = 500, msg = 0,
#   warm = TRUE, X = R_inv, W = R
# )
# R_time <- Sys.time() - R_time
#
# Fortran_time <- Sys.time()
# fortran_ans <- ROptim_to_fortran(
#   S = S, lambda,
#   thr = 1.0e-4, maxIt = 120, maxItLasso = 500,
#   start = "warm", R_init = R, Rinv_init = R_inv
# )
# Fortran_time <- Sys.time() - Fortran_time
#
# R_time
# Fortran_time
#
# as.numeric(R_time) / as.numeric(Fortran_time)
