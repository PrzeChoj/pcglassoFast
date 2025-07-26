## Refactored version of DOptim and helpers with minor readability improvements

#' Diagonal Newton optimization for f(d)
#'
#' @param A symmetric matrix (p x p); A = R * S
#' @param d0 initial vector (defaults to ones)
#' @param tol convergence tolerance
#' @param max_newton_iter max outer Newton iterations
#' @param max_ls_steps max line-search iterations
#' @param alpha regularization parameter
#' @return list with D, iterations, and final value
DOptim <- function(
    A,
    d0 = NULL,
    tol = 1e-4,
    max_newton_iter = 100,
    max_ls_steps = 15,
    alpha = 0,
    diagonal_Newton = TRUE) {
  if (is.null(d0)) {
    d0 <- rep(1, ncol(A))
  }
  stopifnot(all(d0 > 0))

  gradient_line_search(
    d0, A, alpha,
    tol = tol,
    max_iter = max_newton_iter,
    max_ls_steps = max_ls_steps,
    diagonal_Newton = diagonal_Newton
  )
}

#' @importFrom Matrix diag
#' @importFrom rlang warn
gradient_line_search <- function(
    d, A, alpha,
    tol = 1e-4,
    max_iter = 100,
    max_ls_steps = 15,
    diagonal_Newton = TRUE) {
  iter <- 0
  prev_val <- -Inf
  curr_val <- f_d(d, A, alpha)
  diagA <- Matrix::diag(A)

  while ((curr_val - prev_val) > tol && iter < max_iter) {
    prev_val <- curr_val

    g <- gradient_d(d, A, alpha)
    step <- if (diagonal_Newton) {
      H_diag <- -2 * (diagA + (1 - alpha) / (d^2))
      -g / (H_diag - 1e-8)
    } else {
      H <- -2 * (diag((1 - alpha) / d^2, nrow(A)) + A)
      -solve(H - diag(1e-8, nrow(H)), g)
    }

    # Line Search
    step_size <- 1
    success <- FALSE
    for (bt in seq_len(max_ls_steps)) {
      d_new <- d + step_size * step
      # for negative d_new, the f_d() function returns -Infty
      val_new <- f_d(d_new, A, alpha)
      if (val_new >= prev_val) {
        d <- d_new
        curr_val <- val_new
        success <- TRUE
        break
      }
      step_size <- step_size * 0.5
    }
    if (!success) {
      rlang::warn(paste0("Line search failed to improve objective in D after ", max_ls_steps, " steps. This should not occur. Please open issue to let us know."))
      break
    }

    iter <- iter + 1
  }

  if (iter == max_iter) {
    rlang::warn(paste0("Optimization of diagonal D reached the max number of iterations (", max_iter, "). Consider increasing the `max_iter_D_newton` parameter in `pcglassoFast()` function."))
  }

  list(D = d, iter = iter, val = curr_val)
}

#' Evaluate function value at d
#'
#' This function computes the value of f(d) = 2(1-alpha)*sum(log(d)) - d^T (A) d
#' at the provided diagonal vector d.
#'
#' @noRd
#'
#' @param d A numeric vector representing the diagonal (p x 1).
#' @param A A symmetric matrix (p x p).
#' @param alpha A scalar representing the regularization parameter.
#' @return Returns the scalar value of the function f at d.
#' @examples
#' d <- as.matrix(c(1, 2))
#' A <- matrix(c(1, 2, 2, 4), 2, 2)
#' alpha <- 0.5
#' f.d(d, A, alpha) # -24.3
f_d <- function(d, A, alpha) {
  if (any(d <= 0)) {
    return(-Inf)
  }
  2 * (1 - alpha) * sum(log(d)) - sum(d * c(A %*% d))
}

#' Compute the gradient of f at d
#'
#' This function calculates the gradient of the function f(d) at the given diagonal vector d,
#' where f(d) = 2(1-alpha)*sum(log(d)) - d^T (A) d.
#'
#' @noRd
#'
#' @param d A numeric vector representing the diagonal (p x 1).
#' @param A A symmetric matrix (p x p).
#' @param alpha A scalar representing the regularization parameter.
#' @return Returns the gradient vector of f at d.
#' @examples
#' d <- c(1, 2)
#' A <- matrix(c(1, 2, 2, 4), 2, 2)
#' alpha <- 0.5
#' gradient_d(d, A, alpha)
gradient_d <- function(d, A, alpha) {
  2 * ((1 - alpha) / d - c(A %*% d))
}
