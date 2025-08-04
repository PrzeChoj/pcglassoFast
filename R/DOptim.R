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
    step_size <- find_step_size(A, alpha, d, step, prev_val, g, max_ls_steps)
    d <- d + step_size * step

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
#' f_d(d, A, alpha) # -24.3
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

zoom_step_size <- function(
    step_lo, step_hi,
    phi, dphi,
    prev_val, dphi0,
    c1, c2,
    max_zoom_iter) {

  for (iter in seq_len(max_zoom_iter)) {
    step_j <- 0.5 * (step_lo + step_hi)

    phi_j  <- phi(step_j)
    phi_lo <- phi(step_lo)

    if (phi_j < prev_val + c1 * step_j * dphi0 || phi_j <= phi_lo) {
      step_hi <- step_j
    } else {
      dphi_j <- dphi(step_j)
      if (abs(dphi_j) <= c2 * dphi0)
        return(step_j)

      if (dphi_j * (step_hi - step_lo) >= 0)
        step_hi <- step_lo
      step_lo <- step_j
    }
  }
  warning("zoom_step_size(): reached max iterations; returning last candidate.")
  step_j
}

find_step_size <- function(
    A, alpha, d, step, prev_val, g,
    max_ls_steps,
    c1 = 1e-4, c2 = 0.9) {

  ## local closures that see d & step
  phi  <- function(s) f_d(d + s * step, A, alpha)
  dphi <- function(s) sum(gradient_d(d + s * step, A, alpha) * step)

  # directional derivative at s = 0
  dphi0 <- sum(g * step)
  if (dphi0 <= 0)
    stop("`step` must be an *ascent* direction.")

  step_prev  <- 0
  phi_prev   <- prev_val
  step_size  <- 1

  for (i in seq_len(max_ls_steps)) {
    phi_curr  <- phi(step_size)

    if (phi_curr < prev_val + c1 * step_size * dphi0 || (i > 1 && phi_curr <= phi_prev)) {
      return(zoom_step_size(step_prev, step_size, phi, dphi, prev_val, dphi0, c1, c2, max_zoom_iter = max_ls_steps))
    }

    dphi_curr <- dphi(step_size)
    if (abs(dphi_curr) <= c2 * dphi0)
      return(step_size)
    if (dphi_curr <= 0)
      return(zoom_step_size(step_size, step_prev, phi, dphi, prev_val, dphi0, c1, c2, max_zoom_iter = max_ls_steps))

    step_prev <- step_size
    phi_prev  <- phi_curr
    step_size <- step_size * 2
  }

  warning("find_step_size(): hit `max_ls_steps`; returning last tried value.")
  step_size
}
