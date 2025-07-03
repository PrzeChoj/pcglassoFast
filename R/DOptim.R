#' using diagonal approximation fo Hessian
DOptim <- function(
    A,
    D0 = NULL,
    tol = 1e-4,
    max.starting.iter = 500,
    max.outer.iter = 100,
    alpha = 0) {
  if (is.null(D0)) {
    #D0 <- get_good_starting_point(A, max.starting.iter)
    D0 <- rep(1, ncol(A))
  }
  Res <- gradient.line.diagH.search(D0, A, alpha, tol = tol, max_iter = max.outer.iter, max_inner = 15)
  return(Res)
}


# TODO: Check, maybe the previous optimal D would be
# a better starting point for this algorithm than e.
#' @importFrom rlang warn
get_good_starting_point <- function(A, max_iter = 100) {
  myNorm <- function(A, x) {
    e <- rep(1, ncol(A))
    norm(diag(x) %*% A %*% diag(x) %*% e - e, type = "2")
  }

  n <- ncol(A)
  e <- rep(1, n)
  if (myNorm(A, e) < 1) {
    return(as.vector(e))
  }
  if (myNorm(A, 1 / diag(A)) < 1) {
    return(as.vector(1 / diag(A)))
  }

  iter <- 0
  b <- e - A %*% e
  t0 <- 1
  x0 <- e

  x0 <- x0 + newton(A, x0, b, t0)
  hatx <- sqrt(t0) * x0

  while (myNorm(A, hatx) >= 1 && iter < max_iter) {
    t0 <- t0 * (1 - 1 / (4 * sqrt(n) + 1))
    x0 <- x0 + newton(A, x0, b, t0)

    hatx <- sqrt(t0) * x0
    iter <- iter + 1
  }

  if (iter == max_iter) {
    # Meybe the error here would be better, but in my tests this does not work so bad as the theory suggests
    rlang::warn("no good_starting_point found")
  }

  as.vector(hatx)
}

newton <- function(A, x0, b = 0, t0 = 1) {
  rhs <- 1 / x0 - t0 * A %*% x0 - t0 * b
  lhs <- diag(x0^(-2)) + t0 * A
  y <- backsolve(lhs, rhs)

  as.vector(y)
}

#' @importFrom Matrix diag
gradient.line.diagH.search <- function(d, A, alpha, tol = 1e-4, max_iter = 100, max_inner = 15) {
  iter <- 0
  val.old <- -Inf
  val <- f.d(d, A, alpha)
  diag.A <- Matrix::diag(A)
  inner_fail <- FALSE
  while (val - val.old > tol && iter < max_iter && !inner_fail) {
    val.old <- val
    val.star <- -Inf
    grad <- gradient.d(d, A, alpha)
    H <- hessian.diag.d(d, diag.A, alpha)
    if (any(H >= 0)) {
      rlang::warn("Hessian should be negative definite. This should not occur. Please open issue to let us know.")
    }
    eps <- 1e-8
    H.safe <- H - eps
    H.inv.g <- -grad / H.safe
    stepsize <- 1
    inner.iter <- 0
    while (val.star < val.old && inner.iter < max_inner) {
      d.star <- d + stepsize * H.inv.g
      val.star <- f.d(d.star, A, alpha)
      inner.iter <- inner.iter + 1
      stepsize <- stepsize * 0.5
    }
    if (val.star >= val.old) {
      d <- d.star
      val <- val.star
    } else {
      inner_fail <- TRUE
      # Should not happen.
      rlang::warn("Fail of the inner iteration of optimization process for diagonal. This should not occur. Please open issue to let us know.")
    }

    iter <- iter + 1
  }
  return(list(D = c(d), iter = iter, val = val))
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
f.d <- function(d, A, alpha) {
  if (any(d <= 0)) {
    return(-Inf)
  }

  2 * (1 - alpha) * sum(log(d)) - sum(d * (A %*% d))
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
#' gradient.d(d, A, alpha)
gradient.d <- function(d, A, alpha) {
  2 * ((1 - alpha) / d - A %*% d)
}

hessian.diag.d <- function(d, diag.A, alpha) {
  H <- diag.A + (1 - alpha) / (d^2)
  return(-2 * H)
}
