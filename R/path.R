#' Pathwise blockwise optimization for pcglasso
#'
#' @param S (p × p matrix) empirical covariance matrix derived from the data.
#' @param alpha (double, \eqn{\alpha \in \mathbb{R}}) on‑diagonal penalty parameter.
#' @param lambdas (numeric vector) sequence of lambda values to fit. If \code{NULL}, an exponentially spaced grid of length \code{nlambda} from
#'   \eqn{\max_{i≠j}|S_{ij}|} down to \code{lambda.min.ratio * lambda_max} is generated.
#' @param nlambda (integer) number of lambda values when \code{lambdas = NULL}.
#' @param lambda.min.ratio (double) ratio \eqn{\lambda_{\min}/\lambda_{\max}} used to build the grid.
#' @param max.edge.fraction (double in (0,1]) maximum allowed fraction of
#'   nonzero off‑diagonal entries ("edges") in the graph.  Once this is exceeded,
#'   the lambda‑path search stops early.
#' @param R0,R0_inv (p × p matrices) initial precision matrix and its inverse; defaults to \code{diag(p)}.
#' @param D0 (numeric vector of length p) initial diagonal entries; default is \code{rep(1,p)}.
#' @param max.iter,tolerance,R.tol.inner,R.tol.outer,R.max.inner.iter,R.max.outer.iter,D.tol,D.max.starting.iter,D.max.outer.iter
#'         Parameters passed to [pcglassoFast()] for each \code{lambda}.
#'
#' @details
#' This function computes the entire regularization path of the PC‑GLASSO estimator by calling
#' \code{\link{pcglassoFast}()} at each value of \eqn{\lambda} in \code{lambdas}, warm‑starting
#' from the previous solution.  You get back a list of \eqn{(R,D)} pairs (and the resulting
#' covariance \eqn{W = D\,R\,D}) as \eqn{\lambda} decreases, stopping early if the
#' graph density exceeds \code{max.edge.fraction}.
#'
#' @md
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{lambdas}}{The grid of lambda values actually used.}
#'   \item{\code{R_path}}{List of estimated correlation matrices \eqn{R}.}
#'   \item{\code{D_path}}{List of estimated diagonal matrices \eqn{D}.}
#'   \item{\code{W_path}}{List of full covariance estimates \eqn{D\,R\,D}.}
#'   \item{\code{loss}}{Numeric vector of final objective values at each \eqn{\lambda}.}
#' }
#'
#' @seealso \code{\link{pcglassoFast}} for a single‑lambda blockwise optimizer.
#' @export
#'
#' @examples
#' p <- 7
#' R.true <- toeplitz(c(1, -0.5, rep(0, p-2)))
#' D.true <- sqrt(rchisq(p, df = 3))
#' S <- solve(diag(D.true) %*% R.true %*% diag(D.true))
#' alpha <- 4/20
#'
#' # get a small path of 20 lambdas but stop when >40% edges
#' resPath <- pcglassoPath(
#'   S, alpha,
#'   nlambda = 20,
#'   lambda.min.ratio = 0.05,
#'   max.edge.fraction = 0.4
#' )
#' resPath
pcglassoPath <- function(
    S,
    alpha,
    lambdas = NULL,
    nlambda = 50,
    lambda.min.ratio = 0.01,
    max.edge.fraction = 1.0,
    # initial R, D (at the largest lambda)
    R0 = diag(nrow(S)),
    R0_inv = NULL,
    D0 = rep(1, nrow(S)),
    # controls passed **through** to pcglassoFast
    max.iter = 100,
    tolerance = 1e-6,
    R.tol.inner = 1e-2,
    R.tol.outer = 1e-3,
    R.max.inner.iter = 10,
    R.max.outer.iter = 100,
    D.tol = 1e-4,
    D.max.starting.iter = 500,
    D.max.outer.iter = 100) {
  p <- nrow(S)
  if (is.null(R0_inv)) R0_inv <- solve(R0)

  # build lambda‐grid if needed
  if (is.null(lambdas)) {
    lam_max <- max(abs(S - diag(diag(S))))
    lam_min <- lambda.min.ratio * lam_max
    lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))
  }

  # prepare storage
  outR <- vector("list", length(lambdas))
  outD <- vector("list", length(lambdas))
  outW <- vector("list", length(lambdas))
  losses <- numeric(length(lambdas))

  # warm start
  R_curr <- R0
  Rinv_curr <- R0_inv
  D_curr <- D0

  for (k in seq_along(lambdas)) {
    lambda_k <- lambdas[k]

    # run full blockwise optimization at this lambda
    fit <- pcglassoFast(
      S = S,
      lambda = lambda_k,
      alpha = alpha,
      R = R_curr,
      R_inv = Rinv_curr,
      D = D_curr,
      max.iter = max.iter,
      tolerance = tolerance,
      R.tol.inner = R.tol.inner,
      R.tol.outer = R.tol.outer,
      R.max.inner.iter = R.max.inner.iter,
      R.max.outer.iter = R.max.outer.iter,
      D.tol = D.tol,
      D.max.starting.iter = D.max.starting.iter,
      D.max.outer.iter = D.max.outer.iter
    )

    # extract
    R_curr <- fit$R
    Rinv_curr <- fit$R_inv
    D_curr <- fit$D

    # record
    outR[[k]] <- R_curr
    outD[[k]] <- D_curr
    outW[[k]] <- diag(D_curr) %*% R_curr %*% diag(D_curr)
    losses[k] <- tail(fit$loss, 1)

    # compute edge fraction and early stop
    off_nz <- sum(R_curr != 0) - p
    edge_frac <- off_nz / (p * (p - 1))
    if (edge_frac > max.edge.fraction) {
      break
    }
  }

  # trim to actual length
  n_used <- k
  lambdas <- lambdas[1:n_used]
  outR <- outR[1:n_used]
  outD <- outD[1:n_used]
  outW <- outW[1:n_used]
  losses <- losses[1:n_used]

  names(outR) <- names(outD) <- names(outW) <- paste0("lam_", round(lambdas, 4))
  list(
    lambdas = lambdas,
    R_path  = outR,
    D_path  = outD,
    W_path  = outW,
    loss    = losses
  )
}
