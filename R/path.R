#' Pathwise blockwise optimization for pcglasso
#'
#' @param S (p × p matrix) empirical covariance matrix derived from the data.
#' @param alpha (double, \eqn{\alpha \in \mathbb{R}}) on‑diagonal penalty parameter.
#' @param lambdas (numeric vector) sequence of lambda values to fit. If \code{NULL}, an exponentially spaced grid of length \code{nlambda} from
#'   \eqn{\max_{i≠j}|S_{ij}|} down to \code{min_lambda_ratio * lambda_max} is generated.
#' @param nlambda (integer) number of lambda values when \code{lambdas = NULL}.
#' @param min_lambda_ratio (double) ratio \eqn{\lambda_{\min}/\lambda_{\max}} used to build the grid.
#' @param max_edge_fraction (double in (0,1]) maximum allowed fraction of
#'   nonzero off‑diagonal entries ("edges") in the graph.  Once this is exceeded,
#'   the lambda‑path search stops early.
#' @param R0,R0_inv (p × p matrices) initial precision matrix and its inverse; defaults to \code{diag(p)}.
#' @param D0 (numeric vector of length p) initial diagonal entries; default is \code{rep(1,p)}.
#' @param max_iter,tolerance,tol_R,max_iter_R_inner,max_iter_R_outer,tol_D,max_iter_D_newton,max_iter_D_ls
#'         Parameters passed to [pcglassoFast()] for each \code{lambda}.
#'
#' @details
#' This function computes the entire regularization path of the PC‑GLASSO estimator by calling
#' \code{\link{pcglassoFast}()} at each value of \eqn{\lambda} in \code{lambdas}, warm‑starting
#' from the previous solution.  You get back a list of \eqn{(R,D)} pairs (and the resulting
#' covariance \eqn{W = D\,R\,D}) as \eqn{\lambda} decreases, stopping early if the
#' graph density exceeds \code{max_edge_fraction}.
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
#' R.true <- toeplitz(c(1, -0.5, rep(0, p - 2)))
#' D.true <- sqrt(rchisq(p, df = 3))
#' S <- solve(diag(D.true) %*% R.true %*% diag(D.true))
#' alpha <- 4 / 20
#'
#' # get a small path of 20 lambdas but stop when >40% edges
#' resPath <- pcglassoPath(
#'   S, alpha,
#'   nlambda = 20,
#'   min_lambda_ratio = 0.01,
#'   max_edge_fraction = 0.4,
#'   verbose = TRUE
#' )
#' resPath
pcglassoPath <- function(
    S,
    alpha,
    lambdas = NULL,
    nlambda = 50,
    min_lambda_ratio = 0.01,
    max_edge_fraction = 1.0,
    # initial R, D (at the largest lambda)
    R0 = diag(nrow(S)),
    R0_inv = solve(R0),
    D0 = rep(1, nrow(S)),
    # controls passed to pcglassoFast
    max_iter = 100,
    tolerance = 1e-6,
    tol_R = 1e-3,
    max_iter_R_inner = 10,
    max_iter_R_outer = 100,
    tol_D = 1e-4,
    max_iter_D_newton = 500,
    max_iter_D_ls = 100,
    diagonal_Newton = TRUE,
    verbose = FALSE) {
  stopifnot(
    is.matrix(S),
    nrow(S) == ncol(S),
    nrow(R0) == ncol(R0),
    nrow(R0) == nrow(S),
    length(D0) == nrow(S),
    is.numeric(alpha),
    is.null(lambdas) || is.numeric(lambdas),
    min_lambda_ratio >= 0 && min_lambda_ratio <= 1,
    max_edge_fraction >= 0 && max_edge_fraction <= 1,
    length(diagonal_Newton) == 1, is.logical(diagonal_Newton), !is.na(diagonal_Newton)
  )

  p <- nrow(S)
  if (is.null(R0_inv)) R0_inv <- solve(R0)

  # build lambda‐grid if needed
  if (is.null(lambdas)) {
    lam_max <- max(abs(cov2cor(S) - diag(ncol(S)))) + 0.001
    lam_min <- min_lambda_ratio * lam_max
    lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))
  }

  # prepare storage
  K <- length(lambdas)
  outW <- outWi <- outD <- outR <- outRi <- vector("list", K)
  losses <- iters <- numeric(K)

  # warm start
  R_curr <- R0
  Rinv_curr <- R0_inv
  D_curr <- D0

  for (k in 1:K) {
    if (verbose) {
      print(paste0("Path iteration: ", k))
    }
    lambda_k <- lambdas[k]

    # run full blockwise optimization at this lambda
    fit <- pcglassoFast(
      S = S,
      lambda = lambda_k,
      alpha = alpha,
      R = R_curr,
      R_inv = Rinv_curr,
      D = D_curr,
      max_iter = max_iter,
      tolerance = tolerance,
      tol_R = tol_R,
      max_iter_R_inner = max_iter_R_inner,
      max_iter_R_outer = max_iter_R_outer,
      tol_D = tol_D,
      max_iter_D_newton = max_iter_D_newton,
      max_iter_D_ls = max_iter_D_ls,
      diagonal_Newton = diagonal_Newton,
      verbose = verbose
    )

    # extract
    R_curr <- fit$R
    Rinv_curr <- fit$R_inv
    D_curr <- fit$D

    # record
    outR[[k]] <- R_curr
    outRi[[k]] <- Rinv_curr
    outD[[k]] <- D_curr
    outW[[k]] <- R_curr * (D_curr %o% D_curr)
    outWi[[k]] <- Rinv_curr * ((1/D_curr) %o% (1/D_curr))
    iters[k] <- fit$n_iters
    losses[k] <- tail(fit$loss, 1)

    # compute edge fraction and early stop
    edge_frac <- (sum(R_curr != 0) - p) / (p * (p - 1))
    if (edge_frac > max_edge_fraction) {
      break
    }
  }

  # trim to actual length
  used <- 1:k
  lambdas <- lambdas[used]
  outR <- outR[used]
  outRi <- outRi[used]
  outD <- outD[used]
  outW <- outW[used]
  outWi <- outWi[used]
  losses <- losses[used]
  iters <- iters[used]

  names(outR) <- names(outD) <- names(outW) <- paste0("lam_", round(lambdas, 4))
  list(
    lambdas = lambdas,
    R_path = outR,
    Ri_path = outRi,
    D_path = outD,
    W_path = outW,
    Wi_path = outWi,
    loss = losses,
    iters = iters
  )
}


#' Loss evaluation
#'
#' @description
#' computes the loss for the solution of the pcglassoPath, or a list or an array of
#' matrices
#' the BIC_gamma is from Extended Bayesian Information Criteria for Gaussian
#' Graphical Models
#' @param precision_array (p x p x k)
#' @param Sigma (p x p) the covariance matrix
#' @param n     (int)   number of observations
#' @param gamma (1 x 1) scaling parameter for \eqn{BIC_{gamma}}
#' @return list (k x 1) $loglik the average loglikelihood
#'                      $forbenious
#'                      $n_param
#'                      $BIC_gamma the average BIC_gamma
#' @export
evaluate_loss_path <- function(precision_array, Sigma, n, gamma = 0.5) {
  if ("list" %in% is(precision_array)) {
    W <- precision_array$W_path
    precision_array <- simplify2array(W)
  }
  k <- dim(precision_array)[3]
  p <- dim(precision_array)[1]
  loglik <- forbenious <- n_param <- BIC_gamma <- nEdges <- rep(0, k)

  for (i in 1:k) {
    P <- precision_array[, , i]
    L.P <- chol(P)
    loglik[i] <- n * (sum(log(diag(L.P))) - 0.5 * sum(diag(P %*% Sigma)) - 0.5 * p * log(2 * pi))
    forbenious[i] <- sum((solve(P) - Sigma)^2)
    n_param[i] <- (p * p - sum(P == 0) - p) / 2 + p
    nEdges[i] <- (sum(P != 0) - p) / 2
    BIC_gamma[i] <- -2 * loglik[i] + nEdges[i] * (log(n) + 4 * gamma * log(p))
  }

  return(list(
    loglik = loglik,
    forbenious = forbenious,
    n_param = n_param,
    BIC_gamma = BIC_gamma,
    nEdges = nEdges
  ))
}
