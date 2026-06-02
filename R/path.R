#' Regularization Path for Partial Correlation Graphical LASSO
#'
#' Compute the complete regularization path of PCGLASSO over a sequence of sparsity
#' penalties. Leverages warm-starting from previous solutions to efficiently compute
#' solutions across multiple \eqn{\lambda} values while decoupling correlation (R) and
#' variance (D) estimation.
#'
#' @param S (p x p matrix) Empirical covariance matrix.
#' @param alpha (double, \eqn{\alpha < 1}) Penalty on diagonal scaling (variance estimation).
#'   Must be strictly less than 1. Can be negative; typical value: \code{4/n} (n = sample size).
#' @param lambdas (numeric vector) Sequence of sparsity penalties to fit.
#'   If \code{NULL}, an exponentially spaced decreasing grid of length \code{nlambda} from
#'   \eqn{\lambda_{\max} = \max_{i \neq j} |S_{ij}|} down to
#'   \eqn{\lambda_{\min} = \text{min\_lambda\_ratio} \cdot \lambda_{\max}} is generated
#'   automatically.
#' @param nlambda (integer) Number of lambda values when \code{lambdas = NULL}.
#'   Default: 50.
#' @param min_lambda_ratio (double in [0, 1]) Ratio \eqn{\lambda_{\min}/\lambda_{\max}}
#'   used to construct the exponential grid. Default: 0.01.
#' @param max_edge_fraction (double in (0, 1]) Maximum allowed fraction of nonzero
#'   off-diagonal entries ("edges") relative to \eqn{p(p-1)}.
#'   Once exceeded, the path computation stops early. Default: 1.0 (no early stopping).
#'   Values < 1 are only meaningful when \code{lambdas} is decreasing (the default),
#'   where solutions progress from sparse to dense.
#' @param R0,R0_inv (p x p matrix) Initial correlation matrix and its inverse
#'   (unit diagonal). Leave at defaults unless providing a warm-start solution.
#' @param D0 (numeric vector, length p) Initial diagonal scaling vector.
#'   Should be left at default unless providing a warm-start solution.
#' @param max_iter (integer) Maximum iterations of outer blockwise loop per \eqn{\lambda}.
#' @param tolerance (double > 0) Outer loop convergence threshold per \eqn{\lambda}; stops when
#'   objective improvement in a single iteration < tolerance.
#' @param solver_R (character) Optimization method for R-step: \code{"dual"} (Fortran,
#'   default) or \code{"primal"} (C++, alternative).
#' @param tol_R,max_iter_R,max_iter_R_outer,tol_D,max_iter_D_newton,max_iter_D_ls,diagonal_Newton
#'   Tuning parameters passed to \code{\link{pcglassoFast}} for each \eqn{\lambda}.
#' @param verbose (integer, 0--5) Logging detail. 0 = silent; 1 = termination reason;
#'   2 = objective; 3--4 = tolerance adaptation; 5 = detailed R-step diagnostics.
#'
#' @details
#' \strong{Objective Function:}
#'
#' At each \eqn{\lambda}, maximizes \eqn{f(R, D) = \log(\det(R)) + (1 - \alpha) \log(\det(D^2)) - \text{tr}(DSDR) - \lambda ||R||_1}
#'
#' where \eqn{||R||_1} denotes the L1-norm of off-diagonal elements only.
#'
#' \strong{Warm-Starting Strategy:}
#'
#' The algorithm solves for decreasing values of \eqn{\lambda} in sequence.
#' At each \eqn{\lambda_k}, the solution from \eqn{\lambda_{k-1}} is used as the starting point,
#'  warm-start (R, D). This significantly reduces the computational cost compared to
#' fitting each \eqn{\lambda} independently, especially when \eqn{\lambda} values are close.
#'
#' \strong{Early Stopping:}
#'
#' Path computation can terminate before reaching \code{lambdas[end]} if the number of
#' unique edges (nonzero off-diagonal entries in the upper triangle) exceeds
#' \code{max_edge_fraction * p * (p-1) / 2}.
#' The returned path includes only the values actually computed.
#'
#' \strong{Lambda Grid Generation:}
#'
#' When \code{lambdas = NULL}, the grid is constructed as:
#' \eqn{\lambda_k = \exp\left(\log(\lambda_{\max}) + \frac{k-1}{n_\lambda - 1} (\log(\lambda_{\min}) - \log(\lambda_{\max}))\right)}
#'
#' This exponential spacing is more appropriate for graphical lasso-type problems than linear spacing.
#'
#' @return List with elements:
#' \itemize{
#'   \item \code{lambdas}: Numeric vector of \eqn{\lambda} values actually used
#'     (may be shorter than input if early stopping occurs).
#'   \item \code{R_path}: Named list of estimated correlation matrices (unit diagonal).
#'   \item \code{Ri_path}: Named list of estimated inverse correlation matrices.
#'   \item \code{D_path}: Named list of estimated diagonal scaling vectors.
#'   \item \code{W_path}: Named list of covariance estimates \eqn{W = DRD}.
#'   \item \code{Wi_path}: Named list of estimated precision matrices \eqn{W^{-1}}.
#'   \item \code{objective}: Numeric vector of final objective values at each \eqn{\lambda}.
#'   \item \code{iters}: Numeric vector of iterations needed for convergence per \eqn{\lambda}.
#'   \item \code{path_optimization_time}: Total runtime of the entire path computation.
#' }
#'
#' @seealso
#' \code{\link{pcglassoFast}} for a single-lambda optimization.
#' \code{\link{evaluate_objective_path}} for post-hoc evaluation of solutions along the path.
#'
#' @export
#'
#' @importFrom utils tail
#' @importFrom stats cov2cor
#'
#' @examples
#' # Simulate hub network: variable 1 connected to all others
#' p <- 7
#' n <- 40
#' R.true <- diag(1, p, p)
#' R.true[1, 2:p] <- 1/sqrt(p)
#' R.true[2:p, 1] <- 1/sqrt(p)
#' D.true <- sqrt(rchisq(p, 3))
#' S_inv.true <- diag(D.true) %*% R.true %*% diag(D.true)
#' Z <- MASS::mvrnorm(n, rep(0, p), solve(S_inv.true))
#' S <- cov(Z)
#'
#' # Fit PCGLASSO
#' result <- pcglassoPath(S, alpha = 0, verbose = 1, nlambda = 6)
pcglassoPath <- function(
    S,
    alpha,
    lambdas = NULL,
    nlambda = 50,
    min_lambda_ratio = 0.01,
    max_edge_fraction = 1.0,
    # initial R, D (at the largest lambda)
    R0 = .default_R0(S),
    R0_inv = solve(R0),
    D0 = rep(1, nrow(S))/sqrt(diag(S)),
    # controls passed to pcglassoFast
    max_iter = 1000,
    tolerance = 1e-3,
    solver_R = c("dual", "primal"),
    tol_R = 1e-8,
    max_iter_R = 100,
    max_iter_R_outer = 500000,
    tol_D = 1e-8,
    max_iter_D_newton = 5000,
    max_iter_D_ls = 100,
    diagonal_Newton = TRUE,
    verbose = 0) {
  solver_R <- tolower(solver_R)
  solver_R <- match.arg(solver_R)
  stopifnot(
    is.matrix(S),
    nrow(S) == ncol(S),
    all(is.finite(S)),
    is.finite(mean(diag(S))),
    mean(diag(S)) > 0)
  R0 <- tryCatch(
    R0,  # call promise
    error = function(e) {
      stop("Failed to compute R0 (default_R0 / user-supplied). ",
           "Likely: S not invertible even after jitter, or invalid values. ",
           "Original error: ", conditionMessage(e),
           call. = FALSE)
    }
  )
  stopifnot(
    !is.null(R0),
    nrow(R0) == ncol(R0),
    all(diag(R0) == 1),
    nrow(R0) == nrow(S),
    length(D0) == nrow(S),
    all(is.finite(D0)), all(D0 > 0),
    is.numeric(alpha), alpha < 1,
    is.numeric(nlambda), nlambda >= 1,
    is.numeric(tolerance), tolerance > 0,
    max_iter >= 1,
    is.null(lambdas) || is.numeric(lambdas),
    min_lambda_ratio >= 0 && min_lambda_ratio <= 1,
    max_edge_fraction >= 0 && max_edge_fraction <= 1,
    length(diagonal_Newton) == 1, is.logical(diagonal_Newton), !is.na(diagonal_Newton),
    verbose %in% 0:5 # can be TRUE (1) or FALSE (0)
  )

  # TODO: Test, The `max_edge_fraction` can be smaller than 1 only for the decreasing lambdas.

  path_time_start <- Sys.time()

  p <- nrow(S)
  if (is.null(R0_inv)) R0_inv <- solve(R0)

  # build lambda‐grid if needed
  if (is.null(lambdas)) {
    lam_max <- max(abs(stats::cov2cor(S) - diag(ncol(S)))) + 0.001
    lam_min <- min_lambda_ratio * lam_max
    lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))
  }

  # prepare storage
  K <- length(lambdas)
  outW <- outWi <- outD <- outR <- outRi <- vector("list", K)
  objective_values <- iters <- numeric(K)

  # warm start
  R_curr <- R0
  Rinv_curr <- R0_inv
  D_curr <- D0

  for (k in 1:K) {
    lambda_k <- lambdas[k]
    if (verbose != 0) {
      message(paste0("Path iteration: ", k, " of ", K, "; lambda = ", round(lambda_k, 3)))
    }

    # run full blockwise optimization at this lambda
    fit <- pcglassoFast(
      S = S,
      lambda = lambda_k,
      alpha = alpha,
      R0 = R_curr,
      R0_inv = Rinv_curr,
      D0 = D_curr,
      max_iter = max_iter,
      tolerance = tolerance,
      solver_R = solver_R,
      tol_R = tol_R,
      max_iter_R = max_iter_R,
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
    objective_values[k] <- utils::tail(fit$objective, 1)
    # compute edge fraction and early stop
    edge_frac <- (sum(R_curr != 0) - p) / (p * (p - 1))
    if (edge_frac > max_edge_fraction) {
      break
    }
  }

  path_time_end <- Sys.time()
  path_time_full <- path_time_end - path_time_start
  if (verbose >= 1) {
    print(paste0("Path took ", round(path_time_full, 3), " ", attr(path_time_full, "units")))
  }

  # trim to actual length
  used <- 1:k
  lambdas <- lambdas[used]
  outR <- outR[used]
  outRi <- outRi[used]
  outD <- outD[used]
  outW <- outW[used]
  outWi <- outWi[used]
  objective_values <- objective_values[used]
  iters <- iters[used]

  names(outR) <- names(outD) <- names(outW) <- paste0("lam_", round(lambdas, 4))
  list(
    lambdas = lambdas,
    R_path = outR,
    Ri_path = outRi,
    D_path = outD,
    W_path = outW,
    Wi_path = outWi,
    objective = objective_values,
    iters = iters,
    path_optimization_time = path_time_full
  )
}


#' Objective function evaluation
#'
#' @description
#' computes the objective for the solution of the pcglassoPath, or a list or an array of
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
#'
#' @importFrom methods is
#' @export
evaluate_objective_path <- function(precision_array, Sigma, n, gamma = 0.5) {
  if ("list" %in% methods::is(precision_array)) {
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
