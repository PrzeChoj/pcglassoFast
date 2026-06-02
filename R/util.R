#' Invert cov2cor for one or many correlation matrices
#'
#' Reconstruct covariance(s) from correlation(s) plus variances.
#'
#' This function inverts the transformation performed by \code{\link[stats]{cov2cor}}.
#' Given a correlation matrix C (with unit diagonal) and a vector (or matrix) of
#' variances, it reconstructs the corresponding covariance matrix via the transformation
#' \eqn{Sigma = D C D}, where \eqn{D = diag(sqrt(variances))}.
#'
#' @param C A numeric correlation matrix (p x p) or an array of shape p x p x k.
#'   Each slice C[,,i] must be a valid correlation matrix (1 on the diagonal).
#' @param diagSigma Either:
#'   - a numeric vector of length p (variances for each variable; applied to all k slices if C is 3D),
#'   - or a numeric matrix p x k (slice-specific variances; must match dimensions of C).
#'   If \code{C} is 2D, only the length-p vector form is used.
#'
#' @details
#' The transformation formula is:
#' \deqn{Sigma = D C D}
#' where \eqn{D = diag(sqrt(v_1), ..., sqrt(v_p))} is the diagonal matrix of
#' standard deviations derived from variances.
#'
#' @return If \code{C} is 2D, a p x p covariance matrix
#'   \eqn{Sigma = D C D}, where \eqn{D = diag(sqrt(diagSigma))}.
#'   If \code{C} is 3D (p x p x k), returns a p x p x k array of covariance matrices.
#'
#' @seealso
#' \code{\link[stats]{cov2cor}} for the forward transformation.
#' \code{\link{pcglassoFast}} for the PCGLASSO solver (output includes \code{R} and \code{D}).
#'
#' @examples
#' # 2D example: recover covariance from correlation + variances
#' Sigma <- matrix(c(4, 1.2, 1.2, 9), 2, 2)
#' C <- cov2cor(Sigma)
#' Sigma2 <- cov2cor_inv(C, diag(Sigma))
#' all.equal(Sigma, Sigma2)
#'
#' # 3D example: two correlation matrices sharing the same variances
#' C3 <- array(0, c(2, 2, 2))
#' C3[, , 1] <- diag(2)
#' C3[, , 2] <- cov2cor(Sigma)
#' S3 <- cov2cor_inv(C3, diag(Sigma))
#'
#' # 3D with slice-specific variances
#' vars <- cbind(c(4, 9), c(1, 16)) # 2x2: first slice vars (4,9), second (1,16)
#' S4 <- cov2cor_inv(C3, vars)
#'
#' @export
cov2cor_inv <- function(C, diagSigma) {
  # check C
  if (!is.numeric(C) || length(dim(C)) < 2 || dim(C)[1] != dim(C)[2]) {
    stop("`C` must be a numeric p x p matrix or p x p x k array.")
  }
  p <- dim(C)[1]
  k <- if (length(dim(C)) == 3) dim(C)[3] else 1L

  # check diagSigma
  if (is.matrix(diagSigma)) {
    if (!all(dim(diagSigma) == c(p, k))) {
      stop("`diagSigma` must be length p (vector) or a p x k matrix when C has k slices.")
    }
    vars_mat <- diagSigma
  } else {
    if (!is.numeric(diagSigma) || length(diagSigma) != p) {
      stop("`diagSigma` must be a vector of length p or a p x k matrix matching C.")
    }
    # replicate the same variances across k slices
    vars_mat <- matrix(diagSigma, nrow = p, ncol = k)
  }

  # function to build one covariance from one slice
  makeSigma <- function(Cmat, vars) {
    D <- diag(sqrt(vars), nrow = p, ncol = p)
    return(D %*% Cmat %*% D)
  }

  if (k == 1L) {
    # 2D case
    Sigma <- makeSigma(C, vars_mat[, 1])
    dn <- rownames(C)
    if (!is.null(dn)) rownames(Sigma) <- colnames(Sigma) <- dn
    return(Sigma)
  } else {
    # 3D case
    out <- array(0, c(p, p, k))
    dimnames_C <- dimnames(C)
    for (i in seq_len(k)) {
      out[, , i] <- makeSigma(C[, , i], vars_mat[, i])
      if (!is.null(dimnames_C[[1]])) {
        rownames(out[, , i]) <- colnames(out[, , i]) <- dimnames_C[[1]]
      }
    }
    return(out)
  }
}

#' Compare two symmetric matrices via Frobenius norms, RMSE, and error rates
#'
#' Given a true symmetric matrix Q and its estimate Q_est (both p x p),
#' this function computes a comprehensive set of accuracy and sparsity metrics.
#' It is designed to evaluate the quality of estimated precision matrices (or other
#' symmetric matrices) relative to their ground truth.
#'
#' @param Q Numeric p x p symmetric matrix of true values.
#' @param Q_est Numeric p x p symmetric matrix of estimates; must match dims of Q.
#'
#' @details
#' The function computes 12 metrics, grouped by type:
#'
#' \strong{Overall reconstruction error:}
#' \itemize{
#'   \item `frob_norm`: Frobenius norm of the difference matrix.
#'         Formula: \eqn{||Q\_est - Q||_F = sqrt(sum((Q\_est - Q)^2))}{||Q_est - Q||_F}
#'   \item `rmse`: Root-mean-square error per element.
#'         Formula: \eqn{sqrt(sum((Q\_est - Q)^2) / p^2)}{sqrt(sum((Q_est - Q)^2) / p^2)}
#' }
#'
#' \strong{Diagonal accuracy:}
#' \itemize{
#'   \item `frob_diag`: Frobenius norm of diagonal differences only.
#'   \item `rmse_diag`: RMSE on diagonal elements.
#'         Formula: \eqn{sqrt(sum((Q\_est[i,i] - Q[i,i])^2) / p)}{sqrt(sum((Q_est[i,i] - Q[i,i])^2) / p)}
#' }
#'
#' \strong{Off-diagonal zeros (sparsity pattern):}
#' \itemize{
#'   \item `frob_offdiag_zero`: Frobenius norm on off-diagonal entries where Q[i,j] = 0 (true zeros).
#'   \item `rmse_offdiag_zero`: RMSE on true zeros (divided by count of true zeros).
#'   \item `true_non0_rate`: Proportion of true zeros correctly estimated as zero.
#'         Formula: \eqn{P(Q\_est[i,j] = 0 | Q[i,j] = 0)}{P(Q_est[i,j] = 0 | Q[i,j] = 0)}
#'         (Higher is better.)
#'   \item `false_non0_rate`: Proportion of true zeros incorrectly estimated as nonzero.
#'         Formula: \eqn{P(Q\_est[i,j] != 0 | Q[i,j] = 0)}{P(Q_est[i,j] != 0 | Q[i,j] = 0)}
#'         (Lower is better; false positives in the sparsity pattern.)
#' }
#'
#' \strong{Off-diagonal nonzeros (edge recovery):}
#' \itemize{
#'   \item `frob_offdiag_nonzero`: Frobenius norm on off-diagonal entries where Q[i,j] != 0 (true edges).
#'   \item `rmse_offdiag_nonzero`: RMSE on true nonzero entries.
#'   \item `true_0_rate`: Proportion of true nonzeros correctly estimated as nonzero.
#'         Formula: \eqn{P(Q\_est[i,j] != 0 | Q[i,j] != 0)}{P(Q_est[i,j] != 0 | Q[i,j] != 0)}
#'         (Higher is better; sensitivity/power.)
#'   \item `false_0_rate`: Proportion of true nonzeros incorrectly estimated as zero.
#'         Formula: \eqn{P(Q\_est[i,j] = 0 | Q[i,j] != 0)}{P(Q_est[i,j] = 0 | Q[i,j] != 0)}
#'         (Lower is better; false negatives in the sparsity pattern.)
#' }
#'
#' @return A one-row `data.frame` with columns: `frob_norm`, `rmse`, `frob_diag`, `rmse_diag`,
#'   `frob_offdiag_zero`, `rmse_offdiag_zero`, `frob_offdiag_nonzero`, `rmse_offdiag_nonzero`,
#'   `true_non0_rate`, `false_non0_rate`, `true_0_rate`, `false_0_rate`.
#'
#' @seealso
#' \code{\link{irrepPCGLASSO}}, \code{\link{irrepGLASSO}} for theoretical recovery guarantees.
#' \code{\link{pcglassoFast}}, \code{\link{pcglassoPath}} for PCGLASSO estimation.
#' \code{\link{evaluate_objective_path}} for model selection metrics along a regularization path.
#'
#' @examples
#' p <- 5
#' M <- matrix(rnorm(p^2), p, p)
#' Q <- (M + t(M)) / 2
#' diag(Q) <- runif(p, 0, 2)
#' Q[1, 3] <- Q[3, 1] <- 0 # enforce some zeros
#' Q <- (Q + t(Q)) / 2
#' E <- Q + matrix(rnorm(p^2, 0, 0.1), p, p)
#' E <- (E + t(E)) / 2
#' E[2, 4] <- E[4, 2] <- 0 # drop one entry
#' compare_matrices(Q, E)
#' @export
compare_matrices <- function(Q, Q_est) {
  if (!is.matrix(Q) || !is.matrix(Q_est) || any(dim(Q) != dim(Q_est))) {
    stop("Q and Q_est must be numeric matrices of the same dimensions.")
  }
  p <- nrow(Q)
  # ensure symmetry
  if (any(abs(Q - t(Q)) > sqrt(.Machine$double.eps))) stop("Q must be symmetric.")
  if (any(abs(Q_est - t(Q_est)) > sqrt(.Machine$double.eps))) stop("Q_est must be symmetric.")

  # elementwise difference
  D <- Q_est - Q

  # 1) overall Frobenius norm
  frob_norm <- sqrt(sum(D^2))
  # 2) RMSE per element
  rmse <- sqrt(sum(D^2) / (p * p))

  # diagonal errors
  dD <- diag(D)
  frob_diag <- sqrt(sum(dD^2))
  rmse_diag <- sqrt(sum(dD^2) / p)

  # off-diagonal masks
  off <- row(Q) != col(Q)
  mask_zero <- (Q == 0) & off
  mask_nz <- (Q != 0) & off
  # Frobenius and RMSE for off diag zeros
  n_zero <- sum(mask_zero)
  sq_zero <- sum(D[mask_zero]^2)
  frob_offdiag_zero <- sqrt(sq_zero)
  rmse_offdiag_zero <- if (n_zero > 0) sqrt(sq_zero / n_zero) else NA_real_
  # for nonzeros
  n_nz <- sum(mask_nz)
  sq_nz <- sum(D[mask_nz]^2)
  frob_offdiag_nonzero <- sqrt(sq_nz)
  rmse_offdiag_nonzero <- if (n_nz > 0) sqrt(sq_nz / n_nz) else NA_real_

  # false positive/negative rates
  # false positive: Q==0 & Q_est!=0 among Q==0
  n_fp <- sum(mask_zero & (Q_est != 0))
  n_tp <- sum(mask_zero & (Q_est == 0))
  false_non0_rate <- if (n_zero > 0) n_fp / n_zero else NA_real_
  true_non0_rate <- if (n_zero > 0) n_tp / n_zero else NA_real_
  # false negative: Q!=0 & Q_est==0 among Q!=0
  n_fn <- sum(mask_nz & (Q_est == 0))
  n_tn <- sum(mask_nz & (Q_est != 0))
  false_0_rate <- if (n_nz > 0) n_fn / n_nz else NA_real_
  true_0_rate <- if (n_nz > 0) n_tn / n_nz else NA_real_

  data.frame(
    frob_norm,
    rmse,
    frob_diag,
    rmse_diag,
    frob_offdiag_zero,
    rmse_offdiag_zero,
    frob_offdiag_nonzero,
    rmse_offdiag_nonzero,
    true_non0_rate,
    false_non0_rate,
    false_0_rate,
    true_0_rate,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.default_R0 <- function(S) {
  p <- nrow(S)
  eps <- max(1e-8, 1e-8 * mean(diag(S), na.rm=TRUE))
  cov2cor(solve(S + eps * diag(p)))
}
