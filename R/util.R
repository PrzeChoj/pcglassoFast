#' Invert cov2cor for one or many correlation matrices
#'
#' Reconstruct covariance(s) from correlation(s) plus variances.
#'
#' @param C A numeric correlation matrix (p×p) or an array of shape p×p×k.
#'   Each slice C[,,i] must be a valid correlation matrix (1's on the diagonal).
#' @param diagSigma Either
#'   - a numeric vector of length p (variances for each variable),
#'   - or a numeric matrix p×k (variances for each slice).
#'   If \code{C} is 2D, only the length-p option is used.
#'
#' @return If \code{C} is 2D, a p×p covariance matrix
#'   \eqn{\Sigma = D\,C\,D}, where \eqn{D = \mathrm{diag}(\sqrt{\mathrm{diagSigma}})}.
#'   If \code{C} is 3D (p×p×k), returns a p×p×k array of \eqn{\Sigma^{(i)}}.
#'
#' @examples
#' # 2D example
#' Sigma <- matrix(c(4, 1.2, 1.2, 9), 2, 2)
#' C     <- cov2cor(Sigma)
#' Sigma2 <- cov2cor.inv(C, diag(Sigma))
#' all.equal(Sigma, Sigma2)
#'
#' # 3D example: two correlation matrices sharing the same variances
#' C3 <- array(0, c(2,2,2))
#' C3[,,1] <- diag(2)
#' C3[,,2] <- cov2cor(Sigma)
#' S3 <- cov2cor.inv(C3, diag(Sigma))
#' dim(S3)  # 2 2 2
#'
#' # 3D with slice‐specific variances
#' vars <- cbind(c(4,9), c(1,16))  # 2×2: first slice variances (4,9), second (1,16)
#' S4 <- cov2cor.inv(C3, vars)
#' all(dim(S4)==c(2,2,2))
#'
#' @export
cov2cor.inv <- function(C, diagSigma) {
  # check C
  if (!is.numeric(C) || length(dim(C)) < 2 || dim(C)[1] != dim(C)[2]) {
    stop("`C` must be a numeric p×p matrix or p×p×k array.")
  }
  p <- dim(C)[1]
  k <- if (length(dim(C))==3) dim(C)[3] else 1L

  # check diagSigma
  if (is.matrix(diagSigma)) {
    if (!all(dim(diagSigma) == c(p, k))) {
      stop("`diagSigma` must be length p (vector) or a p×k matrix when C has k slices.")
    }
    vars_mat <- diagSigma
  } else {
    if (!is.numeric(diagSigma) || length(diagSigma)!=p) {
      stop("`diagSigma` must be a vector of length p or a p×k matrix matching C.")
    }
    # replicate the same variances across k slices
    vars_mat <- matrix(diagSigma, nrow=p, ncol=k)
  }

  # function to build one covariance from one slice
  makeSigma <- function(Cmat, vars) {
    D <- diag(sqrt(vars), nrow = p, ncol = p)
    return(D %*% Cmat %*% D)
  }

  if (k == 1L) {
    # 2D case
    Sigma <- makeSigma(C, vars_mat[,1])
    dn <- rownames(C)
    if (!is.null(dn)) rownames(Sigma) <- colnames(Sigma) <- dn
    return(Sigma)
  } else {
    # 3D case
    out <- array(0, c(p, p, k))
    dimnames_C <- dimnames(C)
    for (i in seq_len(k)) {
      out[,,i] <- makeSigma(C[,,i], vars_mat[,i])
      if (!is.null(dimnames_C[[1]])) {
        rownames(out[,,i]) <- colnames(out[,,i]) <- dimnames_C[[1]]
      }
    }
    return(out)
  }
}

#' Compare two symmetric matrices via Frobenius norms, RMSE, and error rates
#'
#' Given a true symmetric matrix Q and its estimate Q.est (both p×p),
#' this function computes:
#' 1. `frob_norm`: overall Frobenius norm \(\|Q.est - Q\|_F\).
#' 2. `rmse`: root-mean-square error per element: \(\sqrt{\sum_{i,j}(Q.est - Q)^2 / p^2}\).
#' 3. `frob_diag`: Frobenius norm of diagonal differences.
#' 4. `rmse_diag`: RMSE on diagonal: \(\sqrt{\sum_{i}(Q.est_{ii} - Q_{ii})^2 / p}\).
#' 5. `frob_offdiag_zero`: Frobenius norm on off-diagonal entries where Q_{ij} == 0.
#' 6. `rmse_offdiag_zero`: RMSE on those off-diagonal zero entries: divide by number of such pairs.
#' 7. `frob_offdiag_nonzero`: Frobenius norm on off-diagonal entries where Q_{ij} != 0.
#' 8. `rmse_offdiag_nonzero`: RMSE on those off-diagonal nonzero entries.
#' 9. `false_pos_rate`: proportion of truly-zero entries with Q.est != 0.
#'10.' `false_neg_rate`: proportion of truly-nonzero entries with Q.est == 0.
#'
#' @param Q Numeric p×p symmetric matrix of true values.
#' @param Q.est Numeric p×p symmetric matrix of estimates; must match dims of Q.
#' @return A one-row `data.frame` with the above metrics.
#' @examples
#' set.seed(1)
#' p <- 5
#' M <- matrix(rnorm(p^2), p, p)
#' Q <- (M + t(M)) / 2
#' diag(Q) <- runif(p, 0, 2)
#' Q[1, 3] <- 0  # enforce some zeros
#' Q <- (Q + t(Q)) / 2
#' E <- Q + matrix(rnorm(p^2, 0, 0.1), p, p)
#' E <- (E + t(E)) / 2
#' E[2,4] <- 0     # drop one entry
#' compare_matrices(Q, E)
#' @export
compare_matrices <- function(Q, Q.est) {
  if (!is.matrix(Q) || !is.matrix(Q.est) || any(dim(Q) != dim(Q.est))) {
    stop("Q and Q.est must be numeric matrices of the same dimensions.")
  }
  p <- nrow(Q)
  # ensure symmetry
  if (any(abs(Q - t(Q)) > sqrt(.Machine$double.eps))) stop("Q must be symmetric.")
  if (any(abs(Q.est - t(Q.est)) > sqrt(.Machine$double.eps))) stop("Q.est must be symmetric.")

  # elementwise difference
  D <- Q.est - Q

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
  mask_nz   <- (Q != 0) & off
  # Frobenius and RMSE for off diag zeros
  n_zero <- sum(mask_zero)
  sq_zero <- sum(D[mask_zero]^2)
  frob_offdiag_zero    <- sqrt(sq_zero)
  rmse_offdiag_zero    <- if (n_zero>0) sqrt(sq_zero / n_zero) else NA_real_
  # for nonzeros
  n_nz <- sum(mask_nz)
  sq_nz <- sum(D[mask_nz]^2)
  frob_offdiag_nonzero <- sqrt(sq_nz)
  rmse_offdiag_nonzero <- if (n_nz>0) sqrt(sq_nz / n_nz) else NA_real_

  # false positive/negative rates
  # false positive: Q==0 & Q.est!=0 among Q==0
  n_fp <- sum(mask_zero & (Q.est != 0))
  false_non0_rate <- if (n_zero>0) n_fp / n_zero else NA_real_
  # false negative: Q!=0 & Q.est==0 among Q!=0
  n_fn <- sum(mask_nz & (Q.est == 0))
  false_0_rate <- if (n_nz>0) n_fn / n_nz else NA_real_

  data.frame(
    frob_norm,
    rmse,
    frob_diag,
    rmse_diag,
    frob_offdiag_zero,
    rmse_offdiag_zero,
    frob_offdiag_nonzero,
    rmse_offdiag_nonzero,
    false_non0_rate,
    false_0_rate,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
