#' Compute irrepresentability condition for GLASSO and PCGLASSO
#'
#' Assesses whether the irrepresentability condition (a theoretical requirement for
#' asymptotic recovery guarantees) holds for a given precision matrix.
#'
#' When the irrepresentability condition holds (result < 1), there exists a sequence
#' of regularization parameters lambda_n that decreases to zero with sample size n,
#' such that for this sequence, the probability that the PCGLASSO (or GLASSO) estimator
#' recovers the sparsity pattern (zero/nonzero structure) of the true precision
#' matrix converges to 1 as n approaches infinity.
#'
#' @description
#' Key difference: PCGLASSO's irrepresentability condition depends only on the
#' correlation structure (matrix R), NOT on the variance scales (matrix D).
#' This is a theoretical advantage over GLASSO, which depends on all entries of K.
#' See examples below for a concrete demonstration.
#'
#' @param K True precision matrix (p x p, symmetric, positive definite).
#'
#' @returns A non-negative scalar:
#'   \item{}{< 1: Irrepresentability condition holds. Recovery guarantees apply (asymptotically).}
#'   \item{}{>= 1: Irrepresentability condition violated. No recovery guarantees (but methods may
#'     still work in practice).}
#'
#' @details
#' \strong{PCGLASSO vs. GLASSO comparison:}
#'
#' For a given correlation structure R (unit diagonal) with variance scales D:
#'
#' - PCGLASSO: irrepresentability condition depends only on R. Changing D
#'   does not change the result. This makes PCGLASSO more robust to heterogeneous variances.
#'
#' - GLASSO: irrepresentability condition depends on the full precision matrix
#'   K = DRD. Changing D changes the result, sometimes significantly.
#'
#' See the examples for a concrete demonstration.
#'
#' @seealso
#' \code{\link{pcglassoFast}}, \code{\link{pcglassoPath}} for PCGLASSO estimation.
#' \code{\link{compare_matrices}} for empirical evaluation of estimated vs. true matrices.
#'
#' @export
#'
#' @importFrom pracma kron
#'
#' @examples
#' p <- 7
#' R.true <- diag(1, p, p)
#' R.true[1, 2:p] <- -1/sqrt(p)
#' R.true[2:p, 1] <- -1/sqrt(p)
#'
#' # Scenario 1: Homogeneous variances
#' D.true <- rep(1, p)
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#' irr_pc1 <- irrepPCGLASSO(Q.true)
#' cat("PCGLASSO irrepresentability (homogeneous D):", round(irr_pc1, 4), "\n")
#' # Result: irr_pc1 < 1 (condition holds)
#'
#' # Scenario 2: Heterogeneous variances
#' set.seed(456)
#' D.true <- sqrt(rchisq(p, 3))
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#' irr_pc2 <- irrepPCGLASSO(Q.true)
#' cat("PCGLASSO irrepresentability (heterogeneous D):", round(irr_pc2, 4), "\n")
#' # Result: irr_pc2 == irr_pc1 (VARIANCE-INVARIANT!)
#'
irrepPCGLASSO <- function(K) {
  p <- ncol(K)
  D <- diag(sqrt(diag(K)))
  D_inv <- diag(1 / diag(D))
  R <- D_inv %*% K %*% D_inv
  Rinverse <- solve(R)

  non_zero_indices <- which( abs(K) >  0.00000001 )
  zero_indices     <- which( abs(K) <= 0.00000001 )
  stopifnot(length(zero_indices) > 0) # TODO: serve this edgecase
  diagonal_indices <- which( diag(p) == 1)

  P_diag <- diag(p*p); diag(P_diag)[-diagonal_indices] <- 0
  P_diag_perpendicular <- diag(p*p); diag(P_diag_perpendicular)[diagonal_indices] <- 0
  Gamma_part_1 <- P_diag_perpendicular %*% pracma::kron(Rinverse, Rinverse)
  Gamma_part_2 <- P_diag %*% (pracma::kron(Rinverse, diag(p)) + pracma::kron(diag(p), Rinverse))
  Gamma <- Gamma_part_1 + 0.5 * Gamma_part_2

  GammaSS <- Gamma[non_zero_indices, non_zero_indices]
  GammaSSplus <- solve(GammaSS)
  GammaScS <- Gamma[zero_indices, non_zero_indices]
  Q <- GammaScS %*% GammaSSplus

  my_pi <- as.vector(sign(K - diag(diag(K))))
  my_pi <- my_pi[non_zero_indices]

  max(abs(Q %*% t(t(my_pi))))
}

#' @rdname irrepPCGLASSO
#' @export
#'
#' @examples
#' # Scenario 1 GLASSO: Homogeneous variances
#' D.true <- rep(1, p)
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#' irr_gl1 <- irrepGLASSO(Q.true)
#' cat("GLASSO irrepresentability (homogeneous D):", round(irr_gl1, 4), "\n")
#'
#' # Scenario 2 GLASSO: Heterogeneous variances (same D as in irrepPCGLASSO)
#' set.seed(456)
#' D.true <- sqrt(rchisq(p, 3))
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#' irr_gl2 <- irrepGLASSO(Q.true)
#' cat("GLASSO irrepresentability (heterogeneous D):", round(irr_gl2, 4), "\n")
#' # Key difference: irr_gl2 != irr_gl1 (VARIANCE-DEPENDENT!)
#' # GLASSO is sensitive to variance scales.
irrepGLASSO <- function(K){
  non_zero_indices <- which( abs(K) >  0.00000001 )
  zero_indices     <- which( abs(K) <= 0.00000001 )
  stopifnot(length(zero_indices) > 0) # TODO: serve this edgecase
  Kinverse <- solve(K)

  Gamma <- pracma::kron(Kinverse, Kinverse)
  GammaSS <- Gamma[non_zero_indices, non_zero_indices]
  GammaSSplus <- solve(GammaSS)
  GammaScS <- Gamma[zero_indices, non_zero_indices]
  Q <- GammaScS %*% GammaSSplus

  my_pi <- as.vector(sign(K - diag(diag(K))))
  my_pi <- my_pi[non_zero_indices]

  max(abs(Q %*% t(t(my_pi))))
}
