#' Compute irrepresentability for GLASSO and PCGLASSO
#'
#' When the irrepresentability condition holds, then
#'  for infinite data (\eqn{n = \infty}) there exists a lambda
#'  that the PCGLASSO (or GLASSO) estimator can recover the true
#'  underlying structure of the precision matrix.
#'
#' @param K True precition matrix
#'
#' @returns A non-negative number.
#'  When smaller than 1, it means the irrepresentability condition holds.
#'  When bigger than 1, it means the irrepresentability condition is violated.
#'
#' @export
#'
#' @importFrom pracma kron
#'
#' @examples
#' set.seed(123)
#' p <- 7
#' R.true <- toeplitz(c(1, -0.2, 0, 0, 0, 0, 0))
#' D.true <- sqrt(rchisq(p, 3))
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#'
#' irrepPCGLASSO(Q.true) # `0.3999698`
#'
#' D.true <- sqrt(rchisq(p, 30))
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#'
#' irrepPCGLASSO(Q.true) # `0.3999698`, the same
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
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' p <- 7
#' R.true <- toeplitz(c(1, -0.2, 0, 0, 0, 0, 0))
#' D.true <- sqrt(rchisq(p, 3))
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#'
#' irrepGLASSO(Q.true) # `2.588011`
#'
#' D.true <- sqrt(rchisq(p, 30))
#' Q.true <- diag(D.true) %*% R.true %*% diag(D.true)
#'
#' irrepGLASSO(Q.true) # `0.5570903`, different; depends on `D.true`
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
