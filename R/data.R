#’ Precision matrix from a real‐data PC‐GLASSO path
#’
#’ This dataset contains a single \(124 \times 124\) precision matrix \(Q\)
#’ obtained along the solution path of the PC‐GLASSO algorithm applied
#’ to a real‐world data set.  This matrix reflects the estimated conditional
#’ dependencies among variables at a particular regularization level.
#’
#’ @format A numeric matrix with \code{124} rows and \code{124} columns.
#’   Each diagonal entry is strictly positive, and off‐diagonal entries
#’   represent the estimated negative partial covariances (i.e.\ penalized
#’   precision entries) between variables.
#’
#’
#’ @examples
#’ \dontrun{
#’ # Load the data
#’ data(Q_simulated_pcglasso)
#’
#’ # Inspect sparsity
#’ sum(Q_simulated_pcglasso == 0) / length(Q_simulated_pcglasso)
#’
#’ # Compute implied covariance
#’ Sigma <- solve(Q_simulated_pcglasso)
#’ }
"Q_simulated_pcglasso"


#’ Precision matrix from a simulated GLASSO path
#’
#’ This dataset contains a single \(124\times 124\) precision matrix \(Q\)
#’ estimated by the graphical lasso (GLASSO) on simulated multivariate data.
#’ It illustrates the effect of \(\ell_1\) penalization on the precision
#’ matrix recovery in a controlled setting.
#’
#’ @format A numeric matrix with \code{124} rows and \code{p} columns.
#’   Diagonal entries are positive; off‐diagonals correspond to
#’   penalized precision estimates (zeros indicate conditional independence).
#’
#’
#’ @examples
#’ \dontrun{
#’ # Load the data
#’ data(Q_simulated_glasso)
#’
#’ # Visualize the sparsity pattern
#’ image(Q_simulated_glasso != 0, main = "Nonzero Pattern")
#’
#’ # Reconstruct the covariance
#’ Sigma_glasso <- solve(Q_simulated_glasso)
#’ }
"Q_simulated_glasso"
