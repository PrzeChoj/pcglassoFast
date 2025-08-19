#' glassoFast: a Fast Graphical LASSO
#'
#' A fast implementation of the graphical LASSO (Friedman et al., 2008),
#' using the FORTRAN algorithm by Sustik and Calderhead (2012).
#' Avoids convergence issues found in the original \code{glasso} function.
#'
#' @references
#' Friedman J., Hastie T., Tibshirani R. (2008). Sparse inverse covariance estimation with the graphical lasso. *Biostatistics*, 9(3), 432–441.
#' Sustik M.A., Calderhead B. (2012). GLASSOFAST: An efficient GLASSO implementation. *UTCS Technical Report TR-12-29*.
#'
#' @author
#' Julien Clavel \email{julien.clavel@hotmail.fr}
#'
#' @useDynLib pcglassoFast, .registration = TRUE
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom rlang warn
## usethis namespace: end
NULL
