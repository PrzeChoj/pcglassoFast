# GLASSO algorithm of Friedman et al. 2008 with FORTRAN implementation of Sustik and Calderhead 2012.
# Ported to R by J. Clavel <julien.clavel@hotmail.fr> / <clavel@biologie.ens.fr> - 2017.

#' Fast graphical LASSO
#'
#' A faster alternative to the \code{glasso} function in the \pkg{glasso} package.
#' This implementation wraps the FORTRAN subroutine by Sustik and Calderhead (2012).
#'
#' @param S Covariance matrix (a p by p symmetric matrix).
#' @param rho Regularization parameter (a non-negative value or a p by p matrix).
#' @param thr Threshold for convergence. Default is 1e-4.
#' @param maxIt Maximum number of iterations. Default is 10,000.
#' @param start Type of start: \code{"cold"} or \code{"warm"}.
#' @param w.init Optional starting values for the estimated covariance matrix (p x p). Used only for warm starts.
#' @param wi.init Optional starting values for the inverse covariance matrix (p x p). Used only for warm starts.
#' @param trace Logical. If \code{TRUE}, prints iteration info.
#'
#' @details
#' Estimates a sparse inverse covariance matrix using a lasso (L1) penalty, following the approach of Friedman et al. (2008).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{w}{Estimated covariance matrix}
#'   \item{wi}{Estimated inverse covariance matrix}
#'   \item{errflag}{Memory allocation error flag: 0 = no error}
#'   \item{niter}{Number of iterations}
#' }
#'
#' @author Julien Clavel
#' @references
#' Friedman J., Hastie T., Tibshirani R. (2008). Sparse inverse covariance estimation with the graphical lasso. *Biostatistics*, 9, 432–441.\cr
#' Sustik M.A., Calderhead B. (2012). GLASSOFAST: An efficient GLASSO implementation. *UTCS Technical Report TR-12-29*, 1–3.
#' @seealso \code{\link[glasso]{glasso}}
#' @examples
#' set.seed(100)
#' p <- 5
#' x <- matrix(rnorm(p*p), ncol=p)
#' s <- var(x)
#' glassoFast(s, rho = 0.1)
#'
#' @keywords glasso covariance matrix regularization penalized likelihood
glassoFast <-
function(S, rho, thr=1.0e-4, maxIt=1e4, start=c("cold","warm"), w.init=NULL, wi.init=NULL, trace=FALSE){

  n=nrow(S)           # dimension of S
  if(is.matrix(rho)){
      if(length(rho)!=n*n) stop("The input matrix for \"rho\" must be of size ",n," by ",n)
      L = rho         # matrix of regularization parameters
  }else{
      L = matrix(rho,n,n) # matrix of regularization parameters
  }

  # cold or warm start
  start.type=match.arg(start)
  if(start.type=="cold"){
    is=0
    W=X=matrix(0,nrow=n,ncol=n)
  }
  if(start.type=="warm"){
    is=1
    if(is.null(w.init) | is.null(wi.init)){
      stop("Warm start specified: w.init and wi.init must be non-null")
    }
    W=w.init
    X=wi.init
  }

  Wd = WXj = numeric(n)

  msg=1*trace
  info = 0
  mode(n)="integer"
  mode(S)="double"
  mode(L)="double"
  mode(thr)="double"
  mode(maxIt)="integer"
  mode(msg)="integer"
  mode(is)="integer"
  mode(X)="double"
  mode(W)="double"
  mode(info)="integer"


  LASSO<-.Fortran("glassofast",
                 n,
                 S,
                 L,
                 thr,
                 maxIt,
                 msg,
                 is,
                 X,
                 W,
                 Wd,
                 WXj,
                 info)

  results <- list(w=LASSO[[9]], wi=LASSO[[8]], errflag=LASSO[[12]], niter=LASSO[[5]])
  return(results)
}
