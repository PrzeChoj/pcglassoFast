# Wrapper for the Fortran routine box_qp_f
# The Fortran routine must be registered as part of package initialization

call_box_qp_f <-
  function(Q, u, b, rho, Maxiter = 10^3, tol = 10^-4) {
    Q <- as.matrix(Q)
    u <- as.numeric(u)
    b <- as.numeric(b)
    qq <- nrow(Q)

    if (ncol(Q) != qq) stop("Q should be a square matrix")
    if (length(u) != qq || length(b) != qq) stop("u and b should match Q dimensions")

    mode(Q) <- "double"
    mode(u) <- "double"
    mode(b) <- "double"
    mode(Maxiter) <- "integer"
    mode(tol) <- "double"
    mode(rho) <- "double"
    mode(qq) <- "integer"

    junk <- .Fortran("box_qp_f", Q, uu = u, b, rho, Maxiter, tol, qq, grad_vec = double(qq))

    return(list(grad_vec = junk$grad_vec, u = junk$uu))
  }
