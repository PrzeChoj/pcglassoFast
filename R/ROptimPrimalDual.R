# Note: box_qp_f is a Fortran routine that must be registered at package load
# The wrapper function call_box_qp_f (in BoxQP.R) handles parameter type coercion
# and calls the underlying Fortran routine via .Fortran()

.dp_pcg_objective <- function(R, S, lambda) {
  chol_R <- chol(R)
  logdet_R <- 2 * sum(log(diag(chol_R)))
  off_diag <- row(R) != col(R)

  -logdet_R + sum(S * R) + lambda * sum(abs(R[off_diag]))
}

ROptimPrimalDual <-
  function(S, R = NULL, U = NULL, lambda, outer.Maxiter = 100, outer.tol = 10^-5,
           qp.Maxiter = 1000, qp.tol = 10^-7, obj.seq = FALSE) {
    if (is.null(S)) stop("S is required as input.", call. = FALSE)
    if (!is.matrix(S)) S <- as.matrix(S)
    if (!is.numeric(S)) stop("S should be numeric.", call. = FALSE)

    p <- nrow(S)
    if (p != ncol(S) || !isTRUE(all.equal(S, t(S), tolerance = sqrt(.Machine$double.eps)))) {
      stop("S should be square and symmetric.", call. = FALSE)
    }

    if (missing(lambda) || !is.numeric(lambda) || length(lambda) != 1 ||
      !is.finite(lambda) || lambda < 0) {
      stop("lambda should be a finite, non-negative number.", call. = FALSE)
    }

    if (is.null(R)) {
      R <- diag(p)
    } else {
      if (!is.matrix(R)) R <- as.matrix(R)
      if (!is.numeric(R)) stop("Initial R should be numeric.", call. = FALSE)
      if (nrow(R) != p || ncol(R) != p) stop("Initial R dimensions should match S.", call. = FALSE)
      if (!isTRUE(all.equal(R, t(R), tolerance = sqrt(.Machine$double.eps)))) {
        stop("Initial R should be symmetric.", call. = FALSE)
      }
      if (!isTRUE(all.equal(diag(R), rep(1, p), tolerance = sqrt(.Machine$double.eps)))) {
        stop("Initial R must have diagonal equal to 1.", call. = FALSE)
      }
    }

    diag(R) <- 1
    if (inherits(try(chol(R), silent = TRUE), "try-error")) {
      stop("Initial R must be positive definite.", call. = FALSE)
    }

    if (is.null(U)) {
      U <- matrix(0, p, p)
    } else {
      if (!is.matrix(U)) U <- as.matrix(U)
      if (!is.numeric(U)) stop("Initial U should be numeric.", call. = FALSE)
      if (nrow(U) != p || ncol(U) != p) stop("Initial U dimensions should match S.", call. = FALSE)
      U <- (U + t(U)) / 2
    }
    U[] <- pmin(lambda, pmax(-lambda, U))
    diag(U) <- 0

    if (p == 1) {
      obj.vals <- if (obj.seq) .dp_pcg_objective(R, S, lambda) else NULL
      result <- list(R = R, U = U, rel.err = 0, sparse.nos = 0, time.counter.QP = c(0, 0, 0))
      if (obj.seq) result$obj.vals <- obj.vals
      return(result)
    }

    rel.err <- rep(0, outer.Maxiter)
    obj.vals <- rep(0, outer.Maxiter)
    sparse.nos <- rep(0, outer.Maxiter)
    time.counter.QP <- array(0, dim = c(outer.Maxiter * p, 3))

    ii <- 0
    tol <- Inf

    for (outer.iter in seq_len(outer.Maxiter)) {
      R.old <- R

      for (j in seq_len(p)) {
        ii <- ii + 1
        I <- setdiff(seq_len(p), j)

        A <- R[I, I, drop = FALSE]
        s <- as.numeric(S[I, j])

        U[I, j] <- pmin(lambda, pmax(-lambda, U[I, j]))

        t <- proc.time()
        obj <- call_box_qp_f(
          Q = A,
          u = as.numeric(U[I, j]),
          b = s,
          rho = lambda,
          Maxiter = qp.Maxiter,
          tol = qp.tol
        )
        t <- proc.time() - t
        time.counter.QP[ii, ] <- as.numeric(c(t[1], t[2], t[3]))

        u <- as.numeric(obj$u)
        v <- s + u
        g <- 0.5 * as.numeric(obj$grad_vec)

        tval <- sum(v * g)
        disc <- 1 + 4 * tval
        if (disc <= 0) {
          stop("PCGLASSO block update produced a non-positive discriminant.", call. = FALSE)
        }
        omega <- (1 + sqrt(disc)) / 2
        r <- -g / omega

        R[I, j] <- r
        R[j, I] <- r
        R[j, j] <- 1

        U[I, j] <- u
        U[j, I] <- u
        U[j, j] <- 0
      }

      diag(R) <- 1
      tol <- max(abs(R - R.old)) / max(1, max(abs(R.old)))
      rel.err[outer.iter] <- tol
      sparse.nos[outer.iter] <- sum(abs(R[row(R) != col(R)]) <= 10^-9)
      if (obj.seq) obj.vals[outer.iter] <- .dp_pcg_objective(R, S, lambda)

      if (tol < outer.tol && outer.iter > 1) break
    }

    time.counter.QP <- colSums(time.counter.QP[seq_len(ii), , drop = FALSE])

    diag(R) <- 1
    if (inherits(try(chol(R), silent = TRUE), "try-error")) {
      stop("Final R is not positive definite.", call. = FALSE)
    }

    # Symmetrize R to correct numerical asymmetries
    R_symetric <- (R + t(R)) / 2
    diag(R_symetric) <- 1

    result <- list(
      R = R,
      R_symetric = R_symetric,
      Rinv = NULL,
      dual_box = U,
      outer.count = outer.iter,
      time.counter.QP = time.counter.QP,
      rel.err = rel.err[seq_len(outer.iter)],
      sparse.nos = sparse.nos[seq_len(outer.iter)]
    )
    if (obj.seq) result$obj.vals <- obj.vals[seq_len(outer.iter)]

    return(result)
  }
