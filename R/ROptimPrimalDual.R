# The column sweep is implemented in primalDualSweepCpp(); the R wrapper keeps
# validation, stopping rules, and optional objective tracing in one place.

.dp_pcg_objective <- function(R, S, lambda) {
  chol_R <- chol(R)
  logdet_R <- 2 * sum(log(diag(chol_R)))
  off_diag <- row(R) != col(R)

  -logdet_R + sum(S * R) + lambda * sum(abs(R[off_diag]))
}

ROptimPrimalDual <-
  function(S, R = NULL, U = NULL, lambda, outer.Maxiter = 100, outer.tol = 10^-5,
           qp.Maxiter = 1000, qp.tol = 2*.Machine$double.eps, obj.seq = FALSE,
           stopping.rule = c("hybrid", "max"),
           track.qp.time = FALSE) {
    stopping.rule <- match.arg(stopping.rule)
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
      result <- list(
        R = R, U = U, rel.err = 0, rms.err = 0, sparse.nos = 0,
        time.counter.QP = c(user = 0, system = 0, elapsed = 0)
      )
      if (obj.seq) result$obj.vals <- obj.vals
      return(result)
    }

    if (!obj.seq && stopping.rule == "max" && !track.qp.time) {
      opt <- primalDualOuterCpp(
        S = S,
        R = R,
        U = U,
        lambda = lambda,
        outerMaxIter = outer.Maxiter,
        outerTol = outer.tol,
        qpMaxIter = qp.Maxiter,
        qpTol = qp.tol
      )

      R <- opt$R
      U <- opt$U
      diag(R) <- 1
      if (inherits(try(chol(R), silent = TRUE), "try-error")) {
        stop("Final R is not positive definite.", call. = FALSE)
      }

      R_symetric <- (R + t(R)) / 2
      diag(R_symetric) <- 1

      return(list(
        R = R,
        R_symetric = R_symetric,
        Rinv = NULL,
        dual_box = U,
        outer.count = opt$outer.count,
        time.counter.QP = c(user = 0, system = 0, elapsed = 0),
        rel.err = opt$rel.err,
        rms.err = opt$rms.err,
        sparse.nos = opt$sparse.nos
      ))
    }

    max.err <- numeric(0)
    rms.err <- numeric(0)
    obj.vals <- numeric(0)
    sparse.nos <- numeric(0)
    time.counter.QP <- c(user = 0, system = 0, elapsed = 0)
    off_diag <- row(R) != col(R)
    current_obj <- if (obj.seq || stopping.rule == "hybrid") {
      .dp_pcg_objective(R, S, lambda)
    } else {
      NA_real_
    }

    for (outer.iter in seq_len(outer.Maxiter)) {
      R.old <- R
      old_obj <- current_obj

      if (track.qp.time) {
        t <- proc.time()
        sweep <- primalDualSweepCpp(S, R, U, lambda, qp.Maxiter, qp.tol)
        time.counter.QP <- time.counter.QP + as.numeric((proc.time() - t)[1:3])
      } else {
        sweep <- primalDualSweepCpp(S, R, U, lambda, qp.Maxiter, qp.tol)
      }
      R <- sweep$R
      U <- sweep$U

      diag(R) <- 1
      diff_R <- R - R.old
      max_change <- max(abs(diff_R)) / max(1, max(abs(R.old)))
      rms_change <- sqrt(mean(diff_R[off_diag]^2))
      max.err[outer.iter] <- max_change
      rms.err[outer.iter] <- rms_change
      sparse.nos[outer.iter] <- sum(abs(R[off_diag]) <= 10^-9)

      if (obj.seq || stopping.rule == "hybrid") {
        current_obj <- .dp_pcg_objective(R, S, lambda)
      }
      if (obj.seq) obj.vals[outer.iter] <- current_obj

      if (outer.iter > 1) {
        converged <- max_change < outer.tol
        if (stopping.rule == "hybrid") {
          objective_change <- abs(current_obj - old_obj) / max(1, abs(old_obj))
          objective_not_worse <- current_obj <= old_obj + max(outer.tol, 1e-12)
          max_change_controlled <- max_change < 10 * outer.tol
          converged <- converged ||
            (rms_change < outer.tol && objective_change < outer.tol &&
               objective_not_worse && max_change_controlled)
        }
        if (converged) break
      }
    }

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
      rel.err = max.err,
      rms.err = rms.err,
      sparse.nos = sparse.nos
    )
    if (obj.seq) result$obj.vals <- obj.vals[seq_len(outer.iter)]

    return(result)
  }
