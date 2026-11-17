# Shared numerical support policy for graph summaries; fitted matrices are unchanged.
.validate_edge_tol <- function(edge_tol) {
  if (!is.numeric(edge_tol) || is.complex(edge_tol) || length(edge_tol) != 1L ||
      !is.null(dim(edge_tol)) || !is.finite(edge_tol) || edge_tol < 0) {
    stop("edge_tol must be a finite, nonnegative real numeric scalar.", call. = FALSE)
  }
  invisible(edge_tol)
}

.pcg_edge_mask <- function(M, edge_tol = 1e-8) {
  .validate_edge_tol(edge_tol)
  if (!is.matrix(M) || !is.numeric(M) || is.complex(M) ||
      nrow(M) == 0L || nrow(M) != ncol(M) || any(!is.finite(M))) {
    stop("M must be a finite, nonempty, square real numeric matrix.", call. = FALSE)
  }

  if (edge_tol == 0) {
    # Preserve literal support before scaling, including tiny nonzero entries.
    mask <- M != 0
  } else {
    diagonal <- diag(M)
    if (any(diagonal <= 0)) {
      stop("Positive edge_tol requires positive diagonal entries; ",
           "use edge_tol = 0 for literal support of generic matrices.", call. = FALSE)
    }
    scale <- sqrt(diagonal)
    # Successive divisions avoid forming products of diagonal scaling factors.
    scaled <- sweep(sweep(M, 1L, scale, "/"), 2L, scale, "/")
    mask <- abs(scaled) > edge_tol
  }
  diag(mask) <- FALSE
  mask
}
