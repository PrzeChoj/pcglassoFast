library(parallel)
library(pbmcapply)
library(space)
# simulation_experiment.R
#' Run simulation experiments comparing PC-GLasso, Glasso, and Correlation-Glasso
#'
#' This script defines functions to:
#' 1. Simulate Gaussian data from a known precision Q.
#' 2. Fit PC-GLasso, Glasso, and Correlation-Glasso on training data.
#' 3. Select tuning parameters via BIC (on full data) and cross-validation (log-likelihood on test).
#' 4. Compute loss metrics using compare_matrices().
#' 5. Repeat for multiple sample sizes and simulation replicates.
#'
#' @import pcglassoFast glasso Matrix
#' @importFrom stats cov
#' @export
run_single <- function(Q, n, split_train = 0.7,
                       alpha.grid = sort(unique(c(seq(-0.1, 0.1, length.out = 10), 0))),
                       nlambda = 100, lambda.min.ratio = 0.01) {
  p <- ncol(Q)
  ## 1) Simulate data: x ~ N(0, Q^{-1})
  L <- Cholesky(Matrix(forceSymmetric(Q), sparse = TRUE), LDL = FALSE, perm = TRUE)
  z <- matrix(rnorm(n * p), nrow = p, ncol = n)
  x <- solve(L, solve(L, z, system = "P"), system = "Lt")
  x <- Matrix::solve(L, x, system = "Pt")
  data <- as.matrix(t(x))

  ## covariances
  S_full  <- cov(data)
  # split for CV
  n_train <- floor(split_train * n)
  idx     <- sample.int(n, n_train)
  train   <- data[idx, , drop = FALSE]
  test    <- data[-idx, , drop = FALSE]
  S_train <- cov(train)
  S_test  <- cov(test)
  n_test  <- nrow(test)

  # lambda grid from training
  lam_max <- max(abs(S_train - diag(diag(S_train))))
  lam_min <- lambda.min.ratio * lam_max
  lambdas <- exp(seq(log(lam_max), log(lam_min), length.out = nlambda))

  ### PC-GLasso: BIC on full, CV on train
  best_bic <- list(bic = Inf)
  for (a in alpha.grid) {
    path <- pcglassoPath(S_full, alpha = a,
                         max.edge.fraction = 0.3,
                         lambda.min.ratio = lambda.min.ratio,
                         nlambda = nlambda)
    loss <- loss.evaluation(path, Sigma = S_full, n = n, gamma = 0)
    # find lambda index with min BIC
    i <- which.min(loss$BIC_gamma)
    if (loss$BIC_gamma[i] < best_bic$bic) {
      best_bic <- list(alpha = a,
                       lambda = path$lambda[i],
                       bic   = loss$BIC_gamma[i],
                       W     = path$W_path[[i]])
    }
  }
  Q_pc_bic <- (best_bic$W + t(best_bic$W))/2

  # 3) PC-GLasso: select by CV on S_train/S_test
  best_cv <- list(loglik = -Inf)
  for (a in alpha.grid) {
    path <- pcglassoPath(S_train, alpha = a,
                         max.edge.fraction = 0.3,
                         lambda.min.ratio = lambda.min.ratio,
                         nlambda = nlambda)
    loss <- loss.evaluation(path, Sigma = S_test, n = n_test, gamma = 0)
    j <- which.max(loss$loglik)
    if (loss$loglik[j] > best_cv$loglik) {
      best_cv <- list(alpha = a,
                      lambda = path$lambda[j],
                      loglik = loss$loglik[j],
                      W      = path$W_path[[j]])
    }
  }
  Q_pc_cv <- (best_cv$W + t(best_cv$W))/2


  ### Glasso: BIC on full, CV on train
  gl_full_path <- glasso::glassopath(S_full, rholist = lambdas,
                                     penalize.diagonal = FALSE)
  loss_gl_full <- loss.evaluation(gl_full_path$wi, Sigma = S_full, n = n, gamma = 0.)
  idx_gl_bic   <- which.min(loss_gl_full$BIC)
  Q_gl_bic     <- (gl_full_path$wi[,,idx_gl_bic] +
                     t(gl_full_path$wi[,,idx_gl_bic])) / 2

  gl_tr_path   <- glasso::glassopath(S_train, rholist = lambdas,
                                     penalize.diagonal = FALSE)
  loss_gl_cv   <- loss.evaluation(gl_tr_path$wi, Sigma = S_test, n = n_test, gamma = 0.)
  idx_gl_cv    <- which.max(loss_gl_cv$loglik)
  Q_gl_cv      <- (gl_tr_path$wi[,,idx_gl_cv] +
                     t(gl_tr_path$wi[,,idx_gl_cv])) / 2

  ### Correlation-Glasso: BIC on full, CV on train
  C_full       <- cov2cor(S_full)
  cg_full_path <- glasso::glassopath(C_full, rholist = lambdas,
                                     penalize.diagonal = FALSE)
  vars_full    <- diag(S_full)
  loss_cg_full <- loss.evaluation(cov2cor.inv(cg_full_path$wi, 1/vars_full), Sigma = S_full,
                                  n = n, gamma = 0.)
  idx_cg_bic   <- which.min(loss_cg_full$BIC)
  Theta_cg_bic <- cov2cor.inv(cg_full_path$wi[,,idx_cg_bic], 1/vars_full)
  Q_cg_bic     <- (Theta_cg_bic + t(Theta_cg_bic)) / 2

  C_tr         <- cov2cor(S_train)
  cg_tr_path   <- glasso::glassopath(C_tr, rholist = lambdas,
                                     penalize.diagonal = FALSE)
  vars_tr      <- diag(S_train)
  loss_cg_cv   <- loss.evaluation(cov2cor.inv(cg_tr_path$wi, 1/vars_tr), Sigma = S_test,
                                  n = n_test, gamma = 0.)
  idx_cg_cv    <- which.max(loss_cg_cv$loglik)
  Theta_cg_cv  <- cov2cor.inv(cg_tr_path$wi[,,idx_cg_cv], 1/vars_tr)
  Q_cg_cv      <- (Theta_cg_cv + t(Theta_cg_cv)) / 2


  ## — SPACE (BIC and CV via loss.evaluation) —
  # full-data BIC selection
  l1_full    <- 1/sqrt(n) * qnorm(1 - 1/(2 * p^2))
  scale_full <- seq(0.5, 2, length.out = 2 * nlambda)
  res_space_f <- array(0, dim = c(p, p, length(scale_full)))
  for (i in seq_along(scale_full)) {
    sp <- space.joint(as.matrix(scale(data)),
                      lam1 = n * l1_full * scale_full[i],
                      lam2 = 0, iter = 10)
    Theta <- cov2cor.inv(sp$ParCor ,1/sp$sig.fit)
    res_space_f[,,i] <- (Theta + t(Theta)) / 2
  }
  loss_space_full <- loss.evaluation(res_space_f, Sigma = S_full, n = n, gamma = 0)
  idx_space_bic   <- which.min(loss_space_full$BIC_gamma)
  Q_space_bic     <- res_space_f[,, idx_space_bic]

  # train-test CV selection
  l1_tr    <- 1/sqrt(n_train) * qnorm(1 - 1/(2 * p^2))
  scale_tr <- seq(0.5, 2, length.out = 2 * nlambda)
  res_space_t <- array(0, dim = c(p, p, length(scale_tr)))
  for (i in seq_along(scale_tr)) {
    sp <- space.joint(as.matrix(scale(train)),
                      lam1 = n_train * l1_tr * scale_tr[i],
                      lam2 = 0, iter = 10)
    Theta <- cov2cor.inv(sp$ParCor ,1/sp$sig.fit)
    res_space_t[,,i] <- (Theta + t(Theta)) / 2
  }
  loss_space_cv <- loss.evaluation(res_space_t, Sigma = S_test, n = n_test, gamma = 0)
  idx_space_cv  <- which.max(loss_space_cv$loglik)
  Q_space_cv    <- res_space_t[,, idx_space_cv]
  ## collect metrics
  res_list <- list(
    PcGL_bic = compare_matrices(Q, Q_pc_bic),
    PcGL_cv  = compare_matrices(Q, Q_pc_cv),
    GL_bic = compare_matrices(Q, Q_gl_bic),
    GL_cv  = compare_matrices(Q, Q_gl_cv),
    CorGL_bic = compare_matrices(Q, Q_cg_bic),
    CorGL_cv  = compare_matrices(Q, Q_cg_cv),
    Space_bic  = compare_matrices(Q, Q_space_bic),
    Space_cv   = compare_matrices(Q, Q_space_cv)
  )

  df <- do.call(rbind, res_list)
  df$method <- names(res_list)
  df$n      <- n
  df$alpha <-  c(best_bic$alpha, best_cv$alpha, rep(NA,length(names(res_list))-2))
  df
}



#' Run simulation experiments comparing PC-GLasso, Glasso, and Correlation-Glasso
#'
#' This script defines functions to:
#' 1. Simulate Gaussian data from a known precision matrix \\(Q\\).
#' 2. Fit PC-GLasso (over an alpha.grid), Glasso, and Correlation-Glasso.
#' 3. Select tuning parameters by:
#'    - **BIC**: evaluated on the full sample covariance (no split).
#'    - **Cross-Validation (CV)**: train/test split with test log-likelihood.
#' 4. Compute loss metrics (MSE) via `compare_matrices()`.
#' 5. Repeat over multiple sample sizes (`ns`) and replicates (`sim`) in parallel.
#'
#' @import pcglassoFast glasso Matrix parallel
#' @importFrom stats cov
#' @param Q Numeric \\(p\\times p\\) true precision matrix (symmetric, PD).
#' @param n Integer, sample size per replicate.
#' @param split_train Numeric in (0,1), fraction of data for training in CV (default 0.7).
#' @param nlambda Integer, number of lambda values on the regularization path (default 100).
#' @param lambda.min.ratio Numeric, ratio of minimum to maximum lambda (default 0.01).
#' @param alpha.grid Numeric vector of alpha values for PC-GLasso (includes 0).
#' @param ns Integer vector of sample sizes (for `run_experiments()`).
#' @param sim Integer, number of replicates per `n` (for `run_experiments()`).
#' @param mc_cores Integer, number of cores for parallel execution (default: all cores).
#' @return
run_experiments <- function(Q,
                            ns = c(200,500,1000),
                            sim = 50,
                            mc_cores = parallel::detectCores(),
                            seed=1234, ...) {
   grid <- expand.grid(n = ns, rep = seq_len(sim))
   RNGkind("L'Ecuyer-CMRG")
   # 2) Initialize the master seed
   set.seed(seed)
   results <- pbmclapply(
     seq_len(nrow(grid)),
     function(i) {
       row <- grid[i, ]
       df  <- run_single(Q, n = row$n, ...)
       cbind(n = row$n, rep = row$rep, df)
     },
     mc.cores    = mc_cores,
     mc.set.seed = TRUE
   )

  do.call(rbind, results)
}


#' Summarize & plot RMSE components and error rates across sample sizes
#'
#' @param results Data.frame from run_experiments(), with columns:
#'   n, rep,
#'   frob_norm, rmse,
#'   frob_diag, rmse_diag,
#'   frob_offdiag_zero, rmse_offdiag_zero,
#'   frob_offdiag_nonzero, rmse_offdiag_nonzero,
#'   false_pos_rate, false_neg_rate,
#'   method, alpha
#' @return A list with:
#'   - table: data.frame of mean metrics by n & method
#'   - plots: list(
#'       rmse_overall, rmse_diag,
#'       rmse_offdiag_zero, rmse_offdiag_nonzero,
#'       fp_rate, fn_rate,
#'       rmse_grid, rate_grid
#'     )
#' @import ggplot2 reshape2 patchwork
#' @export
summarize_plot_results <- function(results) {
  required <- c(
    "n","method",
    "rmse","rmse_diag","rmse_offdiag_zero","rmse_offdiag_nonzero",
    "false_pos_rate","false_neg_rate"
  )
  missing_cols <- setdiff(required, names(results))
  if (length(missing_cols))
    stop("Missing columns: ", paste(missing_cols, collapse=", "))

  # Mean over reps
  summary_df <- aggregate(
    cbind(
      rmse, rmse_diag, rmse_offdiag_zero, rmse_offdiag_nonzero,
      false_pos_rate, false_neg_rate
    ) ~ n + method,
    data = results, FUN = mean
  )

  # Table to return
  table <- summary_df

  # Individual RMSE plots
  p_rmse_overall <- ggplot(summary_df, aes(n, rmse, color=method)) +
    geom_line() + geom_point() + labs(y="RMSE Overall") + theme_minimal()
  p_rmse_diag    <- ggplot(summary_df, aes(n, rmse_diag, color=method)) +
    geom_line() + geom_point() + labs(y="RMSE Diagonal") + theme_minimal()
  p_rmse_off0    <- ggplot(summary_df, aes(n, rmse_offdiag_zero, color=method)) +
    geom_line() + geom_point() + labs(y="RMSE Off-diag (true zero)") + theme_minimal()
  p_rmse_offnz   <- ggplot(summary_df, aes(n, rmse_offdiag_nonzero, color=method)) +
    geom_line() + geom_point() + labs(y="RMSE Off-diag (true non-zero)") + theme_minimal()

  # Individual error-rate plots
  p_fp <- ggplot(summary_df, aes(n, false_pos_rate, color=method)) +
    geom_line() + geom_point() + labs(y="False Positive Rate") + theme_minimal()
  p_fn <- ggplot(summary_df, aes(n, false_neg_rate, color=method)) +
    geom_line() + geom_point() + labs(y="False Negative Rate") + theme_minimal()

  # Combined grids via patchwork::wrap_plots
  rmse_grid <- patchwork::wrap_plots(
    p_rmse_overall, p_rmse_diag,
    p_rmse_off0,    p_rmse_offnz,
    ncol = 2
  )
  rate_grid <- patchwork::wrap_plots(
    p_fp, p_fn,
    ncol = 2
  )

  list(
    table      = table,
    plots      = list(
      rmse_overall        = p_rmse_overall,
      rmse_diag           = p_rmse_diag,
      rmse_offdiag_zero   = p_rmse_off0,
      rmse_offdiag_nonzero= p_rmse_offnz,
      fp_rate             = p_fp,
      fn_rate             = p_fn,
      rmse_grid           = rmse_grid,
      rate_grid           = rate_grid
    )
  )
}
