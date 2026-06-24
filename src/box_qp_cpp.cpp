#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

void box_qp_dense(
    const arma::mat& Q,
    arma::vec& u,
    const arma::vec& b,
    const double rho,
    const int maxIter,
    const double tol,
    arma::vec& grad) {
  grad = 2.0 * Q * (u + b);
  double objcur = arma::dot(grad, b + u);
  double objold = objcur;

  for (int outer = 1; outer <= maxIter; ++outer) {
    for (arma::uword col = 0; col < Q.n_rows; ++col) {
      const double uold = u[col];
      const double diag = Q(col, col);
      const double bb = grad[col] - 2.0 * diag * u[col];
      const double tt = std::abs(bb) / (2.0 * diag);
      u[col] = std::copysign(std::min(tt, rho), -bb);

      if (u[col] != uold) {
        grad += 2.0 * (u[col] - uold) * Q.col(col);
      }
    }

    objcur = arma::dot(grad, b + u);
    const double dlx = std::abs(objcur - objold) / (std::abs(objold) + 1e-6);
    objold = objcur;
    if (dlx < tol || outer > maxIter - 1) break;
  }
}

void box_qp_full_without_column(
    const arma::mat& R,
    const arma::subview_col<double>& b,
    arma::vec& u,
    const arma::uword skipped,
    const double rho,
    const int maxIter,
    const double tol,
    arma::vec& grad,
    arma::vec& z) {
  z = u;
  z += b;
  z[skipped] = 0.0;
  u[skipped] = 0.0;

  grad = 2.0 * R * z;
  double objcur = arma::dot(grad, z);
  double objold = objcur;

  for (int outer = 1; outer <= maxIter; ++outer) {
    for (arma::uword col = 0; col < R.n_rows; ++col) {
      if (col == skipped) continue;

      const double uold = u[col];
      const double diag = R(col, col);
      const double bb = grad[col] - 2.0 * diag * u[col];
      const double tt = std::abs(bb) / (2.0 * diag);
      u[col] = std::copysign(std::min(tt, rho), -bb);

      if (u[col] != uold) {
        const double update = u[col] - uold;
        z[col] += update;
        grad += 2.0 * update * R.col(col);
      }
    }

    objcur = arma::dot(grad, z);
    const double dlx = std::abs(objcur - objold) / (std::abs(objold) + 1e-6);
    objold = objcur;
    if (dlx < tol || outer > maxIter - 1) break;
  }
}

void primal_dual_sweep_inplace(
    const arma::mat& S,
    arma::mat& R,
    arma::mat& U,
    const double lambda,
    const int qpMaxIter,
    const double qpTol) {
  const arma::uword p = R.n_rows;
  arma::vec u(p);
  arma::vec grad;
  arma::vec z(p);

  for (arma::uword j = 0; j < p; ++j) {
    u = U.col(j);
    box_qp_full_without_column(
      R, S.col(j), u, j, lambda, qpMaxIter, qpTol, grad, z
    );

    const double tval = arma::dot(z, 0.5 * grad);
    const double disc = 1.0 + 4.0 * tval;
    if (disc <= 0.0) {
      Rcpp::stop("PCGLASSO block update produced a non-positive discriminant.");
    }

    const double omega = (1.0 + std::sqrt(disc)) / 2.0;
    for (arma::uword row = 0; row < p; ++row) {
      if (row == j) continue;
      const double r = -0.5 * grad[row] / omega;
      R(row, j) = r;
      R(j, row) = r;
      U(row, j) = u[row];
      U(j, row) = u[row];
    }

    R(j, j) = 1.0;
    U(j, j) = 0.0;
  }
}

void summarize_outer_change(
    const arma::mat& R,
    const arma::mat& R_old,
    double& max_change,
    double& rms_change,
    double& sparse_nos) {
  const arma::uword p = R.n_rows;
  double max_abs_old = 1.0;
  double max_abs_diff = 0.0;
  double sum_sq_offdiag = 0.0;
  double sparse_count = 0.0;

  for (arma::uword col = 0; col < p; ++col) {
    for (arma::uword row = 0; row < p; ++row) {
      max_abs_old = std::max(max_abs_old, std::abs(R_old(row, col)));
      max_abs_diff = std::max(max_abs_diff, std::abs(R(row, col) - R_old(row, col)));
      if (row != col) {
        const double diff = R(row, col) - R_old(row, col);
        sum_sq_offdiag += diff * diff;
        if (std::abs(R(row, col)) <= 1e-9) {
          sparse_count += 1.0;
        }
      }
    }
  }

  max_change = max_abs_diff / max_abs_old;
  rms_change = std::sqrt(sum_sq_offdiag / static_cast<double>(p * (p - 1)));
  sparse_nos = sparse_count;
}

} // namespace

// [[Rcpp::export]]
Rcpp::List boxQpCpp(
    const arma::mat& Q,
    arma::vec u,
    const arma::vec& b,
    double rho,
    int maxIter,
    double tol) {
  arma::vec grad;
  box_qp_dense(Q, u, b, rho, maxIter, tol, grad);

  return Rcpp::List::create(
    Rcpp::Named("grad_vec") = grad,
    Rcpp::Named("u") = u
  );
}

// [[Rcpp::export]]
Rcpp::List primalDualSweepCpp(
    const arma::mat& S,
    arma::mat R,
    arma::mat U,
    double lambda,
    int qpMaxIter,
    double qpTol) {
  primal_dual_sweep_inplace(S, R, U, lambda, qpMaxIter, qpTol);

  return Rcpp::List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("U") = U
  );
}

// [[Rcpp::export]]
Rcpp::List primalDualOuterCpp(
    const arma::mat& S,
    arma::mat R,
    arma::mat U,
    double lambda,
    int outerMaxIter,
    double outerTol,
    int qpMaxIter,
    double qpTol) {
  arma::vec rel_err(outerMaxIter);
  arma::vec rms_err(outerMaxIter);
  arma::vec sparse_nos(outerMaxIter);
  arma::mat R_old(R.n_rows, R.n_cols);
  int outer_count = outerMaxIter;

  for (int outer = 0; outer < outerMaxIter; ++outer) {
    R_old = R;
    primal_dual_sweep_inplace(S, R, U, lambda, qpMaxIter, qpTol);

    double max_change = 0.0;
    double rms_change = 0.0;
    double sparse_count = 0.0;
    summarize_outer_change(R, R_old, max_change, rms_change, sparse_count);

    rel_err[outer] = max_change;
    rms_err[outer] = rms_change;
    sparse_nos[outer] = sparse_count;

    if (outer > 0 && max_change < outerTol) {
      outer_count = outer + 1;
      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("R") = R,
    Rcpp::Named("U") = U,
    Rcpp::Named("outer.count") = outer_count,
    Rcpp::Named("rel.err") = rel_err.subvec(0, outer_count - 1),
    Rcpp::Named("rms.err") = rms_err.subvec(0, outer_count - 1),
    Rcpp::Named("sparse.nos") = sparse_nos.subvec(0, outer_count - 1)
  );
}
