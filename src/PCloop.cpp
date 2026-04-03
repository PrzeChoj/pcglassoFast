#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
arma::mat Qinv_down_cpp(const arma::mat &Qinv, int i);
void Qinv_up_cpp(arma::mat &Qinv,
                 const arma::mat &Qinv_ii,
                 const arma::vec &Q_i,
                 int i);
void updateBeta(arma::vec& beta, const arma::mat& Qinv_ii, arma::vec& Qinv_beta, const arma::vec& s, double lambda, int p) ;

// [[Rcpp::export]]
double updateLoopCpp(
    const arma::mat S,     // (p x p) matrix
    arma::mat & Q,           // (p x p)
    arma::mat & Qinv,        // (p x p)
    double  loglik_s,
    double lambda,
    double tol_inner,
    int max_inner_iter
) {
  // Dimensions
  int p = Q.n_rows;

  // For each i in 1..p (1-based in R)
  for(int i = 1; i <= p; i++) {
    // Qinv_ii = downdated submatrix of Qinv
    arma::mat Qinv_ii = Qinv_down_cpp(Qinv, i);

    int i0 = i - 1;  // convert to 0-based index

    // Build an index of all except i0
    arma::uvec idx = arma::regspace<arma::uvec>(0, p - 1);
    idx.shed_row(i0);

    // s is the row i0 of S, ignoring column i0
    //  -> a rowvector of length (p-1)
    arma::rowvec s_r = S.row(i0);
    arma::vec s = s_r.cols(idx).t();

    // beta = Q[i, -i], i.e. row i0 of Q, ignoring column i0
    // Q_i0 is 1 × p
    arma::rowvec Q_i0 = Q.row(i0);
    // We want a column vector of length (p-1), so transpose the 1 × (p-1) subview
    arma::vec beta = Q_i0.cols(idx).t();

    // If needed, store beta_prev (not strictly used here, but consistent with your loop)
    arma::vec beta_prev = beta;

    int inner_count = 0;
    bool critin = true;

    while(critin) {
      // Qinv_beta = Qinv_ii %*% beta
      arma::vec Qinv_beta = Qinv_ii * beta;
      double loglik_old_inner = loglik_s;
      // --- Remove terms from loglik_s
      double val = arma::as_scalar(beta.t() * Qinv_beta); // (1×k)(k×1)
      loglik_s -= 0.5 * std::log(1.0 - val);          // remove determinant term
      loglik_s += arma::as_scalar(s.t() * beta);         // remove tr(QS) term
      loglik_s += 2.0 * lambda * arma::sum(arma::abs(beta)); // remove penalty

      updateBeta(beta, Qinv_ii, Qinv_beta, s, 2.0 * lambda, p);

      // Recompute Qinv_beta for the u
      val = arma::as_scalar(beta.t() * Qinv_beta);

      // --- Add terms back in
      loglik_s += 0.5 * std::log(1.0 - val);
      loglik_s -= arma::as_scalar(s.t() * beta);
      loglik_s -= 2.0 * lambda * arma::sum(arma::abs(beta));

      double err_inner = loglik_s - loglik_old_inner;
      inner_count++;

      // Inner loop convergence
      critin = (err_inner > tol_inner) && (inner_count < max_inner_iter);
      beta_prev = beta;
    }

    // Update Q:
    // Q[i, -i] = beta (which is (p-1)×1, so we need to transpose to store in a row)
    arma::rowvec row_i = Q.row(i0);
    row_i.cols(idx) = beta.t();
    Q.row(i0) = row_i;

    // Q[-i, i] = beta
    arma::colvec col_i = Q.col(i0);
    col_i.rows(idx) = beta;
    Q.col(i0) = col_i;

    // Update Qinv
    Qinv_up_cpp(Qinv, Qinv_ii, beta, i);
  }

  return(loglik_s);
}


