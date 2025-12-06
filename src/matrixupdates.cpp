#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;
// Define a custom sign function
double sign(double value) {
  return (value > 0) - (value < 0);
}

// Assuming the updateBtQb_cpp function is defined correctly
void updateBtQb_cpp(double& betaT_Q_beta, const double Q_beta_j, const double beta_j, const double Q_jj, const bool remove = true) {
  if(beta_j == 0)
    return;
  double contribution_j = 2 * beta_j * Q_beta_j - std::pow(beta_j, 2) * Q_jj;
  if(remove) betaT_Q_beta -= contribution_j;
  else betaT_Q_beta += contribution_j;

}


// [[Rcpp::export]]
void updateBeta(arma::vec& beta, const arma::mat& Qinv_ii, arma::vec& Qinv_beta, const arma::vec& s, double lambda, int p) {
  double c = as_scalar(beta.t() * Qinv_beta);
  for(int j = 0; j < (p - 1); ++j) {
    double a = Qinv_ii(j, j);
    double b = Qinv_beta(j) - Qinv_ii(j, j) * beta(j);
    updateBtQb_cpp(c, Qinv_beta(j), beta(j), a);
    double xi;
    if(c>1){
      xi = -b;
    }else{
      xi =  -b / (1 - c) - s(j);
    }
    if(std::abs(xi) > lambda || c > 1) {
      double lambda_s = sign(xi) * lambda;
      double s_l = s(j) + lambda_s;
      double a_tilde = -s_l * a;
      double b_tilde = a - 2 * s_l * b;
      double c_tilde = s_l * (1 - c) + b;
      double b_div_a = b_tilde / (2 * a_tilde);
      double beta_old = beta(j);
      beta(j) = -b_div_a + sign(a_tilde) * std::sqrt(b_div_a*b_div_a - c_tilde / a_tilde);
      Qinv_beta += Qinv_ii.col(j) * (beta(j) - beta_old);
      updateBtQb_cpp(c, Qinv_beta(j), beta(j), a, false);
    } else {
      double beta_old = beta(j);
      if(beta_old != 0) {
        beta(j) = 0;
        Qinv_beta += Qinv_ii.col(j) * (beta(j) - beta_old);
      }
    }
  }
}

// [[Rcpp::export]]
arma::mat Qinv_down_subet(const arma::mat& Qinv,const arma::mat& qtq, int i) {
  // Adjust the index for zero-based C++ indexing
  i -= 1;

  // Define indices for rows and columns to be included
  arma::uvec rows = arma::linspace<arma::uvec>(0, Qinv.n_rows - 1, Qinv.n_rows);

  // Remove the i-th index from rows and columns vectors
  rows.shed_row(i);

  // Extract the submatrix excluding the i-th row and column directly
  arma::mat Qinv_sub = Qinv.submat(rows, rows);



  // Compute the result
  Qinv_sub -= qtq;

  return Qinv_sub;
}

// Function for updating inverse matrix using Schur complement
// [[Rcpp::export]]
void Qinv_up_cpp(arma::mat &Qinv,
             const arma::mat &Qinv_ii,
             const arma::vec &Q_i,
             int i)
{
  // Convert R's 1-based index to C++ 0-based index
  int i0 = i - 1;
  int p  = Qinv.n_rows;

  // Compute Qbeta = Qinv_ii * Q_i
  arma::vec Qbeta = Qinv_ii * Q_i;

  // Compute Schur = 1 / (1 - Q_i^T * Qbeta)
  double Schur = 1.0 / (1.0 - arma::as_scalar(Q_i.t() * Qbeta));

  // SQbeta = Schur * Qbeta
  arma::vec SQbeta = Schur * Qbeta;

  // Create an index vector [0, 1, ..., p-1] and remove i0
  arma::uvec idx = arma::regspace<arma::uvec>(0, p - 1);  // 0..p-1
  idx.shed_row(i0);                                       // remove i0 from the index

  // Update Qinv[-i, -i] = Qinv_ii + (SQbeta * Qbeta^T)
  Qinv.submat(idx, idx) = Qinv_ii + SQbeta * Qbeta.t();

  // Qinv[-i, i] = -SQbeta
  Qinv.submat(idx, arma::uvec({static_cast<unsigned int>(i0)})) = -SQbeta;

  // Qinv[i, -i] = (Qinv[-i, i])^T
  Qinv.submat(arma::uvec({static_cast<unsigned int>(i0)}), idx) = (-SQbeta).t();

  // Qinv[i, i] = Schur
  Qinv(i0, i0) = Schur;

  // Qinv is now updated in-place (no return value needed)
}

// [[Rcpp::export]]
arma::mat Qinv_down_cpp(const arma::mat &Qinv, int i)
{
  // Convert from 1-based R index to 0-based C++ index
  int i0 = i - 1;
  int p = Qinv.n_rows;

  // Create an index vector for everything except i0
  // idx = [0, 1, ..., p-1], then remove i0
  arma::uvec idx = arma::regspace<arma::uvec>(0, p - 1);
  idx.shed_row(i0);

  // Extract Qinv[-i, -i]
  arma::mat Qinv_ii = Qinv.submat(idx, idx);

  // Extract Qinv[-i, i] as a (p-1)-vector
  arma::vec Qinv_i = Qinv.submat(idx, arma::uvec({static_cast<unsigned int>(i0)}));

  // Compute the scalar denominator
  double denom = Qinv(i0, i0);

  // Compute the rank-1 update term: Qinv_i * Qinv_i^T / denom
  arma::mat rank1 = (Qinv_i * Qinv_i.t()) / denom;

  // Return Qinv_ii - rank1
  return (Qinv_ii - rank1);
}

