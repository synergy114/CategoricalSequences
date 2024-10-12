#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec multiply_matrix_vector(const arma::mat& mat, const arma::vec& vec) {
  return mat * vec; // Matrix-vector multiplication
}