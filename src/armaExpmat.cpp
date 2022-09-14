#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

//' armaExpmat
//'
//' Computes the matrix exponential using RcppArmadillo
//' @param m symmetric matrix for which the matrix exponential should be computed
//' @returns matrix
// [[Rcpp::export]]
arma::mat armaExpmat(const arma::mat& m){
  return(arma::expmat(m));
}
