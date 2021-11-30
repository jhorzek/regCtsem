#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the matrix exponential using RcppArmadillo
using namespace Rcpp;
// [[Rcpp::export]]
arma::mat armaExpmat(const arma::mat& m){
  return(arma::expmat(m));
}
