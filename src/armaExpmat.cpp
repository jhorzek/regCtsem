#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

// Computes the matrix exponential using RcppArmadillo
using namespace Rcpp;
// [[Rcpp::export]]
arma::mat armaExpmat(arma::mat m){
  return(arma::expmat(m));
}
