#include <RcppArmadillo.h>
using namespace Rcpp;

// This function creates a deep copy of a cpp object. This is necessary if we want to ensure that references when extracting elements from an OpenMx model,
// we do not simply build a reference to their location

// [[Rcpp::export]]
Rcpp::NumericVector deepCopyNumericVector(Rcpp::NumericVector x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix deepCopyNumericMatrix(Rcpp::NumericMatrix x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::StringVector deepCopyStringVector(Rcpp::StringVector x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::StringMatrix deepCopyStringMatrix(Rcpp::StringMatrix x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::List deepCopyList(Rcpp::List x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::DataFrame deepCopyDataFrame(Rcpp::DataFrame x) {
  return Rcpp::clone(x);
}
