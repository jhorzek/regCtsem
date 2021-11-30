#include <RcppArmadillo.h>
using namespace Rcpp;

// This function creates a deep copy of a cpp object. This is necessary if we want to ensure that references when extracting elements from an OpenMx model,
// we do not simply build a reference to their location

// [[Rcpp::export]]
Rcpp::NumericVector deepCopyNumericVector(const Rcpp::NumericVector x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::NumericMatrix deepCopyNumericMatrix(const Rcpp::NumericMatrix x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::StringVector deepCopyStringVector(const Rcpp::StringVector x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::StringMatrix deepCopyStringMatrix(const Rcpp::StringMatrix x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::List deepCopyList(const Rcpp::List x) {
  return Rcpp::clone(x);
}
// [[Rcpp::export]]
Rcpp::DataFrame deepCopyDataFrame(const Rcpp::DataFrame x) {
  return Rcpp::clone(x);
}
