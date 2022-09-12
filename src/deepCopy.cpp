#include <RcppArmadillo.h>

// This function creates a deep copy of a cpp object. This is necessary if we want to ensure that references when extracting elements from an OpenMx model,
// we do not simply build a reference to their location

//' deepCopyNumericVector
//'
//' deep copy of Rcpp data
//' @param x Rcpp NumericVector
//' @returns copy of Rcpp NumericVector
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector deepCopyNumericVector(const Rcpp::NumericVector x) {
  return Rcpp::clone(x);
}

//' deepCopyNumericMatrix
//'
//' deep copy of Rcpp data
//' @param x Rcpp NumericMatrix
//' @returns copy of Rcpp NumericMatrix
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix deepCopyNumericMatrix(const Rcpp::NumericMatrix x) {
  return Rcpp::clone(x);
}

//' deepCopyStringVector
//'
//' deep copy of Rcpp data
//' @param x Rcpp StringVector
//' @returns copy of Rcpp StringVector
//' @keywords internal
// [[Rcpp::export]]
Rcpp::StringVector deepCopyStringVector(const Rcpp::StringVector x) {
  return Rcpp::clone(x);
}

//' deepCopyStringMatrix
//'
//' deep copy of Rcpp data
//' @param x Rcpp StringMatrix
//' @returns copy of Rcpp StringMatrix
//' @keywords internal
// [[Rcpp::export]]
Rcpp::StringMatrix deepCopyStringMatrix(const Rcpp::StringMatrix x) {
  return Rcpp::clone(x);
}

//' deepCopyList
//'
//' deep copy of Rcpp data
//' @param x Rcpp List
//' @returns copy of Rcpp List
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List deepCopyList(const Rcpp::List x) {
  return Rcpp::clone(x);
}

//' deepCopyDataFrame
//'
//' deep copy of Rcpp data
//' @param x Rcpp DataFrame
//' @returns copy of Rcpp DataFrame
//' @keywords internal
// [[Rcpp::export]]
Rcpp::DataFrame deepCopyDataFrame(const Rcpp::DataFrame x) {
  return Rcpp::clone(x);
}
