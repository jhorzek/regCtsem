#include <RcppArmadillo.h>
#include "computeRAMExpectations.h"
// [[Rcpp :: depends ( RcppArmadillo )]]

//' computeRAMExpectedCovariance
//'
//' Computes the expected covariance in a RAM model give filter matrix F, directed effects matrix A, and covariance matrix S.
//' The implementation closely follows that of Driver, C. C., Oud, J. H. L., & Voelkle, M. C. (2017). Continuous Time Structural Equation Modelling With R Package ctsem. Journal of Statistical Software, 77(5), 1–36. https://doi.org/10.18637/jss.v077.i05
//' @param F filter matrix
//' @param A direct effects matrix
//' @param S undirected effects matrix
//' @returns matrix with model implied covariance
//' @keywords internal
// [[Rcpp::export]]
arma::mat computeRAMExpectedCovariance(const arma::mat& F, const arma::mat& A, const arma::mat& S){

  arma::mat I = arma::eye(arma::size(A));
  arma::mat expectedCovariance = F*arma::inv(I-A)*S*arma::trans(arma::inv(I-A))*arma::trans(F);

  return(expectedCovariance);
}
