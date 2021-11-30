#ifndef fillRAMMatrices_H
#define fillRAMMatrices_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat fillA(arma::mat A,
                const bool hasDiscreteDRIFTUnique,
                const Rcpp::List& discreteDRIFTUnique,
                const bool hasDiscreteTRAITUnique,
                const Rcpp::List& discreteTRAITUnique,
                const arma::mat& LAMBDA,
                const Rcpp::DataFrame& AParameterIndicators);
arma::mat fillS(arma::mat S,
                const arma::mat& T0VAR,
                const arma::mat& MANIFESTVAR,
                const bool hasDiscreteTRAITUnique,
                const arma::mat& TRAITVAR,
                const bool hasDiscreteDIFFUSIONUnique,
                const Rcpp::List& discreteDIFFUSIONUnique,
                const Rcpp::DataFrame& SParameterIndicators);
arma::colvec fillM(arma::colvec M,
                   const arma::colvec& T0MEANS,
                   const arma::colvec& MANIFESTMEANS,
                   const bool hasDiscreteCINTUnique,
                   const Rcpp::List& discreteCINTUnique,
                   const Rcpp::DataFrame& cppMParameterIndicators);

#endif
