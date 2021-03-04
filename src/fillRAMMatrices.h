#ifndef fillRAMMatrices_H
#define fillRAMMatrices_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat fillA(arma::mat A,
                bool hasDiscreteDRIFTUnique, Rcpp::List discreteDRIFTUnique,
                bool hasDiscreteTRAITUnique, Rcpp::List discreteTRAITUnique,
                arma::mat LAMBDA,
                Rcpp::DataFrame AParameterIndicators);
arma::mat fillS(arma::mat S, arma::mat T0VAR, arma::mat MANIFESTVAR,
                bool hasDiscreteTRAITUnique,arma::mat TRAITVAR,
                bool hasDiscreteDIFFUSIONUnique, Rcpp::List discreteDIFFUSIONUnique,
                Rcpp::DataFrame SParameterIndicators);
arma::colvec fillM(arma::colvec M, arma::colvec T0MEANS, arma::colvec MANIFESTMEANS,
                   bool hasDiscreteCINTUnique, Rcpp::List discreteCINTUnique,
                   Rcpp::DataFrame cppMParameterIndicators);

#endif
