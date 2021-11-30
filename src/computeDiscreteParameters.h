#ifndef computeDiscreteParameters_H
#define computeDiscreteParameters_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat computeDRIFTHASH(const arma::mat& DRIFTValues);
Rcpp::List computeDiscreteDRIFTs(const arma::mat& DRIFTValues, Rcpp::List discreteDRIFTUnique);
Rcpp::List computeDiscreteTRAITs(const Rcpp::List& discreteDRIFTUnique, Rcpp::List discreteTRAITUnique);
Rcpp::List computeDiscreteCINTs(Rcpp::List discreteCINTUnique, const arma::mat& DRIFTInverseValues, const Rcpp::List& discreteDRIFTUnique, const arma::colvec& CINTValues);
Rcpp::List computeDRIFTHASHExponentials(const arma::mat& DRIFTHASH, Rcpp::List DRIFTHASHExponentialUnique);
Rcpp::List computeDiscreteDIFFUSIONs(const arma::mat& DRIFTHASHInverse, const arma::mat& DIFFUSION,
                                     const Rcpp::List& DRIFTHASHExponentialUnique, Rcpp::List discreteDIFFUSIONUnique);


#endif
