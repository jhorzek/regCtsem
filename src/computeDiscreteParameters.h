#ifndef computeDiscreteParameters_H
#define computeDiscreteParameters_H

#include <RcppArmadillo.h>
// [[Rcpp :: depends ( RcppArmadillo )]]

arma::mat computeDRIFTHASH(arma::mat DRIFTValues);
Rcpp::List computeDiscreteDRIFTs(arma::mat DRIFTValues, Rcpp::List discreteDRIFTUnique);
Rcpp::List computeDiscreteTRAITs(Rcpp::List discreteDRIFTUnique, Rcpp::List discreteTRAITUnique);
Rcpp::List computeDiscreteCINTs(Rcpp::List discreteCINTUnique,arma::mat DRIFTInverseValues, Rcpp::List discreteDRIFTUnique, arma::colvec CINTValues);
Rcpp::List computeDRIFTHASHExponentials(arma::mat DRIFTHASH, Rcpp::List DRIFTHASHExponentialUnique);
Rcpp::List computeDiscreteDIFFUSIONs(arma::mat DRIFTHASHInverse, arma::mat DIFFUSION,
                                     Rcpp::List DRIFTHASHExponentialUnique, Rcpp::List discreteDIFFUSIONUnique);


#endif
