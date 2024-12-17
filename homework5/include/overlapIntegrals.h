#ifndef OVERLAPINTEGRALS_H
#define OVERLAPINTEGRALS_H

#include "Molecule.h"

double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB);
double overlapIntegral3D(arma::rowvec centers_a, arma::rowvec centers_b, double alpha, double beta, arma::ivec lmn_a, arma::ivec lmn_b);
double calcContractedOverlap(AO AO1, AO AO2);
arma::mat calcOverlapMatrix(Molecule molecule);

#endif // OVERLAPINTEGRALS_H
