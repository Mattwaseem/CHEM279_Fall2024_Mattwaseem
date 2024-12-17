#ifndef FOCKMATRIX_H
#define FOCKMATRIX_H

#include <armadillo>
#include <vector>
#include <unordered_map>
#include "DensityMatrix.h"
#include "OverlapMatrix.h"
#include "GammaCalculator.h"
#include "Constants.h"
#include "CartesianGaussian.h"

class FockMatrix
{
public:
    FockMatrix(const DensityMatrix &densityMatrix,
               const OverlapMatrix &overlapMatrix,
               const arma::mat &H_core,
               const std::vector<int> &atomicNumbersPerBasisFunction,
               const std::vector<double> &alphas,
               const std::vector<double> &d_total,
               const std::vector<CartesianGaussian> &basisFunctions,
               const Constants &constants);

    arma::mat computeFAlpha();
    arma::mat computeFBeta();

private:
    DensityMatrix densityMatrix_;
    OverlapMatrix overlapMatrix_;
    arma::mat H_core_;
    std::vector<int> atomicNumbersPerBasisFunction_;
    std::vector<double> alphas_;
    std::vector<double> d_total_;
    std::vector<CartesianGaussian> basisFunctions_;
    GammaCalculator gammaCalculator_;
    std::unordered_map<size_t, double> I_mu_map_;
    std::unordered_map<size_t, double> A_mu_map_;
    std::unordered_map<size_t, double> beta_A_map_;
    Constants constants_;

    void initializeParameters();
};

#endif
