#include "CartesianGaussian.h"
#include <armadillo>
#include <vector>

CartesianGaussian::CartesianGaussian(const arma::vec &center,
                                     const std::vector<double> &exponents,
                                     const arma::Col<long long> &angularMomentum,
                                     const std::vector<double> &coefficients)
    : center_(center), exponents_(exponents), angularMomentum_(angularMomentum), contractionCoeffs_(coefficients)
{
}

const arma::vec &CartesianGaussian::getCenter() const { return center_; }
const arma::Col<long long> &CartesianGaussian::getAngularMomentum() const { return angularMomentum_; }
const std::vector<double> &CartesianGaussian::getExponents() const { return exponents_; }
const std::vector<double> &CartesianGaussian::getContractionCoeffs() const { return contractionCoeffs_; }
