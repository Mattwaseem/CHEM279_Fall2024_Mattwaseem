
#include "basis_function.hpp"

BasisFunction::BasisFunction(const arma::vec &center, const arma::ivec &angularMomentum, const std::vector<double> &exponents, const std::vector<double> &coefficients)
    : center_(center), angularMomentum_(angularMomentum), exponents_(exponents), coefficients_(coefficients) {}

const arma::vec &BasisFunction::getCenter() const
{
    return center_;
}

const arma::ivec &BasisFunction::getAngularMomentum() const
{
    return angularMomentum_;
}

const std::vector<double> &BasisFunction::getExponents() const
{
    return exponents_;
}

const std::vector<double> &BasisFunction::getCoefficients() const
{
    return coefficients_;
}
