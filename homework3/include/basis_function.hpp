#ifndef BASIS_FUNCTION_HPP
#define BASIS_FUNCTION_HPP

#include <armadillo>
#include <vector>
#include <string>

class BasisFunction
{
public:
    BasisFunction(const arma::vec &center, const arma::ivec &angularMomentum, const std::vector<double> &exponents, const std::vector<double> &coefficients);

    // Getters
    const arma::vec &getCenter() const;
    const arma::ivec &getAngularMomentum() const;
    const std::vector<double> &getExponents() const;
    const std::vector<double> &getCoefficients() const;

private:
    arma::vec center_;
    arma::ivec angularMomentum_;
    std::vector<double> exponents_;
    std::vector<double> coefficients_;
};

#endif
