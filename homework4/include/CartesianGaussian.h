#ifndef CARTESIANGAUSSIAN_H
#define CARTESIANGAUSSIAN_H

#include <armadillo>
#include <vector>

class CartesianGaussian
{
public:
    CartesianGaussian(const arma::vec &center,
                      const std::vector<double> &exponents,
                      const arma::Col<long long> &angularMomentum,
                      const std::vector<double> &coefficients);

    const arma::vec &getCenter() const;
    const arma::Col<long long> &getAngularMomentum() const;
    const std::vector<double> &getExponents() const;
    const std::vector<double> &getContractionCoeffs() const;

private:
    arma::vec center_;
    std::vector<double> exponents_;
    arma::Col<long long> angularMomentum_;
    std::vector<double> contractionCoeffs_;
};

#endif
