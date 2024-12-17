#ifndef GAMMACALCULATOR_H
#define GAMMACALCULATOR_H

#include <armadillo>
#include <vector>
#include "Constants.h"

class GammaCalculator
{
public:
    GammaCalculator(const Constants &constants);
    double calculateGamma(size_t mu, size_t nu, const std::vector<double> &alphas, const std::vector<double> &d_total,
                          int atomicNumber_mu, int atomicNumber_nu, const arma::vec &position_mu, const arma::vec &position_nu);

private:
    Constants constants_;
    double Zero_zero(double distance, int atomicNumber_mu, int atomicNumber_nu);
};

#endif
