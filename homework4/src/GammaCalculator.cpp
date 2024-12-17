#include "GammaCalculator.h"
#include <cmath>
#include <armadillo>
#include <iostream>

GammaCalculator::GammaCalculator(const Constants &constants) : constants_(constants)
{
}

double GammaCalculator::Zero_zero(double distance, int atomicNumber_mu, int atomicNumber_nu)
{
    if (distance < 1e-8)
    {
        double gamma_AA = constants_.getBondingParameter(atomicNumber_mu) * (1.0 / 27.211);
        return gamma_AA;
    }
    else
    {
        double gamma_AB = (constants_.getBondingParameter(atomicNumber_mu) + constants_.getBondingParameter(atomicNumber_nu)) / 2.0 * (1.0 / 27.211);
        return gamma_AB;
    }
}

double GammaCalculator::calculateGamma(size_t mu, size_t nu, const std::vector<double> &alphas, const std::vector<double> &d_total, int atomicNumber_mu, int atomicNumber_nu, const arma::vec &position_mu, const arma::vec &position_nu)
{
    double distance = arma::norm(position_mu - position_nu) * constants_.angstrom_to_bohr;
    double Gamma_AB = Zero_zero(distance, atomicNumber_mu, atomicNumber_nu);
    return Gamma_AB;
}
