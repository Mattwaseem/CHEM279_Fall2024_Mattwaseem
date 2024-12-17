#include "CartesianGaussian.hpp"
#include <cmath>
#include <armadillo>

// Constructor for primitive Gaussians
CartesianGaussian::CartesianGaussian(const arma::vec &center, double exponent, const arma::ivec &angularMomentum)
    : center_(center), angularMomentum_(angularMomentum), exponents_({exponent})
{
    contractionCoeffs_ = {1.0};
    normalizationConstant_ = calculateNormalizationConstant();
}

// Constructor for contracted Gaussians
CartesianGaussian::CartesianGaussian(const arma::vec &center, const std::vector<double> &exponents, const arma::ivec &angularMomentum, const std::vector<double> &contractionCoeffs)
    : center_(center), angularMomentum_(angularMomentum), exponents_(exponents), contractionCoeffs_(contractionCoeffs)
{
    normalizationConstant_ = calculateNormalizationConstant();
}

// Calculate the normalization constant for the Gaussian
double CartesianGaussian::calculateNormalizationConstant() const
{
    double totalNormalization = 0.0;

    for (size_t i = 0; i < exponents_.size(); ++i)
    {
        double alpha = exponents_[i];
        double c = contractionCoeffs_[i];
        int l = angularMomentum_(0);
        int m = angularMomentum_(1);
        int n = angularMomentum_(2);

        double prefactor = std::pow(2 * alpha / M_PI, 0.75);

        double norm_x = std::sqrt(std::pow(2 * alpha, l) / doubleFactorial(2 * l - 1));
        double norm_y = std::sqrt(std::pow(2 * alpha, m) / doubleFactorial(2 * m - 1));
        double norm_z = std::sqrt(std::pow(2 * alpha, n) / doubleFactorial(2 * n - 1));

        double normalization_constant = prefactor * norm_x * norm_y * norm_z;

        totalNormalization += std::pow(c, 2) * normalization_constant;
    }

    return std::sqrt(totalNormalization);
}

// Calculate the product center of the Gaussian
arma::vec CartesianGaussian::productCenter(const CartesianGaussian &other) const
{
    double alpha1 = exponents_[0];
    double alpha2 = other.getExponents()[0];
    arma::vec R1 = center_;
    arma::vec R2 = other.getCenter();

    double gamma = alpha1 + alpha2;
    arma::vec RP = (alpha1 * R1 + alpha2 * R2) / gamma;

    return RP;
}

// Helper function to calculate double factorial
int CartesianGaussian::doubleFactorial(int n) const
{
    if (n <= 0)
        return 1;

    int result = 1;
    for (int i = n; i > 0; i -= 2)
    {
        result *= i;
    }

    return result;
}
