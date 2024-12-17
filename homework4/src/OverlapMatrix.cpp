#include "OverlapMatrix.h"
#include <cmath>
#include <iostream>
#include <armadillo>

OverlapMatrix::OverlapMatrix(const std::vector<CartesianGaussian> &basisFunctions)
    : basisFunctions_(basisFunctions), overlapMatrix_(basisFunctions.size(), basisFunctions.size(), arma::fill::zeros)
{
}

double OverlapMatrix::computeOverlap3DPrimitive(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension, double alpha1, double alpha2)
{
    double x1 = g1.getCenter()[dimension];
    double x2 = g2.getCenter()[dimension];
    int l1 = g1.getAngularMomentum()[dimension];
    int l2 = g2.getAngularMomentum()[dimension];

    double gamma = alpha1 + alpha2;
    double RP = (alpha1 * x1 + alpha2 * x2) / gamma;

    double prefactor = std::exp(-alpha1 * alpha2 * std::pow(x1 - x2, 2) / gamma) * std::sqrt(M_PI / gamma);

    double overlap = 0.0;
    for (int i = 0; i <= l1; ++i)
    {
        for (int j = 0; j <= l2; ++j)
        {
            if ((i + j) % 2 == 0)
            {
                int comb_l1_i = std::tgamma(l1 + 1) / (std::tgamma(i + 1) * std::tgamma(l1 - i + 1));
                int comb_l2_j = std::tgamma(l2 + 1) / (std::tgamma(j + 1) * std::tgamma(l2 - j + 1));
                int double_fact = doubleFactorial(i + j - 1);

                overlap += comb_l1_i * comb_l2_j * double_fact *
                           std::pow(RP - x1, l1 - i) * std::pow(RP - x2, l2 - j) / std::pow(2 * gamma, (i + j) / 2.0);
            }
        }
    }

    return prefactor * overlap;
}

double OverlapMatrix::computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2)
{
    const std::vector<double> &exponents1 = g1.getExponents();
    const std::vector<double> &exponents2 = g2.getExponents();
    const std::vector<double> &coeffs1 = g1.getContractionCoeffs();
    const std::vector<double> &coeffs2 = g2.getContractionCoeffs();

    double totalOverlap = 0.0;

    for (size_t p = 0; p < exponents1.size(); ++p)
    {
        for (size_t q = 0; q < exponents2.size(); ++q)
        {
            double alpha1 = exponents1[p];
            double alpha2 = exponents2[q];

            double Sx = computeOverlap3DPrimitive(g1, g2, 0, alpha1, alpha2);
            double Sy = computeOverlap3DPrimitive(g1, g2, 1, alpha1, alpha2);
            double Sz = computeOverlap3DPrimitive(g1, g2, 2, alpha1, alpha2);

            double primitiveOverlap = Sx * Sy * Sz;

            double norm1_x = std::sqrt(std::pow(2 * alpha1, g1.getAngularMomentum()[0]) / doubleFactorial(2 * g1.getAngularMomentum()[0] - 1));
            double norm1_y = std::sqrt(std::pow(2 * alpha1, g1.getAngularMomentum()[1]) / doubleFactorial(2 * g1.getAngularMomentum()[1] - 1));
            double norm1_z = std::sqrt(std::pow(2 * alpha1, g1.getAngularMomentum()[2]) / doubleFactorial(2 * g1.getAngularMomentum()[2] - 1));

            double norm2_x = std::sqrt(std::pow(2 * alpha2, g2.getAngularMomentum()[0]) / doubleFactorial(2 * g2.getAngularMomentum()[0] - 1));
            double norm2_y = std::sqrt(std::pow(2 * alpha2, g2.getAngularMomentum()[1]) / doubleFactorial(2 * g2.getAngularMomentum()[1] - 1));
            double norm2_z = std::sqrt(std::pow(2 * alpha2, g2.getAngularMomentum()[2]) / doubleFactorial(2 * g2.getAngularMomentum()[2] - 1));

            double norm1 = norm1_x * norm1_y * norm1_z;
            double norm2 = norm2_x * norm2_y * norm2_z;

            totalOverlap += coeffs1[p] * coeffs2[q] * norm1 * norm2 * primitiveOverlap;
        }
    }

    if (std::abs(totalOverlap) < 1e-10)
    {
        totalOverlap = 0.0;
    }

    return totalOverlap;
}

void OverlapMatrix::computeOverlapMatrix()
{
    size_t n = basisFunctions_.size();
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j <= i; ++j)
        {
            double overlap = computeTotalOverlap(basisFunctions_[i], basisFunctions_[j]);

            if (std::abs(overlap) < 1e-10)
            {
                overlap = 0.0;
            }

            overlapMatrix_(i, j) = overlap;
            overlapMatrix_(j, i) = overlap;
        }
    }

    std::cout << "Overlap matrix size: " << overlapMatrix_.n_rows << " x " << overlapMatrix_.n_cols << std::endl;
}

arma::mat OverlapMatrix::getMatrix() const
{
    return overlapMatrix_;
}

void OverlapMatrix::printMatrix() const
{
    std::cout << "Overlap matrix:" << std::endl;
    overlapMatrix_.print();
}

int OverlapMatrix::doubleFactorial(int n)
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
