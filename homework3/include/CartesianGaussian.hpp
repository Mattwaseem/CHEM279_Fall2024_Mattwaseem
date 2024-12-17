#pragma once
#ifndef CARTESIANGAUSSIAN_HPP
#define CARTESIANGAUSSIAN_HPP
#include <armadillo>
#include <vector>

class CartesianGaussian
{
private:
    arma::vec center_;                      // Center of the Gaussian
    arma::ivec angularMomentum_;            // Angular momentum (l, m, n)
    std::vector<double> exponents_;         // Exponents for contracted Gaussians (single primitive Gaussians will have only one value)
    std::vector<double> contractionCoeffs_; // Contraction coefficients for contracted Gaussians
    double normalizationConstant_;          // Normalization constant for the Gaussian
    int doubleFactorial(int n) const;       // Function to calculate double factorial

public:
    // Constructor for primitive Gaussians
    CartesianGaussian(const arma::vec &center, double exponent, const arma::ivec &angularMomentum);

    // Constructor for contracted Gaussians (with contraction coefficients)
    CartesianGaussian(const arma::vec &center, const std::vector<double> &exponents, const arma::ivec &angularMomentum, const std::vector<double> &contractionCoeffs);

    // Getters
    arma::vec getCenter() const { return center_; }
    const std::vector<double> &getExponents() const { return exponents_; } // Return exponents for contracted Gaussians
    arma::ivec getAngularMomentum() const { return angularMomentum_; }
    double getNormalizationConstant() const { return normalizationConstant_; }
    const std::vector<double> &getContractionCoeffs() const { return contractionCoeffs_; } // Return contraction coefficients

    // Function to compute the product center of two Gaussians
    arma::vec productCenter(const CartesianGaussian &other) const;

private:
    // Function to calculate the normalization constant for the Gaussian
    double calculateNormalizationConstant() const;
};

#endif // CARTESIANGAUSSIAN_HPP
