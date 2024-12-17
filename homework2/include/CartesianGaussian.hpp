#pragma once

#include <armadillo>

class CartesianGaussian
{
private:
    arma::vec center_;
    double exponent_;
    arma::ivec angularMomentum_;
    double normalizationConstant_; // normalization
    int doubleFactorial(int n) const;

public:
    // Constructor
    CartesianGaussian(const arma::vec &center, double exponent, const arma::ivec &angularMomentum);

    // Getters
    arma::vec getCenter() const { return center_; }
    double getExponent() const { return exponent_; }
    arma::ivec getAngularMomentum() const { return angularMomentum_; }
    double getNormalizationConstant() const { return normalizationConstant_; }

    // Function to compute the product center of two Gaussians
    arma::vec productCenter(const CartesianGaussian &other) const;

private:
    // Function to calculate the normalization constant
    double calculateNormalizationConstant() const;
};
