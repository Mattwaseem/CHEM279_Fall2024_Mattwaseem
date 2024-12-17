#include "CartesianGaussian.hpp"
#include <cmath>

// Constructor for CartesianGaussian
CartesianGaussian::CartesianGaussian(const arma::vec &center, double exponent, const arma::ivec &angularMomentum)
    : center_(center), exponent_(exponent), angularMomentum_(angularMomentum)
{
    // Debugging print statement for constructor arguments
    // std::cout << "Creating CartesianGaussian with center = " << center_.t()
    //           << "exponent = " << exponent_
    //           << ", angular momentum = " << angularMomentum_.t() << std::endl;

    normalizationConstant_ = calculateNormalizationConstant();

    // Debugging print statement for normalization
    // std::cout << "Normalization constant calculated: " << normalizationConstant_ << std::endl;
}

// Calculate the normalization constant for the Gaussian
double CartesianGaussian::calculateNormalizationConstant() const
{
    double alpha = this->exponent_;
    int l = this->angularMomentum_(0);
    int m = this->angularMomentum_(1);
    int n = this->angularMomentum_(2);

    // Debug print statements to check l, m, and n before calculating double factorial
    // std::cout << "Alpha: " << alpha << std::endl;
    // std::cout << "Angular momentum (l, m, n): " << l << ", " << m << ", " << n << std::endl;

    // Check if l, m, or n are negative
    if (l < 0 || m < 0 || n < 0)
    {
        std::cerr << "Error: Negative angular momentum component detected!" << std::endl;
    }

    double prefactor = std::pow(2 * alpha / M_PI, 0.75);
    double norm_x = std::sqrt(std::pow(2 * alpha, l) / doubleFactorial(2 * l - 1));
    double norm_y = std::sqrt(std::pow(2 * alpha, m) / doubleFactorial(2 * m - 1));
    double norm_z = std::sqrt(std::pow(2 * alpha, n) / doubleFactorial(2 * n - 1));
    // Debugging print statements
    // std::cout << "Alpha: " << alpha << std::endl;
    // std::cout << "Angular momentum (l, m, n): " << l << ", " << m << ", " << n << std::endl;
    // std::cout << "Prefactor: " << prefactor << std::endl;
    // std::cout << "Norm_x: " << norm_x << ", Norm_y: " << norm_y << ", Norm_z: " << norm_z << std::endl;
    // std::cout << "Normalization constant: " << (prefactor * norm_x * norm_y * norm_z) << std::endl;

    return prefactor * norm_x * norm_y * norm_z;
}

// Calculate the product center of the Gaussian
arma::vec CartesianGaussian::productCenter(const CartesianGaussian &other) const
{
    double alpha1 = this->exponent_;
    double alpha2 = other.getExponent();
    arma::vec R1 = this->center_;
    arma::vec R2 = other.getCenter();

    // Center of product Gaussian (RP)
    double gamma = alpha1 + alpha2;
    arma::vec RP = (alpha1 * R1 + alpha2 * R2) / gamma;

    return RP;
}

// Helper function to calculate double factorial
int CartesianGaussian::doubleFactorial(int n) const
{
    // if (n <= 0)
    //     std::cout << "Double factorial of " << n << " is 1 (edge case)" << std::endl; // my Debugging print
    // return 1;                                                                         // edge case not needed since input file not -1
    int result = 1;
    for (int i = n; i > 0; i -= 2)
    {
        result *= i;
        // std::cout << "Calculating double factorial: n = " << n << ", i = " << i << ", intermediate result = " << result << std::endl; // Debugging
    }
    // std::cout << "Double factorial of " << n << " is " << result << std::endl; // Debugging print statement
    return result;
}
