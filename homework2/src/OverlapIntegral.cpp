#include "OverlapIntegral.hpp"
#include <cmath>

// Function to compute overlap (x, y, or z)
double OverlapIntegral::computeOverlap1D(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension)
{
    // Get centers, exponents, and angular momentum
    double x1 = g1.getCenter()[dimension];
    double x2 = g2.getCenter()[dimension];
    double alpha1 = g1.getExponent();
    double alpha2 = g2.getExponent();
    int l1 = g1.getAngularMomentum()[dimension];
    int l2 = g2.getAngularMomentum()[dimension];

    // Calculate product center
    double gamma = alpha1 + alpha2;
    double RP = (alpha1 * x1 + alpha2 * x2) / gamma;

    // Exponential prefactor
    double prefactor = std::exp(-alpha1 * alpha2 * std::pow(x1 - x2, 2) / gamma) * std::sqrt(M_PI / gamma);
    // // Debugging: Print the intermediate values due to not getting the expected result
    // std::cout << "Dimension: " << dimension << std::endl;
    // std::cout << "Centers x1: " << x1 << ", x2: " << x2 << std::endl;
    // std::cout << "Exponents alpha1: " << alpha1 << ", alpha2: " << alpha2 << ", gamma: " << gamma << std::endl;
    // std::cout << "Product center RP: " << RP << std::endl;
    // std::cout << "Exponential prefactor: " << prefactor << std::endl;

    // Double summation to calculate overlap integral
    double overlap = 0.0;
    for (int i = 0; i <= l1; ++i)
    {
        for (int j = 0; j <= l2; ++j)
        {
            // Only considers terms where (i + j) is even
            if ((i + j) % 2 == 0)
            {
                int comb_l1_i = std::tgamma(l1 + 1) / (std::tgamma(i + 1) * std::tgamma(l1 - i + 1));
                int comb_l2_j = std::tgamma(l2 + 1) / (std::tgamma(j + 1) * std::tgamma(l2 - j + 1));
                int double_fact = doubleFactorial(i + j - 1);
                // Debugging statements to print values for each combination
                std::cout << "i: " << i << ", j: " << j << std::endl;
                std::cout << "comb_l1_i: " << comb_l1_i << ", comb_l2_j: " << comb_l2_j << std::endl;
                std::cout << "double_fact: " << double_fact << std::endl;
                std::cout << "RP - x1: " << RP - x1 << ", RP - x2: " << RP - x2 << std::endl;
                std::cout << "Power terms: " << std::pow(RP - x1, l1 - i) << ", " << std::pow(RP - x2, l2 - j) << std::endl;
                std::cout << "Denominator term: " << std::pow(2 * gamma, (i + j) / 2.0) << std::endl;

                overlap += comb_l1_i * comb_l2_j * double_fact *
                           std::pow(RP - x1, l1 - i) * std::pow(RP - x2, l2 - j) / std::pow(2 * gamma, (i + j) / 2.0);
            }
        }
    }

    return prefactor * overlap;
}

// Function to calculate the full overlap integral in 3D
double OverlapIntegral::computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2)
{
    double Sx = computeOverlap1D(g1, g2, 0); // Overlap in x dimension
    double Sy = computeOverlap1D(g1, g2, 1); //  y
    double Sz = computeOverlap1D(g1, g2, 2); // z

    // Debugging the individual overlaps
    std::cout << "Sx: " << Sx << ", Sy: " << Sy << ", Sz: " << Sz << std::endl;

    // Factor in the normalization constants for both Gaussians
    double normalization = g1.getNormalizationConstant() * g2.getNormalizationConstant();

    // Debugging the normalization
    std::cout << "Normalization: " << normalization << std::endl;

    // double totalOverlap = Sx * Sy * Sz * normalization; //  this line was inccorrect causing double normalization and answer to be incorrect.
    double totalOverlap = Sx * Sy * Sz;

    // Debug the final total overlap result
    std::cout << "Final Total Overlap: " << totalOverlap << std::endl;
    return totalOverlap;
}

// Helper function to calculate double factorial
int OverlapIntegral::doubleFactorial(int n)
{
    if (n <= 0)
        return 1; // edge case for n = 0 or negative not really needed. 
    int result = 1;
    for (int i = n; i > 0; i -= 2)
    {
        result *= i;
    }
    return result;
}
