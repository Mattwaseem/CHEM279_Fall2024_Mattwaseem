#ifndef OVERLAPINTEGRAL_HPP
#define OVERLAPINTEGRAL_HPP
#pragma once

#include "CartesianGaussian.hpp"

class OverlapIntegral
{
public:
    // Function to compute overlap in a single dimension (x, y, or z)
    static double computeOverlap1D(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension);

    // Function to calculate the full overlap integral in 3D
    static double computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2);

private:
    // Helper function to calculate double factorial
    static int doubleFactorial(int n);
};

#endif // OVERLAPINTEGRAL_HPP
