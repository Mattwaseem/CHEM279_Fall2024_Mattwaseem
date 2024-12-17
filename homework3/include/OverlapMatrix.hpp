#ifndef OVERLAPMATRIX_HPP
#define OVERLAPMATRIX_HPP
#pragma once

#include "CartesianGaussian.hpp"
#include <armadillo> // Include Armadillo for matrix operations
#include <vector>

class OverlapMatrix
{
public:
    // Constructor to initialize the OverlapMatrix object with basis functions
    OverlapMatrix(const std::vector<CartesianGaussian> &basisFunctions);

    // Function to compute overlap in a single dimension (x, y, or z) for primitive Gaussians with provided exponents
    double computeOverlap3DPrimitive(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension, double alpha1, double alpha2);

    // Function to calculate the full overlap integral in 3D for contracted Gaussians
    double computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2);

    // Function to compute the full overlap matrix
    void computeOverlapMatrix();

    // Function to print the overlap matrix
    void printMatrix() const;

    // Function to get the overlap matrix (now returns an Armadillo matrix)
    const arma::mat &getMatrix() const { return overlapMatrix_; }

private:
    std::vector<CartesianGaussian> basisFunctions_; // Basis functions for the matrix
    arma::mat overlapMatrix_;                       // Overlap matrix stored as an Armadillo matrix

    // Helper function to calculate double factorial
    int doubleFactorial(int n);
};

#endif
