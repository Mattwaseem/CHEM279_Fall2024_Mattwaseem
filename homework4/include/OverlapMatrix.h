#ifndef OVERLAPMATRIX_H
#define OVERLAPMATRIX_H

#include <vector>
#include <armadillo>
#include "CartesianGaussian.h"

class OverlapMatrix
{
public:
    OverlapMatrix(const std::vector<CartesianGaussian> &basisFunctions);
    void computeOverlapMatrix();
    arma::mat getMatrix() const;
    void printMatrix() const;

    double computeOverlap3DPrimitive(const CartesianGaussian &g1, const CartesianGaussian &g2, int dimension, double alpha1, double alpha2);
    double computeTotalOverlap(const CartesianGaussian &g1, const CartesianGaussian &g2);
    int doubleFactorial(int n);

private:
    std::vector<CartesianGaussian> basisFunctions_;
    arma::mat overlapMatrix_;
};

#endif
