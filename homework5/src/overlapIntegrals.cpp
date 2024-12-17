#include "overlapIntegrals.h"
#include <cmath>
#include "factorial.h"

double overlapIntegral1D(double alpha, double beta, double center_a, double center_b, int lA, int lB)
{
    double prefactor = exp(-alpha * beta * pow(center_a - center_b, 2) / (alpha + beta));
    prefactor *= sqrt(M_PI / (alpha + beta));
    double center_product = (alpha * center_a + beta * center_b) / (alpha + beta);
    double sum = 0.0;
    for (int i = 0; i <= lA; i++)
    {
        for (int j = 0; j <= lB; j++)
        {
            if ((i + j) % 2 == 0)
            {
                sum += binomialCoef(lA, i) * binomialCoef(lB, j) *
                       (doubleFactorial(i + j - 1) * pow(center_product - center_a, lA - i) *
                        pow(center_product - center_b, lB - j)) /
                       pow(2 * (alpha + beta), static_cast<double>(i + j) / 2);
            }
        }
    }
    double integral = prefactor * sum;
    return integral;
}

double overlapIntegral3D(arma::rowvec centers_a, arma::rowvec centers_b, double alpha, double beta, arma::ivec lmn_a, arma::ivec lmn_b)
{
    double integral = overlapIntegral1D(alpha, beta, centers_a(0), centers_b(0), lmn_a(0), lmn_b(0)) *
                      overlapIntegral1D(alpha, beta, centers_a(1), centers_b(1), lmn_a(1), lmn_b(1)) *
                      overlapIntegral1D(alpha, beta, centers_a(2), centers_b(2), lmn_a(2), lmn_b(2));
    return integral;
}

double calcContractedOverlap(AO AO1, AO AO2)
{
    double contracted_overlap = 0.0;
    for (int k = 0; k < 3; k++)
    {
        for (int l = 0; l < 3; l++)
        {
            double unnorm_overlap = overlapIntegral3D(AO1.center, AO2.center,
                                                      AO1.exponents(k), AO2.exponents(l),
                                                      AO1.lmn, AO2.lmn);
            contracted_overlap += AO1.contraction_coeffs(k) * AO2.contraction_coeffs(l) *
                                  AO1.norm_constants(k) * AO2.norm_constants(l) * unnorm_overlap;
        }
    }
    return contracted_overlap;
}

arma::mat calcOverlapMatrix(Molecule molecule)
{
    arma::mat overlap_matrix = arma::zeros<arma::mat>(molecule.nBasisFunctions, molecule.nBasisFunctions);
    for (int i = 0; i < molecule.nBasisFunctions; i++)
    {
        for (int j = 0; j < molecule.nBasisFunctions; j++)
        {
            overlap_matrix(i, j) = calcContractedOverlap(molecule.basisFunctionsList[i], molecule.basisFunctionsList[j]);
        }
    }
    return overlap_matrix;
}
