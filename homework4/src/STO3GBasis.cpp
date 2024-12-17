#include "STO3GBasis.h"

void STO3GBasis::addBasisFunction(double exponent, double coefficient)
{
    exponents.push_back(exponent);
    coefficients.push_back(coefficient);
}
