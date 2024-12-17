#ifndef STO3GBASIS_H
#define STO3GBASIS_H

#include <vector>

class STO3GBasis
{
public:
    std::vector<double> exponents;
    std::vector<double> coefficients;

    void addBasisFunction(double exponent, double coefficient);
};

#endif
