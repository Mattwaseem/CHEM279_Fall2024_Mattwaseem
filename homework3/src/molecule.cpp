#include "molecule.hpp"
#include "basis_function.hpp"
#include "CartesianGaussian.hpp"
#include <vector>

Molecule::Molecule() : numBasisFunctions(0), numElectrons(0) {}

void Molecule::computeBasisFunctions()
{
    int numCarbon = 0;
    int numHydrogen = 0;

    basisFunctions.clear();
    numElectrons = 0;

    for (const auto &atom : atoms)
    {
        arma::vec center = {atom.x, atom.y, atom.z};
        arma::ivec angularMomentumS = {0, 0, 0};

        std::vector<double> hydrogenExponents = {3.42525091, 0.62391373, 0.16885540};
        std::vector<double> hydrogenCoefficients = {0.15432897, 0.53532814, 0.44463454};
        std::vector<double> carbonExponents = {71.6168370, 13.0450963, 3.5305122};
        std::vector<double> carbonCoefficients = {0.15432897, 0.53532814, 0.44463454};

        if (atom.element == "H")
        {
            numHydrogen++;
            CartesianGaussian s_orbital(center, hydrogenExponents, angularMomentumS, hydrogenCoefficients);
            basisFunctions.push_back(s_orbital);
            numElectrons += 1;
        }
        else if (atom.element == "C")
        {
            numCarbon++;
            CartesianGaussian s_orbital(center, carbonExponents, angularMomentumS, carbonCoefficients);
            basisFunctions.push_back(s_orbital);

            std::vector<arma::ivec> angularMomentumP = {
                {1, 0, 0},
                {0, 1, 0},
                {0, 0, 1}};

            for (const auto &angularMomentum : angularMomentumP)
            {
                CartesianGaussian p_orbital(center, carbonExponents, angularMomentum, carbonCoefficients);
                basisFunctions.push_back(p_orbital);
            }

            numElectrons += 4;
        }
    }

    numBasisFunctions = basisFunctions.size();

    std::cout << "Number of basis functions: " << numBasisFunctions << std::endl;
    std::cout << "Number of electrons: " << numElectrons << std::endl;

    if (numBasisFunctions % 2 != 0)
    {
        std::cerr << "Error: The number of basis functions is not an even number, which may cause issues." << std::endl;
    }
}
