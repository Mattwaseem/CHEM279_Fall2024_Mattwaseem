#include "Molecule.h"

Molecule::Molecule() : numBasisFunctions(0), numElectrons(0) {}

void Molecule::computeBasisFunctions()
{
    basisFunctions.clear();
    numElectrons = 0;

    for (const auto &atom : atoms)
    {
        arma::vec center = {atom.x, atom.y, atom.z};
        arma::ivec angularMomentumS = {0, 0, 0};

        std::vector<double> exponents;
        std::vector<double> coefficients;

        if (atom.element == "H")
        {
            exponents = {3.42525091, 0.62391373, 0.16885540};
            coefficients = {0.15432897, 0.53532814, 0.44463454};
            numElectrons += 1;

            CartesianGaussian s_orbital(center, exponents, angularMomentumS, coefficients);
            basisFunctions.push_back(s_orbital);
        }
        else if (atom.element == "C" || atom.element == "N" || atom.element == "O" || atom.element == "F")
        {
            exponents = {71.6168370, 13.0450963, 3.5305122};
            coefficients = {0.15432897, 0.53532814, 0.44463454};
            numElectrons += 4;

            CartesianGaussian s_orbital(center, exponents, angularMomentumS, coefficients);
            basisFunctions.push_back(s_orbital);

            std::vector<arma::ivec> angularMomentumP = {arma::ivec{1, 0, 0}, arma::ivec{0, 1, 0}, arma::ivec{0, 0, 1}};
            for (const auto &angMom : angularMomentumP)
            {
                CartesianGaussian p_orbital(center, exponents, angMom, coefficients);
                basisFunctions.push_back(p_orbital);
            }
        }
    }

    numBasisFunctions = basisFunctions.size();
}

const std::vector<CartesianGaussian> &Molecule::getBasisFunctions() const
{
    return basisFunctions;
}

const std::vector<Atom> &Molecule::getAtoms() const
{
    return atoms;
}

int Molecule::getNumBasisFunctions() const
{
    return numBasisFunctions;
}

void Molecule::addAtom(const Atom &atom)
{
    atoms.emplace_back(atom);
}
