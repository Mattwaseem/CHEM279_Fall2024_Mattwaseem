#ifndef MOLECULE_HPP
#define MOLECULE_HPP

#include "basis_function.hpp"
#include "CartesianGaussian.hpp" // Assuming needed for Gaussian functions
#include <vector>
#include <armadillo>

struct Atom
{
    std::string element;
    double x, y, z;
};

class Molecule
{
public:
    Molecule(); // Constructor

    // Function to compute basis functions and electrons based on the molecule's atoms
    void computeBasisFunctions();

    // Getter functions for basis functions, electrons, and atoms
    const std::vector<CartesianGaussian> &getBasisFunctions() const { return basisFunctions; }
    int getNumBasisFunctions() const { return numBasisFunctions; }
    int getNumElectrons() const { return numElectrons; }
    std::vector<Atom> &getAtoms() { return atoms; }

private:
    std::vector<Atom> atoms;
    std::vector<CartesianGaussian> basisFunctions;
    int numBasisFunctions;
    int numElectrons;
};

#endif
