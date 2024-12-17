#ifndef MOLECULE_H
#define MOLECULE_H

#include "CartesianGaussian.h"
#include <vector>
#include <string>

struct Atom
{
    std::string element;
    double x, y, z;
};

class Molecule
{
public:
    Molecule();
    void computeBasisFunctions();
    const std::vector<CartesianGaussian> &getBasisFunctions() const;
    const std::vector<Atom> &getAtoms() const;
    int getNumBasisFunctions() const;

    void addAtom(const Atom &atom);

private:
    std::vector<Atom> atoms;
    std::vector<CartesianGaussian> basisFunctions;
    int numBasisFunctions;
    int numElectrons;
};

#endif // MOLECULE_H
