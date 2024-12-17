#ifndef MOLECULE_H
#define MOLECULE_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <armadillo>
#include <cmath>
#include <cassert>

class AO
{
public:
    std::string AO_type;
    std::string chemSym;
    int valence;
    arma::rowvec center;
    arma::ivec lmn;
    arma::vec exponents;
    arma::vec contraction_coeffs;
    arma::vec norm_constants;

    AO(std::string AO_type, std::string chemSym, int valence, arma::rowvec center, arma::ivec lmn, arma::vec exponents, arma::vec contraction_coeffs);
    void calcNormConstants();
};

class Molecule
{
public:
    Molecule(const char *filename);

    int nAtoms;
    arma::ivec atomicNumbers;
    std::vector<std::string> atomicSymbols;
    arma::ivec atomValences;
    arma::mat coordinates;
    int charge;

    int hydrogenCount;
    int heavyAtomCount;
    int nElectrons;
    int pAlpha;
    int qBeta;
    int nBasisFunctions;
    std::vector<AO> basisFunctionsList;

    std::map<std::string, arma::vec> exponentsMap;
    std::map<std::string, arma::vec> sContractionMap;
    std::map<std::string, arma::vec> pContractionMap;

    void printMoleculeInfo();

private:
    int countBasisFunctions();
    int countElectrons();
    std::vector<AO> buildBasisFunctionsList();
};

#endif // MOLECULE_H
