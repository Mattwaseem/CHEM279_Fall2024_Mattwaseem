#include "Molecule.h"
#include "overlapIntegrals.h"
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cassert>
#include <cmath>

AO::AO(std::string AO_type, std::string chemSym, int valence, arma::rowvec center, arma::ivec lmn, arma::vec exponents, arma::vec contraction_coeffs)
    : AO_type(AO_type), chemSym(chemSym), valence(valence), center(center), lmn(lmn), exponents(exponents),
      contraction_coeffs(contraction_coeffs)
{
    calcNormConstants();
}

void AO::calcNormConstants()
{
    double selfOverlap_1 = overlapIntegral3D(center, center, exponents(0), exponents(0), lmn, lmn);
    double selfOverlap_2 = overlapIntegral3D(center, center, exponents(1), exponents(1), lmn, lmn);
    double selfOverlap_3 = overlapIntegral3D(center, center, exponents(2), exponents(2), lmn, lmn);
    norm_constants = {1.0 / sqrt(selfOverlap_1), 1.0 / sqrt(selfOverlap_2), 1.0 / sqrt(selfOverlap_3)};
}

Molecule::Molecule(const char *filename)
{
    std::ifstream infile(filename);
    assert(infile.good());

    infile >> nAtoms >> charge;

    atomicNumbers.resize(nAtoms);
    atomicSymbols.resize(nAtoms);
    atomValences.resize(nAtoms);
    coordinates.resize(nAtoms, 3);

    std::string line;
    std::getline(infile, line);

    for (int i = 0; i < nAtoms; i++)
    {
        std::getline(infile, line);
        std::istringstream iss(line);
        int atomicNumber;
        iss >> atomicNumber;
        atomicNumbers(i) = atomicNumber;
        for (int j = 0; j < 3; j++)
        {
            iss >> coordinates(i, j);
        }
    }
    infile.close();

    std::map<int, std::string> atomicSymbolsMap = {{1, "H"},
                                                   {6, "C"},
                                                   {7, "N"},
                                                   {8, "O"},
                                                   {9, "F"}};

    std::map<int, int> valenceElectronsMap = {{1, 1},
                                              {6, 4},
                                              {7, 5},
                                              {8, 6},
                                              {9, 7}};

    for (int i = 0; i < nAtoms; i++)
    {
        atomicSymbols[i] = atomicSymbolsMap[atomicNumbers(i)];
        atomValences(i) = valenceElectronsMap[atomicNumbers(i)];
    }

    hydrogenCount = 0;
    heavyAtomCount = 0;
    for (int i = 0; i < nAtoms; i++)
    {
        if (atomicNumbers(i) == 1)
        {
            hydrogenCount++;
        }
        else
        {
            heavyAtomCount++;
        }
    }

    nBasisFunctions = countBasisFunctions();
    nElectrons = countElectrons();
    int remainder = nElectrons % 2;
    pAlpha = nElectrons / 2 + remainder;
    qBeta = nElectrons / 2;

    exponentsMap["H"] = {3.42525091, 0.62391373, 0.16885540};
    sContractionMap["H"] = {0.15432897, 0.53532814, 0.44463454};

    exponentsMap["C"] = {2.94124940, 0.68348310, 0.22228990};
    sContractionMap["C"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["C"] = {0.15591627, 0.60768372, 0.39195739};

    exponentsMap["N"] = {3.78045590, 0.87849660, 0.28571440};
    sContractionMap["N"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["N"] = {0.15591627, 0.60768372, 0.39195739};

    exponentsMap["O"] = {5.03315130, 1.16959610, 0.38038900};
    sContractionMap["O"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["O"] = {0.15591627, 0.60768372, 0.39195739};

    exponentsMap["F"] = {6.46480320, 1.50228120, 0.48858850};
    sContractionMap["F"] = {-0.09996723, 0.39951283, 0.70011547};
    pContractionMap["F"] = {0.15591627, 0.60768372, 0.39195739};

    basisFunctionsList = buildBasisFunctionsList();
}

int Molecule::countBasisFunctions()
{
    return 4 * heavyAtomCount + hydrogenCount;
}

int Molecule::countElectrons()
{
    int nElectronsTotal = 0;
    for (int i = 0; i < nAtoms; i++)
    {
        nElectronsTotal += atomValences(i);
    }
    return nElectronsTotal - charge;
}

std::vector<AO> Molecule::buildBasisFunctionsList()
{
    std::vector<AO> basisFunctionsList;
    for (int i = 0; i < nAtoms; i++)
    {
        std::string atomSym = atomicSymbols[i];
        int valence = atomValences[i];
        std::string orbType;
        if (atomSym == "H")
        {
            orbType = "1s";
        }
        else
        {
            orbType = "2s";
        }
        arma::ivec lmn = {0, 0, 0};
        AO AO_s = AO(orbType, atomSym, valence, coordinates.row(i), lmn, exponentsMap[atomSym], sContractionMap[atomSym]);
        basisFunctionsList.push_back(AO_s);
        if (orbType == "2s")
        {
            arma::ivec lmn_px = {1, 0, 0};
            arma::ivec lmn_py = {0, 1, 0};
            arma::ivec lmn_pz = {0, 0, 1};
            AO AO_2px = AO("2px", atomSym, valence, coordinates.row(i), lmn_px, exponentsMap[atomSym], pContractionMap[atomSym]);
            AO AO_2py = AO("2py", atomSym, valence, coordinates.row(i), lmn_py, exponentsMap[atomSym], pContractionMap[atomSym]);
            AO AO_2pz = AO("2pz", atomSym, valence, coordinates.row(i), lmn_pz, exponentsMap[atomSym], pContractionMap[atomSym]);
            basisFunctionsList.push_back(AO_2px);
            basisFunctionsList.push_back(AO_2py);
            basisFunctionsList.push_back(AO_2pz);
        }
    }
    return basisFunctionsList;
}

void Molecule::printMoleculeInfo()
{
    std::cout << "Number of atoms: " << nAtoms << std::endl;
    std::cout << "Charge: " << charge << std::endl;
    std::cout << "Atomic numbers: " << atomicNumbers.t();
    std::cout << "Atomic symbols: ";
    for (int i = 0; i < nAtoms; i++)
    {
        std::cout << atomicSymbols[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "Coordinates:\n"
              << coordinates << std::endl;
    std::cout << "Number of basis functions: " << nBasisFunctions << std::endl;
    std::cout << "Number of electrons: " << nElectrons << std::endl;
    std::cout << "Number of alpha electrons: " << pAlpha << std::endl;
    std::cout << "Number of beta electrons: " << qBeta << std::endl;
}
