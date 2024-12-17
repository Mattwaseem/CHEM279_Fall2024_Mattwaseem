#include <iostream>
#include "OverlapMatrix.hpp"
#include "Hamiltonian.hpp"
#include <armadillo>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include "molecule.hpp"
#include "input_parser.h"
#include <iomanip>
#include <cmath> // For sqrt

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];
    auto [numAtoms, charge, atomicNumbers, xCoords, yCoords, zCoords] = InputParser::parseMoleculeInput(inputFile);

    if (numAtoms == 0)
    {
        std::cerr << "Failed to parse the input file: " << inputFile << std::endl;
        return 1;
    }

    std::cout << "Number of atoms: " << numAtoms << ", Charge: " << charge << std::endl;

    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::cout << "Atomic Number: " << atomicNumbers[i] << ", X: " << xCoords[i] << ", Y: " << yCoords[i] << ", Z: " << zCoords[i] << std::endl;
    }

    Molecule molecule;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        std::string element = (atomicNumbers[i] == 1) ? "H" : "C";
        molecule.getAtoms().push_back({element, xCoords[i], yCoords[i], zCoords[i]});
    }

    molecule.computeBasisFunctions();

    std::cout << "Number of basis functions: " << molecule.getNumBasisFunctions() << std::endl;
    std::cout << "Number of electrons: " << molecule.getNumElectrons() << std::endl;

    auto basisFunctions = molecule.getBasisFunctions();

    // Debug statements to print basis function details
    for (const auto &basis : basisFunctions)
    {
        std::cout << "Basis function center: " << basis.getCenter().t() << std::endl;
        std::cout << "Angular momentum: " << basis.getAngularMomentum().t() << std::endl;
        std::cout << "Exponents: ";
        for (const auto &exp : basis.getExponents())
            std::cout << exp << " ";
        std::cout << std::endl;
        std::cout << "Coefficients: ";
        for (const auto &coef : basis.getContractionCoeffs())
            std::cout << coef << " ";
        std::cout << std::endl;
    }

    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();

    std::cout << std::fixed << std::setprecision(4);
    overlapMatrix.printMatrix();

    std::vector<double> diagEnergies;
    int numElectrons = molecule.getNumElectrons() - charge;
    for (size_t i = 0; i < atomicNumbers.size(); ++i)
    {
        if (atomicNumbers[i] == 1)
        {
            diagEnergies.push_back(-13.6);
        }
        else if (atomicNumbers[i] == 6)
        {
            diagEnergies.push_back(-21.4);
            diagEnergies.push_back(-11.4);
            diagEnergies.push_back(-11.4);
            diagEnergies.push_back(-11.4);
        }
        else
        {
            std::cerr << "Unsupported element with atomic number: " << atomicNumbers[i] << std::endl;
            return 1;
        }
    }

    if (diagEnergies.size() != static_cast<size_t>(molecule.getNumBasisFunctions())) // Cast to match types
    {
        std::cerr << "Error: The number of diagonal energies (" << diagEnergies.size()
                  << ") does not match the number of basis functions ("
                  << molecule.getNumBasisFunctions() << ")!" << std::endl;
        return 1;
    }

    Hamiltonian hamiltonian(overlapMatrix.getMatrix(), diagEnergies);
    hamiltonian.computeHamiltonianMatrix();

    arma::mat H = hamiltonian.getMatrix();
    arma::mat S = overlapMatrix.getMatrix();

    arma::vec eigval;
    arma::mat eigvec;
    arma::eig_sym(eigval, eigvec, S);
    arma::mat S_inv_sqrt = eigvec * arma::diagmat(1.0 / arma::sqrt(eigval)) * eigvec.t();

    arma::mat H_ortho = S_inv_sqrt.t() * H * S_inv_sqrt;

    arma::eig_sym(eigval, eigvec, H_ortho);

    arma::mat C = S_inv_sqrt * eigvec;

    std::cout << "Hamiltonian matrix" << std::endl;
    H.print();

    std::cout << "X_mat (S^(-1/2))" << std::endl;
    S_inv_sqrt.print();

    std::cout << "MO coefficients (C Matrix): " << std::endl;
    C.print();

    arma::mat MO_overlap = C.t() * S * C;
    std::cout << "MO overlap matrix:" << std::endl;
    MO_overlap.print();

    // Compute nuclear repulsion energy
    double nuclearRepulsionEnergy = 0.0;
    for (int i = 0; i < numAtoms; ++i) // Adjusted to int for matching types
    {
        for (int j = i + 1; j < numAtoms; ++j) // Adjusted to int for matching types
        {
            double distance = std::sqrt(
                std::pow(xCoords[i] - xCoords[j], 2) +
                std::pow(yCoords[i] - yCoords[j], 2) +
                std::pow(zCoords[i] - zCoords[j], 2));
            nuclearRepulsionEnergy += (atomicNumbers[i] * atomicNumbers[j]) / distance;
        }
    }

    // nuclear repulsion not added -- update to my total energy
    std::cout << "Nuclear repulsion energy: " << nuclearRepulsionEnergy << std::endl;

    // Add nuclear repulsion with total energy
    double total_energy = arma::sum(eigval) + nuclearRepulsionEnergy;
    std::cout << "The molecule in file " << inputFile << " has energy " << std::fixed << std::setprecision(6) << total_energy << std::endl;

    std::cout << "Eigenvalues of the overlap matrix S:" << std::endl;
    eigval.print();

    return 0;
}
