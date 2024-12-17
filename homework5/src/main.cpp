// main.cpp
#include <iostream>
#include <fstream>
#include "Molecule.h"
#include "FileInputParser.h"
#include "CNDO.h"
#include "overlapIntegrals.h"

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input file>" << std::endl;
        return 1;
    }

    // Parse the input molecule
    FileInputParser parser;
    Molecule molecule = parser.parseMolecule(argv[1]);

    // Calculate the overlap matrix
    arma::mat overlapMatrix = calcOverlapMatrix(molecule);

    // Initialize CNDO with the molecule and overlap matrix
    CNDO cndo(molecule, overlapMatrix);

    // Perform the SCF cycle
    cndo.scfCycle();

    // Retrieve energies
    double nuclearRepulsionEnergy = cndo.calcNuclearRepulsionEnergy();
    double totalEnergy = cndo.calcTotalEnergy();
    double electronEnergy = totalEnergy - nuclearRepulsionEnergy;

    // Convert energies from Hartrees to eV
    double nuclearRepulsionEnergy_eV = nuclearRepulsionEnergy * 27.2114;
    double totalEnergy_eV = totalEnergy * 27.2114;
    double electronEnergy_eV = electronEnergy * 27.2114;

    // Print the energies
    std::cout << "Nuclear Repulsion Energy: " << nuclearRepulsionEnergy_eV << " eV." << std::endl;
    std::cout << "Electron Energy: " << electronEnergy_eV << " eV." << std::endl;
    std::cout << "Total Energy: " << totalEnergy_eV << " eV." << std::endl;

    // Calculate and print the energy gradient
    arma::mat gradient = cndo.calcTotalEnergyGradMat();

    std::cout << "Gradient:" << std::endl
              << gradient << std::endl;

    return 0;
}
