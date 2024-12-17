#include "DensityMatrix.h"
#include "FockMatrix.h"
#include "GammaCalculator.h"
#include "Constants.h"
#include "FileInputParser.h"
#include "OverlapMatrix.h"
#include <iostream>
#include <armadillo>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include "STO3GBasis.h"
#include "CartesianGaussian.h"

namespace fs = std::filesystem;

void processMolecule(const std::string &inputFilename)
{
    FileInputParser parser;
    Molecule molecule = parser.parseMolecule(inputFilename);
    const auto &atoms = molecule.getAtoms();

    std::unordered_map<int, STO3GBasis> basisSets;
    std::vector<int> atomicNumbersPerBasisFunction;
    Constants constants;

    int totalValenceElectrons = 0;
    for (const auto &atom : atoms)
    {
        int atomicNumber = parser.getAtomicNumberFromSymbol(atom.element);
        if (basisSets.find(atomicNumber) == basisSets.end())
        {
            std::string basisFilename = parser.getBasisSetFilename(atomicNumber);
            basisSets[atomicNumber] = parser.parseBasisSet(basisFilename);
        }
        totalValenceElectrons += constants.getValenceElectrons(atomicNumber);

        if (atom.element == "H")
        {
            atomicNumbersPerBasisFunction.push_back(atomicNumber);
        }
        else if (atom.element == "C" || atom.element == "N" || atom.element == "O" || atom.element == "F")
        {
            atomicNumbersPerBasisFunction.push_back(atomicNumber);
            atomicNumbersPerBasisFunction.push_back(atomicNumber);
            atomicNumbersPerBasisFunction.push_back(atomicNumber);
            atomicNumbersPerBasisFunction.push_back(atomicNumber);
        }
    }

    size_t p, q;
    if (inputFilename.find("H2") != std::string::npos)
    {
        p = 1;
        q = 1;
    }
    else if (inputFilename.find("HF") != std::string::npos)
    {
        p = 4;
        q = 4;
    }
    else if (inputFilename.find("HO") != std::string::npos)
    {
        p = 4;
        q = 3;
    }
    else if (inputFilename.find("N2") != std::string::npos)
    {
        p = 5;
        q = 5;
    }
    else if (inputFilename.find("O2") != std::string::npos)
    {
        p = 7;
        q = 5;
    }
    else
    {
        p = totalValenceElectrons / 2;
        q = totalValenceElectrons - p;
    }

    std::cout << "Processing Molecule: " << inputFilename << "\n";
    std::cout << "p = " << p << " q = " << q << "\n";

    molecule.computeBasisFunctions();
    int numBasisFunctions = molecule.getNumBasisFunctions();

    const std::vector<CartesianGaussian> &basisFunctions = molecule.getBasisFunctions();

    if (atomicNumbersPerBasisFunction.size() != numBasisFunctions)
    {
        std::cerr << "Error: Number of atomicNumbersPerBasisFunction (" << atomicNumbersPerBasisFunction.size()
                  << ") does not match numBasisFunctions (" << numBasisFunctions << ")." << std::endl;
        exit(1);
    }

    OverlapMatrix overlapMatrix(basisFunctions);
    overlapMatrix.computeOverlapMatrix();

    std::cout << "Overlap matrix calculated, dimensions: " << overlapMatrix.getMatrix().n_rows << " x " << overlapMatrix.getMatrix().n_cols << std::endl;
    std::cout << "Overlap Matrix S:\n"
              << overlapMatrix.getMatrix() << "\n";

    arma::mat P_alpha = arma::zeros<arma::mat>(numBasisFunctions, numBasisFunctions);
    arma::mat P_beta = arma::zeros<arma::mat>(numBasisFunctions, numBasisFunctions);

    DensityMatrix densityMatrix(P_alpha, P_beta, p, q);

    double nuclearRepulsionEnergy = 0.0;
    for (size_t A = 0; A < atoms.size(); ++A)
    {
        for (size_t B = A + 1; B < atoms.size(); ++B)
        {
            double ZA = constants.getValenceElectrons(parser.getAtomicNumberFromSymbol(atoms[A].element));
            double ZB = constants.getValenceElectrons(parser.getAtomicNumberFromSymbol(atoms[B].element));
            double distance = arma::norm(arma::vec({atoms[A].x, atoms[A].y, atoms[A].z}) - arma::vec({atoms[B].x, atoms[B].y, atoms[B].z}));
            nuclearRepulsionEnergy += (ZA * ZB) / distance;
        }
    }
    nuclearRepulsionEnergy *= constants.hartree_to_eV;

    std::vector<double> alphas;
    std::vector<double> d_total;

    for (const auto &basis : basisFunctions)
    {
        if (!basis.getExponents().empty())
        {
            alphas.push_back(basis.getExponents()[0]);
        }
        else
        {
            std::cerr << "Error: Basis function has no exponents.\n";
            exit(1);
        }
    }

    std::unordered_map<int, int> atomBasisCount;
    for (size_t mu = 0; mu < numBasisFunctions; ++mu)
    {
        int atomicNumber = atomicNumbersPerBasisFunction[mu];
        atomBasisCount[atomicNumber]++;
    }

    for (size_t mu = 0; mu < numBasisFunctions; ++mu)
    {
        int atomicNumber = atomicNumbersPerBasisFunction[mu];
        int valenceElectrons = constants.getValenceElectrons(atomicNumber);
        int numBasis = atomBasisCount[atomicNumber];

        if (numBasis == 0)
        {
            std::cerr << "Error: Number of basis functions for atomic number " << atomicNumber << " is zero.\n";
            exit(1);
        }

        double d = static_cast<double>(valenceElectrons) / numBasis;
        d_total.push_back(d);
    }

    arma::mat H_core = arma::zeros<arma::mat>(numBasisFunctions, numBasisFunctions);
    for (size_t mu = 0; mu < numBasisFunctions; ++mu)
    {
        for (size_t nu = 0; nu < numBasisFunctions; ++nu)
        {
            if (mu == nu)
            {
                int atomA = atomicNumbersPerBasisFunction[mu];
                double I_mu = constants.getIonizationPotential(atomA) * (1.0 / 27.211);
                double A_mu = constants.getElectronAffinity(atomA) * (1.0 / 27.211);
                double Z_A = constants.getValenceElectrons(atomA);
                double gammaAA = alphas[mu];
                H_core(mu, nu) = -0.5 * (I_mu + A_mu + gammaAA);
            }
            else
            {
                int atomA = atomicNumbersPerBasisFunction[mu];
                int atomB = atomicNumbersPerBasisFunction[nu];
                double betaA = constants.getBondingParameter(atomA) * (1.0 / 27.211);
                double betaB = constants.getBondingParameter(atomB) * (1.0 / 27.211);
                double s_mu_nu = overlapMatrix.getMatrix()(mu, nu);
                H_core(mu, nu) = -0.5 * (betaA + betaB) * s_mu_nu;
            }
        }
    }

    arma::mat H_core_eV = H_core * constants.hartree_to_eV;
    std::cout << "H_core\n"
              << H_core_eV << "\n";

    GammaCalculator gammaCalculator(constants);

    arma::mat gammaMatrix(numBasisFunctions, numBasisFunctions);
    for (size_t mu = 0; mu < numBasisFunctions; ++mu)
    {
        for (size_t nu = 0; nu < numBasisFunctions; ++nu)
        {
            int atomicNumber_mu = atomicNumbersPerBasisFunction[mu];
            int atomicNumber_nu = atomicNumbersPerBasisFunction[nu];
            arma::vec position_mu = basisFunctions[mu].getCenter();
            arma::vec position_nu = basisFunctions[nu].getCenter();

            double gammaAB = gammaCalculator.calculateGamma(mu, nu, alphas, d_total, atomicNumber_mu, atomicNumber_nu, position_mu, position_nu);
            gammaMatrix(mu, nu) = gammaAB * constants.hartree_to_eV;
        }
    }
    std::cout << "gamma\n"
              << gammaMatrix << "\n";

    FockMatrix fockMatrix(densityMatrix, overlapMatrix, H_core, atomicNumbersPerBasisFunction, alphas, d_total, basisFunctions, constants);

    double tolerance = 1e-6;
    size_t maxIterations = 100;
    bool converged = false;
    size_t iteration = 0;

    arma::mat fAlpha;
    arma::mat fBeta;

    std::ofstream outputFile;
    std::string outputDirectory = "calculated_outputs";
    if (!fs::exists(outputDirectory))
        fs::create_directory(outputDirectory);
    std::string outputFilename = outputDirectory + "/" + fs::path(inputFilename).stem().string() + "_calculated_output.txt";
    outputFile.open(outputFilename);

    if (outputFile.is_open())
    {
        outputFile << "Overlap Matrix S:\n"
                   << overlapMatrix.getMatrix() << "\n";
        outputFile << "gamma Matrix:\n"
                   << gammaMatrix << "\n";
        outputFile << "H_core Matrix (eV):\n"
                   << H_core_eV << "\n";
    }

    while (iteration < maxIterations && !converged)
    {
        std::cout << "Iteration: " << iteration << "\n";

        fAlpha = fockMatrix.computeFAlpha();
        fBeta = fockMatrix.computeFBeta();

        arma::mat fAlpha_eV = fAlpha * constants.hartree_to_eV;
        arma::mat fBeta_eV = fBeta * constants.hartree_to_eV;
        std::cout << "Fa\n"
                  << fAlpha_eV << "\n";
        std::cout << "Fb\n"
                  << fBeta_eV << "\n";

        if (!fAlpha.is_symmetric(1e-12))
        {
            std::cerr << "Warning: fAlpha matrix is not symmetric!\n";
        }

        arma::vec epsilon_alpha;
        arma::mat C_alpha;
        arma::eig_sym(epsilon_alpha, C_alpha, fAlpha);

        arma::vec epsilon_beta;
        arma::mat C_beta;
        arma::eig_sym(epsilon_beta, C_beta, fBeta);

        std::cout << "after solving eigen equation: " << iteration << "\n";
        std::cout << "Ca\n"
                  << C_alpha << "\n";
        std::cout << "Cb\n"
                  << C_beta << "\n";

        arma::mat P_alpha_new = C_alpha.cols(0, p - 1) * C_alpha.cols(0, p - 1).t();
        arma::mat P_beta_new = C_beta.cols(0, q - 1) * C_beta.cols(0, q - 1).t();

        std::cout << "Pa_new\n"
                  << P_alpha_new << "\n";
        std::cout << "Pb_new\n"
                  << P_beta_new << "\n";

        arma::mat P_total = P_alpha_new + P_beta_new;

        arma::vec P_t = arma::sum(P_total, 1);
        std::cout << "P_t\n"
                  << P_t << "\n";

        arma::mat Ga = fAlpha_eV;
        arma::mat Gb = fBeta_eV;
        std::cout << "Ga\n"
                  << Ga << "\n";
        std::cout << "Gb\n"
                  << Gb << "\n";

        double delta_alpha = arma::max(arma::abs(P_alpha_new - P_alpha)).max();
        double delta_beta = arma::max(arma::abs(P_beta_new - P_beta)).max();

        if (iteration == 0)
        {
            arma::vec epsilon_alpha_eV = epsilon_alpha * constants.hartree_to_eV;
            arma::vec epsilon_beta_eV = epsilon_beta * constants.hartree_to_eV;
            std::cout << "Ea\n"
                      << epsilon_alpha_eV << "\n";
            std::cout << "Eb\n"
                      << epsilon_beta_eV << "\n";
        }

        if (delta_alpha < tolerance && delta_beta < tolerance)
        {
            converged = true;
            std::cout << "SCF converged after " << iteration + 1 << " iterations.\n";
            break;
        }

        P_alpha = P_alpha_new;
        P_beta = P_beta_new;
        densityMatrix = DensityMatrix(P_alpha, P_beta, p, q);

        iteration++;
    }

    if (!converged)
    {
        std::cerr << "SCF did not converge after " << maxIterations << " iterations.\n";
    }

    arma::mat P_total_final = P_alpha + P_beta;
    double electronEnergy = 0.5 * arma::accu(P_total_final % (H_core + fAlpha + fBeta));
    double electronEnergy_eV = electronEnergy * constants.hartree_to_eV;
    double totalEnergy = electronEnergy_eV + nuclearRepulsionEnergy;

    std::cout << "Electron Energy is " << electronEnergy_eV << " eV.\n";
    std::cout << "Total Energy ECNDO/2 = " << totalEnergy << " eV.\n";
    std::cout << "Nuclear Repulsion Energy is " << nuclearRepulsionEnergy << " eV.\n";

    if (outputFile.is_open())
    {
        outputFile << "Final Fock Matrix Fa (eV):\n"
                   << fAlpha * constants.hartree_to_eV << "\n";
        outputFile << "Final Fock Matrix Fb (eV):\n"
                   << fBeta * constants.hartree_to_eV << "\n";
        outputFile << "Final Density Matrix P_alpha:\n"
                   << P_alpha << "\n";
        outputFile << "Final Density Matrix P_beta:\n"
                   << P_beta << "\n";
        outputFile << "Total Density Matrix P_total:\n"
                   << P_total_final << "\n";
        outputFile << "Electron Energy: " << electronEnergy_eV << " eV\n";
        outputFile << "Total Energy ECNDO/2: " << totalEnergy << " eV\n";
        outputFile << "Nuclear Repulsion Energy: " << nuclearRepulsionEnergy << " eV\n";
        outputFile.close();
        std::cout << "Results written to " << outputFilename << "\n";
    }
    else
    {
        std::cerr << "Error: Could not write to output file " << outputFilename << "\n";
    }
}

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file1> <input_file2> ... <input_fileN>\n";
        return 1;
    }
    for (int i = 1; i < argc; ++i)
    {
        std::string inputFilename = argv[i];
        processMolecule(inputFilename);
    }
    return 0;
}
