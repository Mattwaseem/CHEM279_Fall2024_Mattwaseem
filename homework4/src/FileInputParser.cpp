#include "FileInputParser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "STO3GBasis.h"

Molecule FileInputParser::parseMolecule(const std::string &filename)
{
    Molecule molecule;
    std::ifstream file(filename);
    if (file.is_open())
    {
        int numAtoms, charge;
        int atomicNumber;
        double x, y, z;
        std::string line;
        std::getline(file, line);
        std::istringstream headerStream(line);
        headerStream >> numAtoms >> charge;
        std::unordered_map<int, std::string> atomicNumberToSymbol = {
            {1, "H"},
            {6, "C"},
            {7, "N"},
            {8, "O"},
            {9, "F"}};

        for (int i = 0; i < numAtoms; ++i)
        {
            if (!std::getline(file, line))
            {
                std::cerr << "Error: Unexpected end of file while reading atom data." << std::endl;
                exit(1);
            }

            std::istringstream iss(line);
            if (iss >> atomicNumber >> x >> y >> z)
            {
                if (atomicNumber == 0)
                {
                    std::cerr << "Error: Atomic number is 0, which is invalid. Check input file format." << std::endl;
                    exit(1);
                }

                if (atomicNumberToSymbol.find(atomicNumber) != atomicNumberToSymbol.end())
                {
                    // Use the addAtom method instead of emplace_back
                    Atom atom = {atomicNumberToSymbol[atomicNumber], x, y, z};
                    molecule.addAtom(atom);
                }
                else
                {
                    std::cerr << "Error: Unsupported atomic number " << atomicNumber << " in input file." << std::endl;
                    exit(1);
                }
            }
            else
            {
                std::cerr << "Error: Failed to parse atomic data in line: " << line << std::endl;
                exit(1);
            }
        }
    }
    else
    {
        std::cerr << "Error: Could not open molecule file " << filename << std::endl;
        exit(1);
    }
    return molecule;
}

STO3GBasis FileInputParser::parseBasisSet(const std::string &filename)
{
    STO3GBasis basis;
    std::ifstream file(filename);
    if (file.is_open())
    {
        double exponent, coefficient;
        while (file >> exponent >> coefficient)
        {
            basis.exponents.push_back(exponent);
            basis.coefficients.push_back(coefficient);
        }
    }
    else
    {
        std::cerr << "Error: Could not open basis set file " << filename << std::endl;
        exit(1);
    }
    return basis;
}

std::string FileInputParser::getBasisSetFilename(int atomicNumber)
{
    switch (atomicNumber)
    {
    case 1:
        return "STO3G_info/H_STO3G.txt";
    case 6:
        return "STO3G_info/C_STO3G.txt";
    case 7:
        return "STO3G_info/N_STO3G.txt";
    case 8:
        return "STO3G_info/O_STO3G.txt";
    case 9:
        return "STO3G_info/F_STO3G.txt";
    default:
        std::cerr << "Error: Unsupported atomic number " << atomicNumber << std::endl;
        exit(1);
    }
}

int FileInputParser::getAtomicNumberFromSymbol(const std::string &element)
{
    static const std::unordered_map<std::string, int> elementToAtomicNumber = {
        {"H", 1},
        {"C", 6},
        {"N", 7},
        {"O", 8},
        {"F", 9}};

    auto it = elementToAtomicNumber.find(element);
    if (it != elementToAtomicNumber.end())
    {
        return it->second;
    }
    else
    {
        std::cerr << "Error: Unsupported element symbol \"" << element << "\" encountered." << std::endl;
        return 0;
    }
}
