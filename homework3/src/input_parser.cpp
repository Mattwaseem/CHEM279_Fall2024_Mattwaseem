#include "input_parser.h"
#include <fstream>
#include <sstream>
#include <iostream>

std::tuple<int,
           int,
           std::vector<int>,
           std::vector<double>,
           std::vector<double>,
           std::vector<double>>
InputParser::parseMoleculeInput(const std::string &filePath)
{
    std::ifstream inputFile(filePath);
    if (!inputFile.is_open())
    {
        std::cerr << "Error opening input files: " << filePath << std::endl;
        return {};
    }

    int numAtoms, charge;
    inputFile >> numAtoms >> charge;

    std::vector<int> atomicNumbers;
    std::vector<double> xCoords;
    std::vector<double> yCoords;
    std::vector<double> zCoords;

    int atomicNumber;
    double x, y, z;

    for (int i = 0; i < numAtoms; ++i)
    {
        inputFile >> atomicNumber >> x >> y >> z;
        atomicNumbers.push_back(atomicNumber);
        xCoords.push_back(x);
        yCoords.push_back(y);
        zCoords.push_back(z);
    }

    inputFile.close();

    return std::make_tuple(numAtoms, charge, atomicNumbers, xCoords, yCoords, zCoords); // returning tuple
}
