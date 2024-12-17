#include "FileInputParser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include "overlapIntegrals.h"

Molecule FileInputParser::parseMolecule(const std::string &filename)
{
    Molecule molecule(filename.c_str());
    return molecule;
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
