#ifndef FILEINPUTPARSER_H
#define FILEINPUTPARSER_H

#include <string>
#include "Molecule.h"
#include "STO3GBasis.h"

class FileInputParser
{
public:
    Molecule parseMolecule(const std::string &filename);
    STO3GBasis parseBasisSet(const std::string &filename);
    std::string getBasisSetFilename(int atomicNumber);
    int getAtomicNumberFromSymbol(const std::string &element);

private:
};

#endif
