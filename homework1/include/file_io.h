#ifndef FILE_IO_H
#define FILE_IO_H

#include <vector>
#include <string>

struct Atom
{
    int atomic_number;
    double x, y, z;
};

// Fxn declaration to read atoms from a file
std::vector<Atom> read_atoms_from_file(const std::string &filename);

#endif // FILE_IO_H
