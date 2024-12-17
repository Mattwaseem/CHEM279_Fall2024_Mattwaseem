#include "file_io.h"
#include <fstream>
#include <stdexcept>

std::vector<Atom> read_atoms_from_file(const std::string &filename)
{
    std::ifstream infile(filename);
    if (!infile)
    {
        throw std::runtime_error("Error opening file.");
    }

    std::vector<Atom> atoms;
    int num_atoms;
    int atomic_number;
    double x, y, z;

    infile >> num_atoms; // This ignores the first line which is the number of atoms in the file based on the formatting of the txt files

    while (infile >> atomic_number >> x >> y >> z)
    {
        if (atomic_number != 79)
        {
            throw std::runtime_error("Error: only Au atoms are allowed.");
        }
        atoms.push_back(Atom{atomic_number, x, y, z});
    }
    return atoms;
}
