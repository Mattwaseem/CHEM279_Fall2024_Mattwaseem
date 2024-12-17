#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include <filesystem>

struct Atom
{
    int atomic_number;
    double x, y, z;
};

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

    infile >> num_atoms; // ignores the first line which is the number of atoms in the system

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

double calculate_distance(const Atom &a1, const Atom &a2)
{
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

const double epsilon_au = 5.29; // kcal/mol
const double sigma_au = 2.951;  // Angstrom

double calculate_lj_energy(double distance)
{
    double ratio = sigma_au / distance;
    double ratio_6 = std::pow(ratio, 6);
    double ratio_12 = ratio_6 * ratio_6;
    return 4 * epsilon_au * (ratio_12 - ratio_6);
}

double calculate_total_energy(const std::vector<Atom> &atoms)
{
    double total_energy = 0.0;
    int n = atoms.size();
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double distance = calculate_distance(atoms[i], atoms[j]);
            total_energy += calculate_lj_energy(distance);
        }
    }
    return total_energy;
}

int main()
{
    try
    {
        std::string directory_path = "sample_input/Energy";

        // Loop through all files in the specified directory
        for (const auto &entry : std::filesystem::directory_iterator(directory_path))
        {
            std::string filename = entry.path().string();
            std::vector<Atom> atoms = read_atoms_from_file(filename);
            std::cout << "Processing file: " << filename << "\n";
            std::cout << "Atoms in the cluster:\n";

            for (const auto &atom : atoms)
            {
                std::cout << atom.atomic_number << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
            }

            double total_energy = calculate_total_energy(atoms);
            std::cout << "Total Lennard-Jones energy: " << total_energy << " kcal/mol\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}
