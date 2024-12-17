#include <iostream>
#include <filesystem>
#include "file_io.h"
#include "LJ_PE.h"
#include <iomanip> // precision improvement
#include <fstream>
#include <cmath>

int main()
{
    // Open file to write truncation error results
    std::ofstream outfile("truncation_error_Problem2.txt");
    outfile << "h,Forward_Error,Central_Error\n"; // My CSV header

    try
    {
        std::string directory_path = "sample_input/Energy";
        double h_values[] = {0.1, 0.01, 0.001, 0.0001}; // Different step sizes for truncation error calculation based on hw PDF

        // Loop through all files in the specified directory
        for (const auto &entry : std::filesystem::directory_iterator(directory_path))
        {
            std::string filename = entry.path().string();
            std::vector<Atom> atoms = read_atoms_from_file(filename);
            std::cout << "Processing file: " << filename << "\n";
            std::cout << "Atoms in the cluster (Atomic Number and Coordinates in Matrix Format):\n";

            // Print coordinates in matrix format to match sample output for easy comparison of my answers
            for (const auto &atom : atoms)
            {
                std::cout << std::fixed << std::setprecision(4);
                std::cout << std::setw(10) << atom.atomic_number << std::setw(10) << atom.x << std::setw(10) << atom.y << std::setw(10) << atom.z << "\n";
            }

            // Calculate total Lennard-Jones energy
            double total_energy = calculate_total_energy(atoms);
            std::cout << "\nE_LJ = " << std::fixed << std::setprecision(5) << total_energy << " kcal/mol\n";

            // Loop over each step size to calculate truncation errors
            for (double h : h_values)
            {
                double forward_error_sum = 0.0;
                double central_error_sum = 0.0;
                double total_error_sum = 0.0; // for total error

                for (int i = 0; i < atoms.size(); ++i)
                {
                    // Analytical force
                    std::vector<double> analytical_force = calculate_lj_force(atoms[0], atoms[i]);

                    // Forward difference force for atom i
                    double forward_force = forward_difference_force(atoms, i, h);

                    // Central difference force for atom i
                    double central_force = central_difference_force(atoms, i, h);

                    // Calculate second derivative for truncation error calculation
                    double second_derivative_val = second_derivative(atoms, i, h); // Ensure this line is here

                    // Calculate truncation errors (sum of squared errors for x-component of force)
                    forward_error_sum += std::pow(analytical_force[0] - forward_force, 2);
                    central_error_sum += std::pow(analytical_force[0] - central_force, 2);

                    // Calculate total error
                    total_error_sum += total_error(h, calculate_total_energy(atoms), second_derivative_val); // Now passes total_energy and second_derivative_val
                }

                // Calculate average errors
                double forward_error = std::sqrt(forward_error_sum / atoms.size());
                double central_error = std::sqrt(central_error_sum / atoms.size());
                double total_error_avg = std::sqrt(total_error_sum / atoms.size());

                // Output the step size and errors
                outfile << h << "," << forward_error << "," << central_error << "," << total_error_avg << "\n";

                std::cout << "\nFor step size h = " << h << ":\n";
                std::cout << "  Forward difference error: " << forward_error << "\n";
                std::cout << "  Central difference error: " << central_error << "\n";
                std::cout << "  Total error: " << total_error_avg << "\n";
            }
            // Perform optimization using steepest descent
            std::cout << "\nPerforming optimization using steepest descent...\n";
            std::vector<Atom> optimized_atoms = steepest_descent_optimization(atoms);

            //  need to add line search also

            // Output optimized results
            std::cout << "\nOptimized Lennard-Jones energy: " << calculate_total_energy(optimized_atoms) << " kcal/mol\n";
            std::cout << "Optimized atomic positions:\n";
            for (const auto &atom : optimized_atoms)
            {
                std::cout << "Atom " << atom.atomic_number << " at (" << atom.x << ", " << atom.y << ", " << atom.z << ")\n";
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    outfile.close();
    return 0;
}
