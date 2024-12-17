#include "LJ_PE.h"
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include <iomanip> // so my terminal output matches the txt.output files

// Fxn to calculate the distance between two atoms
double calculate_distance(const Atom &a1, const Atom &a2)
{
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to calculate the Lennard-Jones energy between two atoms
std::tuple<double, double, double> calculate_lj_energy_components(double distance)
{
    if (distance == 0)
        return {0.0, 0.0, 0.0}; // Prevent division by zero
    double ratio = sigma_au / distance;
    double ratio_6 = std::pow(ratio, 6);
    double ratio_12 = ratio_6 * ratio_6;

    double energy = epsilon_au * (ratio_12 - 2 * ratio_6);
    return {energy, ratio_6, ratio_12};
}

// Function to calculate the total Lennard-Jones energy for a cluster of atoms
double calculate_total_energy(const std::vector<Atom> &atoms)
{
    double total_energy = 0.0;
    int n = atoms.size();
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double distance = calculate_distance(atoms[i], atoms[j]);
            auto [lj_energy, sigma_r_6, sigma_r_12] = calculate_lj_energy_components(distance);
            total_energy += lj_energy;
        }
    }
    return total_energy;
}
// -------------------------------------------------------- Problem 2 ------------------------------------------------------------------
// Fxn to calculate the force between two atoms
std::vector<double> calculate_lj_force(const Atom &a1, const Atom &a2)
{
    std::vector<double> force(3, 0.0);     // Force vector
    double r = calculate_distance(a1, a2); //  between the atoms

    if (r == 0)
        return force; // Prevent division by zero (was getting Nan values, maybe do this differently)

    // Reuses the energy fxn
    auto [energy, sigma_r_6, sigma_r_12] = calculate_lj_energy_components(r);

    // Calculates the magnitude of the force
    double force_magnitude = 24 * epsilon_au * (2 * sigma_r_12 / (r * r) - sigma_r_6 / (r * r));

    // Apply this force along the direction of r_ij (normalized vector - per Dr.Mayank's office hours)
    force[0] = force_magnitude * ((a1.x - a2.x) / r); // Fx
    force[1] = force_magnitude * ((a1.y - a2.y) / r); // Fy
    force[2] = force_magnitude * ((a1.z - a2.z) / r); // Fz

    return force;
}
// Numerical approx. for second derivative using central difference method for atom displacement
double second_derivative(const std::vector<Atom> &atoms, int atom_index, double h) // *****
{
    std::vector<Atom> atoms_plus_h = atoms;
    std::vector<Atom> atoms_minus_h = atoms;

    // chang in the atom's position by +/- h in the x-direction
    atoms_plus_h[atom_index].x += h;
    atoms_minus_h[atom_index].x -= h;

    // Calculate total energy for both perturbed configurations
    double energy_plus_h = calculate_total_energy(atoms_plus_h);
    double energy_minus_h = calculate_total_energy(atoms_minus_h);
    double energy_original = calculate_total_energy(atoms);

    // Use the central difference approximation for the second derivative
    return (energy_plus_h - 2 * energy_original + energy_minus_h) / (h * h);
}

// Fxn to calculate truncation error using
double truncation_error(double h, double second_derivative) // *****
{
    return 0.5 * h * std::abs(second_derivative); // *****
}

// Function to calculate round-off error
double round_off_error(double h, double function_value) // *****
{
    const double epsilon = 1e-16;                  // Machine precision for double precision // *****
    return epsilon * std::abs(function_value) / h; // *****
}

// Function to calculate total error
double total_error(double h, double function_value, double second_derivative) // *****
{
    double trunc_error = truncation_error(h, second_derivative); // *****
    double round_error = round_off_error(h, function_value);     // *****
    return trunc_error + round_error;                            // *****
}

// Numerical approximation for second derivative using central difference method
double second_derivative(double (*LJ_potential)(double), double x, double h) // *****
{
    return (LJ_potential(x + h) - 2 * LJ_potential(x) + LJ_potential(x - h)) / (h * h); // *****
}
//  implementing the forward, central diff force equations from lecture and lab slides

// Forward difference approximation for the force
double forward_difference_force(const std::vector<Atom> &atoms, int atom_index, double h)
{
    Atom original_atom = atoms[atom_index];
    std::vector<Atom> atoms_plus_h = atoms;
    atoms_plus_h[atom_index].x += h;
    double energy_plus_h = calculate_total_energy(atoms_plus_h);
    double energy_original = calculate_total_energy(atoms);
    return -(energy_plus_h - energy_original) / h;
}

// Central difference approximation for the force
double central_difference_force(const std::vector<Atom> &atoms, int atom_index, double h)
{
    Atom original_atom = atoms[atom_index];
    std::vector<Atom> atoms_plus_h = atoms;
    atoms_plus_h[atom_index].x += h;
    std::vector<Atom> atoms_minus_h = atoms;
    atoms_minus_h[atom_index].x -= h;
    double energy_plus_h = calculate_total_energy(atoms_plus_h);
    double energy_minus_h = calculate_total_energy(atoms_minus_h);

    std::cout << "Energy at x+h: " << energy_plus_h << "\n";
    std::cout << "Energy at x-h: " << energy_minus_h << "\n";
    std::cout << std::setprecision(15) << "Energy at x+h: " << energy_plus_h << "\n";  // debugging comment
    std::cout << std::setprecision(15) << "Energy at x-h: " << energy_minus_h << "\n"; // debugging comment

    return -(energy_plus_h - energy_minus_h) / (2 * h);
}

// Function to print forces in matrix format to verify against the expected output txt file
void print_forces_matrix(const std::vector<std::vector<double>> &forces)
{
    for (const auto &row : forces)
    {
        for (const auto &val : row)
        {
            std::cout << std::setw(10) << val << " ";
        }
        std::cout << "\n";
    }
}
// -------------------------------------------------------- Problem 3 -------------------------------------------------------------

// Line search fxn for steepest descent optimization
double line_search(const std::vector<Atom> &atoms, const std::vector<std::vector<double>> &forces)
{
    double best_alpha = STEP_SIZE;
    double best_energy = calculate_total_energy(atoms);
    for (double alpha = 0.001; alpha < 1.0; alpha += 0.01)
    {
        std::vector<Atom> test_positions = atoms;
        for (int i = 0; i < atoms.size(); ++i)
        {
            test_positions[i].x += alpha * forces[i][0];
            test_positions[i].y += alpha * forces[i][1];
            test_positions[i].z += alpha * forces[i][2];
        }
        double test_energy = calculate_total_energy(test_positions);
        if (test_energy < best_energy)
        {
            best_energy = test_energy;
            best_alpha = alpha;
        }
    }
    return best_alpha;
}

// Steepest descent optimization fxn

const double ENERGY_TOLERANCE = 1e-5; // Convergence tolerance for energy change
const int MAX_ITERATIONS = 1000;      // Maximum number of iterations for the optimizer

std::vector<Atom> steepest_descent_optimization(std::vector<Atom> &atoms)
{
    double prev_energy = calculate_total_energy(atoms);
    std::vector<Atom> new_positions = atoms;

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter)
    {
        std::cout << "Iteration " << iter + 1 << ":\n";

        // Calculate forces for each atom pair
        std::vector<std::vector<double>> forces(atoms.size(), std::vector<double>(3, 0.0));
        for (int i = 0; i < atoms.size(); ++i)
        {
            for (int j = 0; j < atoms.size(); ++j)
            {
                if (i == j)
                // to prevent double counting of distance for the same atom, if on the same atom skip to the next
                {
                    continue;
                }
                std::vector<double> force = calculate_lj_force(atoms[i], atoms[j]);
                for (int k = 0; k < 3; ++k)
                {
                    forces[i][k] += force[k];
                    // forces[j][k] -= force[k];
                    // don't need to account for negative since for loop condition change where the force application
                }
            }
        }

        // Perform line search to determine the best step size
        double alpha = line_search(atoms, forces);

        // Update atomic positions
        for (int i = 0; i < atoms.size(); ++i)
        {
            new_positions[i].x += alpha * forces[i][0];
            new_positions[i].y += alpha * forces[i][1];
            new_positions[i].z += alpha * forces[i][2];
        }

        double new_energy = calculate_total_energy(new_positions);

        std::cout << "Energy after iteration " << iter + 1 << ": " << new_energy << " kcal/mol\n";

        // Convergence check
        if (std::fabs(new_energy - prev_energy) < ENERGY_TOLERANCE)
        {
            std::cout << "Converged after " << iter + 1 << " iterations.\n";
            return new_positions;
        }

        prev_energy = new_energy;
        atoms = new_positions;
    }
    std::cout << "Reached maximum iterations without convergence.\n";
    return new_positions;
}
