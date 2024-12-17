#ifndef LJ_PE_H
#define LJ_PE_H

#include <vector>
#include "file_io.h" // this includes my Atom structure

// Constants for gold (Au) atoms
const double epsilon_au = 5.29; // kcal/mol
const double sigma_au = 2.951;  // Angstrom

// Fxn declarations for Lennard-Jones energy and distance calculation
double calculate_distance(const Atom &a1, const Atom &a2);
std::tuple<double, double, double> calculate_lj_energy_components(double distance);
double calculate_total_energy(const std::vector<Atom> &atoms);

// Function to declareand calculate force
std::vector<double> calculate_lj_force(const Atom &a1, const Atom &a2);

// forward diff and central diff eqxn fxns
double forward_difference_force(const std::vector<Atom> &atoms, int atom_index, double h);
double central_difference_force(const std::vector<Atom> &atoms, int atom_index, double h);

// Function to calculate truncation error
double truncation_error(double h, double second_derivative);

// Function to calculate round-off error
double round_off_error(double h, double function_value);

// Function to calculate total error (sum of truncation and round-off errors)
double total_error(double h, double function_value, double second_derivative);

// Numerical approximation for second derivative (central difference method)
double second_derivative(const std::vector<Atom> &atoms, int atom_index, double h);

// Steepest descent optimizer
const double STEP_SIZE = 0.01;

std::vector<Atom> steepest_descent_optimization(std::vector<Atom> &atoms);
double line_search(const std::vector<Atom> &atoms, const std::vector<std::vector<double>> &forces);

#endif
