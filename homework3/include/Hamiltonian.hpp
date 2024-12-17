#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP
#pragma once

#include <vector>
#include <armadillo>

class Hamiltonian
{
public:
    // Constructor accepting the overlap matrix and diagonal energies
    Hamiltonian(const arma::mat &overlapMatrix, const std::vector<double> &diagEnergies);

    // Function to compute the Hamiltonian matrix
    void computeHamiltonianMatrix();

    // Getter function to return the Hamiltonian matrix (as a constant reference)
    const arma::mat &getMatrix() const;

private:
    arma::mat overlapMatrix_;          // Overlap matrix
    std::vector<double> diagEnergies_; // Diagonal energies for the Hamiltonian
    arma::mat hamiltonianMatrix_;      // Hamiltonian matrix
};

#endif
