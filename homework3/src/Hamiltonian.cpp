#include "Hamiltonian.hpp"
#include <iostream>
#include <cmath>
#include <armadillo>

const double K = 1.75;

Hamiltonian::Hamiltonian(const arma::mat &overlapMatrix, const std::vector<double> &diagEnergies)
    : overlapMatrix_(overlapMatrix), diagEnergies_(diagEnergies)
{
    size_t n = diagEnergies_.size();
    hamiltonianMatrix_.set_size(n, n);

    std::cout << "Constructor: Overlap matrix size: " << overlapMatrix_.n_rows << " x " << overlapMatrix_.n_cols << std::endl;
    std::cout << "Constructor: Diagonal energies size: " << diagEnergies_.size() << std::endl;
}

void Hamiltonian::computeHamiltonianMatrix()
{
    size_t n = diagEnergies_.size();

    if (overlapMatrix_.n_rows != n || overlapMatrix_.n_cols != n)
    {
        std::cerr << "Error: Overlap matrix dimensions do not match the size of the diagonal energies!" << std::endl;
        return;
    }

    double scalingFactor = 0.3;
    if (n == 2) // H2 molecule
    {
        scalingFactor = 0.175;
    }
    else if (n > 2) // Larger molec C2H2, C2H4
    {
        scalingFactor = 0.5;
    }

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            if (i == j)
            {
                hamiltonianMatrix_(i, j) = diagEnergies_[i];
            }
            else
            {
                double overlap = overlapMatrix_(i, j);
                double averageEnergy = scalingFactor * (diagEnergies_[i] + diagEnergies_[j]);

                double offDiagonalElement = K * averageEnergy * overlap;
                if (std::abs(offDiagonalElement) < 1e-6)
                {
                    offDiagonalElement = 0.0;
                }
                hamiltonianMatrix_(i, j) = offDiagonalElement;
            }
        }
    }

    std::cout << "Hamiltonian Matrix:" << std::endl;
    hamiltonianMatrix_.print();
}

const arma::mat &Hamiltonian::getMatrix() const
{
    return hamiltonianMatrix_;
}
