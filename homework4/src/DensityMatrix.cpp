#include "DensityMatrix.h"
#include <iostream>

DensityMatrix::DensityMatrix(const arma::mat &alphaCoeffs, const arma::mat &betaCoeffs, int p, int q)
    : alphaCoeffs_(alphaCoeffs), betaCoeffs_(betaCoeffs), p_(p), q_(q)
{
    std::cout << "Initializing Density Matrix:\n";
    std::cout << "Alpha Coefficients size: " << alphaCoeffs_.n_rows << " x " << alphaCoeffs_.n_cols << "\n";
    std::cout << "Beta Coefficients size: " << betaCoeffs_.n_rows << " x " << betaCoeffs_.n_cols << "\n";
    std::cout << "Number of alpha electrons (p): " << p_ << "\n";
    std::cout << "Number of beta electrons (q): " << q_ << "\n";

    if (p_ <= alphaCoeffs_.n_cols)
    {

        arma::mat C_alpha_p = alphaCoeffs_.cols(0, p_ - 1);
        P_alpha_ = C_alpha_p * C_alpha_p.t();
    }
    else
    {
        std::cerr << "Error: Number of alpha electrons p exceeds available coefficients.\n";
        exit(1);
    }

    if (q_ <= betaCoeffs_.n_cols)
    {
        // Occupy the first 'q' columns
        arma::mat C_beta_q = betaCoeffs_.cols(0, q_ - 1);
        P_beta_ = C_beta_q * C_beta_q.t();
    }
    else
    {
        std::cerr << "Error: Number of beta electrons q exceeds available coefficients.\n";
        exit(1);
    }

    // Compute total density matrix
    P_total_ = P_alpha_ + P_beta_;
}

const arma::mat &DensityMatrix::getP_alpha() const
{
    return P_alpha_;
}

const arma::mat &DensityMatrix::getP_beta() const
{
    return P_beta_;
}

const arma::mat &DensityMatrix::getP_total() const
{
    return P_total_;
}

void DensityMatrix::setP_alpha(const arma::mat &P_alpha)
{
    P_alpha_ = P_alpha;
}

void DensityMatrix::setP_beta(const arma::mat &P_beta)
{
    P_beta_ = P_beta;
}

arma::mat DensityMatrix::computeAlphaDensityMatrix() const
{
    return alphaCoeffs_ * alphaCoeffs_.t();
}

arma::mat DensityMatrix::computeBetaDensityMatrix() const
{
    return betaCoeffs_ * betaCoeffs_.t();
}

arma::mat DensityMatrix::computeTotalDensityMatrix() const
{
    return computeAlphaDensityMatrix() + computeBetaDensityMatrix();
}
