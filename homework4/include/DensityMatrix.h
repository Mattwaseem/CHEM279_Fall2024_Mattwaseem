#ifndef DENSITYMATRIX_H
#define DENSITYMATRIX_H

#include <armadillo>

class DensityMatrix
{
public:
    DensityMatrix(const arma::mat &alphaCoeffs, const arma::mat &betaCoeffs, int p, int q);

    // Getter functions
    const arma::mat &getP_alpha() const;
    const arma::mat &getP_beta() const;
    const arma::mat &getP_total() const;

    // Setter functions (optional, if needed)
    void setP_alpha(const arma::mat &P_alpha);
    void setP_beta(const arma::mat &P_beta);

    arma::mat computeAlphaDensityMatrix() const;
    arma::mat computeBetaDensityMatrix() const;
    arma::mat computeTotalDensityMatrix() const;

private:
    arma::mat alphaCoeffs_;
    arma::mat betaCoeffs_;
    int p_; // Number of alpha electrons
    int q_; // Number of beta electrons

    arma::mat P_alpha_; // Alpha density matrix
    arma::mat P_beta_;  // Beta density matrix
    arma::mat P_total_;
};

#endif
