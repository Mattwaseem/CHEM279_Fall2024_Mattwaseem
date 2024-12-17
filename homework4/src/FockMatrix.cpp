#include "FockMatrix.h"
#include <iostream>

FockMatrix::FockMatrix(const DensityMatrix &densityMatrix,
                       const OverlapMatrix &overlapMatrix,
                       const arma::mat &H_core,
                       const std::vector<int> &atomicNumbersPerBasisFunction,
                       const std::vector<double> &alphas,
                       const std::vector<double> &d_total,
                       const std::vector<CartesianGaussian> &basisFunctions,
                       const Constants &constants)
    : densityMatrix_(densityMatrix),
      overlapMatrix_(overlapMatrix),
      H_core_(H_core),
      atomicNumbersPerBasisFunction_(atomicNumbersPerBasisFunction),
      alphas_(alphas),
      d_total_(d_total),
      basisFunctions_(basisFunctions),
      gammaCalculator_(constants),
      constants_(constants)
{
    initializeParameters();
}

void FockMatrix::initializeParameters()
{
    for (size_t mu = 0; mu < atomicNumbersPerBasisFunction_.size(); ++mu)
    {
        int atom = atomicNumbersPerBasisFunction_[mu];
        I_mu_map_[mu] = constants_.getIonizationPotential(atom) * (1.0 / 27.211);
        A_mu_map_[mu] = constants_.getElectronAffinity(atom) * (1.0 / 27.211);
        beta_A_map_[mu] = constants_.getBondingParameter(atom) * (1.0 / 27.211);
    }
}

arma::mat FockMatrix::computeFAlpha()
{
    arma::mat Fa = H_core_;
    arma::mat S = overlapMatrix_.getMatrix();
    arma::mat P_alpha = densityMatrix_.getP_alpha();
    arma::mat P_total = densityMatrix_.getP_total();

    size_t n = basisFunctions_.size();

    for (size_t mu = 0; mu < n; ++mu)
    {
        for (size_t nu = 0; nu < n; ++nu)
        {
            int atomicNumber_mu = atomicNumbersPerBasisFunction_[mu];
            int atomicNumber_nu = atomicNumbersPerBasisFunction_[nu];

            double gammaAB = gammaCalculator_.calculateGamma(mu, nu, alphas_, d_total_, atomicNumber_mu, atomicNumber_nu, basisFunctions_[mu].getCenter(), basisFunctions_[nu].getCenter());

            if (mu == nu)
            {
                double I_mu = I_mu_map_[mu];
                double A_mu = A_mu_map_[mu];
                double Z_A = constants_.getValenceElectrons(atomicNumber_mu);
                double p_AA = 0.0;
                for (size_t k = 0; k < n; ++k)
                {
                    if (atomicNumbersPerBasisFunction_[k] == atomicNumber_mu)
                        p_AA += P_total(k, k);
                }

                Fa(mu, mu) += -0.5 * (I_mu + A_mu) + (p_AA - Z_A) * gammaAB - P_alpha(mu, mu) * gammaAB - 0.5 * gammaAB * gammaAB;
            }
            else
            {
                double beta_A = beta_A_map_[mu];
                double beta_B = beta_A_map_[nu];
                Fa(mu, nu) += 0.5 * (-beta_A - beta_B) * S(mu, nu) - P_alpha(mu, nu) * gammaAB;
            }
        }
    }

    return Fa;
}

arma::mat FockMatrix::computeFBeta()
{
    arma::mat Fb = H_core_;
    arma::mat S = overlapMatrix_.getMatrix();
    arma::mat P_beta = densityMatrix_.getP_beta();
    arma::mat P_total = densityMatrix_.getP_total();

    size_t n = basisFunctions_.size();

    for (size_t mu = 0; mu < n; ++mu)
    {
        for (size_t nu = 0; nu < n; ++nu)
        {
            int atomicNumber_mu = atomicNumbersPerBasisFunction_[mu];
            int atomicNumber_nu = atomicNumbersPerBasisFunction_[nu];

            double gammaAB = gammaCalculator_.calculateGamma(mu, nu, alphas_, d_total_, atomicNumber_mu, atomicNumber_nu, basisFunctions_[mu].getCenter(), basisFunctions_[nu].getCenter());

            if (mu == nu)
            {
                double I_mu = I_mu_map_[mu];
                double A_mu = A_mu_map_[mu];
                double Z_A = constants_.getValenceElectrons(atomicNumber_mu);
                double p_AA = 0.0;
                for (size_t k = 0; k < n; ++k)
                {
                    if (atomicNumbersPerBasisFunction_[k] == atomicNumber_mu)
                        p_AA += P_total(k, k);
                }

                Fb(mu, mu) += -0.5 * (I_mu + A_mu) + (p_AA - Z_A) * gammaAB - P_beta(mu, mu) * gammaAB - 0.5 * gammaAB * gammaAB;
            }
            else
            {
                double beta_A = beta_A_map_[mu];
                double beta_B = beta_A_map_[nu];
                Fb(mu, nu) += 0.5 * (-beta_A - beta_B) * S(mu, nu) - P_beta(mu, nu) * gammaAB;
            }
        }
    }

    return Fb;
}
