#include "CNDO.h"
#include <cassert>
#include <cmath>
#include <iostream>

const double AU_TO_EV = 27.2114;

CNDO::CNDO(Molecule molecule, arma::mat overlapMatrix)
    : molecule(molecule), overlapMatrix(overlapMatrix)
{
    int index = 0;
    for (int A = 0; A < molecule.nAtoms; A++)
    {
        int numAOs = calcNumAOs(molecule.atomicSymbols[A]);
        for (int i = 0; i < numAOs; i++)
        {
            aoIndexToAtom[index] = A;
            index++;
        }
    }

    diagCNDOPara["H"]["1s"] = 7.176;
    diagCNDOPara["C"]["2s"] = 14.051;
    diagCNDOPara["C"]["2px"] = 5.572;
    diagCNDOPara["C"]["2py"] = 5.572;
    diagCNDOPara["C"]["2pz"] = 5.572;
    diagCNDOPara["N"]["2s"] = 19.316;
    diagCNDOPara["N"]["2px"] = 7.275;
    diagCNDOPara["N"]["2py"] = 7.275;
    diagCNDOPara["N"]["2pz"] = 7.275;
    diagCNDOPara["O"]["2s"] = 25.390;
    diagCNDOPara["O"]["2px"] = 9.111;
    diagCNDOPara["O"]["2py"] = 9.111;
    diagCNDOPara["O"]["2pz"] = 9.111;
    diagCNDOPara["F"]["2s"] = 32.272;
    diagCNDOPara["F"]["2px"] = 11.080;
    diagCNDOPara["F"]["2py"] = 11.080;
    diagCNDOPara["F"]["2pz"] = 11.080;

    offDiagCNDOPara["H"] = 9.0;
    offDiagCNDOPara["C"] = 21.0;
    offDiagCNDOPara["N"] = 25.0;
    offDiagCNDOPara["O"] = 31.0;
    offDiagCNDOPara["F"] = 39.0;

    alphaCoeffMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaCoeffMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    alphaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    totalDensity = arma::zeros(molecule.nAtoms);

    alphaEnergy.set_size(molecule.nBasisFunctions);
    betaEnergy.set_size(molecule.nBasisFunctions);

    gammaMatrix = calcGammaMatrix();
    hCoreMat = calcHCoreMat();
    nuclearRepulsionEnergy = calcNuclearRepulsionEnergy();
}

int CNDO::calcNumAOs(std::string chemSym)
{
    if (chemSym == "H")
    {
        return 1;
    }
    else
    {
        return 4;
    }
}

double CNDO::calc2eIntegral(AO AO1, AO AO2)
{
    if (!(arma::accu(AO1.lmn) == 0 && arma::accu(AO2.lmn) == 0))
    {
        std::cerr << "Error: 2e integrals only implemented for s orbitals" << std::endl;
        return 0.0;
    }

    arma::vec dprime_a = AO1.contraction_coeffs % AO1.norm_constants;
    arma::vec dprime_b = AO2.contraction_coeffs % AO2.norm_constants;
    int len = AO1.exponents.n_elem;
    double gamma = 0.0;
    for (int k1 = 0; k1 < len; k1++)
    {
        for (int k2 = 0; k2 < len; k2++)
        {
            double sigmaA = 1.0 / (AO1.exponents(k1) + AO1.exponents(k2));
            for (int l1 = 0; l1 < len; l1++)
            {
                for (int l2 = 0; l2 < len; l2++)
                {
                    double sigmaB = 1.0 / (AO2.exponents(l1) + AO2.exponents(l2));
                    double I2e = pg2eIntegral(AO1.center, AO2.center, sigmaA, sigmaB);
                    gamma += dprime_a(k1) * dprime_a(k2) * dprime_b(l1) * dprime_b(l2) * I2e;
                }
            }
        }
    }
    return gamma;
}

double CNDO::pg2eIntegral(arma::rowvec center_a, arma::rowvec center_b, double sigmaA, double sigmaB)
{
    double U = pow(M_PI * sigmaA, 1.5) * pow(M_PI * sigmaB, 1.5);
    double V2 = 1.0 / (sigmaA + sigmaB);
    double distance = arma::norm(center_a - center_b, 2);
    if (distance == 0.0)
    {
        return U * sqrt(2 * V2) * sqrt(2 / M_PI);
    }
    double sqrtT = sqrt(V2) * distance;
    double result = U / distance * std::erf(sqrtT);
    return result;
}

arma::mat CNDO::calcGammaMatrix()
{
    std::vector<AO> sBasisFunctionsList;
    for (int i = 0; i < molecule.nBasisFunctions; i++)
    {
        if (molecule.basisFunctionsList[i].AO_type == "1s" || molecule.basisFunctionsList[i].AO_type == "2s")
        {
            sBasisFunctionsList.push_back(molecule.basisFunctionsList[i]);
        }
    }

    arma::mat gamma_matrix = arma::zeros<arma::mat>(molecule.nAtoms, molecule.nAtoms);
    for (size_t i = 0; i < sBasisFunctionsList.size(); i++)
    {
        for (size_t j = 0; j < sBasisFunctionsList.size(); j++)
        {
            gamma_matrix(i, j) = calc2eIntegral(sBasisFunctionsList[i], sBasisFunctionsList[j]);
        }
    }
    return gamma_matrix;
}

double CNDO::calcNuclearRepulsionEnergy()
{
    double nuclearRepulsionEnergy = 0.0;
    for (int A = 0; A < molecule.nAtoms; A++)
    {
        for (int B = 0; B < A; B++)
        {
            double distance = arma::norm(molecule.coordinates.row(A) - molecule.coordinates.row(B), 2);
            nuclearRepulsionEnergy += molecule.atomValences[A] * molecule.atomValences[B] / distance;
        }
    }
    return nuclearRepulsionEnergy;
}

arma::vec CNDO::calcTotalDensity()
{
    arma::vec totalDensity = arma::zeros(molecule.nAtoms);
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++)
    {
        int A = aoIndexToAtom[mu];
        totalDensity[A] += alphaDensityMat(mu, mu) + betaDensityMat(mu, mu);
    }
    return totalDensity;
}

arma::mat CNDO::calcFockMat(arma::mat densityMat)
{
    arma::mat fockMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++)
    {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++)
        {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            std::string chemSymA = molecule.atomicSymbols[A];
            std::string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double gammaAB = gammaMatrix(A, B);
            double pAA = totalDensity[A];
            double ZA = molecule.atomValences[A];
            if (mu == nu)
            {
                std::string AO_type = molecule.basisFunctionsList[mu].AO_type;
                fockMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] + ((pAA - ZA) - (densityMat(mu, mu) - 0.5)) * gammaAA;
                for (int C = 0; C < molecule.nAtoms; C++)
                {
                    if (A != C)
                    {
                        double pCC = totalDensity[C];
                        double ZC = molecule.atomValences[C];
                        double gammaAC = gammaMatrix(A, C);
                        fockMat(mu, nu) += (pCC - ZC) * gammaAC;
                    }
                }
            }
            else
            {
                fockMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) / 2.0 * overlapMatrix(mu, nu) - (densityMat(mu, nu) * gammaAB);
            }
        }
    }
    return fockMat;
}

arma::mat CNDO::calcHCoreMat()
{
    arma::mat hCoreMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++)
    {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++)
        {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            std::string chemSymA = molecule.atomicSymbols[A];
            std::string chemSymB = molecule.atomicSymbols[B];
            double gammaAA = gammaMatrix(A, A);
            double ZA = molecule.atomValences[A];
            if (mu == nu)
            {
                std::string AO_type = molecule.basisFunctionsList[mu].AO_type;
                hCoreMat(mu, nu) = -diagCNDOPara[chemSymA][AO_type] - (ZA - 0.5) * gammaAA;
                for (int C = 0; C < molecule.nAtoms; C++)
                {
                    if (A != C)
                    {
                        double ZC = molecule.atomValences[C];
                        double gammaAC = gammaMatrix(A, C);
                        hCoreMat(mu, nu) -= ZC * gammaAC;
                    }
                }
            }
            else
            {
                hCoreMat(mu, nu) = (-offDiagCNDOPara[chemSymA] - offDiagCNDOPara[chemSymB]) / 2.0 * overlapMatrix(mu, nu);
            }
        }
    }
    return hCoreMat;
}

arma::mat CNDO::calcDensityMat(arma::mat coeffMatA, std::string type)
{
    arma::mat densityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    int numElectrons = (type == "alpha") ? molecule.pAlpha : molecule.qBeta;
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++)
    {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++)
        {
            for (int i = 0; i < numElectrons; i++)
            {
                densityMat(mu, nu) += coeffMatA(mu, i) * coeffMatA(nu, i);
            }
        }
    }
    return densityMat;
}

double CNDO::calcTotalEnergy()
{
    double totalEnergy = 0.0;
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++)
    {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++)
        {
            totalEnergy += (alphaDensityMat(mu, nu) * (hCoreMat(mu, nu) + alphaFockMat(mu, nu)) +
                            betaDensityMat(mu, nu) * (hCoreMat(mu, nu) + betaFockMat(mu, nu)));
        }
    }
    totalEnergy /= 2.0;
    totalEnergy += nuclearRepulsionEnergy;
    return totalEnergy;
}

void CNDO::scfCycle()
{
    alphaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    betaDensityMat = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    bool converged = false;

    while (!converged)
    {
        totalDensity = calcTotalDensity();
        alphaFockMat = calcFockMat(alphaDensityMat);
        betaFockMat = calcFockMat(betaDensityMat);
        arma::eig_sym(alphaEnergy, alphaCoeffMat, alphaFockMat);
        arma::eig_sym(betaEnergy, betaCoeffMat, betaFockMat);
        arma::mat oldAlphaDensityMat = alphaDensityMat;
        arma::mat oldBetaDensityMat = betaDensityMat;
        alphaDensityMat = calcDensityMat(alphaCoeffMat, "alpha");
        betaDensityMat = calcDensityMat(betaCoeffMat, "beta");
        if (arma::abs(alphaDensityMat - oldAlphaDensityMat).max() < 1e-6 &&
            arma::abs(betaDensityMat - oldBetaDensityMat).max() < 1e-6)
        {
            converged = true;
            totalEnergy = calcTotalEnergy();
        }
    }
}

arma::mat CNDO::calcXMatrix()
{
    arma::mat X = arma::zeros(molecule.nBasisFunctions, molecule.nBasisFunctions);
    for (int mu = 0; mu < molecule.nBasisFunctions; mu++)
    {
        for (int nu = 0; nu < molecule.nBasisFunctions; nu++)
        {
            int A = aoIndexToAtom[mu];
            int B = aoIndexToAtom[nu];
            double betaA = offDiagCNDOPara[molecule.atomicSymbols[A]];
            double betaB = offDiagCNDOPara[molecule.atomicSymbols[B]];
            X(mu, nu) = 0.5 * (alphaDensityMat(mu, nu) + betaDensityMat(mu, nu)) * (betaA + betaB);
        }
    }
    return X;
}

arma::mat CNDO::calcYMatrix()
{
    arma::mat Y = arma::zeros(molecule.nAtoms, molecule.nAtoms);
    for (int A = 0; A < molecule.nAtoms; A++)
    {
        for (int B = 0; B < molecule.nAtoms; B++)
        {
            Y(A, B) = 0.5 * totalDensity[A] * totalDensity[B];
        }
    }
    return Y;
}

arma::mat CNDO::calcOverlapGradientMat()
{
    arma::mat overlapGradientMat;
    return overlapGradientMat;
}

arma::mat CNDO::calcGammaGradientMat()
{
    arma::mat gammaGradientMat;
    return gammaGradientMat;
}

arma::mat CNDO::calcNuclearEnergyGradMat()
{
    arma::mat dE_nuclear = arma::zeros(molecule.nAtoms, 3);
    for (int A = 0; A < molecule.nAtoms; A++)
    {
        for (int B = 0; B < A; B++)
        {
            arma::rowvec RA = molecule.coordinates.row(A);
            arma::rowvec RB = molecule.coordinates.row(B);
            arma::rowvec diff = RA - RB;
            double distance = arma::norm(diff, 2);
            if (distance > 1e-8)
            {
                arma::rowvec dE_dRA = (molecule.atomValences[A] * molecule.atomValences[B]) * (diff / (distance * distance * distance));
                arma::rowvec dE_dRB = -dE_dRA;
                dE_nuclear.row(A) += dE_dRA;
                dE_nuclear.row(B) += dE_dRB;
            }
        }
    }
    return dE_nuclear;
}

arma::mat CNDO::calcTotalEnergyGradMat()
{
    arma::mat dE_nuclear = calcNuclearEnergyGradMat();
    arma::mat dS = calcOverlapGradientMat();
    arma::mat dGamma = calcGammaGradientMat();
    arma::mat X = calcXMatrix();
    arma::mat Y = calcYMatrix();

    arma::mat gradient = arma::zeros(molecule.nAtoms, 3);

    return gradient;
}
