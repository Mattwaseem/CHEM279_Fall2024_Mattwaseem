#ifndef CNDO_H
#define CNDO_H

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <armadillo>
#include "Molecule.h"

class CNDO
{
public:
    CNDO(Molecule molecule, arma::mat overlapMatrix);

    void scfCycle();
    arma::mat calcTotalEnergyGradMat();
    double calcNuclearRepulsionEnergy();
    double calcTotalEnergy();

private:
    Molecule molecule;
    arma::mat overlapMatrix;
    std::map<int, int> aoIndexToAtom;
    std::map<std::string, std::map<std::string, double>> diagCNDOPara;
    std::map<std::string, double> offDiagCNDOPara;
    arma::mat alphaCoeffMat;
    arma::mat betaCoeffMat;
    arma::mat alphaDensityMat;
    arma::mat betaDensityMat;
    arma::vec totalDensity;
    arma::vec alphaEnergy;
    arma::vec betaEnergy;
    arma::mat gammaMatrix;
    arma::mat alphaFockMat;
    arma::mat betaFockMat;
    arma::mat hCoreMat;
    double nuclearRepulsionEnergy;
    double totalEnergy;

    int calcNumAOs(std::string chemSym);
    double calc2eIntegral(AO AO1, AO AO2);
    double pg2eIntegral(arma::rowvec center_a, arma::rowvec center_b, double sigmaA, double sigmaB);
    arma::mat calcGammaMatrix();
    arma::vec calcTotalDensity();
    arma::mat calcFockMat(arma::mat densityMat);
    arma::mat calcHCoreMat();
    arma::mat calcDensityMat(arma::mat coeffMatA, std::string type);
    arma::mat calcXMatrix();
    arma::mat calcYMatrix();

    arma::mat calcOverlapGradientMat();
    arma::mat calcGammaGradientMat();
    arma::mat calcNuclearEnergyGradMat();
};

#endif // CNDO_H
