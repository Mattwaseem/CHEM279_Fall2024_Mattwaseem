#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <unordered_map>
#include <vector>
#include <string>

struct AtomicParams
{
    double I_s;
    double I_p;
    double beta;
};

class Constants
{
public:
    Constants();
    const double hartree_to_eV = 27.211;
    const double angstrom_to_bohr = 1.88973;

    double getIonizationPotential(int atomicNumber) const;
    double getElectronAffinity(int atomicNumber) const;
    double getValenceElectrons(int atomicNumber) const;
    double getBondingParameter(int atomicNumber) const;

private:
    std::unordered_map<int, double> electronAffinities;
    std::unordered_map<int, double> ionizationPotentials;
    std::unordered_map<int, double> bondingParameters;
    std::unordered_map<int, double> valenceElectrons;
};

#endif
