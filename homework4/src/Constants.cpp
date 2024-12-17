#include "Constants.h"
#include <unordered_map>
#include <iostream>
#include <stdexcept>

Constants::Constants()
{
    ionizationPotentials = {{1, 7.176}, {6, 14.051}, {7, 19.316}, {8, 25.390}, {9, 32.272}};
    electronAffinities = {{1, 5.572}, {6, 7.275}, {7, 9.111}, {8, 11.080}, {9, 3.401}};
    bondingParameters = {{1, 9.0}, {6, 21.0}, {7, 25.0}, {8, 31.0}, {9, 39.0}};
    valenceElectrons = {{1, 1}, {6, 4}, {7, 5}, {8, 6}, {9, 7}};
}

double Constants::getIonizationPotential(int atomicNumber) const
{
    auto it = ionizationPotentials.find(atomicNumber);
    if (it == ionizationPotentials.end())
    {
        std::cerr << "Ionization potential not found for atomic number: " << atomicNumber << std::endl;
        throw std::out_of_range("Ionization potential not found for atomic number: " + std::to_string(atomicNumber));
    }
    return it->second;
}

double Constants::getElectronAffinity(int atomicNumber) const
{
    auto it = electronAffinities.find(atomicNumber);
    if (it == electronAffinities.end())
    {
        std::cerr << "Electron affinity not found for atomic number: " << atomicNumber << std::endl;
        throw std::out_of_range("Electron affinity not found for atomic number: " + std::to_string(atomicNumber));
    }
    return it->second;
}

double Constants::getBondingParameter(int atomicNumber) const
{
    auto it = bondingParameters.find(atomicNumber);
    if (it == bondingParameters.end())
    {
        std::cerr << "Bonding parameter not found for atomic number: " << atomicNumber << std::endl;
        throw std::out_of_range("Bonding parameter not found for atomic number: " + std::to_string(atomicNumber));
    }
    return it->second;
}

double Constants::getValenceElectrons(int atomicNumber) const
{
    auto it = valenceElectrons.find(atomicNumber);
    if (it == valenceElectrons.end())
    {
        std::cerr << "Valence electrons not defined for atomic number: " << atomicNumber << std::endl;
        throw std::out_of_range("Valence electrons not defined for atomic number: " + std::to_string(atomicNumber));
    }
    return it->second;
}
