#pragma once
#include <string>
#include <vector>
#include <tuple>

class InputParser
{
public:
    static std::tuple<int,
                      int,
                      std::vector<int>,
                      std::vector<double>,
                      std::vector<double>,
                      std::vector<double>>
    parseMoleculeInput(const std::string &filePath);
};
