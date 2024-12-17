#pragma once
#include <string>
#include <vector>
#include <tuple>

class InputParser
{
public:
    // This will help me parse through the input files for appropriate numerical data points
    static std::tuple<double, double, int, double, double, int, int, int> parse1Dinput(const std::string &filePath);

    // This will parse for 3d overlap will be used for problem 2 data points
    // static std::tuple<std::vector<double>, std::vector<double>, double,
    // double,std::vector<int>, std::vector<int>> parse3dInput(const std::string &filePath);
    //---------------- corrected --------------
    // Function to parse the 3D input file for problem 2
    static std::tuple<std::vector<double>, double, std::vector<int>,
                      std::vector<double>, double, std::vector<int>>
    parse3DInput(const std::string &filePath);
};
