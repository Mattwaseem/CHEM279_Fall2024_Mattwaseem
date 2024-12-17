#include "InputParser.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
// Fxn to parse the input file for 1D numerical overlap integrals for question 1
std::tuple<double, double, int, double, double, int, int, int> InputParser::parse1Dinput(const std::string &filePath)
{
    std::ifstream file(filePath);
    double x1, exp1, x2, exp2;
    int l1, l2, ang1, ang2; // Corrected l1 and l2 to be int

    if (file.is_open())
    {
        // Reads the first line for Gaussian 1
        file >> x1 >> exp1 >> l1;

        // second line
        file >> x2 >> exp2 >> l2;

        file.close();
    }
    else
    {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        exit(1);
    }

    return std::make_tuple(x1, exp1, l1, x2, exp2, l2, 0.0, 0.0); // ang1, ang2); // Corrected return tuple to match variable types
}

// Fxn to parse the input file for 3D analytical overlap integrals for question 2
std::tuple<std::vector<double>, double, std::vector<int>,
           std::vector<double>, double, std::vector<int>>
InputParser::parse3DInput(const std::string &filePath)
{
    std::ifstream file(filePath);
    std::vector<double> center1(3), center2(3);
    double exp1, exp2;
    std::vector<int> ang1(3), ang2(3);

    if (file.is_open())
    {
        std::string line;

        // Read the first line for shell 1
        if (std::getline(file, line))
        {
            std::istringstream iss(line);
            iss >> center1[0] >> center1[1] >> center1[2] >> exp1 >> ang1[0] >> ang1[1] >> ang1[2];
            std::cout << "Parsed Shell 1: Center = (" << center1[0] << ", " << center1[1] << ", " << center1[2]
                      << "), Exponent = " << exp1 << ", Angular Momentum = (" << ang1[0] << ", " << ang1[1]
                      << ", " << ang1[2] << ")" << std::endl;
        }

        // Read the second line for shell 2
        if (std::getline(file, line))
        {
            std::istringstream iss(line);
            iss >> center2[0] >> center2[1] >> center2[2] >> exp2 >> ang2[0] >> ang2[1] >> ang2[2];
            std::cout << "Parsed Shell 2: Center = (" << center2[0] << ", " << center2[1] << ", " << center2[2]
                      << "), Exponent = " << exp2 << ", Angular Momentum = (" << ang2[0] << ", " << ang2[1]
                      << ", " << ang2[2] << ")" << std::endl;
        }

        file.close(); // Properly close the file after parsing
    }
    else
    {
        std::cerr << "Failed to open file: " << filePath << std::endl;
        exit(1);
    }

    // Return the tuple with values, not references
    return std::make_tuple(
        std::vector<double>(center1),
        exp1,
        std::vector<int>(ang1),
        std::vector<double>(center2),
        exp2,
        std::vector<int>(ang2));
}
