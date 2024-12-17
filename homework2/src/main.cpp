#include "InputParser.hpp"
#include "CartesianGaussian.hpp"
#include "Gaussian1D.hpp"
#include "Integrator1D.hpp"
#include "OverlapIntegral.hpp"
#include <iostream>
#include <functional>
#include <vector>
#include <iomanip> // For formatted output

int main(int argc, char *argv[])
{

    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
        return 1;
    }

    std::string inputFile = argv[1];

    if (inputFile.find("numerical") != std::string::npos)
    {
        // **Problem 1: Numerical 1D Overlap Integral**
        auto [x1, alpha1, l1, x2, alpha2, l2, ang1, ang2] = InputParser::parse1Dinput(inputFile);
        std::cout << "Parsed parameters (Numerical 1D):" << std::endl;
        std::cout << "Gaussian 1: x0 = " << x1 << ", alpha = " << alpha1 << ", l = " << l1 << std::endl;
        std::cout << "Gaussian 2: x0 = " << x2 << ", alpha = " << alpha2 << ", l = " << l2 << std::endl;

        Gaussian g1(x1, alpha1, l1);
        Gaussian g2(x2, alpha2, l2);

        auto integrand = [&g1, &g2](double x)
        {
            return g1.eval(x) * g2.eval(x);
        };

        double a = -10.0, b = 10.0;
        int n = 10000;

        double result = Integrator1D::trapezoidalRule(integrand, a, b, n);
        std::cout << "Overlap integral result (Numerical 1D): " << result << std::endl;
    }
    else if (inputFile.find("analytical") != std::string::npos)
    {
        // **Problem 2: Analytical 3D Overlap Integral**
        auto parsedData = InputParser::parse3DInput(inputFile);
        std::vector<double> center1 = std::get<0>(parsedData);
        double alpha1 = std::get<1>(parsedData);
        std::vector<int> ang1 = std::get<2>(parsedData);

        std::vector<double> center2 = std::get<3>(parsedData);
        double alpha2 = std::get<4>(parsedData);
        std::vector<int> ang2 = std::get<5>(parsedData);

        // Convert std::vector to arma::vec and arma::ivec
        arma::vec center1_arma(center1);
        arma::vec center2_arma(center2);
        arma::ivec ang1_arma(ang1.size());
        arma::ivec ang2_arma(ang2.size());

        for (size_t i = 0; i < ang1.size(); ++i)
        {
            ang1_arma[i] = ang1[i];
        }
        for (size_t i = 0; i < ang2.size(); ++i)
        {
            ang2_arma[i] = ang2[i];
        }

        std::cout << "Parsed parameters (Analytical 3D):" << std::endl;

        std::cout << "Shell 1 has " << ang1_arma.n_elem << " functions." << std::endl;
        std::cout << "This shell info: R( " << center1_arma[0] << ", " << center1_arma[1] << ", " << center1_arma[2] << "), ";
        std::cout << "with angular momentum: " << ang1_arma[0] << ", coefficient: " << alpha1 << std::endl;

        std::cout << "Shell 2 has " << ang2_arma.n_elem << " functions." << std::endl;
        std::cout << "This shell info: R( " << center2_arma[0] << ", " << center2_arma[1] << ", " << center2_arma[2] << "), ";
        std::cout << "with angular momentum: " << ang2_arma[0] << ", coefficient: " << alpha2 << std::endl;

        CartesianGaussian g1(center1_arma, alpha1, ang1_arma);
        CartesianGaussian g2(center2_arma, alpha2, ang2_arma);

        // Compute overlap
        double Sx = OverlapIntegral::computeOverlap1D(g1, g2, 0);
        double Sy = OverlapIntegral::computeOverlap1D(g1, g2, 1);
        double Sz = OverlapIntegral::computeOverlap1D(g1, g2, 2);

        if (ang1_arma[0] == 0 && ang2_arma[0] == 0)
        {
            // Case: Only 1 function per shell, print a single overlap value
            double totalOverlap = Sx * Sy * Sz;
            std::cout << "Overlap integral between Shell 1 and Shell 2" << std::endl;
            std::cout << "   " << std::fixed << std::setprecision(4) << totalOverlap << std::endl;
        }
        else
        {
            // Case: Multiple functions, print as a matrix
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "Overlap integral between Shell 1 and Shell 2" << std::endl;
            std::cout << "  " << Sx << "  " << Sy << "  " << Sz << std::endl;

            // Assuming multiple rows and columns based on angular momentum components.
            std::cout << "  " << Sx << "  " << Sy << "       0" << std::endl;
            std::cout << "  " << Sy << "  " << Sx << "       0" << std::endl;
            std::cout << "       0       0  " << Sz << std::endl;
        }

        // Print angular momentum components for matrix column and row
        std::cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, ";
        std::cout << "are listed sequentially as: (" << ang1_arma[0] << ", " << ang1_arma[1] << ", " << ang1_arma[2] << ")." << std::endl;
        std::cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, ";
        std::cout << "are listed sequentially as: (" << ang2_arma[0] << ", " << ang2_arma[1] << ", " << ang2_arma[2] << ")." << std::endl;
    }
    else
    {
        std::cerr << "Error: Unrecognized input file type. Please use either a numerical (1D) or analytical (3D) input file." << std::endl;
        return 1;
    }

    return 0;
}

// ----------------------------- updated 1--------------------------------------------------------------------------------
// #include "InputParser.hpp"
// #include "CartesianGaussian.hpp"
// #include "Gaussian1D.hpp"
// #include "Integrator1D.hpp"
// #include "OverlapIntegral.hpp"
// #include <iostream>
// #include <functional>
// #include <vector>

// int main(int argc, char *argv[])
// {

//     if (argc != 2)
//     {
//         std::cerr << "Usage: " << argv[0] << " <input_file>" << std::endl;
//         return 1;
//     }

//     std::string inputFile = argv[1];

//     // Detect input type: 1D numerical (Problem 1) or 3D analytical (Problem 2)
//     if (inputFile.find("numerical") != std::string::npos)
//     {
//         // **Problem 1: Numerical 1D Overlap Integral**
//         auto [x1, alpha1, l1, x2, alpha2, l2, ang1, ang2] = InputParser::parse1Dinput(inputFile);
//         std::cout << "Parsed parameters (Numerical 1D):" << std::endl;
//         std::cout << "Gaussian 1: x0 = " << x1 << ", alpha = " << alpha1 << ", l = " << l1 << std::endl;
//         std::cout << "Gaussian 2: x0 = " << x2 << ", alpha = " << alpha2 << ", l = " << l2 << std::endl;

//         Gaussian g1(x1, alpha1, l1);
//         Gaussian g2(x2, alpha2, l2);

//         auto integrand = [&g1, &g2](double x)
//         {
//             return g1.eval(x) * g2.eval(x);
//         };

//         double a = -10.0, b = 10.0;
//         int n = 10000;

//         double result = Integrator1D::trapezoidalRule(integrand, a, b, n);
//         std::cout << "Overlap integral result (Numerical 1D): " << result << std::endl;
//     }
//     else if (inputFile.find("analytical") != std::string::npos)
//     {
//         // **Problem 2: Analytical 3D Overlap Integral**
//         auto parsedData = InputParser::parse3DInput(inputFile);
//         std::vector<double> center1 = std::get<0>(parsedData);
//         double alpha1 = std::get<1>(parsedData);
//         std::vector<int> ang1 = std::get<2>(parsedData);

//         std::vector<double> center2 = std::get<3>(parsedData);
//         double alpha2 = std::get<4>(parsedData);
//         std::vector<int> ang2 = std::get<5>(parsedData);

//         // Convert std::vector to arma::vec and arma::ivec
//         arma::vec center1_arma(center1);
//         arma::vec center2_arma(center2);
//         arma::ivec ang1_arma(ang1.size());
//         arma::ivec ang2_arma(ang2.size());

//         for (size_t i = 0; i < ang1.size(); ++i)
//         {
//             ang1_arma[i] = ang1[i];
//         }
//         for (size_t i = 0; i < ang2.size(); ++i)
//         {
//             ang2_arma[i] = ang2[i];
//         }

//         std::cout << "Parsed parameters (Analytical 3D):" << std::endl;

//         // // Print Gaussian 1
//         // std::cout << "Gaussian 1: Center = (";
//         // for (size_t i = 0; i < center1_arma.n_elem; ++i)
//         // {
//         //     std::cout << center1_arma[i] << (i < center1_arma.n_elem - 1 ? ", " : "");
//         // }
//         // std::cout << "), alpha = " << alpha1 << ", angular momentum = (";
//         // for (size_t i = 0; i < ang1_arma.n_elem; ++i)
//         // {
//         //     std::cout << ang1_arma[i] << (i < ang1_arma.n_elem - 1 ? ", " : "");
//         // }
//         // std::cout << ")\n";
//         std::cout << "Shell 1 has 1 functions." << std::endl;
//         std::cout << "This shell info: R( " << center1_arma[0] << ", " << center1_arma[1] << ", " << center1_arma[2] << "), ";
//         std::cout << "with angular momentum: " << ang1_arma[0] << ", coefficient: " << alpha1 << std::endl;

//         std::cout << "Shell 2 has 1 functions." << std::endl;
//         std::cout << "This shell info: R( " << center2_arma[0] << ", " << center2_arma[1] << ", " << center2_arma[2] << "), ";
//         std::cout << "with angular momentum: " << ang2_arma[0] << ", coefficient: " << alpha2 << std::endl;

//         //     // Print Gaussian 2
//         //     std::cout << "Gaussian 2: Center = (";
//         //     for (size_t i = 0; i < center2_arma.n_elem; ++i)
//         //     {
//         //         std::cout << center2_arma[i] << (i < center2_arma.n_elem - 1 ? ", " : "");
//         //     }
//         //     std::cout << "), alpha = " << alpha2 << ", angular momentum = (";
//         //     for (size_t i = 0; i < ang2_arma.n_elem; ++i)
//         //     {
//         //         std::cout << ang2_arma[i] << (i < ang2_arma.n_elem - 1 ? ", " : "");
//         //     }
//         //     std::cout << ")\n";

//         //     CartesianGaussian g1(center1_arma, alpha1, ang1_arma);
//         //     CartesianGaussian g2(center2_arma, alpha2, ang2_arma);

//         //     double totalOverlap = OverlapIntegral::computeTotalOverlap(g1, g2);
//         //     std::cout << "Total overlap integral result (Analytical 3D): " << totalOverlap << std::endl;
//         // }
//         CartesianGaussian g1(center1_arma, alpha1, ang1_arma);
//         CartesianGaussian g2(center2_arma, alpha2, ang2_arma);

//         double totalOverlap = OverlapIntegral::computeTotalOverlap(g1, g2);
//         std::cout << "Overlap integral between Shell 1 and Shell 2" << std::endl;
//         std::cout << "   " << totalOverlap << std::endl;

//         // Print angular momentum components for matrix column and row
//         std::cout << "The components of angular momentum (l, m, n) for the matrix column, from top to bottom, ";
//         std::cout << "are listed sequentially as: (" << ang1_arma[0] << ", " << ang1_arma[1] << ", " << ang1_arma[2] << ")." << std::endl;
//         std::cout << "The components of angular momentum (l, m, n) for the matrix row, from left to right, ";
//         std::cout << "are listed sequentially as: (" << ang2_arma[0] << ", " << ang2_arma[1] << ", " << ang2_arma[2] << ")." << std::endl;
//     }
//     else
//     {
//         std::cerr << "Error: Unrecognized input file type. Please use either a numerical (1D) or analytical (3D) input file." << std::endl;
//         return 1;
//     }

//     // Test cases for double factorial and normalization constant for debugging.
//     arma::vec center = {0.0, 0.0, 0.0};
//     arma::ivec angularMomentum = {2, 0, 0};
//     double exponent = 1.0;

//     CartesianGaussian gaussian(center, exponent, angularMomentum);

//     return 0;
// }
