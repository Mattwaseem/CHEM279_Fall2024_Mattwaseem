#include <iostream>
#include "Gaussian1D.hpp" // Include the header file for the class declarations
#include <cmath>          // For mathematical operations

// Constructor
Gaussian::Gaussian(double x0_in, double alpha_in, int l_in)
    : x0_(x0_in), alpha_(alpha_in), l_(l_in) {}

// Evaluates the Gaussian function at point x
// double Gaussian::eval(double x)
// {
//     return std::pow(x - x0_, l_) * std::exp(-alpha_ * std::pow(x - x0_, 2));
// }
double Gaussian::eval(double x)
{
    double diff = x - x0_;
    double pow_term = std::pow(diff, l_);
    double exp_term = std::exp(-alpha_ * diff * diff);

    // Check for NaN or Inf
    if (std::isnan(pow_term) || std::isnan(exp_term) || std::isinf(pow_term) || std::isinf(exp_term))
    {
        std::cerr << "Error: NaN or Inf encountered in eval function!" << std::endl;
        std::cerr << "Parameters - x: " << x << ", x0_: " << x0_ << ", alpha_: " << alpha_ << ", l_: " << l_ << std::endl;
        std::cerr << "diff: " << diff << ", pow_term: " << pow_term << ", exp_term: " << exp_term << std::endl;
    }

    return pow_term * exp_term;
}

// Getters for Gaussian parameters
double Gaussian::get_x0() const { return x0_; }
double Gaussian::get_alpha() const { return alpha_; }
int Gaussian::get_l() const { return l_; }

// Sets parameters from parsed input
void Gaussian::setParameters(double x0_in, double alpha_in, int l_in)
{
    x0_ = x0_in;
    alpha_ = alpha_in;
    l_ = l_in;
}
