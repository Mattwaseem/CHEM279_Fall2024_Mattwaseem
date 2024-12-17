#ifndef INTEGRATOR1D_HPP
#define INTEGRATOR1D_HPP

#include <functional>

// Class for performing numerical integration
class Integrator1D
{
public:
    // Uses the trapezoidal rule for numerical integration
    static double trapezoidalRule(const std::function<double(double)> &func, double a, double b, int n);
};

#endif // INTEGRATOR1D_HPP
