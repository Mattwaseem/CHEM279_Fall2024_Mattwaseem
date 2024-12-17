#include "Integrator1D.hpp"

//  area of a trapezoid is A = ((a + b) /(2)) * h  or (1/2 * a + 1/2 * b) * h
double Integrator1D::trapezoidalRule(const std::function<double(double)> &func, double a, double b, int n)
{
    double h = (b - a) / n;                      // Step size (width of each subinterval)
    double integral = 0.5 * (func(a) + func(b)); // Initialize with half the first and last points

    // Sum up the function values at the interior points
    for (int i = 1; i < n; ++i)
    {
        double x = a + i * h; // Current x value
        integral += func(x);  // Add the function value at this point
    }

    integral *= h; // Multiply by the step size to get the total integral value
    return integral;
}
