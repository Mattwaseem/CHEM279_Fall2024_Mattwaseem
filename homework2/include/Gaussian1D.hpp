#ifndef GAUSSIAN1D_HPP
#define GAUSSIAN1D_HPP

#include <cmath>

// Interfaces integrand fxn
class Integrand
{
public:
    virtual double eval(double x) = 0;    // Pure virtual fxn to evaluate at point x
    virtual double get_x0() const = 0;    // Gets center param.
    virtual double get_alpha() const = 0; // Gets the alpa of the guassian curve
    virtual int get_l() const = 0;        // Gets the angular momentum component
};

// Class Gaussian fxn derived from Integrand above
class Gaussian : public Integrand
{
private:
    double x0_;    // Center of the Gaussian
    double alpha_; // Alpha controls the spread of the guassian curve
    int l_;        // Angular momentum component

public:
    // Constructor with default values
    Gaussian(double x0_in = 0.0, double alpha_in = 1.0, int l_in = 0);

    // Override the eval function to calculate the value of the Gaussian at x
    double eval(double x) override;

    // Getters for Gaussian parameters
    double get_x0() const override;
    double get_alpha() const override;
    int get_l() const override;

    // Setter for parameters from parsed input
    void setParameters(double x0_in, double alpha_in, int l_in);
};

#endif // GAUSSIAN1D_HPP
