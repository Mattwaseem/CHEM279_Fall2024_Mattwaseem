#include "factorial.h"

double factorial(int n)
{
    int result = 1;

    for (int i = n; i > 0; i--)
    {
        result *= i;
    }
    return result;
}

double binomialCoef(int m, int n)
{

    if (m < n)
    {
        return 0;
    }
    return factorial(m) / (factorial(n) * factorial(m - n));
}

double doubleFactorial(int n)
{
    int result = 1;

    // Multiply every other integer from n to 1
    for (int i = n; i > 0; i -= 2)
    {
        result *= i;
    }
    return result;
}
