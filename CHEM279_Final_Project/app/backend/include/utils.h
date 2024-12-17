#pragma once
#include <math.h>
#include <armadillo>
#include <iostream>
#include <omp.h>

using namespace std;

/// @brief Boltzman constant. Units kJ/K
static const float boltz_k = 1.380649e-23;

/// @brief 18 grams per mole / avagadros number * conv to kg. Units kg
static const float h2o_m = 2.99e-26;

// Define a custom reduction operator for arma::fvec
#pragma omp declare reduction(+ : arma::fvec : omp_out += omp_in) \
    initializer(omp_priv = arma::zeros<arma::fvec>(3))

float sign(float x)
{
    if (x > 0.0f)
        return 1.0f;
    else if (x < 0.0f)
        return -1.0f;
    else
        return 0.0f;
}

/// @brief Largely converted from my own personal code.
/// @param w First quaternian value.
/// @param x Second quaternian value.
/// @param y Third quaternian value.
/// @param z Fourth quaternian value.
/// @return 4x4 Arma matrix.
arma::fmat QuatToMatrix(float w, float x, float y, float z)
{
    // Create a 4x4 matrix
    arma::fmat rotMatrix = arma::zeros<arma::fmat>(4, 4);

    // Compute the elements of the rotation matrix
    rotMatrix(0, 0) = 1 - 2 * (y * y + z * z);
    rotMatrix(0, 1) = 2 * (x * y - z * w);
    rotMatrix(0, 2) = 2 * (x * z + y * w);

    rotMatrix(1, 0) = 2 * (x * y + z * w);
    rotMatrix(1, 1) = 1 - 2 * (x * x + z * z);
    rotMatrix(1, 2) = 2 * (y * z - x * w);

    rotMatrix(2, 0) = 2 * (x * z - y * w);
    rotMatrix(2, 1) = 2 * (y * z + x * w);
    rotMatrix(2, 2) = 1 - 2 * (x * x + y * y);

    // Set the last row and column for homogeneous coordinates
    rotMatrix(3, 3) = 1;

    return rotMatrix;
}

// 3d conversion: done
float SmoothingKernelPoly6(float dst, float radius)
{
    if (dst < radius)
    {
        float scale = 315.0f / (64.0f * M_PI * pow(abs(radius), 9.0));
        float v = radius * radius - dst * dst;
        return v * v * v * scale;
    }
    return 0;
}

// 3d conversion: done
float SM_PIkyKernelPow3(float dst, float radius)
{
    if (dst < radius)
    {
        float scale = 15.0f / (M_PI * pow(radius, 6));
        float v = radius - dst;
        return v * v * v * scale;
    }
    return 0;
}

// 3d conversion: done
// Integrate[(h-r)^2 r^2 Sin[θ], {r, 0, h}, {θ, 0, π}, {φ, 0, 2*π}]
float SM_PIkyKernelPow2(float dst, float radius)
{
    if (dst < radius)
    {
        float scale = 15.0f / (2.0f * M_PI * pow(radius, 5.0f));
        float v = radius - dst;
        return v * v * scale;
    }
    return 0;
}

// 3d conversion: done
float DerivativeSM_PIkyPow3(float dst, float radius)
{
    if (dst <= radius)
    {
        float scale = 45.0f / (pow(radius, 6) * M_PI);
        float v = radius - dst;
        return -v * v * scale;
    }
    return 0.0f;
}

// 3d conversion: done
float DerivativeSM_PIkyPow2(float dst, float radius)
{
    if (dst <= radius)
    {
        float scale = 15.0f / (pow(radius, 5) * M_PI);
        float v = radius - dst;
        return -v * scale;
    }
    return 0.0f;
}

float DensityKernel(float dst, float radius)
{
    // return SmoothingKernelPoly6(dst, radius);
    return SM_PIkyKernelPow2(dst, radius);
}

float NearDensityKernel(float dst, float radius)
{
    return SM_PIkyKernelPow3(dst, radius);
}

float DensityDerivative(float dst, float radius)
{
    return DerivativeSM_PIkyPow2(dst, radius);
}

float NearDensityDerivative(float dst, float radius)
{
    return DerivativeSM_PIkyPow3(dst, radius);
}

float TempToVelocity(float temp, size_t num_part){
    // newtonian kinetic energy. in meters/sec multiplied by 1e10 to convert to angstroms.
    return std::sqrt((3.0f * temp * boltz_k) / h2o_m) * 1e10;
}

// Quaternion multiplication using arma::fvec
arma::fvec quat_mult(const arma::fvec &q1, const arma::fvec &q2)
{
    // q1 = (x1, y1, z1, w1)
    // q2 = (x2, y2, z2, w2)
    float x1 = q1(0), y1 = q1(1), z1 = q1(2), w1 = q1(3);
    float x2 = q2(0), y2 = q2(1), z2 = q2(2), w2 = q2(3);

    // Perform quaternion multiplication
    float x3 = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2;
    float y3 = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2;
    float z3 = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2;
    float w3 = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2;

    return arma::fvec({x3, y3, z3, w3});
}

// Quaternion conjugate
arma::fvec quat_conj(const arma::fvec &q)
{
    // Conjugate of quaternion (x, y, z, w) is (-x, -y, -z, w)
    return arma::fvec({-q(0), -q(1), -q(2), q(3)});
}

/// @brief Defined as q x r x q* where q* is the conjugate of the quaternian.
/// @param r Vector to be rotated.
/// @param q Quaternian rotating vector.
/// @return Rotated vector.
arma::fvec rot_vec_by_quat(arma::fvec r, arma::fvec q)
{
    arma::fmat r_q = arma::fvec({r(0), r(1), r(2), 0.0f});
    arma::fvec rot_q = quat_mult(q, r_q);
    rot_q = quat_mult(rot_q, quat_conj(q));
    return arma::fvec({rot_q(0), rot_q(1), rot_q(2)});
}