#pragma once
#include <memory>
#include <random>
#include <vector>
#include <chrono>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>

// Function to generate a cubic grid of coordinates
std::vector<std::tuple<int, int, int>> generate_grid(int num_points)
{
    // Calculate the side length of the cube
    int root = std::cbrt(num_points);
    int side_length = root;
    if (root * root * root != num_points)
        side_length += 1;

    // Generate all potential grid coordinates
    std::vector<std::tuple<int, int, int>> all_positions;
    for (int x = 0; x < side_length; ++x)
    {
        for (int y = 0; y < side_length; ++y)
        {
            for (int z = 0; z < side_length; ++z)
            {
                all_positions.push_back(std::make_tuple(x, y, z));
            }
        }
    }

    // Shuffle the coordinates randomly
    std::srand(static_cast<unsigned int>(std::time(0))); // Use current time as seed
    std::random_shuffle(all_positions.begin(), all_positions.end());

    // Trim the list to the desired number of points (num_points)
    if (all_positions.size() > num_points)
    {
        all_positions.resize(num_points);
    }

    return all_positions;
}

struct SpawnData
{
    /// @brief Array of XYZ coords. Size is 3xN where N is number of particles.
    std::vector<float> points;

    /// @brief Array of XYZ velocities. Size is 3xN where N is number of particles.
    std::vector<float> velocities;

    std::vector<float> rotations;

    std::vector<float> bond_angles;

    SpawnData(
        std::vector<float> &inPoints,
        std::vector<float> &inVels,
        std::vector<float> &inRots,
        std::vector<float> &inAng) : points(inPoints),
                                     velocities(inVels),
                                     rotations(inRots),
                                     bond_angles(inAng) {}
};

class Spawner3D
{
public:
    size_t num_parts;
    float size;
    std::vector<float> center; // Size 3
    float jitterStrength;
    float init_temp;
    float part_diff = 4.0f;
    float pos_scale = 10.0f;
    float vel_scale = 1e10f;

    /// @brief
    /// @param cen
    /// @param vel
    /// @param nP
    /// @param sz
    /// @param jS
    Spawner3D(const float (&cen)[3],
              size_t nP,
              float sz,
              float jS,
              float in_temp) : center(cen, cen + 3),
                               num_parts(nP), 
                               size(sz), 
                               jitterStrength(jS),
                               init_temp(in_temp) {}

    SpawnData GetSpawnData()
    {
        std::random_device rd;
        auto current_time = std::chrono::system_clock::now().time_since_epoch();
        auto seed = std::chrono::duration_cast<std::chrono::milliseconds>(current_time).count();
        std::mt19937 gen(seed);
        std::uniform_real_distribution<> dis(-size, size);
        std::vector<float> points = std::vector<float>(3 * num_parts);
        std::vector<float> velocities = std::vector<float>(3 * num_parts);
        std::vector<float> rotations = std::vector<float>(4 * num_parts);
        std::vector<float> bond_angles = std::vector<float>(num_parts);

        float vel = TempToVelocity(init_temp, num_parts) * vel_scale;
        std::vector<std::tuple<int, int, int>> grid_combos = generate_grid(num_parts);
        int root = std::cbrt(num_parts);
        int side_length = root;
        if (root * root * root != num_parts)
            side_length += 1;
        float world_size = (side_length + 1) * part_diff;
        float center = world_size *0.5f - part_diff;

        for (size_t i = 0; i < num_parts; i++)
        {
            std::tuple<int, int, int> pos_offset = grid_combos[i];
            points[i * 3] = std::get<0>(pos_offset)*part_diff + dis(gen) * jitterStrength - center;
            points[i * 3 + 1] = std::get<1>(pos_offset)*part_diff + dis(gen) * jitterStrength - center;
            points[i * 3 + 2] = std::get<2>(pos_offset)*part_diff + dis(gen) * jitterStrength - center;
            arma::fvec dir = arma::randu<arma::fvec>(3);
            dir = arma::normalise(2 * dir - 1.0f);
            velocities[i * 3] = dir(0) * vel;
            velocities[i * 3 + 1] = dir(1) * vel;
            velocities[i * 3 + 2] = dir(2) * vel;
            arma::fvec rot = arma::randu<arma::fvec>(4);
            rot = arma::normalise(2 * rot - 1.0f);
            rotations[i * 4] = rot(0);
            rotations[i * 4 + 1] = rot(1);
            rotations[i * 4 + 2] = rot(2);
            rotations[i * 4 + 3] = rot(3);
            bond_angles[i] = 109.5;
        }

        return SpawnData(points, velocities, rotations, bond_angles);
    }
};