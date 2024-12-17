#include <napi.h>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <iostream>
#include <armadillo>
#include "utils.h"
#include "water_sim_3d.h"
#include "spawner3d.h"

// Store the global particles.
size_t num_particles;

/// @brief Float vector of size 3xN where N is number of particles and this holds XYZ.
static std::vector<float> particles;

/// @brief Float vector of size 3xN where N is number of particles and this holds XYZ.
static std::vector<float> velocities;

/// @brief Float vector of size 4xN where N is number of particles and this holds quaternions.
static std::vector<float> rotations;

/// @brief Float vector of size N where N is number of particles and this holds the angles.
static std::vector<float> bond_angles;

static float prevTime;
static float deltaTime;
static Sim3d water_sim;
float sim_scale;

/// @brief Represents a quaternian.
static arma::fvec rot_of_bounds = {0.0f, 0.0f, 0.0f, 1.0f};
// static arma::fvec rot_of_bounds = {0.92388f, 0.22068f, 0.22068f, 0.22068f};

/// @brief XYZ bounds of our bounding box.
static arma::fvec scale_of_bounds = {2.0f, 2.0f, 2.0f};

/// @brief XYZ position of the bounds
static arma::fvec pos_of_bounds = {0.0f, 0.0f, 0.0f};
arma::fmat local_to_world;
arma::fmat world_to_local;

bool ready_for_next_step = true;

Napi::Buffer<char> pack_data(Napi::Env &env)
{
    size_t pos_buff_size = num_particles * 3 * sizeof(float);
    size_t vel_buff_size = num_particles * 3 * sizeof(float);
    size_t rot_buff_size = num_particles * 4 * sizeof(float);
    size_t ang_buff_size = num_particles * 1 * sizeof(float);

    // Create a new combined buffer
    size_t total_size = pos_buff_size + vel_buff_size + rot_buff_size + ang_buff_size;
    Napi::Buffer<char> combined_buffer = Napi::Buffer<char>::New(env, total_size);

    // Copy data from each individual buffer into the combined buffer
    std::memcpy(combined_buffer.Data(), particles.data(), pos_buff_size);
    std::memcpy(combined_buffer.Data() + pos_buff_size, velocities.data(), vel_buff_size);
    std::memcpy(combined_buffer.Data() + pos_buff_size + vel_buff_size, rotations.data(), rot_buff_size);
    std::memcpy(combined_buffer.Data() + pos_buff_size + vel_buff_size + rot_buff_size, bond_angles.data(), ang_buff_size);

    return combined_buffer;
}

Napi::Buffer<char> SetupSimulation(const Napi::CallbackInfo &info)
{
    particles.clear();
    velocities.clear();
    rotations.clear();
    bond_angles.clear();
    Napi::Env env = info.Env();

    // Check that we have exactly five arguments
    if (info.Length() < 6 ||
        !info[0].IsNumber() || // num_part
        !info[1].IsNumber() || // sim_scale
        !info[2].IsArray() ||  // center (array of 3 floats)
        !info[3].IsNumber() || // initTemp (array of 3 floats)
        !info[4].IsNumber() || // size
        !info[5].IsNumber())   // jitterStrength
    {
        Napi::TypeError::New(env, "Expected five arguments: num_part, center, init_temp, size, jitterStrength")
            .ThrowAsJavaScriptException();
        return Napi::Buffer<char>::New(env, 0); // Return an empty array on error
    }

    // Extract num_P_per_axisPerAxis
    num_particles = info[0].As<Napi::Number>().Uint32Value();
    sim_scale = info[1].As<Napi::Number>().FloatValue();

    // Extract the 'center' array and ensure it has 3 elements
    Napi::Array centerArray = info[2].As<Napi::Array>();
    if (centerArray.Length() != 3)
    {
        Napi::TypeError::New(env, "center array must have exactly 3 elements").ThrowAsJavaScriptException();
        return Napi::Buffer<char>::New(env, 0); // Return empty array on error
    }

    // Extract the 'init_temp'
    float init_temp = info[3].As<Napi::Number>().Uint32Value();

    // Extract size and jitterStrength
    float size = info[4].As<Napi::Number>().FloatValue();
    float jitterStrength = info[5].As<Napi::Number>().FloatValue();

    float center[3];
    for (size_t i = 0; i < 3; i++)
    {
        center[i] = centerArray.Get(i).As<Napi::Number>().FloatValue();
    }

    // Create the Spawner3D object
    Spawner3D part_spawner(center, num_particles, size, jitterStrength, init_temp);
    SpawnData spawnData = part_spawner.GetSpawnData();

    // Store the particle coordinates in the particles vector
    for (size_t i = 0; i < num_particles; i++)
    {
        // Storing the x, y, z values directly in the particles vector
        particles.push_back(spawnData.points[i * 3]);     // x
        particles.push_back(spawnData.points[i * 3 + 1]); // y
        particles.push_back(spawnData.points[i * 3 + 2]); // z

        velocities.push_back(spawnData.velocities[i * 3]);     // x
        velocities.push_back(spawnData.velocities[i * 3 + 1]); // y
        velocities.push_back(spawnData.velocities[i * 3 + 2]); // z

        rotations.push_back(spawnData.rotations[i * 4]);     // x
        rotations.push_back(spawnData.rotations[i * 4 + 1]); // y
        rotations.push_back(spawnData.rotations[i * 4 + 2]); // z
        rotations.push_back(spawnData.rotations[i * 4 + 3]); // w

        bond_angles.push_back(spawnData.bond_angles[i]);
    }
    // std::cout << velocities[0] << " " << velocities[1] << " " << velocities[2] << std::endl;
    // Send points over as a buffer.
    water_sim = Sim3d(particles, velocities, rotations, bond_angles, num_particles);
    Napi::Buffer<char> buffer = pack_data(env);
    return buffer; // Return the array of points
}

/// @brief Thank G for my history with video game programming or this would be completely
/// impossible. It's essentially a direct rip from my own code in my personal 3D game engine
/// made in C++.
/// @param local_to_world
/// @param world_to_local
void UpdateModelMats()
{
    arma::fmat trans_mat = arma::eye<arma::fmat>(4, 4);
    for (int i = 0; i < 3; i++)
        trans_mat(i, 3) = pos_of_bounds(i);

    arma::fmat rot_mat = QuatToMatrix(rot_of_bounds[0], rot_of_bounds[1], rot_of_bounds[2], rot_of_bounds[3]);
    arma::fmat scale_mat = arma::eye<arma::fmat>(4, 4);
    for (int i = 0; i < 3; i++)
        scale_mat(i, i) = scale_of_bounds(i);

    local_to_world = trans_mat * rot_mat * scale_mat;
    world_to_local = arma::inv(local_to_world);
}

static std::chrono::high_resolution_clock::time_point last_time = std::chrono::high_resolution_clock::now(); // Initialize with the current time
Napi::Buffer<char> StepSim(const Napi::CallbackInfo &info)
{
    Napi::Env env = info.Env();
    Napi::Buffer<char> buffer;
    if (!ready_for_next_step)
    {
        buffer = Napi::Buffer<char>::New(env, 0);
        return buffer;
    }
    ready_for_next_step = false;

    // Get current time
    auto current_time = std::chrono::high_resolution_clock::now();
    std::chrono::nanoseconds nano_now = std::chrono::duration_cast<std::chrono::nanoseconds>(current_time - last_time);
    float time_diff = std::chrono::duration<float>(nano_now).count();
    // std::cout << time_diff << std::endl;
    water_sim.RunSimulationFrame(time_diff, local_to_world, world_to_local);
    last_time = current_time;
    particles = water_sim.positions;
    velocities = water_sim.velocities;
    rotations = water_sim.rotations;
    bond_angles = water_sim.bond_angles;

    buffer = pack_data(env);

    ready_for_next_step = true;
    return buffer;
}

Napi::Boolean SetBoxParams(const Napi::CallbackInfo &info)
{
    Napi::Env env = info.Env();

    // Check that we have exactly three arguments
    if (info.Length() < 3 ||
        !info[0].IsNumber() || // Width
        !info[1].IsNumber() || // Height
        !info[2].IsNumber())   // Depth
    {
        Napi::TypeError::New(env, "Expected five arguments: height, width, depth")
            .ThrowAsJavaScriptException();
        return Napi::Boolean::New(env, false); // Return false on error
    }
    float width = info[0].As<Napi::Number>().FloatValue();
    float height = info[1].As<Napi::Number>().FloatValue();
    float depth = info[2].As<Napi::Number>().FloatValue();
    scale_of_bounds = {width, height, depth};
    UpdateModelMats();
    return Napi::Boolean::New(env, true);
}

// Initialize the addon and expose the function to JavaScript
Napi::Object Init(Napi::Env env, Napi::Object exports)
{
    exports.Set("SetupSimulation", Napi::Function::New(env, SetupSimulation));
    exports.Set("StepSim", Napi::Function::New(env, StepSim));
    exports.Set("SetBoxParams", Napi::Function::New(env, SetBoxParams));
    // exports.Set("generateCoordinates", Napi::Function::New(env, GenerateCoordinates));
    // exports.Set("updateParticlePositions", Napi::Function::New(env, UpdateParticlePositions));
    return exports;
}

NODE_API_MODULE(addon, Init)
