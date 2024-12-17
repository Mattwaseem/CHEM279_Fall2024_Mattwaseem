#pragma once
#include <math.h>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <armadillo>
#include "utils.h"
#include "spawner3d.h"
#include <omp.h>
#include <ctime>

const double k_e = 8.99e23;                   // Coulomb's constant in kJ*A^2*C^-2
// const double k_e = 8.99e29;                   // Coulomb's constant in kJ*A^2*C^-2
static const float sigma_OO = 3.166f;         // Angstroms
static const float epsilon_OO_kcal = 0.6502;  // kcal/mol
static const float epsilon_OO_J = 4.52e-21;   // J
static const float e_OO = 0.00673938f;        // Converted of above for one molecule at eV
static const float q_H2O = 12.816e-18f;       // 10 total electrons? [Matt review]
static const float e_charge = 1.602e-19;
static const float q_O = -0.834 * e_charge;   // Partial charge on Oxygen
static const float q_H = 0.417 * e_charge;    // Partial charge on Hydrogen
static const float hr_104 = 0.957;            // Radial arm length of H to O when at 104.52 bond angle.
static const float hr_109 = 1.0;              // Radial arm length of H to O when at 109.5 bond angle.
static const float mass_H = 1.66e-24f;
static const float time_delta = 2e-15;        // femtosecond
static std::ofstream outFile("../../frontend/err_logs/log.txt", std::ios::app);
static time_t now;

/// @brief Quaternion representing the X axis forward.
static const arma::fvec forward_q = arma::fvec({1.0f, 0.0f, 0.0f, 0.0f});

class Sim3d
{
private:
    float collision_damping = 0.3f;
    std::vector<float> predicted_positions; // 3 x N where N is number of particles.
    int num_particles;                      // N where N is number of particles.

public:
    std::vector<float> positions;   // 3 x N where N is number of particles.
    std::vector<float> velocities;  // 3 x N where N is number of particles.
    std::vector<float> rotations;   // 3 x N where N is number of particles.
    std::vector<float> ang_vels;    // 3 x N where N is number of particles.
    std::vector<float> bond_angles; // 3 x N where N is number of particles.
    std::vector<float> energies;      // N where N is number of particles.
    std::vector<float> forces;      // 3 x N where N is number of particles.
    std::vector<float> forces_coul; // 3 x N where N is number of particles.
    std::vector<float> torques;     // 3 x N where N is number of particles.
    std::vector<float> densities;   // 2 x N where N is the number of particles.
    // event System.Action SimulationStepCompleted; <- Will be interesting
    float time_scale = 1e-10f;
    int iter_per_frame = 1;
    // float gravity = -9.81;
    float gravity = -0.0f;
    void set_collision_damping(float val)
    {
        if (val < 0.0f | val > 1.0f)
            throw std::out_of_range("Tried setting collision damping out of range [0, 1]");
        collision_damping = val;
    }
    // float smoothing_rad = 1.2f;
    float smoothing_rad = 20.0f;
    float target_density = 2.3f;
    float viscosityStrength = 0.01f;
    float LJ_strength = 1.0f;
    // GPUSort gpuSort;
    bool isPaused;
    bool pauseNextFrame;
    float pressureMultiplier = 3.0f;
    float nearPressureMultiplier = 1.5f;

    float world_scale = 10.0f;
    float inv_world_scale = 1.0f / world_scale;

    float PressureFromDensity(float density)
    {
        return (density - target_density) * pressureMultiplier;
    }

    float NearPressureFromDensity(float nearDensity)
    {
        return nearDensity * nearPressureMultiplier;
    }

    /// @brief 4x4 Matrix to describe local space to world space.
    arma::fmat localToWorld;

    /// @brief 4x4 Matrix to describe world space to local space.
    arma::fmat worldToLocal;

    Sim3d() {}

    Sim3d(
        std::vector<float> ps,
        std::vector<float> vs,
        std::vector<float> rs,
        std::vector<float> bas,
        int nP) : positions(ps),
                  velocities(vs),
                  rotations(rs),
                  bond_angles(bas),
                  num_particles(nP)
    {
        predicted_positions.resize(nP * 3, 0.0f);
        densities.resize(nP * 2, 0.0f);
        forces.resize(nP * 3, 0.0f);
        energies.resize(nP, 0.0f);
        forces_coul.resize(nP * 3, 0.0f);
        ang_vels.resize(nP * 3, 0.0f);
        torques.resize(nP * 3, 0.0f);
        localToWorld = arma::eye<arma::fmat>(4, 4);
        worldToLocal = arma::eye<arma::fmat>(4, 4);
    }

    /// @brief
    /// @param [in] frameTime
    /// @return The updated positions after the simulation has ran.
    void RunSimulationFrame(float frame_time, arma::fmat &l2w, arma::fmat &w2l)
    {
        localToWorld = l2w;
        worldToLocal = w2l;
        now = time(0);
        // std::cout << frame_time << std::endl;
        if (!isPaused)
        {
            float time_step = time_scale / (frame_time * 1e-12);

            // Remember, updates shaders, should be handled in node
            // UpdateSettings(timeStep);

            for (int i = 0; i < iter_per_frame; i++)
            {
                std::fill(energies.begin(), energies.end(), 0.0f);
                std::fill(forces.begin(), forces.end(), 0.0f);
                std::fill(forces_coul.begin(), forces_coul.end(), 0.0f);
                std::fill(torques.begin(), torques.end(), 0.0f);
                RunSimulationStep(time_step);
            }
        }
    }

    /// @brief Came from Sebastian Lague's Compute shader.
    /// @param delta_time
    /// @param pos
    /// @param vels
    /// @param pred_pos
    /// @param N
    void CalcExternalForces(float &delta_time)
    {
        for (size_t i = 0; i < num_particles; i++)
        {
            velocities[i * 3 + 1] += gravity * delta_time;
            // pred_pos[i] = positions[i] + velocities[i] * (1.0f / 120.0f);             // Why 1/120?
            // pred_pos[i + 1] = positions[i + 1] + velocities[i + 1] * (1.0f / 120.0f); // Why 1/120?
            // pred_pos[i + 2] = positions[i + 2] + velocities[i + 2] * (1.0f / 120.0f); // Why 1/120?
        }
    }

    void Calculate_LJ_Energy()
    {
        for (size_t i = 0; i < num_particles; i++)
        {
            size_t pos_si = i * 3;
            arma::fvec arma_pos = arma::fvec({positions[pos_si], positions[pos_si + 1], positions[pos_si + 2]});
            for (size_t j = i + 1; j < num_particles; j++){
                size_t pos_si2 = j * 3;
                arma::fvec arma_neighbor = arma::fvec({positions[pos_si2], positions[pos_si2 + 1], positions[pos_si2 + 2]});
                float R = arma::norm(arma_pos - arma_neighbor);
                if (R > smoothing_rad)
                    continue;
                float sr = sigma_OO / R;
                float sr_6 = std::pow(sr, 6);
                float sr_12 = sr_6 * sr_6;
                float E_LJ = 4.0f * epsilon_OO_J * (sr_12 - sr_6);

                energies[i] += E_LJ;
                energies[j] += E_LJ;
            }
        }
    }

    std::vector<arma::fvec> get_hydrogens_pos(size_t i, arma::fvec O_pos)
    {
        size_t rot_si = i * 4;
        std::vector<arma::fvec> H_poss;
        // Get bond angle for hydrogens
        float ba = bond_angles[i];
        float ba_r = (ba/2.0f - 90.f) * M_PI / 180.0f;
        arma::fvec h1_temp_pos = arma::fvec({-hr_109 * std::cos(ba_r),
                                             hr_109 * std::sin(ba_r),
                                             0.0f});
        arma::fvec h2_temp_pos = arma::fvec({hr_109 * std::cos(ba_r),
                                             hr_109 * std::sin(ba_r),
                                             0.0f});

        arma::fvec quat_mole = arma::fvec({rotations[rot_si],
                                           rotations[rot_si + 1],
                                           rotations[rot_si + 2],
                                           rotations[rot_si + 3]});
        arma::fvec quat_m_conj = quat_conj(quat_mole);
        arma::fvec forward_to_quat = quat_mult(quat_mult(quat_mole, quat_conj(forward_q)), quat_m_conj);
        H_poss.push_back(rot_vec_by_quat(h1_temp_pos, forward_to_quat) + O_pos);
        H_poss.push_back(rot_vec_by_quat(h2_temp_pos, forward_to_quat) + O_pos);
        return H_poss;
    }

    /// @brief Calculates the translational forces that the hyrdogens and oxygens apply onto each other.
    /// @param O_pos Position of Oxygen.
    /// @param H_pos Position of Hydrogen.
    /// @return Translation force from hydrogen to oxygen.
    arma::fvec Calc_O_H_Force(arma::fvec O_pos, arma::fvec H_pos)
    {
        arma::fvec OH_dir = arma::normalise(O_pos - H_pos);
        float R2 = arma::dot(OH_dir, OH_dir);
        float F_C_OH = k_e * q_O * q_H / R2;
        return OH_dir * F_C_OH;
    }

    /// @brief Gets the translational and torque components of this H to H interaction.
    /// @param O_pos Oxygen that holds first passed in H.
    /// @param H1_pos H that belongs to the Oxygen.
    /// @param H2_pos H that belongs to the other molecule.
    /// @return Two vectors, first being translational component, second being torque.
    std::vector<arma::fvec> Calc_H_H_Coul_Forces(arma::fvec O_pos, arma::fvec H1_pos, arma::fvec H2_pos)
    {
        // Setup our spatial math
        std::vector<arma::fvec> forces;
        arma::fvec HH_diff = H1_pos - H2_pos;
        float R2 = arma::dot(HH_diff, HH_diff);
        float F_C_HH = k_e * q_H * q_H / R2;
        arma::fvec F_v = arma::normalise(HH_diff) * F_C_HH;

        // Split into our translational and torque components
        arma::fvec HO_r = O_pos - H1_pos;
        arma::fvec HO_dir = arma::normalise(HO_r);
        arma::fvec trans_comp = arma::dot(F_v, HO_dir) * HO_dir;

        arma::fvec tor_comp = arma::cross(HO_r, F_v);

        // Calculate moment of inertia and scale torque by inverse. I = m*r^2
        // Technically this should conver this to rotational velocity for us
        // to easier apply to the molecule later.
        float I = arma::dot(HO_r, HO_r) * mass_H;
        tor_comp /= I;

        forces.push_back(trans_comp);
        forces.push_back(tor_comp);
        return forces;
    }

    std::vector<arma::fvec> Calculate_Coulombic_Forces(size_t i, size_t j)
    {
        // First force is translational, second rotational.
        std::vector<arma::fvec> coul_forces;
        arma::fvec trans_forces_O1 = arma::fvec({0.0f, 0.0f, 0.0f});
        arma::fvec trans_forces_O2 = arma::fvec({0.0f, 0.0f, 0.0f});
        arma::fvec rot_forces_O1 = arma::fvec({0.0f, 0.0f, 0.0f});
        arma::fvec rot_forces_O2 = arma::fvec({0.0f, 0.0f, 0.0f});
        size_t pos_si = i * 3;
        arma::fvec arma_pos = arma::fvec({positions[pos_si],
                                          positions[pos_si + 1],
                                          positions[pos_si + 2]});
        std::vector<arma::fvec> h_poss_1 = get_hydrogens_pos(i, arma_pos);
        arma::fvec h1_pos_1 = h_poss_1[0];
        arma::fvec h2_pos_1 = h_poss_1[1];

        size_t pos_si2 = j * 3;
        arma::fvec arma_neighbor = arma::fvec({positions[pos_si2],
                                               positions[pos_si2 + 1],
                                               positions[pos_si2 + 2]});
        std::vector<arma::fvec> h_poss_2 = get_hydrogens_pos(j, arma_pos);
        arma::fvec h1_pos_2 = h_poss_2[0];
        arma::fvec h2_pos_2 = h_poss_2[1];

        arma::fvec diff = arma_pos - arma_neighbor;
        arma::fvec OO_dir = arma::normalise(diff);
        float R2 = arma::dot(diff, diff);

        // Oxygen_1 - Oxygen_2
        float F_C_OO = k_e * q_O * q_O / R2;
        trans_forces_O1 += F_C_OO * OO_dir;
        trans_forces_O2 -= F_C_OO * OO_dir;

        // Calculate forces for each OH combination
        arma::fvec F_O1_H1_2 = Calc_O_H_Force(arma_pos, h1_pos_2);
        arma::fvec F_O1_H2_2 = Calc_O_H_Force(arma_pos, h2_pos_2);
        arma::fvec F_O2_H1_1 = Calc_O_H_Force(arma_neighbor, h1_pos_1);
        arma::fvec F_O2_H2_1 = Calc_O_H_Force(arma_neighbor, h2_pos_1);

        // First mole's H Force components
        std::vector<arma::fvec> O1_H1_1_H1_2_comps = Calc_H_H_Coul_Forces(arma_pos, h1_pos_1, h1_pos_2);
        std::vector<arma::fvec> O1_H1_1_H2_2_comps = Calc_H_H_Coul_Forces(arma_pos, h1_pos_1, h2_pos_2);
        arma::fvec O1_H1_tor = O1_H1_1_H1_2_comps[1] + O1_H1_1_H2_2_comps[1];
        std::vector<arma::fvec> O1_H2_1_H1_2_comps = Calc_H_H_Coul_Forces(arma_pos, h2_pos_1, h1_pos_2);
        std::vector<arma::fvec> O1_H2_1_H2_2_comps = Calc_H_H_Coul_Forces(arma_pos, h2_pos_1, h2_pos_2);
        arma::fvec O1_H2_tor = O1_H2_1_H1_2_comps[1] + O1_H2_1_H2_2_comps[1];

        // Second mole's H Force components
        std::vector<arma::fvec> O2_H1_2_H1_1_comps = Calc_H_H_Coul_Forces(arma_neighbor, h1_pos_2, h1_pos_1);
        std::vector<arma::fvec> O2_H1_2_H2_1_comps = Calc_H_H_Coul_Forces(arma_neighbor, h1_pos_2, h2_pos_1);
        arma::fvec O2_H1_tor = O2_H1_2_H1_1_comps[1] + O2_H1_2_H2_1_comps[1];
        std::vector<arma::fvec> O2_H2_2_H1_1_comps = Calc_H_H_Coul_Forces(arma_neighbor, h2_pos_2, h1_pos_1);
        std::vector<arma::fvec> O2_H2_2_H2_1_comps = Calc_H_H_Coul_Forces(arma_neighbor, h2_pos_2, h2_pos_1);
        arma::fvec O2_H2_tor = O2_H2_2_H1_1_comps[1] + O2_H2_2_H2_1_comps[1];

        // Sum translation forces acting on Oxygens:
        trans_forces_O1 += F_O1_H1_2 + F_O1_H2_2;
        trans_forces_O2 += F_O2_H1_1 + F_O2_H2_1;
        trans_forces_O1 += O1_H1_1_H1_2_comps[0] + O1_H1_1_H2_2_comps[0] + O1_H2_1_H1_2_comps[0] + O1_H2_1_H2_2_comps[0];
        trans_forces_O2 += O2_H1_2_H1_1_comps[0] + O2_H1_2_H2_1_comps[0] + O2_H2_2_H1_1_comps[0] + O2_H2_2_H2_1_comps[0];

        coul_forces.push_back(trans_forces_O1);
        coul_forces.push_back(trans_forces_O2);
        coul_forces.push_back(O1_H1_tor + O1_H2_tor);
        coul_forces.push_back(O2_H1_tor + O2_H2_tor);
        return coul_forces;
    }

    void Calculate_LJ_Force()
    {
        float sqr_rad = smoothing_rad * smoothing_rad;
        for (size_t i = 0; i < num_particles; i++)
        {
            size_t pos_si = i * 3;
            arma::fvec arma_pos = arma::fvec({positions[pos_si],
                                              positions[pos_si + 1],
                                              positions[pos_si + 2]});

            for (size_t j = i + 1; j < num_particles; j++)
            {
                size_t pos_si2 = j * 3;
                if (i == j)
                    continue;
                arma::fvec arma_neighbor = arma::fvec({positions[pos_si2],
                                                       positions[pos_si2 + 1],
                                                       positions[pos_si2 + 2]});
                arma::fvec diff = arma_pos - arma_neighbor;
                float R2 = arma::dot(diff, diff);
                if (R2 > sqr_rad)
                    continue;
                if (R2 == 0.0f)
                    continue;
                float R = std::sqrt(R2);

                double sr = sigma_OO/R;
                double sr_7 = pow(sr, 7);
                double sr_13 = pow(sr, 13);
                float f_tot = 0.0;

                // LJ
                float F_LJ = 24.0f * epsilon_OO_J * (2.0f*sr_13 - sr_7);
                f_tot += F_LJ;

                // Coulombic Forces
                std::vector<arma::fvec> coul_forces = Calculate_Coulombic_Forces(i, j);

                if (std::isnan(f_tot))
                    continue;

                arma::fvec f_part = arma::fvec({0.0f, 0.0f, 0.0f});
                f_part += f_tot * arma::normalise(diff);
                arma::fvec r_O1_part = coul_forces[2];
                // std::cout << "R_o1: " << r_O1_part.t();
                arma::fvec r_O2_part = coul_forces[3];
                // std::cout << "R_o2: " << r_O2_part.t();

                forces[pos_si] += f_part(0);
                forces[pos_si + 1] += f_part(1);
                forces[pos_si + 2] += f_part(2);
                forces[pos_si2] -= f_part(0);
                forces[pos_si2 + 1] -= f_part(1);
                forces[pos_si2 + 2] -= f_part(2);

                forces_coul[pos_si] += coul_forces[0](0);
                forces_coul[pos_si + 1] += coul_forces[0](1);
                forces_coul[pos_si + 2] += coul_forces[0](2);
                forces_coul[pos_si2] += coul_forces[1](0);
                forces_coul[pos_si2 + 1] += coul_forces[1](1);
                forces_coul[pos_si2 + 2] += coul_forces[1](2);



                torques[pos_si] += r_O1_part(0);
                torques[pos_si + 1] += r_O1_part(1);
                torques[pos_si + 2] += r_O1_part(2);
                torques[pos_si2] += r_O2_part(0);
                torques[pos_si2 + 1] += r_O2_part(1);
                torques[pos_si2 + 2] += r_O2_part(2);
            }
        }
    }
    void UpdatePositions(float delta_time)
    {
        // outFile.open("./err_logs/log.txt", std::ios::app);
        for (size_t i = 0; i < num_particles; i++)
        {
            arma::fvec t_v = arma::fvec({torques[i * 3], torques[i * 3 + 1], torques[i * 3 + 2]});

            // in Angstroms
            arma::fvec pos = arma::fvec({positions[i * 3], positions[i * 3 + 1], positions[i * 3 + 2]});

            // in Angstroms/s
            arma::fvec vel = arma::fvec({velocities[i * 3], velocities[i * 3 + 1], velocities[i * 3 + 2]});

            // Forces, in J [kg m^2/s^2] need to multiply by 1e10 to get into angstroms and divide by mass
            arma::fvec f_v = arma::fvec({forces[i * 3], forces[i * 3 + 1], forces[i * 3 + 2]});
            // Forces, in N [kg m/s^2] need to multiply by 1e10 to get into angstroms and divide by mass
            arma::fvec fc_v = arma::fvec({forces_coul[i * 3], forces_coul[i * 3 + 1], forces_coul[i * 3 + 2]});
            arma::fvec f_tot = f_v * 1e10f*1e10f + fc_v * 1e10f;
            f_tot /= h2o_m;

            pos += 0.5f * vel * time_delta + f_tot * time_delta * time_delta;
            vel += f_tot * time_delta;

            // Update positions based on velocities and forces
            // outFile << now << " " << arma::norm(f_v);
            // outFile << " " << arma::norm(fc_v);
            // outFile << " " << arma::norm(t_v);
            // outFile << " " << energies[i] << std::endl;
            positions[i * 3] = pos(0);//0.5f * velocities[i * 3] * time_delta + forces[i * 3] * time_delta * time_delta;
            positions[i * 3 + 1] = pos(1);//0.5f * velocities[i * 3 + 1] * time_delta + forces[i * 3 + 1] * time_delta * time_delta;
            positions[i * 3 + 2] = pos(2);//0.5f * velocities[i * 3 + 2] * time_delta + forces[i * 3 + 2] * time_delta * time_delta;

            velocities[i * 3] = vel(0);//forces[i * 3] * time_delta;
            velocities[i * 3 + 1] = vel(1);//forces[i * 3 + 1] * time_delta;
            velocities[i * 3 + 2] = vel(2);//forces[i * 3 + 2] * time_delta;

            arma::fvec curr_rot = arma::fvec({rotations[i * 4],
                                              rotations[i * 4 + 1],
                                              rotations[i * 4 + 2],
                                              rotations[i * 4 + 3]});
            arma::fvec curr_rot_conj = quat_conj(curr_rot);

            // Get the current angular velocity as a quaternion
            arma::fvec ang_vel = arma::fvec({ang_vels[i * 3],
                                             ang_vels[i * 3 + 1],
                                             ang_vels[i * 3 + 2]});

            // Get the current torque as a quaternion
            arma::fvec tor = arma::fvec({torques[i * 3], torques[i * 3 + 1], torques[i * 3 + 2]});
            ang_vel += tor * 1e20 * time_delta;
            // arma::fvec torq = arma::fvec({tor(0), tor(1), tor(2), 0.0f});
            arma::fvec ang_velq;
            if (arma::norm(ang_vel) < 1e-6)
                ang_velq = arma::fvec({0.0f, 0.0f, 0.0f, 1.0f});
            else
                ang_velq = arma::normalise(arma::fvec({ang_vel(0), ang_vel(1), ang_vel(2), 0.0f}));

            arma::fvec delta_rot = 0.5f * ang_velq * 1e10 * time_delta * time_delta;
            curr_rot = arma::normalise(quat_mult(quat_mult(delta_rot, curr_rot_conj), curr_rot));

            // Update the angular velocity and rotation quaternion
            ang_vels[i * 3] = ang_vel(0);
            ang_vels[i * 3 + 1] = ang_vel(1);
            ang_vels[i * 3 + 2] = ang_vel(2);

            rotations[i * 4] = curr_rot(0);
            rotations[i * 4 + 1] = curr_rot(1);
            rotations[i * 4 + 2] = curr_rot(2);
            rotations[i * 4 + 3] = curr_rot(3);

            ResolveCollisions(i);
        }
        // outFile.close();
    }

    void ResolveCollisions(size_t si)
    {
        si *= 3;
        // Transform position/velocity to the local space of the bounding box (scale not included)
        std::vector<float> pos = {positions[si] * inv_world_scale,
                                  positions[si + 1] * inv_world_scale,
                                  positions[si + 2] * inv_world_scale};
        arma::fvec arma_pos = arma::join_cols(arma::fvec(pos), arma::fvec({1.0f}));
        arma::fvec pos_local = worldToLocal * arma_pos;

        std::vector<float> vel = {velocities[si], velocities[si + 1], velocities[si + 2]};
        arma::fvec arma_vel = arma::join_cols(arma::fvec(vel), arma::fvec({1.0f}));
        arma::fvec vel_local = worldToLocal * arma_vel;

        // Half the size of the bounding box
        const arma::fvec halfSize = {0.5f, 0.5f, 0.5f};

        // Apply periodic boundary conditions (wrap around the edges)
        for (int axis = 0; axis < 3; ++axis)
        {
            float edgeDst = halfSize[axis] - abs(pos_local[axis]);
            // Check if position is outside the bounds on the negative side
            if (edgeDst <= 0.0f)
            {
                pos_local[axis] = halfSize[axis] * sign(pos_local[axis]); // Wrap to the other side
                vel_local[axis] *= -1 * collision_damping;                // Apply damping to velocity
            }
        }

        // Transform resolved position/velocity back to world space
        arma::fvec corr_pos = localToWorld * pos_local;
        arma::fvec corr_vel = localToWorld * vel_local;

        positions[si] = corr_pos[0] * world_scale;
        positions[si + 1] = corr_pos[1] * world_scale;
        positions[si + 2] = corr_pos[2] * world_scale;
        velocities[si] = corr_vel[0];
        velocities[si + 1] = corr_vel[1];
        velocities[si + 2] = corr_vel[2];
    }

    void RunSimulationStep(float time_step)
    {
        Calculate_LJ_Energy();
        Calculate_LJ_Force();
        UpdatePositions(time_step);
    }
};