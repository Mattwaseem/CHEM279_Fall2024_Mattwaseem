
    void CalculateViscosity(float &delta_time, std::vector<float> &pred_pos, std::vector<float> &vels, size_t si)
    {
        size_t pos_si = si * 3;
        arma::fvec arma_pos = arma::fvec({pred_pos[si],
                                          pred_pos[si + 1],
                                          pred_pos[si + 2]});
        // int3 originCell = GetCell3D(pos, smoothing_rad);
        float sqrRadius = smoothing_rad * smoothing_rad;
        arma::fvec viscosityForce = arma::fvec({0.0f, 0.0f, 0.0f});

        arma::fvec arma_vel = arma::fvec({vels[si],
                                          vels[si + 1],
                                          vels[si + 2]});

// gonna do some fancy omp here
#pragma omp parallel for reduction(+ : viscosityForce) shared(arma_pos, arma_vel, sqrRadius, pred_pos, vels, viscosityStrength)
        for (size_t j = 0; j < num_particles; j++)
        {
            size_t si2 = j * 3;
            if (si == si2)
                continue;

            arma::fvec arma_neighbor = arma::fvec({pred_pos[si2],
                                                   pred_pos[si2 + 1],
                                                   pred_pos[si2 + 2]});
            arma::fvec neighbor_diff = arma_neighbor - arma_pos;
            float sqr_dist_neighbor = arma::dot(neighbor_diff, neighbor_diff);

            // Skip if not within radius
            if (sqr_dist_neighbor > sqrRadius)
                continue;

            // Calculate viscosity
            float dst = sqrt(sqr_dist_neighbor);
            if (dst == 0)
            {
                viscosityForce += ((arma::fvec)arma::normalise(arma::randn<arma::fvec>(3))) * SmoothingKernelPoly6(1e-6, smoothing_rad);
                continue;
            }
            arma::fvec arma_neigh_vel = arma::fvec({vels[si2],
                                                    vels[si2 + 1],
                                                    vels[si2 + 2]});
            viscosityForce += (arma_neigh_vel - arma_vel) * SmoothingKernelPoly6(dst, smoothing_rad);
        }
        if (std::isnan(viscosityForce[0]) || std::isnan(viscosityForce[1]) || std::isnan(viscosityForce[2]))
            viscosityForce = arma::fvec({0.0f, 0.0f, 0.0f});
        // std::cout << "iV: " << viscosityStrength << " " << viscosityForce.t();
        vels[si] += viscosityForce[0] * viscosityStrength * delta_time;
        vels[si + 1] += viscosityForce[1] * viscosityStrength * delta_time;
        vels[si + 2] += viscosityForce[2] * viscosityStrength * delta_time;
    }

    void CalculatePressureForce(float delta_time,
                                std::vector<float> &dens,
                                std::vector<float> &pred_pos,
                                std::vector<float> &vels,
                                size_t si)
    {
        size_t pos_si = si * 3;
        size_t den_si = si * 2;

        // Calculate pressure
        float density = dens[den_si];
        float densityNear = dens[den_si + 1];
        float pressure = PressureFromDensity(density);
        float nearPressure = NearPressureFromDensity(densityNear);
        arma::fvec pressureForce = arma::fvec({0.0f, 0.0f, 0.0f});
        arma::fvec arma_pos = arma::fvec({pred_pos[pos_si],
                                          pred_pos[pos_si + 1],
                                          pred_pos[pos_si + 2]});
        // int3 originCell = GetCell3D(pos, smoothingRadius);
        float sqrRadius = smoothing_rad * smoothing_rad;

        // _neighbor search
#pragma omp parallel for reduction(+ : pressureForce) shared(arma_pos, pressure, nearPressure, sqrRadius)
        for (size_t j = 0; j < num_particles; j++)
        {
            size_t pos_si2 = j * 3;
            size_t den_si2 = j * 2;
            if (pos_si == pos_si2)
            {
                pressureForce += arma::fvec({0.0f, 0.0f, 0.0f});
                continue;
            }

            arma::fvec arma_neigh_pos = arma::fvec({pred_pos[pos_si2],
                                                    pred_pos[pos_si2 + 1],
                                                    pred_pos[pos_si2 + 2]});
            arma::fvec offsetTo_neighbor = arma_neigh_pos - arma_pos;
            float sqrDstTo_neighbor = arma::dot(offsetTo_neighbor, offsetTo_neighbor);

            // Skip if not within radius
            if (sqrDstTo_neighbor > sqrRadius)
                continue;

            // Calculate pressure force
            float density_neighbor = dens[den_si2];
            float nearDensity_neighbor = dens[den_si2 + 1];
            float neighborPressure = PressureFromDensity(density_neighbor);
            float neighborPressureNear = NearPressureFromDensity(nearDensity_neighbor);

            float sharedPressure = (pressure + neighborPressure) / 2;
            float sharedNearPressure = (nearPressure + neighborPressureNear) / 2;

            float dst = sqrt(sqrDstTo_neighbor);
            arma::fvec dir = dst > 0 ? offsetTo_neighbor / dst : ((arma::fvec)arma::normalise(arma::randn<arma::fvec>(3)));
            arma::fvec neigh_addition = dir * DensityDerivative(dst, smoothing_rad) * sharedPressure / density_neighbor;
            arma::fvec neigh_near_addition = dir * DensityDerivative(dst, smoothing_rad) * sharedPressure / nearDensity_neighbor;

            if (density_neighbor == 0)
                neigh_addition = arma::fvec({0.0f, 0.0f, 0.0f});
            if (nearDensity_neighbor == 0)
                neigh_near_addition = arma::fvec({0.0f, 0.0f, 0.0f});

            pressureForce += neigh_addition;
            pressureForce += neigh_near_addition;
        }
        // }
        arma::fvec acceleration = pressureForce / density;
        if (density == 0 || std::isnan(pressureForce[0]) || std::isnan(pressureForce[1]) || std::isnan(pressureForce[2]))
            acceleration = arma::fvec({0.0f, 0.0f, 0.0f});
        // std::cout << "ia: " << si << " " << density << " " << pressureForce.t() << std::endl;
        vels[pos_si] += acceleration[0] * delta_time;
        vels[pos_si + 1] += acceleration[1] * delta_time;
        vels[pos_si + 2] += acceleration[2] * delta_time;
    }

    void CalculateDensities(size_t si, std::vector<float> &pred_pos, std::vector<float> &dens)
    {
        size_t pos_si = si * 3;
        arma::fvec arma_pos = arma::fvec(std::vector<float>({pred_pos[pos_si],
                                                             pred_pos[pos_si + 1],
                                                             pred_pos[pos_si + 2]}));
        // int3 originCell = GetCell3D(pos, smoothingRadius);
        float sqrRadius = smoothing_rad * smoothing_rad;
        float density = 0;
        float nearDensity = 0;

// _neighbor search
#pragma omp parallel for reduction(+ : density, nearDensity) shared(arma_pos, sqrRadius)
        for (size_t j = 0; j < num_particles; j++)
        {
            size_t si2 = j * 3;
            if (si == si2)
                continue;

            arma::fvec arma_neigh_pos = arma::fvec(std::vector<float>({pred_pos[si2],
                                                                       pred_pos[si2 + 1],
                                                                       pred_pos[si2 + 2]}));
            arma::fvec offsetTo_neighbor = arma_neigh_pos - arma_pos;
            float sqrDstTo_neighbor = arma::dot(offsetTo_neighbor, offsetTo_neighbor);

            // Skip if not within radius
            if (sqrDstTo_neighbor > sqrRadius)
                continue;
            // Calculate density and near density
            float dst = sqrt(sqrDstTo_neighbor);
            density += DensityKernel(dst, smoothing_rad);
            nearDensity += NearDensityKernel(dst, smoothing_rad);
            // }
        }
        size_t den_si = si * 2;
        dens[den_si] = density;
        dens[den_si + 1] = nearDensity;
    }