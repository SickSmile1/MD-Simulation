//
// Created by ilia on 19/04/24.
//

#include "verlet.h"

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  Eigen::Array3Xd &forces, const double timestep, const double mass) {
    velocities += 0.5 * forces * timestep /mass;
    positions += velocities * timestep;
}

void verlet_step2(Eigen::Array3Xd &velocities, Eigen::Array3Xd &forces, const double timestep,
                  const double mass) {
    velocities += (0.5 * forces * timestep / mass);
}