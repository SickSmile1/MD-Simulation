//
// Created by ilia on 19/04/24.
//

#ifndef __VERLET_H
#define __VERLET_H
#include "Atoms.h"

void verlet_step1(Eigen::Array3Xd &positions, Eigen::Array3Xd &velocities,
                  Eigen::Array3Xd &forces, const double timestep, const double mass = 1.);
void verlet_step2(Eigen::Array3Xd &velocities, Eigen::Array3Xd &forces,
                  const double timestep, const double mass = 1.);

#endif  // __VERLET_H