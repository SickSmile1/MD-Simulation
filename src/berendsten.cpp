//
// Created by ilia on 17/06/24.
//
#include "iostream"
#include "Atoms.h"

void berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time, const double mass){
    // double ekin = 0, T = 0, lambda = 0, diff_t = 0;
    // pi = vi/mi
    // ekin = sum_i 1/2 (pi^2 / mi)
    // ekin = (atoms.velocities.colwise().squaredNorm()*0.5).sum();
    // T = ekin / (1.5 * atoms.velocities.cols() * 1);
    const double diff_t = temperature/ ( ( mass *
            atoms.velocities.colwise().squaredNorm()*0.5).sum() /
            (1.5 * atoms.velocities.cols() * kB) );
    // T;
    // lambda = sqrt(1+(diff_t-1)*(timestep/relaxation_time));
    atoms.velocities.array() *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time)); // lambda
}

void berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time){
    // double ekin = 0, T = 0, lambda = 0, diff_t = 0;
    // pi = vi/mi
    // ekin = sum_i 1/2 (pi^2 / mi)
    // ekin = (atoms.velocities.colwise().squaredNorm()*0.5).sum();
    // T = ekin / (1.5 * atoms.velocities.cols() * 1);
    const double diff_t = temperature/ ( (atoms.velocities.colwise().squaredNorm()*0.5).sum() /
                                         (1.5 * atoms.velocities.cols()) );
    // T;
    // lambda = sqrt(1+(diff_t-1)*(timestep/relaxation_time));
    //std::cout << 1.5 * atoms.velocities.cols() << std::endl;
    //std::cout << (atoms.velocities.colwise().squaredNorm()*0.5).sum() << std::endl;
    atoms.velocities *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time)); // lambda
}