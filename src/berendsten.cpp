//
// Created by ilia on 17/06/24.
//
#include "iostream"
#include "Atoms.h"
#include "domain.h"

void berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time, const double mass){
    // double ekin = 0, T = 0, lambda = 0, diff_t = 0;
    // pi = vi/mi
    const double diff_t = temperature/ ( ( mass *
            atoms.velocities.colwise().squaredNorm()*0.5).sum() /
            (1.5 * atoms.velocities.cols() * kB) );
    atoms.velocities.array() *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time)); // lambda
}

void berendsen_thermostat(Atoms &atoms, Domain &domain, const double temperature, const double timestep,
                          const double relaxation_time, const double mass){
    // double ekin = 0, T = 0, lambda = 0, diff_t = 0;
    // pi = vi/mi
    double v_squared_norm{MPI::allreduce(
            (mass*atoms.velocities.block(0,0,3,domain.nb_local()).colwise().squaredNorm()*0.5).sum(),
                                     MPI_SUM, MPI_COMM_WORLD)};
    int n_atoms{MPI::allreduce(domain.nb_local(),
                               MPI_SUM, MPI_COMM_WORLD)};
    const double diff_t = temperature/ ( v_squared_norm /
                                         (1.5 * n_atoms * kB) );
    atoms.velocities.array() *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time)); // lambda
}

void berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time){
    // pi = vi/mi
    // ekin = sum_i 1/2 (pi^2 / mi)
    const double diff_t = temperature/ ( (atoms.velocities.colwise().squaredNorm()*0.5).sum() /
                                         (1.5 * atoms.velocities.cols()) );
    atoms.velocities *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time)); // lambda
}