//
// Created by ilia on 17/06/24.
//
#include "Atoms.h"
#include "domain.h"

double berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time, const double mass){
  // version with normal units
  const double diff_t = temperature/ ( ( mass *
          atoms.velocities.colwise().squaredNorm()*0.5).sum() /
          (1.5 * atoms.velocities.cols() * kB) );
  atoms.velocities.array() *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time));
  return sqrt(1+(diff_t-1)*(timestep/relaxation_time));
}

double berendsen_thermostat(Atoms &atoms, Domain &domain, const double temperature, const double timestep,
                          const double relaxation_time, const double mass){
  // get velocities without ghosts from all cores
  double v_squared_norm{MPI::allreduce(
          (mass*atoms.velocities.block(0,0,3,domain.nb_local()).colwise().squaredNorm()*0.5).sum(),
                                   MPI_SUM, MPI_COMM_WORLD)};
  // get overall number of atoms on all domains w/o ghosts
  int n_atoms{MPI::allreduce(domain.nb_local(),
                             MPI_SUM, MPI_COMM_WORLD)};
  // T0/T, T0 = set Temperature, T = current Temperature 
  const double diff_t = temperature/ ( v_squared_norm /
                                       (1.5 * n_atoms * kB) );
  // apply Lambda on velocities
  atoms.velocities.array() *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time));
  // return velocities for tests
  return sqrt(1+(diff_t-1)*(timestep/relaxation_time));
}

void berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time){
  // version for LJ units with kB = 1
  const double diff_t = temperature/ ( (atoms.velocities.colwise().squaredNorm()*0.5).sum() /
                                       (1.5 * atoms.velocities.cols()) );
  atoms.velocities *=  sqrt(1+(diff_t-1)*(timestep/relaxation_time));
}
