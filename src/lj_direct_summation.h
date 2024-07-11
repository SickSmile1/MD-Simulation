//
// Created by ilia on 05/06/24.
//

#ifndef MD_SIMULATION_LJ_DIRECT_SUMMATION_H
#define MD_SIMULATION_LJ_DIRECT_SUMMATION_H
#include "Atoms.h"
#include "neighbors.h"

double lj_direct_summation(Atoms &atoms, const NeighborList, const double epsilon = 1.0,
                           const double sigma = 1.0, const double cutoff = 1.0, double mass = 1.,
                           const double c_energy = 0.);
#endif //MD_SIMULATION_LJ_DIRECT_SUMMATION_H
