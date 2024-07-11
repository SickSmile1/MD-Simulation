//
// Created by ilia on 05/06/24.
//

#ifndef MD_SIMULATION_LJ_H
#define MD_SIMULATION_LJ_H
#include "Atoms.h"

double lj_direct_summation(Atoms &atoms, const double epsilon = 1.0, const double sigma = 1.0, const double c_energy = 0.);

#endif //MD_SIMULATION_LJ_H
