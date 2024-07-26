//
// Created by ilia on 17/06/24.
//

#ifndef MD_SIMULATION_BERENDSTEN_H
#define MD_SIMULATION_BERENDSTEN_H
#include "Atoms.h"
#include "domain.h"

double berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time, const double mass);

double berendsen_thermostat(Atoms &atoms, Domain &domain, const double temperature, const double timestep,
                          const double relaxation_time, const double mass);

void berendsen_thermostat(Atoms &atoms, const double temperature, const double timestep,
                          const double relaxation_time);

#endif //MD_SIMULATION_BERENDSTEN_H
