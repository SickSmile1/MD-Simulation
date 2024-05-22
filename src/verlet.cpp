//
// Created by ilia on 19/04/24.
//

#include "verlet.h"

void verlet_step1(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass; 
    atoms.positions += atoms.velocities * timestep; 
}

void verlet_step2(Atoms &atoms, double timestep, double mass) {
    atoms.velocities += 0.5 * atoms.forces * timestep / mass;

}