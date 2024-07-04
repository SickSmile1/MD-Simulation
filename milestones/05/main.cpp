//
// Created by ilia on 19/06/24.
//

#include "berendsten.h"
#include "xyz.h"

int main() {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions, velocities);
    // usual value for relaxation: 1^-12 -> 10^-12
    berendsen_thermostat(atoms,1.8, 1e-14, 1e-12);
}