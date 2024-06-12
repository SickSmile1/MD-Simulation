//
// Created by ilia on 05/06/24.
//

#include "Atoms.h"
#include "lj_direct_summation.h"
#include "xyz.h"
auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};

int main() {
    Atoms atoms(positions, velocities);
    lj_direct_summation(atoms,1.,1.);

    return 1;
}