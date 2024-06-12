//
// Created by ilia on 05/06/24.
//

#include "Atoms.h"
#include "verlet.h"
#include "lj_direct_summation.h"
#include "xyz.h"

int main() {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions, velocities);
    for (int i=0;i<101;i++) {
        lj_direct_summation(atoms,1.,1.);
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, .1e-3, 1);
        verlet_step2(atoms.velocities, atoms.forces, 1e-3);
    }



    return 1;
}