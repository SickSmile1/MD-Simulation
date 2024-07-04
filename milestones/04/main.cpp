//
// Created by ilia on 05/06/24.
//

#include "Atoms.h"
#include "verlet.h"
#include "lj.h"
#include "xyz.h"
#include "ostream"

int main() {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions, velocities);
    std::ofstream myfile;
    // myfile.open ("energy"); // +std::to_string(j));
    double epot = 0, ekin = 0;
    for (int i=0;i<2;i++) {
        epot = lj_direct_summation(atoms,1.,1.);
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, .1e-3, 1);
        verlet_step2(atoms.velocities, atoms.forces, 1e-3);
        // ekin = 1*(atoms.velocities.colwise().squaredNorm()*0.5).sum();
        // myfile << ekin+epot << "\t";
    }
    // myfile.close();
    return 1;
}