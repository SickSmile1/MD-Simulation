//
// Created by ilia on 03/07/24.
//

#include "Atoms.h"
#include "verlet.h"
#include "lj.h"
#include "xyz.h"
#include "ostream"
#include "iostream"

int main() {
   std::ofstream myfile;
    for (int j = 1; j < 11; j++) {
        auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
        Atoms atoms(positions, velocities);
        myfile.open("energy_"+std::to_string(j));
        double epot = 0, ekin = 0;
        int runs = std::ceil(10000/j);
        for (int i = 0; i < (runs); i++) {
            epot = lj_direct_summation(atoms, 1., 1.);
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, j, 1);
            verlet_step2(atoms.velocities, atoms.forces, j);
            ekin = 1 * (atoms.velocities.colwise().squaredNorm() * 0.5).sum();
            int dp = std::ceil(runs/100);
            if(i%(dp == 0)) {myfile << ekin + epot << "\n";}
        }
        std::cout << "finished run " << j << std::endl;
        myfile.close();
    }
    return 1;
}