//
// Created by ilia on 03/07/24.
//

#include "Atoms.h"
#include "verlet.h"
#include "lj.h"
#include "xyz.h"
#include "ostream"
#include "iostream"
#include "math.h"

int main() {
   std::ofstream myfile;
   double arr[4] = {.001,.0001,.00001,.000001};
    for (int j = 0; j < 4; j++) {
        auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
        Atoms atoms(positions, velocities);
        myfile.open("energy_"+std::to_string(j));
        double epot = 0, ekin = 0;
        int runs = 100000000 * arr[j];//std::ceil(.0001*std::pow(10,j));
        std::cout << runs << ":" << std::pow(10,-(j+1) )<<std::endl;
        int dp = runs/100;
        for (int i = 0; i < (runs); i++) {
            epot = lj_direct_summation(atoms, 1., 1.);
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces, arr[j]);
            verlet_step2(atoms.velocities, atoms.forces, arr[j]);
            ekin = 1 * (atoms.velocities.colwise().squaredNorm() * 0.5).sum();

            if(i%dp == 0) {myfile << ekin + epot << "\n"; }
        }
        std::cout << "finished run " << j << std::endl;
        myfile.close();
    }
    return 1;
}