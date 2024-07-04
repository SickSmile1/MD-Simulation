//
// Created by ilia on 28/06/24.
//

#include "berendsten.h"
#include "verlet.h"
#include "lj_direct_summation.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "xyz.h"
#include "iostream"

int main() {
    auto [names, positions, velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms at(positions);
    at.velocities.setRandom();
    at.velocities *= .1e-5;
    std::ofstream traj("traj2.xyz");
    NeighborList nl;
    nl.update(at,5.5);
    const double timestep = 1; double T = 0;
    const double fixed_mass = 196.96657 / 0.009649;
    for(int i = 0; i < 5000; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        nl.update(at,5.5);
        // lj_direct_summation(at, nl, 1, 1, 5.5, fixed_mass);
        at.forces.setZero();
        double epot = ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
        if (i < 2000) { berendsen_thermostat(at, 650, timestep, 200, fixed_mass); }
        if(i%10==0) {
            write_xyz(traj, at);
            double ekin = fixed_mass*(at.velocities.colwise().squaredNorm()*0.5).sum();
            T = ( ekin / (1.5 * at.positions.cols() * kB) );
            // std::cout << "Temp is: " << T << ", energy: " << ekin+epot << std::endl;
        }
        // if(i%100==0) {std::cout << std::to_string(i)+" steps finished." << std::endl;}
    }
    traj.close();
}