//
// Created by ilia on 01/07/24.
//

#include <gtest/gtest.h>
#include "iostream"

#include "Atoms.h"
#include "neighbors.h"
#include "xyz.h"
#include "ducastelle.h"
#include "berendsten.h"
#include "verlet.h"

TEST(DucastelleTest2, energy) {
    GTEST_SKIP();
    auto [names, positions, velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    NeighborList nl;
    double old_energy = 0, ekin = 0;
    std::ofstream traj("traj2.xyz");
    for(int j = 0; j < 20; j+=2){
        Atoms at(positions);
        at.velocities.setRandom();
        at.velocities *= .1e-5;
        nl.update(at,5.5);
        const double timestep = j; double T = 0;
        const double fixed_mass = 196.96657 / 0.009649;
        for(int i = 0; i < 10000; i++) {
            verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
            nl.update(at,5.5);
            at.forces.setZero();
            double epot = ducastelle(at, nl);
            verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
            if (i < 2000) { berendsen_thermostat(at, 750, timestep, 1000, fixed_mass); }
            if(i%100==0) {
                write_xyz(traj, at);
                ekin = fixed_mass*(at.velocities.colwise().squaredNorm()*0.5).sum();
                T = ( ekin / (1.5 * at.positions.cols() * kB) );
                // std::cout << "Temp is: " << T << ", energy: " << ekin+epot << std::endl;
            }
            if (i == 1490) { old_energy = ekin+epot; std::cout << "step: " << j << std::endl;}
            if (i == 1999) { EXPECT_NEAR(old_energy, ekin+epot, 15); }
        }
    }
    traj.close();
}

TEST(DucastelleTest3, temp_to_file) {
    GTEST_SKIP();
    auto [names, positions, velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    NeighborList nl;
    double old_energy = 0, ekin = 0;
    for(int j = 1; j < 32; j+=5){
        std::ofstream myfile;
        std::ofstream myfile2;
        myfile.open ("temp_"+std::to_string(j));
        myfile2.open ("e_all_"+std::to_string(j));
        Atoms at(positions);
        at.velocities.setRandom();
        at.velocities *= .1e-5;
        nl.update(at,5.5);
        const double timestep = j; double T = 0;
        const double fixed_mass = 196.96657 / 0.009649;
        int runs = 10000/j;
        int dp = std::ceil(runs/100);
        for(int i = 0; i < std::ceil(runs); i++) {
            verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
            nl.update(at,5.5);
            at.forces.setZero();
            double epot = ducastelle(at, nl);
            verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
            if (i < 2000) { berendsen_thermostat(at, 750, timestep, 1000, fixed_mass); }
            if(i> 150 && i%dp==0) {
                ekin = fixed_mass*(at.velocities.colwise().squaredNorm()*0.5).sum();
                T = ( ekin / (1.5 * at.positions.cols() * kB) );
                myfile << T << "\t";
                myfile2 << ekin+epot << "\t";
            }

        }
        myfile.close();
        myfile2.close();
    }
}