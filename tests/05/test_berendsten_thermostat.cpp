//
// Created by ilia on 19/06/24.
//

#include <gtest/gtest.h>
#include "verlet.h"
#include "Atoms.h"
#include "berendsten.h"
#include "lj.h"
#include "cubicLatice.h"
#include "xyz.h"
#include "math.h"
#include "Eigen/Dense"
#include "iostream"
#include "fstream"
#include "cassert"

TEST(Berendsten, temperature) {
    GTEST_SKIP();
    Atoms atoms(2);
    atoms.velocities << 2,2,2,2,2,2;
    atoms.positions << 1,0,2,2,0,2;
    double ekin = 0, T = 0, lambda = 0, diff_t = 0;

    ekin = (atoms.velocities.colwise().squaredNorm()*0.5).sum();
    T = ekin / (1.5 * atoms.velocities.cols() * 1);
    // std::cout << "T: " << T << std::endl;
    // diff_t = temperature/T;
    // lambda = sqrt(1+(diff_t-1)*(timestep/relaxation_time));
}

TEST(Berendsten, equilibrium) {
    //std::ofstream myfile;
    //myfile.open ("temperatures");
    Atoms at(64);
    at.velocities.setRandom();
    at.velocities *= 1e-5;
    const double timestep = 1e-8;
    double T;
    createCubicLatice(at,4);
    for(int i = 0; i < 1000; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, 1);
        // std::cout << at.forces << std::endl;
        double epot = lj_direct_summation(at);
        verlet_step2(at.velocities, at.forces, timestep);
        berendsen_thermostat(at, 0.3, timestep, 1e-7);
        // T = 1 * (at.velocities.colwise().squaredNorm()*0.5).sum()/ (1.5 * at.velocities.cols() * 1);
        //myfile << T << "\t";
    }
    // myfile.close();
}

TEST(Berendsten2, equilibrium_xyz) {
    GTEST_SKIP();
    std::ofstream traj("traj.xyz");
    //myfile.open ("temperatures");
    Atoms at(125);
    const double timestep = 1e-4;
    createCubicLatice(at,5);
    for(int i = 0; i < 2500; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, 1);
        lj_direct_summation(at, 1, 1);
        verlet_step2(at.velocities, at.forces, timestep);
        berendsen_thermostat(at, 0.3, timestep, 1e-3);
        if(i%100==0) { write_xyz(traj, at); }
    }
    traj.close();
}
