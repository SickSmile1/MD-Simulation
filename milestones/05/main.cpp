//
// Created by ilia on 19/06/24.
//

#include "Atoms.h"
#include "verlet.h"
#include "lj.h"
#include "xyz.h"
#include "ostream"
#include "cubicLatice.h"
#include <fmt/core.h>
#include "iostream"
#include "berendsten.h"

void equilibrate(Atoms &atoms, std::string ener, std::string traj, double sigma) {
    std::ofstream myfile;
    std::ofstream myfile1;
    myfile.open ("plot/energy_"+ener);
    myfile1.open ("plot/tray_"+traj);
    double epot = 0, ekin = 0;
    epot = lj_direct_summation(atoms,1.,1.);
    for (int i=0;i<1000;i++) {
        // lj_direct_summation(atoms,1.,1.);
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, .1e-3, 1);
        epot = lj_direct_summation(atoms,1.,1.);
        verlet_step2(atoms.velocities, atoms.forces, 1e-3);
        berendsen_thermostat(atoms,.1, 1e-3, 1e-1);
        ekin = 1*(atoms.velocities.colwise().squaredNorm()*0.5).sum();
        myfile << ekin+epot << "\t";
        if(i%100==0) write_xyz(myfile1,atoms);
    }
    write_xyz(myfile1,atoms);
    myfile.close();
    myfile1.close();
}

int main() {
    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions, velocities);
    equilibrate(atoms,"ij","ij",1.0);
    Atoms at1(64);
    createCubicLatice(at1,4, 1.);
    equilibrate(at1,"cubic_small","cubic_small",1.0);
    // createCubicLatice(at2,10);
    // equilibrate(at2,"cubic_big_1.1", "cubic_big_1.1",1.1);
    for (double i = .7; i<1.5; i+=.1) {
        Atoms at2(1000);
        createCubicLatice(at2,10, i);
        std::string i_s = fmt::format("{:.2f}",i); //std::to_string(i);
        equilibrate(at2,"cubic_big_"+i_s, "cubic_big_"+i_s,1);
    }
    return 0;
}
