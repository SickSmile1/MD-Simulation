//
// Created by ilia on 05/06/24.
//

#include "Atoms.h"
#include "verlet.h"
#include "lj.h"
#include "xyz.h"
#include "ostream"
#include <fmt/core.h>
#include <string>
#include "iostream"
// obsolete
void find_timestep(Atoms &atoms, std::string ener, std::string traj, double ts) {
    std::ofstream myfile;
    std::ofstream myfile1;
    myfile.open ("energy_"+ener);
    myfile1.open ("tray_"+traj);
    double epot = 0, ekin = 0;
    epot = lj_direct_summation(atoms,1.,1.);
    for (int i=0;i<10000;i++) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, ts, 1);
        epot = lj_direct_summation(atoms,1.,1.);
        verlet_step2(atoms.velocities, atoms.forces, ts);
        ekin = 1*(atoms.velocities.colwise().squaredNorm()*0.5).sum();
        myfile << ekin+epot << "\t";
        if(i%100==0) write_xyz(myfile1,atoms);
    }
    write_xyz(myfile1,atoms);
    myfile.close();
    myfile1.close();
}

int main() {
    for (double i = .000005; i<=501; i*=10) {
      auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
      Atoms atoms(positions, velocities);
      std::string i_s = fmt::format("{:.6f}",i); //std::to_string(i);
      find_timestep(atoms,"lj_"+i_s, "lj_"+i_s,i);
    }
    return 0;
}
