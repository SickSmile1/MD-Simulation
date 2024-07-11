//
// Created by ilia on 28/06/24.
//

#include "vector"
#include "berendsten.h"
#include "verlet.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "iostream"
#include "xyz.h"

double ramp_up(Atoms &at, const double timestep, double &T, const double fixed_mass,
             NeighborList &nl, std::ofstream& traj) {
    double epot;
    std::ofstream myfile;
    myfile.open ("energy");
    std::ofstream myfile1;
    myfile1.open ("temp");
    for(int i = 0; i < 10000; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        nl.update(at,5.5);
        // lj_direct_summation(at, nl, 1, 1, 5.5, fixed_mass);
        at.forces.setZero();
        epot = ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
        if (i < 1000) {berendsen_thermostat(at, 650, timestep, 400, fixed_mass); }
        if(i%10==0) {
            write_xyz(traj, at);
            // std::cout << "Temp is: " << T << ", energy: " << ekin+epot << std::endl;
        }
        if(i%100 == 0) {
            double ekin = fixed_mass*(at.velocities.colwise().squaredNorm()*0.5).sum();
            T = ( ekin / (1.5 * at.positions.cols() * kB) );
            myfile << ekin+epot << "\n";
            myfile1 << T << "\n";
        }
        // if(i%100==0) {std::cout << std::to_string(i)+" steps finished." << std::endl;}
    }
    myfile.close();
    myfile1.close();
    return epot;
}

bool equi_energy(std::vector<double> ta) {
    double res1 = std::abs(ta[ta.size()-1])-std::abs(ta[ta.size()-2]);
    double res2 = std::abs(ta[ta.size()-2])-std::abs(ta[ta.size()-3]);
    std::cout << "ta: " << ta[ta.size()-1] << " ta2: " << ta[ta.size()-2] << "ta3" << ta[ta.size()-3] << std::endl;
    std::cout << "res1: " << res1 << " res2: " << res2 << std::endl;
    if (res1 < 20. && res2 < 20.) {return 1;}
    else { return 0;}
}

void melt(Atoms &at, const double timestep, double &T, const double fixed_mass,
             NeighborList &nl, std::ofstream& traj, std::vector<double> ta) {
    std::ofstream myfile;
    myfile.open ("energy1");
    std::ofstream myfile1;
    myfile1.open ("temp1");
    for(int i = 0; i < 130000; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        nl.update(at,5.5);
        // lj_direct_summation(at, nl, 1, 1, 5.5, fixed_mass);
        at.forces.setZero();
        double epot = ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
        if(i%1000==0) {
            write_xyz(traj, at);
            T = at.get_temp(fixed_mass,at.velocities);
            double ekin = at.get_ekin(fixed_mass,at.velocities);
            ta.push_back(ekin+epot);
            if (equi_energy(ta)) at.velocities *= std::sqrt(1+(20/T));
            if(T >= 1200.) {std::cout << i << std::endl; i = 130001;}
            std::cout << "Temp is: " << T << ", energy: " << ekin+epot << ", velocity change: "<<
                std::sqrt(1+(10/T)) << std::endl;
        }
        if(i%100==0) {
            double ekin = at.get_ekin(fixed_mass,at.velocities);
            myfile << ekin+epot << "\n";
            myfile1 << T << "\n";}
    }
    myfile.close();
    myfile1.close();
}


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
    double epot = ramp_up(at,timestep, T, fixed_mass, nl, traj);
    std::vector<double> temp_array;
    double ekin = at.get_ekin(fixed_mass, at.velocities);
    temp_array.push_back(ekin+epot);
    temp_array.push_back(ekin+epot);
    temp_array.push_back(ekin+epot);
    melt(at,timestep, T, fixed_mass, nl, traj, temp_array);

    traj.close();
}
