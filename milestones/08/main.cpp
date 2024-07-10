//
// Created by ilia on 28/06/24.
//

#include "vector"
#include "berendsten.h"
#include "verlet.h"
#include "domain.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "iostream"
#include "xyz.h"
#include "mpi.h"
/*
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
        if (i < 1000) { berendsen_thermostat(at, 650, timestep, 400, fixed_mass); }
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
            T = at.get_temp(fixed_mass,at.velocities,at.positions);
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
}*/

double get_max_pos(Eigen::Array3Xd &p) {
    std::vector<double> max;
    max.push_back(p.row(0).maxCoeff());
    max.push_back(p.row(1).maxCoeff());
    max.push_back(p.row(2).maxCoeff());
    return *std::max_element(std::begin(max), std::end(max));
}

int main(int argc, char** argv) {
    auto [names, positions, velocities]{read_xyz_with_velocities("cluster_923.xyz")};
    Atoms at(positions);
    double maxPos = get_max_pos(at.positions);
    at.velocities.setRandom();
    at.velocities *= .1e-5;
    const double timestep = 1;
    const double fixed_mass = 196.96657 / 0.009649;

    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  {maxPos*1.1,maxPos*1.1,maxPos*1.1},
                  {2,3,1},
                  {0, 0, 0});

    NeighborList nl;
    nl.update(at,5.5);
    ducastelle(at, nl);

    domain.enable(at);
    domain.update_ghosts(at, 11);
    std::cout << at.nb_atoms() << std::endl;

    for(int i = 0; i < 100000; i++) {
        at.velocities.colwise() -= at.velocities.rowwise().mean();
        // get rid of ghosts??
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        domain.exchange_atoms(at);
        domain.update_ghosts(at,11.);
        nl.update(at,5.5);
        at.forces.setZero();
        double epot = ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
        if(i%1000==0) {
            //double ekin = at.get_ekin(fixed_mass,at.velocities);
            double global_pot = MPI::allreduce(at.energies.sum(), MPI_SUM, MPI_COMM_WORLD);
            //double global_kin = MPI::allreduce(ekin, MPI_SUM, MPI_COMM_WORLD);
            // std::cout << global_kin+global_pot << std::endl;
            domain.disable(at);
            nl.update(at, 5.5);
            double epot = ducastelle(at, nl);
            double ekin = at.get_ekin(fixed_mass, at.velocities);
            if (rank ==0) {
                std::cout << "reduce: " << global_pot+ekin << std::endl;
                std::cout << epot + ekin << std::endl;
            }
            domain.enable(at);
            domain.exchange_atoms(at);
            domain.update_ghosts(at,11.);
            nl.update(at,5.5);
            //write_xyz(traj, at);
        }
    }
    MPI_Finalize();
    /*std::ofstream traj("traj2.xyz");
    NeighborList nl;
    nl.update(at,5.5);
    const double timestep = 1; double T = 0;
    const double fixed_mass = 196.96657 / 0.009649;
    // double epot = ramp_up(at,timestep, T, fixed_mass, nl, traj);
    std::vector<double> temp_array;
    double ekin = at.get_ekin(fixed_mass, at.velocities);
    // temp_array.push_back(ekin+epot);
    // temp_array.push_back(ekin+epot);
    // temp_array.push_back(ekin+epot);
    // melt(at,timestep, T, fixed_mass, nl, traj, temp_array);

    traj.close();*/
}
