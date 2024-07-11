//
// Created by ilia on 11/07/24.
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

double get_max_pos(Eigen::Array3Xd &p) {
    std::vector<double> max;
    max.push_back(p.row(0).maxCoeff());
    max.push_back(p.row(1).maxCoeff());
    max.push_back(p.row(2).maxCoeff());
    return *std::max_element(std::begin(max), std::end(max));
}

int main(int argc, char** argv) {
    auto [names, positions, velocities]{read_xyz_with_velocities("whisker_small.xyz")};
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