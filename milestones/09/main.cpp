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

std::tuple<std::vector<double>,std::vector<double>> get_max_pos(Eigen::Array3Xd &p) {
    std::vector<double> max, min;
    max.push_back(p.row(0).maxCoeff());
    max.push_back(p.row(1).maxCoeff());
    max.push_back(p.row(2).maxCoeff());
    min.push_back(p.row(0).minCoeff());
    min.push_back(p.row(1).minCoeff());
    return std::tuple(max,min);
}

double ramp_up(Atoms &at, const double timestep, double &T, const double fixed_mass,
               NeighborList &nl, std::ofstream& traj, Domain domain, int rank) {
    double epot;
    std::cout << "<<starting ramp up>>: " << std::endl;
    std::ofstream myfile;
    myfile.open ("energy");
    std::ofstream myfile1;
    myfile1.open ("temp");
    for(int i = 0; i < 10000; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        nl.update(at,5.5);
        domain.exchange_atoms(at);
        domain.update_ghosts(at,11.);
        nl.update(at,5.5);
        epot = ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
        if (domain.nb_local() == 0 && i > 1) {
            // domain.disable(at);
            auto [maxPos, minPos] = get_max_pos(at.positions);
            Eigen::Array3d a(3,1);
            a << (minPos[0]+maxPos[0])*1.1, (minPos[0]+maxPos[1])*1.1, maxPos[2]*1.1;
            domain.scale(at, a);
            /*domain.enable(at);
            nl.update(at,5.5);
            domain.exchange_atoms(at);
            domain.update_ghosts(at,11.);*/
            std::cout << "changed domain" << std::endl;
        }
        if (i < 1000) { std::cout << "running berendsen" << std::endl; berendsen_thermostat(at, 750, timestep, 400, fixed_mass); }
        // std::cout << "berendsen finished" << std::endl;

        if(i%100 == 0 && rank == 0) {
            epot = ducastelle(at, nl);
            domain.disable(at);
            std::cout << std::to_string(i)+" steps finished." << std::endl;
            double ekin = at.get_ekin(fixed_mass, at.velocities);
            T = ( ekin / (1.5 * at.positions.cols() * kB) );
            std::cout << "epot: " << epot << " ekin: " << ekin << " all: "<<epot + ekin << std::endl;
            write_xyz(traj, at);
            myfile << ekin+epot << "\n";
            myfile1 << T << "\n";
            at.velocities.colwise() -= at.velocities.rowwise().mean();
            // std::cout << "Temp is: " << T << ", energy: " << ekin+epot << std::endl;
            domain.enable(at);
            nl.update(at,5.5);
            domain.exchange_atoms(at);
            domain.update_ghosts(at,11.);

        }
        // if(i%100==0) {}
    }
    myfile.close();
    myfile1.close();
    return epot;
}

int main(int argc, char** argv) {
    auto [names, positions, velocities]{read_xyz_with_velocities("whisker_small.xyz")};
    Atoms at(positions);
    auto [maxPos, minPos] = get_max_pos(at.positions);
    std::cout << maxPos[2] << std::endl;
    at.velocities.setRandom();
    at.velocities *= .1e-6;
    const double timestep = 4; double T;
    const double fixed_mass = 196.96657 / 0.009649;

    double const strain_per_step = 0.000001;
    double const max_strain = 10.;
    std::ofstream traj("traj2.xyz");
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  {
                   (minPos[0]+maxPos[0])*1.5,
                   (minPos[0]+maxPos[1])*1.5,
                   144.24978336},
                  {1,1,6},
                  {0, 0, 0});
    domain.enable(at);
    domain.update_ghosts(at,11.);
    NeighborList nl;
    nl.update(at,5.5);
    std::cout << "new values1: " << domain.nb_local() << std::endl;

    // double epot = ramp_up(at,timestep, T, fixed_mass, nl, traj, domain, rank);

    for(int i = 0; i < 50000; i++) {
        // at.velocities.colwise() -= at.velocities.rowwise().mean();
        // get rid of ghosts??
        if (domain.nb_local() == 0 && i > 1) {
            auto [maxPos, minPos] = get_max_pos(at.positions);
            Eigen::Array3d a(3,1);
            a << (minPos[0]+maxPos[0])*1.1, (minPos[0]+maxPos[1])*1.1, maxPos[2]*1.1;
            domain.scale(at, a);
            std::cout << "changed domain" << std::endl;
        }
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        domain.exchange_atoms(at);
        domain.update_ghosts(at,11.);
        nl.update(at,5.5);
        double epot = ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);
        if (i < 5002) { berendsen_thermostat(at, 750, timestep, 1000, fixed_mass); }
        if(i%1000==0 && i > 1) {
            epot = ducastelle(at, nl);
            domain.disable(at);
            write_xyz(traj, at);
            nl.update(at, 5.5);
            double ekin_total{MPI::allreduce(at.get_ekin(fixed_mass, at.velocities, domain.nb_local()), MPI_SUM, MPI_COMM_WORLD)};
            //double ekin1 = at.get_ekin(fixed_mass, at.velocities);
            //std::cout << "nblocal: " << ekin << " no local: " << ekin1 << std::endl;
            T = ( ekin_total / (1.5 * at.positions.cols() * kB) );
            if (rank ==0) {
                // std::cout << "globpot" << global_pot << "reduce: " << global_pot+ekin << std::endl;
                std::cout << "epot: " << epot << " ekin: " << ekin_total << " T: " << T << std::endl;
            }
            domain.enable(at);
            nl.update(at,5.5);
            domain.exchange_atoms(at);
            domain.update_ghosts(at,11.);
        }
    }
    MPI_Finalize();
    //std::vector<double> temp_array;
    //double ekin = at.get_ekin(fixed_mass, at.velocities);
    // temp_array.push_back(ekin+epot);
    // temp_array.push_back(ekin+epot);
    // temp_array.push_back(ekin+epot);
    // melt(at,timestep, T, fixed_mass, nl, traj, temp_array);
    traj.close();
}