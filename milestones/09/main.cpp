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

void center(Atoms &at, Eigen::Array3d dom_length) {
    Eigen::Array2d center;
    center << (dom_length[0]/2) , (dom_length[1]/2);
    at.positions.block(0,0,2,at.nb_atoms()).colwise() +=
            center -at.positions.block(0,0,2,at.nb_atoms()).rowwise().mean();
    std::cout << center -at.positions.block(0,0,2,at.nb_atoms()).rowwise().mean() << std::endl;
}

int main(int argc, char** argv) {
    auto [names, positions, velocities]{read_xyz_with_velocities("whisker_large.xyz")};
    Atoms at(positions);
    auto [maxPos, minPos] = get_max_pos(at.positions);
    Eigen::Array3d a(3,1);
    a << maxPos[0]*2, maxPos[1]*2, std::ceil(maxPos[2]+1);// 160, 160, 144.249;
    std::cout << "doin big af whhisker with "<<a[2]<<" Lz length" << std::endl;
    // at.velocities.setRandom();
    // at.velocities *= .1e-6;
    const double timestep = 8; double T;
    const double fixed_mass = 196.96657 / 0.009649;

    double const strain_per_step = 0.000001;
    double const max_strain = 10.;
    center(at,a);
    std::ofstream traj("traj_big_whisk.xyz");
    std::ofstream temp("temp");
    std::ofstream ener("ener");
    write_xyz(traj, at);
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  { a[0], a[1], a[2]},
                  {1,1,12},
                  {0, 0, 1});
    domain.enable(at);
    domain.update_ghosts(at,10.);
    NeighborList nl;
    nl.update(at,5);
    std::cout << "new values1: " << domain.nb_local() << std::endl;
    ducastelle(at, nl);

    for(int i = 1; i < 50001; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
        domain.exchange_atoms(at);
        domain.update_ghosts(at,10.);
        nl.update(at,5);
        ducastelle(at, nl);
        verlet_step2(at.velocities, at.forces, timestep, fixed_mass);

        auto [maxPos, minPos] = get_max_pos(at.positions);
        Eigen::Array3d dom_length{ maxPos[0]-minPos[0],maxPos[1]-minPos[1],a[2]};
        if(i%1000==0 && i > 1) {
            double epot_tot{MPI::allreduce(ducastelle(at,nl, domain, dom_length.prod(),5.),
                                           MPI_SUM, MPI_COMM_WORLD)};
            double ekin_total{MPI::allreduce(at.get_ekin(fixed_mass,
                                             at.velocities,domain.nb_local()),
                                             MPI_SUM, MPI_COMM_WORLD)};
            int n_atoms{MPI::allreduce(domain.nb_local(),
                                          MPI_SUM, MPI_COMM_WORLD)};

            domain.disable(at);
            if (rank ==0) {
                std::cout << a[2] << std::endl;
                T = ( ekin_total / (1.5 * n_atoms * kB) );
                write_xyz(traj, at);
                temp << T << "\n";
                ener << ekin_total+epot_tot << "\n";
                std::cout << "epot: " << epot_tot << " ekin: " << ekin_total << " T: " << T << std::endl;
            }
            domain.enable(at);
            domain.update_ghosts(at,10.);
            nl.update(at,5);
            domain.exchange_atoms(at);
        }
        // std::cout << /*at.velocities.colwise() -= */at.velocities.rowwise().mean() << std::endl;
        // at.velocities.colwise() -= at.velocities.rowwise().mean();
        berendsen_thermostat(at, domain,20, timestep, 1000, fixed_mass);
        a[2] += strain_per_step * a[2];
        domain.scale(at, a);
        domain.exchange_atoms(at);
        domain.update_ghosts(at,10.);
        nl.update(at, 5.);
        ducastelle(at, nl);
    }
    MPI_Finalize();
    temp.close();
    ener.close();
    traj.close();
}