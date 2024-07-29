//
// Created by ilia on 26/07/24.
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

int main(int argc, char** argv) {
    auto [names, positions]{read_xyz("cluster_3871.xyz")};

    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  {30, 30, 30},
                  {1, 1, 10},
                  {1, 1, 1});
    for (int i = 2; i < 21; i+=2){
        Atoms at(positions);
        const double timestep = i;
        double T;
        const double mass = 196.96657 / 0.009649;
        //at.velocities.setRandom();
        //at.velocities.colwise().normalize();
        //at.velocities *= std::sqrt((3 * kB * 100) / mass);
        // center
        Eigen::Array3d center;
        // center << 30 , 30, 30;
        at.positions.block(0,0,3,at.nb_atoms()).colwise() +=
                center -at.positions.block(0,0,3,at.nb_atoms()).rowwise().mean();
        // std::ofstream traj("traj_big_cluster.xyz");
        std::ofstream temp("temp_"+std::to_string(i));
        std::ofstream ener("ener_"+std::to_string(i));
        // write_xyz(traj, at);

        domain.enable(at);
        domain.update_ghosts(at, 11.);
        NeighborList nl;
        nl.update(at, 5.5);
        std::cout << "Domein length of Rank " << rank << " is: " << domain.nb_local() << std::endl;
        ducastelle(at, nl);

        for (int i = 1; i < 10001; i++) {
            verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
            domain.exchange_atoms(at);
            domain.update_ghosts(at, 11.);
            nl.update(at, 5.5);
            ducastelle(at, nl,11.);
            verlet_step2(at.velocities, at.forces, timestep, mass);

            if (i % 10 == 0 && i > 1) {
                double epot_tot{MPI::allreduce(ducastelle(at, nl, domain,0,11.),
                                               MPI_SUM, MPI_COMM_WORLD)};
                double ekin_total{MPI::allreduce(at.get_ekin(mass,
                                                             at.velocities, domain.nb_local()),
                                                 MPI_SUM, MPI_COMM_WORLD)};
                int n_atoms{MPI::allreduce(domain.nb_local(),
                                           MPI_SUM, MPI_COMM_WORLD)};

                domain.disable(at);
                if (rank == 0) {
                    // std::cout << "System has "<< n_atoms << " atoms." << std::endl;
                    T = (ekin_total / (1.5 * n_atoms * kB));
                    // write_xyz(traj, at);
                    temp << T << "\n";
                    ener << ekin_total+epot_tot << "\n";
                    std::cout <<  "T: " << T << " total Energy: "<< epot_tot + ekin_total << std::endl;
                }
                domain.enable(at);
                domain.update_ghosts(at, 11.);
                nl.update(at, 5.5);
                domain.exchange_atoms(at);
            }
            // std::cout << /*at.velocities.colwise() -= */at.velocities.rowwise().mean() << std::endl;
            // at.velocities.colwise() -= at.velocities.rowwise().mean();
            // berendsen_thermostat(at, domain,20, timestep, 1000, mass);
        }
        temp.close();
        ener.close();
        domain.disable(at);
    }
    MPI_Finalize();
}
