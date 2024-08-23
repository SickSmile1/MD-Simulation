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
#include <string>

void createClusters(Domain domain,int rank) {
    auto [names, positions]
        {read_xyz("cluster_3871.xyz")};
    std::cout << "doing heat up of periodic boundary system until 600 and 1100" << std::endl;
    Atoms at(positions);
    const double timestep = 2, mass = 196.96657 / 0.009649;
    at.velocities.setRandom();
    at.velocities.colwise().normalize();
    at.velocities *= 1e-5;

    std::ofstream traj("equilibrated_1000k.xyz");
    std::ofstream traj2("equilibrated_600k.xyz");
    domain.enable(at);
    domain.update_ghosts(at, 20.);
    NeighborList nl;
    nl.update(at, 10);
    ducastelle(at, nl, 10);

    for (int i = 1; i < 5001; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
        domain.exchange_atoms(at);
        domain.update_ghosts(at, 20.);
        nl.update(at, 10.);
        ducastelle(at, nl,10.);
        verlet_step2(at.velocities, at.forces, timestep, mass);
        if(i<3000)berendsen_thermostat(at, domain, 600, timestep, timestep*100, mass);
    }
    domain.disable(at);
    if (rank == 0) {write_xyz(traj2,at); }
    traj2.close();
    
    // run up to 1200k
    domain.enable(at);
    domain.update_ghosts(at, 20.);
    nl.update(at, 10);
    ducastelle(at, nl, 10);

    for (int i = 1; i < 5001; i++) {
        verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
        domain.exchange_atoms(at);
        domain.update_ghosts(at, 20.);
        nl.update(at, 10.);
        ducastelle(at, nl,10.);
        verlet_step2(at.velocities, at.forces, timestep, mass);
        if(i<3000)berendsen_thermostat(at, domain, 1000, timestep, timestep*100, mass);
    }
    domain.disable(at);
    if (rank == 0) {write_xyz(traj,at); }
    traj.close();
}

void startStopWrite(Atoms &at, NeighborList &nl, Domain &domain, double mass, std::ofstream &traj,
                    std::ofstream &ener, std::ofstream &temp, int &rank){  
  double epot_tot{MPI::allreduce(ducastelle(at, nl, domain,0,10.),
                                 MPI_SUM, MPI_COMM_WORLD)};
  double ekin_total{MPI::allreduce(at.get_ekin(mass,
                                               at.velocities, domain.nb_local()),
                                   MPI_SUM, MPI_COMM_WORLD)};
  int n_atoms{MPI::allreduce(domain.nb_local(),
                             MPI_SUM, MPI_COMM_WORLD)};

  domain.disable(at);
  if (rank == 0) {
    // std::cout << "System has "<< n_atoms << " atoms." << std::endl;
    double T = (ekin_total / (1.5 * n_atoms * kB));
    // write_xyz(traj, at);
    temp << T << "\n";
    ener << ekin_total+epot_tot << "\n";
    std::cout <<  "T: " << T << " total Energy: "<< epot_tot + ekin_total << std::endl;
  }
  domain.enable(at);
  domain.update_ghosts(at, 20.);
  nl.update(at, 10);
  domain.exchange_atoms(at);
}

void checkTimesteps(Domain domain, int rank, std::string temp, int temps) {
  auto [names, positions]{read_xyz(temp)};
  std::cout << "checking reasonable timesteps for "+temp << std::endl;
  for (int j = 8; j < 48; j+=8){
    Atoms at(positions);
    const double timestep = j, mass = 196.96657 / 0.009649;
    std::ofstream traj("plot/traj_big_"+std::to_string(temps)+"_"+std::to_string(j)+".xyz");
    write_xyz(traj,at);
    std::ofstream temp("plot/temp_s_"+std::to_string(temps)+"_"+std::to_string(j));
    std::ofstream ener("plot/ener_s_"+std::to_string(temps)+"_"+std::to_string(j));

    domain.enable(at);
    domain.update_ghosts(at, 20.);
    NeighborList nl;
    nl.update(at, 10);
    std::cout << "Domein length of Rank " << rank << " is: " << domain.nb_local() << std::endl;
    ducastelle(at, nl, 10);

    for (int i = 1; i < 5001; i++) {
      verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
      domain.exchange_atoms(at);
      domain.update_ghosts(at, 20.);
      nl.update(at, 10.);
      ducastelle(at, nl,10.);
      verlet_step2(at.velocities, at.forces, timestep, mass);
      if (i % 100 == 0 && i > 1) {
        startStopWrite(at,nl,domain,mass,traj,ener,temp, rank);
      }
      if(i<2000)berendsen_thermostat(at, domain, temps, timestep, 800, mass);
    }
    domain.disable(at);
    if (rank == 0) {write_xyz(traj,at); }
    temp.close();
    ener.close();
    traj.close();
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Domain domain(MPI_COMM_WORLD,
                {46.68, 46.68, 46.68},
                {1, 1, 5},
                {1,1,1});
  createClusters(domain, rank);
  checkTimesteps(domain,rank, "equilibrated_600k.xyz", 600);
  checkTimesteps(domain,rank, "equilibrated_1000k.xyz", 1200);
  MPI_Finalize();
}

/*void old() {
    auto [names, positions]
        {read_xyz("cluster_3871.xyz")};
    std::cout << "doing heat up of periodic boundary system until 1100" << std::endl;
    MPI_Init(&argc, &argv);
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Domain domain(MPI_COMM_WORLD,
                  {46.68, 46.68, 46.68},
                  {1, 1, 5},
                  {1,1,1});
    // for (int i = 2; i < 41; i+=2){
        Atoms at(positions);
        const double timestep = 4;
        double T;
        const double mass = 196.96657 / 0.009649;
        at.velocities.setRandom();
        at.velocities.colwise().normalize();
        at.velocities *= 1e-5;
        //at.velocities *= std::sqrt((3 * kB * 100) / mass);
        // center
        // std::ofstream traj("that_try_"+std::to_string(i)+".xyz");
        // write_xyz(traj,at);
        // Eigen::Array3d center;
        // center << 30,30,30;
        // at.positions.block(0,0,3,at.nb_atoms()).colwise() +=
        //         center -at.positions.block(0,0,3,at.nb_atoms()).rowwise().mean();
        // std::ofstream traj("traj_small_cluster"+std::to_string(i)+".xyz");
        // write_xyz(traj,at);
        // std::ofstream temp("temp_s_"+std::to_string(i));
        // std::ofstream ener("ener_s_"+std::to_string(i));

        std::ofstream traj("equilibrated_1100k.xyz");
        domain.enable(at);
        domain.update_ghosts(at, 20.);
        NeighborList nl;
        nl.update(at, 10);
        std::cout << "Domein length of Rank " << rank << " is: " << domain.nb_local() << std::endl;
        ducastelle(at, nl, 10);

        for (int i = 1; i < 5001; i++) {
            verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
            domain.exchange_atoms(at);
            domain.update_ghosts(at, 20.);
            nl.update(at, 10.);
            ducastelle(at, nl,10.);
            verlet_step2(at.velocities, at.forces, timestep, mass);
            if (i % 100 == 0 && i > 1) {
                double epot_tot{MPI::allreduce(ducastelle(at, nl, domain,0,10.),
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
                    // temp << T << "\n";
                    // ener << ekin_total+epot_tot << "\n";
                    std::cout <<  "T: " << T << " total Energy: "<< epot_tot + ekin_total << std::endl;
                }
                domain.enable(at);
                domain.update_ghosts(at, 20.);
                nl.update(at, 10);
                domain.exchange_atoms(at);
            }
            if(i<2000)berendsen_thermostat(at, domain,1100, timestep, 800, mass);
        }
        // temp.close();
        // ener.close();
        domain.disable(at);
        if (rank == 0) {write_xyz(traj,at); }
        traj.close();
    // }
    MPI_Finalize();
}*/

/*if (i % 500 == 0 && i > 1) {
    double epot_tot{MPI::allreduce(ducastelle(at, nl, domain,0,10.),
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
        //write_xyz(traj, at);
        // temp << T << "\n";
        // ener << ekin_total+epot_tot << "\n";
        std::cout <<  "T: " << T << " total Energy: "<< epot_tot + ekin_total << std::endl;
    }
    domain.enable(at);
    domain.update_ghosts(at, 20.);
    nl.update(at, 10);
    domain.exchange_atoms(at);
}*/
// std::cout << /*at.velocities.colwise() -= */at.velocities.rowwise().mean() << std::endl;
// at.velocities.colwise() -= at.velocities.rowwise().mean();
