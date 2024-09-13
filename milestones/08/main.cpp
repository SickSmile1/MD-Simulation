//
// Created by ilia on 26/07/24.
//


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
    double T = (ekin_total / (1.5 * n_atoms * kB));
    temp << T << "\n";
    ener << ekin_total+epot_tot << "\n";
  }
  domain.enable(at);
  domain.update_ghosts(at, 20.);
  nl.update(at, 10);
  domain.exchange_atoms(at);
}

void checkTimesteps(Domain domain, int rank, std::string temp, int temps) {
  auto [names, positions, velocities]{read_xyz_with_velocities(temp)};
  std::cout << "checking reasonable timesteps for "+temp << std::endl;
  for (int j = 2; j < 33; j+=6){
    Atoms at(positions, velocities);
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

    for (int i = 0; i < 5001; i++) {
      verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
      domain.exchange_atoms(at);
      domain.update_ghosts(at, 20.);
      nl.update(at, 10.);
      ducastelle(at, nl,10.);
      verlet_step2(at.velocities, at.forces, timestep, mass);
      if (i % 100 == 0 && i > 1) {
        startStopWrite(at,nl,domain,mass,traj,ener,temp, rank);
      }
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
                {66.68, 66.68, 66.68},
                {2, 2, 1},
                {1,1,1});
  createClusters(domain, rank);
  checkTimesteps(domain,rank, "equilibrated_600k.xyz", 600);
  checkTimesteps(domain,rank, "equilibrated_1000k.xyz", 1000);
  MPI_Finalize();
}

