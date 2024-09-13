//
// Created by ilia on 26/07/24.
//

#include "verlet.h"
#include "domain.h"
#include "ducastelle.h"
#include "neighbors.h"
#include "berendsten.h"
#include "iostream"
#include "xyz.h"
#include "mpi.h"

void melt(Domain &domain, Eigen::Array3Xd positions, std::string name, int rank){
  // initiate file write operations
  std::ofstream temp("temp_melt"+name);
  std::ofstream traj("traj_melt"+name);
  std::ofstream ener("ener_melt"+name);
  ener << "total energy\t" << "pot_energy\t" << "kin_energy\t" << "atoms\n";

  // initialization of atoms, timestep, mass,
  // centering the mass,
  // itroducing minimal velocity in x,y,z directions
  Atoms at(positions);
  const double timestep = 5;
  double T;
  const double mass = 196.96657 / 0.009649;
  Eigen::Array3d center;
  center << 40,40,40;
  at.positions.block(0,0,3,at.nb_atoms()).colwise() +=
          center -at.positions.block(0,0,3,at.nb_atoms()).rowwise().mean();
  write_xyz(traj, at);
  at.velocities.setRandom();
  at.velocities.colwise().normalize();
  at.velocities *= 1e-5;
  at.velocities *= std::sqrt((3 * kB * 1) / mass);

  NeighborList nl;
  nl.update(at, 10);
  
  // equilibrate freshly created mockay files
  for (int i = 1; i < 1001; i++) {
    verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
    ducastelle(at, nl, 10.);
    verlet_step2(at.velocities, at.forces, timestep, mass);
    ducastelle(at, nl, 10.);
    berendsen_thermostat(at, 100, timestep, 100, mass);
    // if (i%10==0) write_xyz(traj,at);
  }
  
  domain.enable(at);
  domain.update_ghosts(at, 20.);
  nl.update(at, 10);
  std::cout << "Domein length of Rank " << rank << " is: " << domain.nb_local() << std::endl;
  ducastelle(at, nl, 10);

  for (int i = 1; i < 400001; i++) {
    // main loop
    verlet_step1(at.positions, at.velocities, at.forces, timestep, mass);
    domain.exchange_atoms(at);
    domain.update_ghosts(at, 20.);
    nl.update(at, 10.);
    ducastelle(at, nl, 10.);
    verlet_step2(at.velocities, at.forces, timestep, mass);
    if (i % 100 == 1) {
      // compute energy and atoms over all cores for file write operations
      double epot_tot{MPI::allreduce(ducastelle(at, nl, domain, 0, 10.),
                                     MPI_SUM, MPI_COMM_WORLD)};
      double ekin_total{MPI::allreduce(at.get_ekin(mass,
                                                   at.velocities, domain.nb_local()),
                                       MPI_SUM, MPI_COMM_WORLD)};
      int n_atoms{MPI::allreduce(domain.nb_local(),
                                 MPI_SUM, MPI_COMM_WORLD)};
      T = (ekin_total / (1.5 * n_atoms * kB));
      domain.disable(at);
      if (rank == 0) {
        // write_xyz(traj, at);
        // printouts and file write operations
        temp << T << "\n";
        ener << ekin_total + epot_tot << "\t" << epot_tot << "\t" << ekin_total << "\t" << n_atoms <<"\n";
        std::cout <<"i:"<<i<< " T: " << T << " total Energy: " << epot_tot + ekin_total << std::endl;
        std::cout << std::sqrt(1 + (10 / T)) << std::endl;
      }
      // add energy by changing velocity
      if (T < 1400 && (i % 500) == 1) { at.velocities *= std::sqrt(1 + (20 / T)); } else if(T>1400){ i = 400000;}
      domain.enable(at);
      domain.exchange_atoms(at);
      domain.update_ghosts(at,10.);
      nl.update(at, 5.);
      ducastelle(at, nl);
    }
  }
  domain.disable(at);
  // traj.close();
  temp.close();
  ener.close();
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Domain domain(MPI_COMM_WORLD,
                {80, 80, 80},
                {2, 2, 1},
                {0, 0, 0});

  // run simulations for different sizes
  auto [names, positions]
          {read_xyz("cluster_1415.xyz")};
  melt(domain, positions, "1415", rank);
  auto [names1, positions1]
          {read_xyz("cluster_2057.xyz")};
  melt(domain, positions1, "2057", rank);
  auto [names2, positions2]
          {read_xyz("cluster_2869.xyz")};
  melt(domain, positions2, "2869", rank);
  auto [names3, positions3]
          {read_xyz("cluster_3871.xyz")};
  melt(domain, positions3, "3871", rank);
  MPI_Finalize();
}

