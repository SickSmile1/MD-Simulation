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
#include <string>

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
    // std::cout << center -at.positions.block(0,0,2,at.nb_atoms()).rowwise().mean() << std::endl;
}

void stretch(std::string name, int rank, double temp, double max_strain) {
    auto [names, positions]{read_xyz(name+".xyz")};
    Atoms at(positions);
    int atom_size = at.nb_atoms();
    auto [maxPos, minPos] = get_max_pos(at.positions);
    Eigen::Array3d a(3,1);
    a << maxPos[0]*2, maxPos[1]*2, std::ceil(maxPos[2]-1);// 160, 160, 144.249;
    center(at, a);
    // std::cout << "doin big whhisker with "<<a[2]<<" Lz length" << std::endl;
    // at.velocities.setRandom();
    // at.velocities *= .1e-6;
    int steps = 100000;
    const double timestep = 8; double T;
    const double fixed_mass = 196.96657 / 0.009649;

    double strain = max_strain/steps;
    double accumulated_strain = 0;
    center(at,a);
    std::string s_temp = std::to_string(int(temp)), s_strain = std::to_string(int(max_strain));
    std::ofstream traj(name+"traj_"+s_temp+"_"+s_strain+".xyz");
    std::ofstream temps(name+"temp_"+s_temp+"_"+s_strain);
    std::ofstream ener(name+"ener_"+s_temp+"_"+s_strain);
    std::ofstream stress(name+"stress_"+s_temp+"_"+s_strain);
    // write_xyz(traj, at);
    Domain domain(MPI_COMM_WORLD,
                  { a[0], a[1], a[2]},
                  {1,1,5},
                  {0, 0, 1});

    domain.enable(at);
    domain.update_ghosts(at,10.);
    NeighborList nl;
    nl.update(at,5);
    // std::cout << "rank "<<rank<<" has domain length: " << domain.nb_local() << std::endl;
    ducastelle(at, nl);

    for(int i = 1; i < steps; i++) {
      verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
      domain.exchange_atoms(at);
      domain.update_ghosts(at,10.);
      nl.update(at,5);
      ducastelle(at, nl);
      verlet_step2(at.velocities, at.forces, timestep, fixed_mass);

      auto [maxPos, minPos] = get_max_pos(at.positions); // TODO global at.pos
      Eigen::Array3d dom_length{ maxPos[0]-minPos[0],maxPos[1]-minPos[1],a[2]};
      if(i%1000==0 && i > 1) {
        double epot_tot{MPI::allreduce(ducastelle(at,nl, domain, dom_length.prod()),
                                       MPI_SUM, MPI_COMM_WORLD)};
        double ekin_total{MPI::allreduce(at.get_ekin(fixed_mass,
                                                     at.velocities,domain.nb_local()),
                                         MPI_SUM, MPI_COMM_WORLD)};
        int n_atoms{MPI::allreduce(domain.nb_local(),
                                   MPI_SUM, MPI_COMM_WORLD)};

        domain.disable(at);
        if (rank ==0) {
          // at.stresses = at.stresses*(1.6/10e-1);
          // std::cout << "stresses: "<< at.stresses << std::endl;
          // std::cout << a[2] << std::endl;
          T = ( ekin_total / (1.5 * n_atoms * kB) );
          write_xyz(traj, at);
          temps << T << "\n";
          ener << ekin_total+epot_tot << "\n";
          stress << at.stresses(2,2)<<"\t" << accumulated_strain << "\n";
          // std::cout << "epot: " << epot_tot << " ekin: " << ekin_total << " T: " << T << std::endl;
        }
        domain.enable(at);
        domain.update_ghosts(at,10.);
        nl.update(at,5);
        domain.exchange_atoms(at);
      }
      berendsen_thermostat(at, domain,temp, timestep, temp, fixed_mass);
      if (rank == 0) {
        accumulated_strain = strain * a[2];
        // std::cout << strain*a[2] << std::endl;
        a[2] += strain * a[2];
      }
      // if (accumulated_strain >= max_strain) i = steps+1;
      domain.scale(at, a);
      domain.exchange_atoms(at);
      domain.update_ghosts(at,10.);
      nl.update(at, 5.);
      dom_length[2] = a[2];
      ducastelle(at, nl, domain,dom_length.prod());
      // if (rank == 0 && i%10==0) std::cout << "stresses: "<<at.stresses(2,2) << std::endl;
    }
    domain.disable(at);
    if (rank==0)write_xyz(traj,at);
    temps.close();
    ener.close();
    traj.close();
    stress.close();
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  stretch("whisker_r20", rank, 20., 20.);
  stretch("whisker_r20", rank, 20.,10.);

  stretch("whisker_r25", rank, 20., 20.);
  stretch("whisker_r25", rank, 20.,10.);
  
  //stretch("whisker_large", rank, 20., 20.);
  //stretch("whisker_large", rank, 200, 20.);


  stretch("whisker_small", rank, 20, 20.);
  // stretch("whisker_small", rank, 200, 20.);
  stretch("whisker_small", rank, 20, 10.);
  // stretch("whisker_small", rank, 200, 10.);

  MPI_Finalize();
}
