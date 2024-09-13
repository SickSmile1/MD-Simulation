//
// Created by ilia on 11/07/24.
//
#include "mpi_support.h"
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
    auto [maxPos, minPos] = get_max_pos(at.positions);
    Eigen::Array3d a(3,1);
    a << maxPos[0]*2, maxPos[1]*2, std::ceil(maxPos[2]);// 160, 160, 144.249;
    center(at, a);
    // std::cout << "doin big whhisker with "<<a[2]<<" Lz length" << std::endl;
    at.velocities.setRandom();
    at.velocities *= .1e-6;
    int steps = 50000;
    const double timestep = 5; double T;
    const double fixed_mass = 196.96657 / 0.009649;

    double strain = max_strain/50000;
    double accumulated_strain = 0;
    center(at,a);
    std::string s_temp = std::to_string(int(temp)), s_strain = std::to_string(int(max_strain*10));
    // std::ofstream traj(name+"traj_"+s_temp+"_"+s_strain+".xyz");
    std::ofstream temps(name+"temp_"+s_temp+"_"+s_strain);
    std::ofstream ener(name+"ener_"+s_temp+"_"+s_strain);
    std::ofstream stress(name+"stress_"+s_temp+"_"+s_strain);
    // write_xyz(traj, at);
    Domain domain(MPI_COMM_WORLD,
                  { a[0], a[1], a[2]},
                  {1,1,20},
                  {0, 0, 1});

    domain.enable(at);
    domain.update_ghosts(at,20.);
    NeighborList nl;
    nl.update(at,10.);
    std::cout << "rank "<<rank<<" has domain length: " << domain.nb_local() << std::endl;
    ducastelle(at, nl);
    int brk = 0, brk1 = 0;

    for(int i = 1; i < steps; i++) {
      verlet_step1(at.positions, at.velocities, at.forces, timestep, fixed_mass);
      domain.exchange_atoms(at);
      domain.update_ghosts(at,20.);
      nl.update(at,10);
      ducastelle(at, nl);
      if (domain.nb_local()<=5){ 
        std::cout <<"<< 0finished computation, local nb_atoms: " << domain.nb_local() << " on rank " 
                  << rank << " >>" << std::endl;
        brk=1;
      }
      brk1 = MPI::allreduce(brk,MPI_SUM,MPI_COMM_WORLD);
      if (brk1 >= 1) i = steps;
      verlet_step2(at.velocities, at.forces, timestep, fixed_mass);

      if(i%1000==0 && i > 1) {
        double epot_tot{MPI::allreduce(ducastelle(at,nl, domain, a.prod()),
                                       MPI_SUM, MPI_COMM_WORLD)};
        double ekin_total{MPI::allreduce(at.get_ekin(fixed_mass,
                                                     at.velocities,domain.nb_local()),
                                         MPI_SUM, MPI_COMM_WORLD)};
        int n_atoms{MPI::allreduce(domain.nb_local(),
                                   MPI_SUM, MPI_COMM_WORLD)};

        domain.disable(at);
        if (rank ==0) {
          std::cout << "<< steps finished: " << i << " >>"<<std::endl;
          T = ( ekin_total / (1.5 * n_atoms * kB) );
          // write_xyz(traj, at);
          temps << T << "\n";
          ener << ekin_total+epot_tot << "\n";
          stress << at.stresses(2,2)<<"\t" << accumulated_strain << "\n";
        }
        domain.enable(at);
        domain.update_ghosts(at,20.);
        nl.update(at,10.);
        // domain.exchange_atoms(at);
        ducastelle(at, nl);
      }
      berendsen_thermostat(at, domain,temp, timestep, 1000, fixed_mass);
      accumulated_strain += strain * a[2];
      a[2] += strain * a[2];
      domain.scale(at, a);
      domain.exchange_atoms(at);
      domain.update_ghosts(at,20.);
      nl.update(at, 10.);
      ducastelle(at, nl, domain, a.prod());
    }
    domain.disable(at);
    // if (rank==0)write_xyz(traj,at);
    temps.close();
    ener.close();
    // traj.close();
    stress.close();
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // stretch("whisker_r20", rank, 20., 20.);
  //stretch("whisker_r20", rank, 20.,10.);

  //stretch("whisker_r25", rank, 20., 20.);
  //stretch("whisker_r25", rank, 20.,10.);
  
  //stretch("whisker_small", rank, 20., 20.);
  //stretch("whisker_small", rank, 200, 20.);

  //stretch("whisker_r20", rank, .00001, .5);
  //stretch("whisker_r20", rank, .00001, .3);
  //stretch("whisker_r20", rank, .00001, .1);
  //stretch("whisker_r20", rank, 100, .5);
  //stretch("whisker_r20", rank, 100, .3);
  //stretch("whisker_r20", rank, 100, .1);

  stretch("whisker_r25", rank, .00001, .5);
  stretch("whisker_r25", rank, .00001, .3);
  stretch("whisker_r25", rank, .00001, .1);
  stretch("whisker_r25", rank, 100, .5);
  stretch("whisker_r25", rank, 100, .3);
  stretch("whisker_r25", rank, 100, .1);

  // stretch("whisker_small", rank, 20, 2.);
  // stretch("whisker_small", rank, 20, 1.5);
  // stretch("whisker_small", rank, 20, 1.);
  // stretch("whisker_small", rank, 100, 2.);
  // stretch("whisker_small", rank, 100, 1.5);
  // stretch("whisker_small", rank, 100, 1.);

  MPI_Finalize();
}
