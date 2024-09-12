//
// Created by ilia on 05/06/24.
//
#include "Atoms.h"
#include "neighbors.h"
#include "iostream"
#include "Eigen/Dense"

double lj_direct_summation(Atoms &atoms,const NeighborList neighbors, const double epsilon = 1.0,
                           const double sigma = 1.0, const double cutoff = 1.0, double mass = 1.,
                           const double c_energy = 0.) {
  double e_pot_all = 0., e_pot = 0., r_norm = 0., pauli = 0., london = 0., force = 0.;
  Eigen::Vector3d d_vec;
  // set f to zero before calculating new f
  atoms.forces.setZero();
  for (auto[i,j]: neighbors){
    if(i<j) {
      // set f to zero before calculating new f
      d_vec = atoms.positions.col(i) - atoms.positions.col(j);
      r_norm = d_vec.norm();
      // if atom is in cutoff range
      if(r_norm <= cutoff) {
        // set f to zero before calculating new f
        london = std::pow((sigma / r_norm), 6);
        pauli = std::pow(london, 2);
        // set f to zero before calculating new f
        e_pot = 4 * epsilon * (pauli - london) - c_energy;
        force = (24 * epsilon) / r_norm * (2 * pauli - london);
        // add f to atoms container
        atoms.forces.col(i).array() += force * (d_vec.array() / r_norm);
        atoms.forces.col(j).array() -= force * (d_vec.array() / r_norm);
        e_pot_all += e_pot;
      }
    }
  }
  return e_pot_all;
}
