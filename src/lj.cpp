//
// Created by ilia on 05/06/24.
//
#include "Atoms.h"
#include "iostream"
#include "Eigen/Dense"

double lj_direct_summation(Atoms &atoms, const double epsilon = 1.0, const double sigma = 1.0,
                           const double c_energy = 0.) {
  double e_pot_all = 0., e_pot = 0., r_norm = 0., pauli = 0., london = 0., force = 0.;
  Eigen::Vector3d d_vec;
  // set f to zero before calculating new f
  atoms.forces.setZero();
  for(long i = 1; i < atoms.positions.cols(); i++){
    for(long j = 0; j < i; j++) {
      // distance vector and vector norm
      d_vec = atoms.positions.col(i) - atoms.positions.col(j);
      r_norm = d_vec.norm();
    // calculate london and pauli
      london = std::pow((sigma/r_norm),6);
      pauli = std::pow(london,2);
    // Epot and force from the analytic expression
      e_pot = 4*epsilon*(pauli - london) - c_energy;
      force = (24*epsilon)/r_norm*(2*pauli-london);
      // add f to atoms container
      atoms.forces.col(i).array() += force*(d_vec.array() / r_norm);
      atoms.forces.col(j).array() -= force*(d_vec.array() / r_norm);
      e_pot_all += e_pot;
    }
  }
  return e_pot_all;
}
