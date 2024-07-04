//
// Created by ilia on 05/06/24.
//
#include "Atoms.h"
#include "iostream"
#include "Eigen/Dense"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0) {
    // Eigen::RowVector2d d_vec(atoms.nb_atoms());
    double e_pot_all = 0., e_pot = 0., r_norm = 0., pauli = 0., london = 0., force = 0.;
    Eigen::Vector3d d_vec;
    atoms.forces.setZero();
    for(std::size_t i = 1; i < atoms.positions.cols(); i++){
        for(std::size_t j = 0; j < i; j++) {
            d_vec = atoms.positions.col(i) - atoms.positions.col(j);
            r_norm = d_vec.norm();
            london = std::pow((sigma/r_norm),6);
            pauli = std::pow(london,2);
            e_pot = 4*epsilon*(pauli - london);
            force = (24*epsilon)/r_norm*(2*pauli-london);
            atoms.forces.col(i).array() += force*(d_vec.array() / r_norm);
            atoms.forces.col(j).array() -= force*(d_vec.array() / r_norm);
            e_pot_all += e_pot;
        }
    }
    return e_pot_all;
}
