//
// Created by ilia on 05/06/24.
//
#include "Atoms.h"
#include "iostream"
#include "Eigen/Dense"

void lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0) {
    // Eigen::RowVector2d d_vec(atoms.nb_atoms());
    double e_pot_all = 0., e_pot = 0.;

    for(std::size_t i = 0; i < atoms.positions.cols(); i++){
        for(std::size_t j = 0; j < i; j++) {
            // double epot = 0.;
            Eigen::Vector3d temp_d = atoms.positions.col(i)-atoms.positions.col(j);
            double r_norm = temp_d.norm();
            double pauli = std::pow((12/r_norm),12);
            double london = std::pow((6/r_norm),6);
            e_pot = 4*sigma*(pauli - london);
            double force = (e_pot/r_norm);
            e_pot_all += e_pot;
            // std::cout << "force:" << force << "e_pot" << e_pot_all << std::endl;
        }
    }
    std::cout << e_pot_all << std::endl;
    // return 1.;
    // atoms.forc = ((24*epsilon)/ atoms.pos.square() ) * (sigma/);
}
