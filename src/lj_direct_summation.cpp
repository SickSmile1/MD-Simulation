//
// Created by ilia on 05/06/24.
//
#include "Atoms.h"
#include "Eigen/Dense"

double lj_direct_summation(Atoms &atoms, double epsilon = 1.0, double sigma = 1.0) {
    // Eigen::RowVector2d d_vec(atoms.nb_atoms());

    /*for(std::size_t i = 0; i < atoms.pos.cols(); i++){
        for(std::size_t j = 0; j < i; j++) {
            double epot = 0.;
            auto temp_d = atoms.pos.col(i)-atoms.pos.col(j);
            temp_d = temp_d.square();
            std::cout << temp_d;
        }
    }*/
    return 1.;
    // atoms.forc = ((24*epsilon)/ atoms.pos.square() ) * (sigma/);
}
