//
// Created by ilia on 19/04/24.
//

#ifndef __ATOMS_H
#define __ATOMS_H
#include "cassert"
#include "Eigen/Dense"

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
const double kB = 8.61733e-5;

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(const Positions_t &p): positions{p}, velocities(3, p.cols()), forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const std::vector<std::basic_string<char>> names, const Positions_t &p): positions{p}, velocities(3, p.cols()), forces{3, p.cols()} {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const int t): positions{3,t}, velocities(3, t), forces{3, t} {
        velocities.setZero();
        forces.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()} {
        assert(p.cols() == v.cols());
        forces.setZero();
    }

   int nb_atoms() const {
        return positions.cols();
    }
};

#endif // __ATOMS_H
