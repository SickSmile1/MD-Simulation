//
// Created by ilia on 19/04/24.
//

#ifndef __ATOMS_H
#define __ATOMS_H

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    Atoms(Positions_t &p): positions{p},velocities(3,p.cols()), forces{3, p.cols()} {
        velocities.setZero()
        forces.setZero();
    }
};

#endif // __ATOMS_H
