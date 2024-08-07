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
using Masses = Eigen::Array3Xd;
using Stresses = Eigen::MatrixXd;
const double kB = 8.61733e-5;
const double toGPa = 160.2176621;

struct Atoms {
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;
    Masses masses;
    Stresses stresses;
    Eigen::ArrayXd energies;

    Atoms(const Positions_t &p): positions{p}, velocities(3, p.cols()), forces{3, p.cols()},
            masses{3, p.cols()}, stresses{3,3} {
        masses.setOnes();
        velocities.setZero();
        forces.setZero();
        stresses.setZero();
    }

    Atoms(const std::vector<std::basic_string<char>> names, const Positions_t &p): positions{p}, velocities(3, p.cols()),
            forces{3, p.cols()}, masses{3, p.cols()}, stresses{3,3}{
        masses.setOnes();
        velocities.setZero();
        forces.setZero();
        stresses.setZero();
    }

    Atoms(const int t): positions{3,t}, velocities(3, t), forces{3, t},
            masses{3, t}, stresses{3,3}{
        masses.setOnes();
        velocities.setZero();
        forces.setZero();
        stresses.setZero();
    }

    Atoms(const Positions_t &p, const Velocities_t &v) :
            positions{p}, velocities{v}, forces{3, p.cols()},
            masses{3, p.cols()}, stresses{3,3} {
        assert(p.cols() == v.cols());
        masses.setOnes();
        forces.setZero();
        stresses.setZero();
    }

    void set_masses (const double m) {
        masses.setOnes();
        masses *= m;
    }

    void resize(int local) {
        positions.conservativeResize(3, local);
        velocities.conservativeResize(3, local);
        forces.conservativeResize(3, local);
        masses.conservativeResize(3, local);
    }

    double e_shift(Eigen::Vector3d at1, Eigen::Vector3d at2, double sigma, double epsilon) const {
        Eigen::Vector3d d_vec = at1 - at2;
        double r_norm = d_vec.norm();
        double london = std::pow((sigma / r_norm), 6);
        double pauli = std::pow(london, 2);
        return (4 * epsilon * (pauli - london));
    }

    double get_temp (const double mass, const Velocities_t &v) {
        return ( mass * (v.colwise().squaredNorm()*0.5).sum()/
                (1.5 * v.cols() * kB) );
    }

    double get_ekin (const double mass, Velocities_t &v) {
        return mass * (v.colwise().squaredNorm()*0.5).sum();
    }

    double get_ekin (const double mass, Velocities_t &v,size_t N_local) {
        auto ekin = mass * (v.colwise().squaredNorm().array()*0.5);
        return ekin(Eigen::seqN(0,N_local)).sum();
    }

    int nb_atoms() const {
        return positions.cols();
    }

    void init_energies(int num) {
        energies.conservativeResize(num);
    }
};

#endif // __ATOMS_H
