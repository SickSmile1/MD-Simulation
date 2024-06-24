//
// Created by ilia on 24/04/24.
//

#include "iostream"
#include "Atoms.h"
#include "verlet.h"
#include "Eigen/Dense"

// quickref eigen: https://www.quantstart.com/articles/Eigen-Library-for-Matrix-Algebra-in-C/
int main() {
    Eigen::Array3Xd p{ 3,10};
    // p.setZero();
    p.col(0) << 10. , 1., 0.;
    Atoms at(p);
    at.positions.col(0) << 0.,-1., 0.;
    std::cout << at.positions << std::endl;
    double mass = 1.59;
    int nb_steps = 10;
    for (int i = 0; i < nb_steps; ++i) {
        // std::cout << "Step: " << i << std::endl;
        verlet_step1(at.positions, at.velocities, at.forces, 1e-3, mass);
        // ... compute forces here ... //
        verlet_step2(at.velocities, at.forces, 1e-3);
    }
    return 1;
}