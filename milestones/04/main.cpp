//
// Created by ilia on 05/06/24.
//

#include "Atoms.h"
// #include "lj_direct_summation.h"
//  #include "xyz.h"
#include "iostream"
// auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};

int main() {
    std::cout << "A fixed-size array:\n";
    Eigen::Array33f a1 = Eigen::Array33f::Zero();
    std::cout << a1 << "\n\n";
    //Eigen::ArrayXd my1(3);
    //
    //Eigen::Array3d my{2,2,2};
    //Eigen::Vector3d that{my-my1};
    //std::cout << that << std::endl;
    return 1;
}