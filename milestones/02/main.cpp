//
// Created by ilia on 24/04/24.
//

#include "iostream"
#include "verlet.h"

int main() {
    double x = 10., y = 1., z = 0.;
    double vx = 0., vy = -1., vz = 0.;
    double fx = 0., fy = 0., fz = 0;
    double mass = 1.59e-8;
    int nb_steps = 10;
    for (int i = 0; i < nb_steps; ++i) {
        std::cout << "Step: " << i << std::endl;
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, 1., mass);
        // ... compute forces here ... //
        verlet_step2(vx, vy, vz, fx, fy, fz, 1., mass);
    }
    return 1;
}