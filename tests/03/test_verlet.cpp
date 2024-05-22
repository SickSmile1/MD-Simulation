//
// Created by ilia on 19/04/24.
//
#include <gtest/gtest.h>
#include "verlet.h"
#include "math.h"

TEST(VerletTest, VerletAtomClass) {
    Eigen::Array3Xd p(3, 10);
    // p.setZero();
    p.col(0) << 10. , 1., 0.;
    Atoms at(p);
    at.velocities.col(0) << 0.,-1., 0.;
    at.forces.col(0) << 1.,0.,0.;
    double mass = 1.59;
    int nb_steps = 10;
    double res = at.forces.col(0)[0]/(2*mass) * std::pow(nb_steps,2) +
          at.velocities.col(0)[0] * nb_steps + at.positions.col(0)[0];
    for (int i = 0; i < nb_steps; ++i) {
        std::cout << "Step: " << i << std::endl;
        verlet_step1(at, 1e-3, mass);
        // ... compute forces here ... //
        verlet_step2(at, 1e-3, mass);
    }
    ASSERT_NEAR(res, at.positions.col(0)[0], 1.5e-3);
}