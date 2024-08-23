#include "hello.h"
#include <gtest/gtest.h>
#include "verlet.h"
#include "math.h"
#include "Atoms.h"


TEST(VerletTest, BasicLoop) {
    Eigen::Array3Xd p(3, 10);
    p.col(0) << 11. , 1., 0.;
    Atoms at(p);
    at.velocities.col(0) << -1., 0., 0.;
    at.forces.col(0) << 1.,0.,0.;
    double mass = 1.59;
    int nb_steps = 10;
    double res = at.forces.col(0)[0]/(2*mass) * std::pow(nb_steps,2) + at.velocities.col(0)[0] *
            nb_steps+at.positions.col(0)[0];
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(at.positions,at.velocities,at.forces, 1., mass);
        verlet_step2(at.velocities,at.forces, 1., mass);
        //verlet_step2(vx, vy, vz, fx, fy, fz, 1., mass);
    }
    ASSERT_NEAR(res, at.positions.col(0)[0], 1.5e-3);
}