// Include a library file to make sure proper includes are set
#include "hello.h"
#include <gtest/gtest.h>
#include "verlet.h"
#include "math.h"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
    // Testing if we can call a function from our MD library
    hello_eigen();
}

TEST(VerletTest, BasicLoop) {
    double x = 11., y = 1., z = 0.;
    double vx = -1., vy = 0., vz = 0.;
    double fx = 1., fy = 0., fz = 0;
    double mass = 1.59;
    int nb_steps = 10;
    double res = fx/(2*mass) * std::pow(nb_steps,2) + vx * nb_steps+x;
    for (int i = 0; i < nb_steps; ++i) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, 1., mass);
        verlet_step2(vx, vy, vz, fx, fy, fz, 1., mass);
    }
    ASSERT_NEAR(res, x, 1.5e-3);
}