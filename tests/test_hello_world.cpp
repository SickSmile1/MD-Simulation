// Include a library file to make sure proper includes are set
#include "hello.h"
#include "verlet.h"
#include "math.h"
#include <gtest/gtest.h>

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
    // Testing if we can call a function from our MD library
    hello_eigen();
}

TEST(Verlet, ms1) {
    // s = F/(2m)*t^2+v0*t+s0
    double mass = 1.994e-13; // carbon-12
    int nb_steps = 10;
    double v0 = -1., s0 = 10.5;
    
    double x = s0, y = 1., z = 0.;
    double vx = v0, vy = 0., vz = 0.;
    double fx = 1., fy = 1., fz = 1;
    
    double res = 1/2*mass* std::pow(nb_steps,2)+v0*nb_steps+s0;
    for (int i = 1; i <= nb_steps;i++) {
        verlet_step1(x, y, z, vx, vy, vz, fx, fy, fz, 1., mass);
        // ... compute forces here ... //
        verlet_step2(vx, vy, vz, fx, fy, fz, 1., mass);
    }
    ASSERT_NEAR(res, x, 1e-13);
}
