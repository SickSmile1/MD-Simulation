//
// Created by ilia on 25/06/24.
//

#include <gtest/gtest.h>
#include "verlet.h"
#include "Atoms.h"
#include "cubicLatice.h"
#include "Eigen/Dense"
#include "iostream"

TEST(CubicLatice, stucture) {
    GTEST_SKIP();
    // 8 27 64 125 216 343 512 729 1000
    Atoms at(27);
    at.positions.setZero();
    //std::cout << at.positions.size() << std::endl;
    createCubicLatice(at,3);
    //std::cout << at.positions << std::endl;
}
