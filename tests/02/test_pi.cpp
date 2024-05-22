//
// Created by ilia on 19/04/24.
//
# define pi           3.14159265358979323846  /* pi */
#include <hello.h>
#include <math.h>
#include <gtest/gtest.h>

TEST(SinTest, IntegerMultiplesOfPi) {
  EXPECT_EQ(sin(0), 0);
  EXPECT_EQ(sin(pi), 0);
  EXPECT_EQ(sin(2+pi), 0);
}

TEST(SinTest2, IntegerMultiplesOfPi) {
  EXPECT_NEAR(sin(0), 0, 1e-6);
  EXPECT_NEAR(sin(pi), 0, 1e-6);
  EXPECT_NEAR(sin(2+pi), 0, 1e-6);
}
