//
// Created by ilia on 19/04/24.
//
#include <math.h>
#include <gtest/gtest.h>

TEST(SinTest, IntegerMultiplesOfPi) {
  ASSERT_NEAR(sin(0), 0, 1e-9);
  ASSERT_NEAR(sin(M_PI), 0,1e-9);
  ASSERT_NEAR(sin(2*M_PI), 0,1e-9);
}

TEST(SinTest2, IntegerMultiplesOfPi) {
  EXPECT_NEAR(sin(0), 0, 1e-6);
  EXPECT_NEAR(sin(M_PI), 0, 1e-6);
  EXPECT_NEAR(sin(2*M_PI), 0, 1e-6);
}
