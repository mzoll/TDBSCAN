//
// Created by netsu on 29/08/2026.
//

#include <gtest/gtest.h>

#include "common_defs.h"


// Demonstrate some basic assertions.
TEST(CommonDefs, TimeAdheres) {
  EXPECT_EQ(double(MyTime_t(42)), 42);

  EXPECT_EQ(MyTime_t(42.), MyTime_t(42.) );

  EXPECT_TRUE(MyTime_t::min() < MyTime_t::max() );

  EXPECT_EQ(MyTime_t(42.)- MyTime_t(42.),0. );
  EXPECT_EQ(MyTime_t(42.)- MyTime_t(0.),42. );
}

TEST(CommonDefs, Position3dAdheres) {
  EXPECT_EQ(Position3d(42,42,42), Position3d(42,42,42));

  EXPECT_EQ(Position3d(0,0,0).magnitude(), 0);
  EXPECT_EQ(Position3d(1,1,1).magnitude(), sqrt(3) );

  EXPECT_EQ(Position3d(0,0,0).distance(Position3d(1,1,1)), sqrt(3) );

  EXPECT_TRUE(Position3d(0,0,0) == Position3d(0,0,0));
  EXPECT_TRUE(Position3d(1,2,3) == Position3d(1,2,3));
  EXPECT_FALSE(Position3d(42,2,1) == Position3d(1,2,3));

  EXPECT_FALSE(Position3d(0,0,0) < Position3d(0,0,0));
  EXPECT_TRUE(Position3d(0,0,0) < Position3d(1,0,0));
  EXPECT_FALSE(Position3d(1,0,0) < Position3d(0,0,0));
}

// Demonstrate some basic assertions.