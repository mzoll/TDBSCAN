//
// Created by mzoll on 29/08/2026.
//

#include <gtest/gtest.h>
#include "tdbscan_algo/common_blibs.h"

using namespace tdbscan;

TEST(CommonBlibs, ScalarBlibAdheres) {
  const ScalarBlib b0({0.}, {0.});
  const ScalarBlib b1({42.}, {1.});

  EXPECT_EQ(b0.timeTo(b0), 0.);
  EXPECT_EQ(b1.timeTo(b1), 0.);
  EXPECT_EQ(b0.timeTo(b1), 1.);
  EXPECT_EQ(b1.timeTo(b0), -1.);

  EXPECT_EQ(b0.distanceTo(b0), 0.);
  EXPECT_EQ(b1.distanceTo(b1), 0.);
  EXPECT_EQ(b0.distanceTo(b1), 42.);
  EXPECT_EQ(b1.distanceTo(b0), -42.); // this actually has sign
}

TEST(CommonBlibs, Blib3dAdheres) {
  const Blib3d b0({0., 0., 0.}, {0.});
  const Blib3d b1({1., 1., 1.}, {42.});

  EXPECT_EQ(b0.timeTo(b0), 0.);
  EXPECT_EQ(b1.timeTo(b1), 0.);
  EXPECT_EQ(b0.timeTo(b1), 42.);
  EXPECT_EQ(b1.timeTo(b0), -42.);

  EXPECT_EQ(b0.distanceTo(b0), 0.);
  EXPECT_EQ(b1.distanceTo(b1), 0.);
  EXPECT_EQ(b0.distanceTo(b1), sqrt(3.));
  EXPECT_EQ(b1.distanceTo(b0), sqrt(3.));
}
