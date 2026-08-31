//
// Created by netsu on 27/08/2026.
//

#include <gtest/gtest.h>

#include "tdbscan_algo/common_defs.h"


using namespace tdbscan;


TEST(CommonDefs, Blib4dAheres) {

  Blib4d({1,2,0}, 0);
  EXPECT_EQ(Blib4d({1,2,3}, 42), Blib4d({1,2,3}, 42));
  EXPECT_NE(Blib4d({1,2,3}, 42), Blib4d({3,2,1}, 42));

  //the lesser operation is
  EXPECT_FALSE(Blib4d({1,2,3}, 42) < Blib4d({1,2,3}, 42));
  EXPECT_TRUE(Blib4d({1,2,3}, 0) < Blib4d({1,2,3}, 42));
}


TEST(CommonDefs, ScalarBlib) {

  ScalarBlib({1.}, 0);
  EXPECT_EQ(ScalarBlib({42},42.), ScalarBlib({42.}, 42));
  EXPECT_NE(ScalarBlib({0}, 42), ScalarBlib({42.}, 42));

  //the lesser operation is
  EXPECT_FALSE(ScalarBlib({42}, 42) < ScalarBlib({42}, 42));
  EXPECT_TRUE(ScalarBlib({1}, 0) < ScalarBlib({2}, 42));

  EXPECT_EQ(ScalarBlib({42},42.).timeDiff(ScalarBlib({42.}, 42)), 0);
  EXPECT_EQ(ScalarBlib({42},42.).timeDiff(ScalarBlib({42.}, 0)), 42);
  EXPECT_EQ(ScalarBlib({42},42.).timeDiff(ScalarBlib({42.}, 42)), 0);
  EXPECT_EQ(ScalarBlib({42},42.).timeDiff(ScalarBlib({42.}, 0)), 42);
}