//
// Created by netsu on 31/08/2026.
//


#include <gtest/gtest.h>
#include "tdbscan_algo/auxil.h"


// Demonstrate some basic assertions.
TEST(Test_Aux, ZeroOne) {
  ZeroOne(0.);
  ZeroOne(1.);
  ZeroOne(.5);

  ASSERT_ANY_THROW(ZeroOne(-0.1));
  ASSERT_ANY_THROW(ZeroOne(1.1));
  ASSERT_EQ(ZeroOne(1.), ZeroOne(1.));
}