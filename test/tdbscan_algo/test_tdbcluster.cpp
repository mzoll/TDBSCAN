//
// Created by netsu on 23/09/2026.
//



#include <gtest/gtest.h>

#include "tdbscan_algo/common_defs.h"
#include "tdbscan_algo/connector.h"

#include "tdbscan_algo/tdbcluster.h"

using namespace tdbscan;

// Demonstrate some basic assertions.
TEST(ClusterTest, ArtificiallyConstructedCase) {

  CausalCluster<ScalarBlib> c1;
  CausalCluster<ScalarBlib> c2;

  ScalarBlib b0(Position1d{0}, {1});
  ScalarBlib b1(Position1d{1}, {1});
  ScalarBlib b2(Position1d{2}, {1});
  ScalarBlib b3(Position1d{3}, {1});
  ScalarBlib b4(Position1d{4}, {1});

  EXPECT_TRUE(c1.isConcruent(c2));
  EXPECT_TRUE(c2.isConcruent(c1));
  EXPECT_TRUE(c1.isSupersetOf(c2));
  EXPECT_TRUE(c2.isSupersetOf(c1));
  EXPECT_TRUE(c1.isSubsetOf(c2));
  EXPECT_TRUE(c2.isSubsetOf(c1));

  EXPECT_EQ(c1.nOverlap(c2), 0);
  EXPECT_EQ(c2.nOverlap(c1), 0);

  c1.insertBlib(b0);

  EXPECT_FALSE(c1.isConcruent(c2));
  EXPECT_FALSE(c2.isConcruent(c1));
  EXPECT_TRUE(c1.isSupersetOf(c2));
  EXPECT_FALSE(c2.isSupersetOf(c1));
  EXPECT_FALSE(c1.isSubsetOf(c2));
  EXPECT_TRUE(c2.isSubsetOf(c1));

  EXPECT_EQ(c1.nOverlap(c2), 0);
  EXPECT_EQ(c2.nOverlap(c1), 0);

  c2.insertBlib(b0);

  EXPECT_TRUE(c1.isConcruent(c2));
  EXPECT_TRUE(c2.isConcruent(c1));
  EXPECT_TRUE(c1.isSupersetOf(c2));
  EXPECT_TRUE(c2.isSupersetOf(c1));
  EXPECT_TRUE(c1.isSubsetOf(c2));
  EXPECT_TRUE(c2.isSubsetOf(c1));

  EXPECT_EQ(c1.nOverlap(c2), 1);
  EXPECT_EQ(c2.nOverlap(c1), 1);

   c1.insertBlib(b1);

   EXPECT_FALSE(c1.isConcruent(c2));
   EXPECT_FALSE(c2.isConcruent(c1));
   EXPECT_TRUE(c1.isSupersetOf(c2));
   EXPECT_FALSE(c2.isSupersetOf(c1));
   EXPECT_FALSE(c1.isSubsetOf(c2));
   EXPECT_TRUE(c2.isSubsetOf(c1));

   EXPECT_EQ(c1.nOverlap(c2), 1);
   EXPECT_EQ(c2.nOverlap(c1), 1);

   c2.insertBlib(b2);

   EXPECT_FALSE(c1.isConcruent(c2));
   EXPECT_FALSE(c2.isConcruent(c1));
   EXPECT_FALSE(c1.isSupersetOf(c2));
   EXPECT_FALSE(c2.isSupersetOf(c1));
   EXPECT_FALSE(c1.isSubsetOf(c2));
   EXPECT_FALSE(c2.isSubsetOf(c1));

   EXPECT_EQ(c1.nOverlap(c2), 1);
   EXPECT_EQ(c2.nOverlap(c1), 1);
}
