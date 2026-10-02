//
// Created by mzoll on 30/08/2026.
//


#include <gtest/gtest.h>
#include "tdbscan_algo/connector.h"

#include "tdbscan_algo/common_blibs.h"

using namespace tdbscan;

class Connector_TRUE : public ConnectorSingle<Blib3d> {
public:
  bool eval(const Blib3d &h1, const Blib3d &h2) const override { return true; };
  Connector_TRUE() : ConnectorSingle<Blib3d>("True") {};
};

class Connector_FALSE : public ConnectorSingle<Blib3d> {
public:
  bool eval(const Blib3d &h1, const Blib3d &h2) const override { return false; };
  Connector_FALSE() : ConnectorSingle<Blib3d>("False") {};
};


TEST(ConnectorTest, Connector_Logic) {
  Blib3d b({0, 0, 0}, 0);

  auto con_true = new Connector_TRUE();
  auto con_false = new Connector_FALSE();

  EXPECT_TRUE(con_true->eval(b, b));
  EXPECT_FALSE(con_false->eval(b, b));

  auto t0 = new ConnectorAssembly_AND<Blib3d>();
  t0->addConnector(con_true);
  t0->addConnector(con_true);
  auto t1 = new ConnectorAssembly_AND<Blib3d>();
  t1->addConnector(con_true);
  t1->addConnector(con_false);
  auto t2 = new ConnectorAssembly_AND<Blib3d>();
  t2->addConnector(con_false);
  t2->addConnector(con_true);
  auto t3 = new ConnectorAssembly_AND<Blib3d>();
  t3->addConnector(con_false);
  t3->addConnector(con_false);

  EXPECT_TRUE(t0->eval(b, b));
  EXPECT_FALSE(t1->eval(b, b));
  EXPECT_FALSE(t2->eval(b, b));
  EXPECT_FALSE(t3->eval(b, b));

  auto f0 = new ConnectorAssembly_OR<Blib3d>();
  f0->addConnector(con_true);
  f0->addConnector(con_true);
  auto f1 = new ConnectorAssembly_OR<Blib3d>();
  f1->addConnector(con_true);
  f1->addConnector(con_false);
  auto f2 = new ConnectorAssembly_OR<Blib3d>();
  f2->addConnector(con_false);
  f2->addConnector(con_true);
  auto f3 = new ConnectorAssembly_OR<Blib3d>();
  f3->addConnector(con_false);
  f3->addConnector(con_false);

  EXPECT_TRUE(f0->eval(b, b));
  EXPECT_TRUE(f1->eval(b, b));
  EXPECT_TRUE(f2->eval(b, b));
  EXPECT_FALSE(f3->eval(b, b));
}



