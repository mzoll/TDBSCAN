/**
 * \file Connection.cxx
 *
 * (c) 2013 the IceCube Collaboration
 *
 * $Id: OMKeyHashTest.cxx 148312 2016-07-11 17:41:28Z mzoll $
 * \version $Revision: 148312 $
 * \date $Date: 2016-07-11 19:41:28 +0200 (mån, 11 jul 2016) $
 * \author mzoll <marcel.zoll@fysik.su.se>
 *
 * Unit test to test the 
 */

#include <I3Test.h>

#include "IceHiveZ/internals/Connector.h"

#include "ToolZ/IC86Topology.h"

#include <sstream>

#include "TestHelpers.h"

TEST_GROUP(Connector);

//prepare some faky input objects
const I3GeometryConstPtr geometry = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
const I3CalibrationConstPtr calibration(new I3Calibration());
const I3DetectorStatusConstPtr status(new I3DetectorStatus());

const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(geometry->omgeo);


#include "ToolZ/Interfaces.h"
TEST(Connector){
  ConnectorBlockPtr connectorBlock = boost::make_shared<ConnectorBlock>(hashedGeo);
  connectorBlock->AddConnector(boost::make_shared<Connector>("ConnectNone",
                                                            hashedGeo,
                                                            boost::make_shared<BoolConnection>(hashedGeo, false),
                                                            boost::make_shared<Relation>(hashedGeo->GetHashService(), false)));
  
  //get me two dummy hits
  AbsHit one(2,0.);
  AbsHit two(2,0.);
  
  //probe evaluation
  connectorBlock->Connected(one, two);
    
  //retrieve connectors
  ConnectorPtr con_cum = connectorBlock->GetConnector(-1); //the connector containing the accumulative connector
  ConnectorPtr con_0 =connectorBlock->GetConnector(0);

  ENSURE_EQUAL(con_cum->Connected(one, two), con_0->Connected(one, two));
};


#if SERIALIZATION_ENABLED
TEST(Connector_Serialize_raw_ptr){
  Connector* con_save = new Connector("ConnectNone",
                                      hashedGeo,
                                      boost::make_shared<BoolConnection>(hashedGeo, false),
                                      boost::make_shared<Relation>(hashedGeo->GetHashService(), false));
  Connector* con_load = NULL;
  
  serialize_object(con_save, con_load);
  delete con_save;
  delete con_load;
};

TEST(Connector_Serialize_boost_shared_ptr){
  ConnectorPtr con_save= boost::make_shared<Connector>("ConnectNone",
                                                        hashedGeo,
                                                        boost::make_shared<BoolConnection>(hashedGeo, false),
                                                        boost::make_shared<Relation>(hashedGeo->GetHashService(), false));
  ConnectorPtr con_load;

  serialize_object(con_save, con_load);
};

TEST(ConnectorBlock_Serialize_raw_ptr){
  ConnectorBlock* connectorBlock_save = new ConnectorBlock(hashedGeo);
  connectorBlock_save->AddConnector(boost::make_shared<Connector>("ConnectNone",
                                                                          hashedGeo,
                                                                          boost::make_shared<BoolConnection>(hashedGeo, false),
                                                                          boost::make_shared<Relation>(hashedGeo->GetHashService(), false)));
  ConnectorBlock* connectorBlock_load = NULL;
  
  serialize_object(connectorBlock_save, connectorBlock_load);
  delete connectorBlock_save;
  delete connectorBlock_load;
};

TEST(ConnectorBlock_Serialize_boost_shared_ptr){
  ConnectorBlockPtr connectorBlock_save = boost::make_shared<ConnectorBlock>(hashedGeo);
  connectorBlock_save->AddConnector(boost::make_shared<Connector>("ConnectNone",
                                                            hashedGeo,
                                                            boost::make_shared<BoolConnection>(hashedGeo, false),
                                                            boost::make_shared<Relation>(hashedGeo->GetHashService(), false)));
  ConnectorBlockPtr connectorBlock_load;
  
  serialize_object(connectorBlock_save, connectorBlock_load);
};

#endif //SERIALIZATION_ENABLED
