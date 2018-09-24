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
#include "ToolZ/Interfaces.h"

#include "dataclasses/I3Constants.h"
#include <sstream>

#include "TestHelpers.h"

TEST_GROUP(ConnectorOLD);

TEST(BLOEDEL) {
  //prepare some faky input objects
  const I3GeometryConstPtr geometry = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
  const I3CalibrationConstPtr calibration(new I3Calibration());
  const I3DetectorStatusConstPtr status(new I3DetectorStatus());
  
  const HashedGeometry hashedGeo(geometry->omgeo);
  
  std::vector<OMKey> omvec;
  for (int str=1; str<=86; str++) {
    for (int om=1; om<=64; om++) {
      omvec.push_back(OMKey(str,om));
    }
  }
  assert( omvec.size() == (86*64) );

  indexmatrix::AsymmetricIndexMatrix_Bool* idxm;

  fileinterfaces::ReadFromFile<indexmatrix::AsymmetricIndexMatrix_Bool>(idxm, "/home/mzoll/saved_HiveSplitter_legacy.conf");
  RelationPtr rel = boost::make_shared<Relation>(hashedGeo.hashService_, idxm);
  
  DynamicConnectionPtr con_A = boost::make_shared<DynamicConnection>(hashedGeo.distService_);
  con_A->speed_ = I3Constants::c;
  con_A->tresidual_early_ = 300.;
  con_A->tresidual_late_ = 300.;
  
  DynamicConnectionPtr con_B = boost::make_shared<DynamicConnection>(hashedGeo.distService_);
  con_B->speed_ = I3Constants::c_ice;
  con_B->tresidual_early_ = 200.;
  con_B->tresidual_late_ = 200.;
  
  DeltaTimeConnectionPtr con_C = boost::make_shared<DeltaTimeConnection>(hashedGeo.distService_);
  con_C->tresidual_early_ = 200.;
  con_C->tresidual_late_ = 200.;
  
  
  ConnectorBlock* connectorBlock_save = new ConnectorBlock(hashedGeo);
  connectorBlock_save->AddConnector(boost::make_shared<Connector>("ParticleCausal",
                                                                          hashedGeo,
                                                                          (ConnectionPtr)con_A,
                                                                          rel));
  connectorBlock_save->AddConnector(boost::make_shared<Connector>("PhotonCausal",
                                                                          hashedGeo,
                                                                          (ConnectionPtr)con_B,
                                                                          rel));
  connectorBlock_save->AddConnector(boost::make_shared<Connector>("StaticCausal",
                                                                          hashedGeo,
                                                                          (ConnectionPtr)con_C,
                                                                          rel));
  
//   connectorBlock_save->ResetHelperServices(distService);
  
  fileinterfaces::WriteToFile(connectorBlock_save, "/home/mzoll/saved_HiveSplitter_legacy_CB.conf");
};