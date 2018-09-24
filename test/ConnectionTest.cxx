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
#include "ToolZ/Hitclasses.h"

#include "TestHelpers.h"

TEST_GROUP(Connection);

//prepare some faky input objects
const I3GeometryConstPtr geometry = boost::make_shared<const I3Geometry>(IC86Topology::Build_IC86_Geometry());
const I3CalibrationConstPtr calibration(new I3Calibration());
const I3DetectorStatusConstPtr status(new I3DetectorStatus());

const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(geometry->omgeo);

TEST(BoolConnection) {
  BoolConnection bc(hashedGeo);
  
  //some dummy hits
  bc.connect_everything_ = false;
  
  ENSURE(! bc.AreConnected(AbsHit(0, 0.), AbsHit(1, INFINITY)));
  
  bc.connect_everything_ = true;
  
  ENSURE( bc.AreConnected(AbsHit(0, 0.), AbsHit(1, INFINITY)));
}


TEST(DynamicConnection) {
  DynamicConnection dc(hashedGeo);

  dc.speed_=0.;
  dc.tresidual_early_=0.;
  dc.tresidual_late_=0.;
  
  //hits at the same time distance 0
  ENSURE( dc.AreConnected(AbsHit(0, 0.), AbsHit(0, 0.)));

  //late hits can not be included by tresidual_early
  dc.tresidual_early_=1.;
  dc.tresidual_late_=0.;
  ENSURE( ! dc.AreConnected(AbsHit(0, 0.), AbsHit(0, 1.)));
  
  dc.tresidual_early_=0.;
  dc.tresidual_late_=1.;
  ENSURE( dc.AreConnected(AbsHit(0, 0.), AbsHit(0, 1.)));
  
  //let it move
  dc.speed_ = 1.;
  dc.tresidual_early_=0.;
  dc.tresidual_late_=0.;
  
  const double dist = hashedGeo->GetDistService()->GetDistance(0,1);
  const double time = dist/dc.speed_;
  ENSURE( dc.AreConnected(AbsHit(0, 0.), AbsHit(1, time)));
  
  dc.tresidual_early_=1.;
  dc.tresidual_late_=0.;
  ENSURE( dc.AreConnected(AbsHit(0, 0.), AbsHit(1, time-1)));
  
  dc.tresidual_early_=0.;
  dc.tresidual_late_=1.;
  ENSURE( dc.AreConnected(AbsHit(0, 0.), AbsHit(1, time+1)));
}

TEST(DeltaTimeConnection) {
  DeltaTimeConnection dtc(hashedGeo);
  
  dtc.tresidual_early_ = 0.;
  dtc.tresidual_late_ = 0.;
  ENSURE( dtc.AreConnected(AbsHit(0, 0.), AbsHit(0, 0.)));
  ENSURE( dtc.AreConnected(AbsHit(0, 0.), AbsHit(1, 0.)));
  
  dtc.tresidual_early_ = 1.;
  ENSURE( ! dtc.AreConnected(AbsHit(0, 0.), AbsHit(0, 1.)));
  ENSURE( ! dtc.AreConnected(AbsHit(0, 0.), AbsHit(1, 1.)));
  
  dtc.tresidual_late_ = 1.;
  ENSURE( dtc.AreConnected(AbsHit(0, 0.), AbsHit(0, 1.)));
  ENSURE( dtc.AreConnected(AbsHit(0, 0.), AbsHit(1, 1.)));
}

#include "dataclasses/I3Constants.h"

TEST(PhotonDiffusionConnection) {
  PhotonDiffusionConnection pdc(hashedGeo);
  
  pdc.tresidual_early_=0.;
  pdc.tresidual_late_=0.;
  pdc.lower_cont_quantile_=0.0;
  pdc.upper_cont_quantile_=1.;
  pdc.min_pdfvalue_=0.;
  
  //just make a photon hit
  pdc.tresidual_early_ = 0.;
  pdc.tresidual_late_ = 0.;
  ENSURE( pdc.AreConnected(AbsHit(0, 0.), AbsHit(0, 0.)));
  
  const double speed = I3Constants::c_ice;
  const double dist = hashedGeo->GetDistService()->GetDistance(0,1);
  const double time = dist/speed;
  ENSURE( pdc.AreConnected(AbsHit(0, 0.), AbsHit(1, time)));
  
  //things cannot go faster than photon pase-nominal speed
  ENSURE( ! pdc.AreConnected(AbsHit(0, 0.), AbsHit(1, time-1)));

  pdc.lower_cont_quantile_=0.01;
  pdc.upper_cont_quantile_=0.99;
  //things are no longer connected, if not in the time-residual interval [lower-upper], even when physical causal
  ENSURE( ! pdc.AreConnected(AbsHit(0, 0.), AbsHit(1, time)));

  
  
//   ///param: permitted time residual early
//   double tresidual_early_;
//   ///param: permitted time residual late
//   double tresidual_late_;
//   /// lower containment quantile [0.05]
//   double lower_cont_quantile_;
//   /// upper containment quantile [0.9]
//   double upper_cont_quantile_;
//   /// minimal required pdf-value ATM DISABLED
//   double min_pdfvalue_;
};

//serialize by org pointer
#if SERIALIZATION_ENABLED
TEST(Serialize_raw_ptr_BoolConnection){
  BoolConnection* con_save= new BoolConnection(hashedGeo, false);
  BoolConnection* con_load = nullptr;

  serialize_object(con_save, con_load);
  delete con_load;
};

TEST(Serialize_boost_shared_ptr_BoolConnection){
  BoolConnectionPtr con_save= boost::make_shared<BoolConnection>(hashedGeo, false);
  BoolConnectionPtr con_load;

  serialize_object(con_save, con_load);
};
#endif //SERIALIZATION_ENABLED 

//serialize by base-pointer
#if SERIALIZATION_ENABLED
TEST(Serialize_raw_ptr_trough_base_ptr){
  Connection* con_save= new BoolConnection(hashedGeo, false);
  Connection* con_load = nullptr;
  
  serialize_object(con_save, con_load);
  delete con_load;
};

TEST(Serialize_boost_shared_ptr_trough_base_ptr){
  ConnectionPtr con_save= boost::make_shared<BoolConnection>(hashedGeo, false);
  ConnectionPtr con_load;
  
  serialize_object(con_save, con_load);
};
#endif //SERIALIZATION_ENABLED
