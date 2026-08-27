/**
 * \file Connection.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: DistanceMapService.cxx 129048 2015-02-13 14:42:33Z mzoll $
 * \version $Revision: 129048 $
 * \date $Date: 2015-02-13 15:42:33 +0100 (fre, 13 feb 2015) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/internals/Connection.h"

#include <boost/foreach.hpp>
#include <cmath>

using namespace std;

//============ CLASS Connection =================

Connection::Connection()
{};

Connection::Connection(
  const HashedGeometryConstPtr& hashedGeo)
: hashedGeo_(hashedGeo)
{};

void Connection::Configure(const HashedGeometryConstPtr& hashedGeo)
{ hashedGeo_ = hashedGeo; };

Connection::~Connection() {};


//============== CLASS BoolConnection =====================

template<> const Connection::SpeedRating ConnectionBase<BoolConnection>::evalSpeedRating_(Connection::FAST);

BoolConnection::BoolConnection()
: ConnectionBase(),
  connect_everything_(false) 
{};

BoolConnection::BoolConnection(
  const bool connect_everything)
: ConnectionBase(),
  connect_everything_(connect_everything) 
{};

BoolConnection::BoolConnection(
  const HashedGeometryConstPtr& hashedGeo)
: ConnectionBase(hashedGeo),
  connect_everything_(false)
{};

BoolConnection::BoolConnection(
  const HashedGeometryConstPtr& hashedGeo,
  const bool connect_everything)
: ConnectionBase(hashedGeo),
  connect_everything_(connect_everything)
{};

// implements this function
bool BoolConnection::AreConnected (
  const AbsHit& h1,
  const AbsHit& h2) const
{
  log_trace_stream("Hits are "<<(connect_everything_ ? "" :"NOT ")<<"CONNECTED; because of connection");
  return (connect_everything_);
};

///Are two hits causally connected (DAQ precision)
///\param h1 the one hit
///\param h2 the other hit
bool BoolConnection::AreConnected (
  const AbsDAQHit& h1,
  const AbsDAQHit& h2) const
{
  log_trace_stream("Hits are "<<(connect_everything_ ? "" :"NOT ")<<"CONNECTED; because of connection");
  return (connect_everything_);
};    

bool BoolConnection::CorrectlyConfigured() const
{ return true; };

void BoolConnection::Configure(const HashedGeometryConstPtr& hashedGeo)
{};


//============ CLASS DeltaTimeConnection =================

template<> const Connection::SpeedRating ConnectionBase<DeltaTimeConnection>::evalSpeedRating_(Connection::MEDIUM_FAST);

DeltaTimeConnection::DeltaTimeConnection()
: DTConnection<DeltaTimeConnection>(),
  tresidual_early_(NAN),
  tresidual_late_(NAN)
{};

DeltaTimeConnection::DeltaTimeConnection(
  const HashedGeometryConstPtr& hashedGeo)
: DTConnection<DeltaTimeConnection>(hashedGeo),
  tresidual_early_(NAN),
  tresidual_late_(NAN)
{};

DeltaTimeConnection::DeltaTimeConnection(
  const double tresidual_early,
  const double tresidual_late)
: DTConnection<DeltaTimeConnection>(),
  tresidual_early_(tresidual_early),
  tresidual_late_(tresidual_late)
{};

DeltaTimeConnection::DeltaTimeConnection(
  const HashedGeometryConstPtr& hashedGeo,
  const double tresidual_early,
  const double tresidual_late)
: DTConnection<DeltaTimeConnection>(hashedGeo),
  tresidual_early_(tresidual_early),
  tresidual_late_(tresidual_late)
{};

inline
bool DeltaTimeConnection::Causal(const double dr, const double dt) const
{
  const bool in_time = (-tresidual_early_<=dt) && (dt<=tresidual_late_);
  
  log_debug_stream("Hits are "<<(!in_time ? "NOT " :"")<<"CONNECTED; because of connection");
  return (in_time);
};

bool DeltaTimeConnection::CorrectlyConfigured() const
{ 
  return (! (std::isnan(tresidual_early_) 
            || std::isnan(tresidual_late_))
          && tresidual_early_>=0. 
          && tresidual_late_>=0.);
};

//   if (isnan(tresidual_early_))
//     log_error("tresidual_early is NAN; DeltaTimeConnection might not function as intended");
//   if (isnan(tresidual_late_))
//     log_error("tresidual_late is NAN; DeltaTimeConnection might not function as intended");
//   if (tresidual_early_<0. || tresidual_late_<0.)
//     log_error("tresidual_early and tresidual_late have to be positive numbers");
//   if (tresidual_early_ > tresidual_late_)
//     log_error("tresidual_early is greater than tresidual_late; DeltaTimeConnection might not function as intended");

//============== CLASS DynamicConnection ========================

template<> const Connection::SpeedRating ConnectionBase<DynamicConnection>::evalSpeedRating_(Connection::MEDIUM);

DynamicConnection::DynamicConnection(
  const HashedGeometryConstPtr& hashedGeo)
: DTConnection<DynamicConnection>(hashedGeo),
  speed_(NAN),
  tresidual_early_(NAN),
  tresidual_late_(NAN)
{};

DynamicConnection::DynamicConnection()
: DTConnection<DynamicConnection>(),
  speed_(NAN),
  tresidual_early_(NAN),
  tresidual_late_(NAN)
{};

bool DynamicConnection::Causal(const double dr, const double dt) const
{
  double time_residual;
  if (speed_)
    time_residual = std::abs(dt)-dr/speed_;
  else
    time_residual = std::abs(dt);
  const bool in_time = (-tresidual_early_<=time_residual) && (time_residual<= tresidual_late_);
  
  log_trace_stream("Hits are "<<(!in_time ? "NOT " :"")<<"CONNECTED; because of connection");
  return (in_time);
};

bool DynamicConnection::CorrectlyConfigured() const
{
  return (! (std::isnan(speed_) 
              || std::isnan(tresidual_early_)
              || std::isnan(tresidual_late_))
            && speed_>=0.);
};

//=========== CLASS PhotonDiffusionConnection ===========

template<> const Connection::SpeedRating ConnectionBase<PhotonDiffusionConnection>::evalSpeedRating_(Connection::MEDIUM_SLOW);

#include <boost/math/special_functions/gamma.hpp>
#include <math.h>
#include "dataclasses/I3Constants.h"

const double PhotonDiffusionConnection::c_ice_ = I3Constants::c_ice;

const double PhotonDiffusionConnection::tau_ = 557.E-9; //sec
const double PhotonDiffusionConnection::lambda_s_ = 98.; //m
const double PhotonDiffusionConnection::lambda_a_ = 33.3; //m

const double PhotonDiffusionConnection::const_z_ 
  = 1./PhotonDiffusionConnection::tau_ + PhotonDiffusionConnection::c_ice_/PhotonDiffusionConnection::lambda_a_;

PhotonDiffusionConnection::PhotonDiffusionConnection()
: DTConnection<PhotonDiffusionConnection>(),
  tresidual_early_(0.),
  tresidual_late_(0.),
  lower_cont_quantile_(0.01),
  upper_cont_quantile_(0.9),
  min_pdfvalue_(0.)
{};

PhotonDiffusionConnection::PhotonDiffusionConnection(
  const HashedGeometryConstPtr& hashedGeo)
: DTConnection<PhotonDiffusionConnection>(hashedGeo),
  tresidual_early_(0.),
  tresidual_late_(0.),
  lower_cont_quantile_(0.01),
  upper_cont_quantile_(0.9),
  min_pdfvalue_(0.)
{};

inline
bool PhotonDiffusionConnection::Causal(const double dr, const double dt) const
{
  if (dr==0. && dt==0.)
    return true;
  
  //negative values are physically acausal; the greater the value, the more delayed the hit is
  const double t_res = std::abs(dt)-dr/I3Constants::c_ice;
  
  //if things too acausal make a quick decision
  if ( t_res + tresidual_early_ < 0)
    return false;

  const double tres_PandelPDF_lowerQuantile = 
    (upper_cont_quantile_==1.) ? 0 : intPandelPDF_Quantile_inv(dr, lower_cont_quantile_);
  
  if (t_res + tresidual_early_ < tres_PandelPDF_lowerQuantile) {
    log_debug("hit too early");
    return false;
  }
    
  const double tres_PandelPDF_upperQuantile = 
    (lower_cont_quantile_==0.) ? INFINITY : intPandelPDF_Quantile_inv(dr, upper_cont_quantile_);
  
  if (t_res - tresidual_late_ > tres_PandelPDF_upperQuantile) {
    log_debug("hit too late");
    return false;
  }  
  
  if (min_pdfvalue_ && t_res>=0) {
    if ( PandelPDF(dr, std::abs(dt))< min_pdfvalue_) {
      log_debug("absolute hit probability too low");
      return false;
    }
  }
  
  log_debug("connected");
  return true;
};


double PhotonDiffusionConnection::PandelPDF(const double r, const double tres) {
  assert(r>=0);
  if (tres<0) //no acausal propagation
    return 0;
  
  const double rls = r/lambda_s_;
  const double rla = r/lambda_a_;
  double p = 1./boost::math::tgamma(rls);
  p *= pow(tau_*tres, -rls) / tres;
  p *= exp(-rla -const_z_*tres);
  
  return p;
};

double PhotonDiffusionConnection::intPandelPDF_0_inf(const double r) {
  const double rls = r/lambda_s_;
  const double rla = r/lambda_a_;
  return exp(-rla) * pow(1 + I3Constants::c_ice*tau_/lambda_a_, -rls); //second term= pow(const_z_/tau_,-a);
};

double PhotonDiffusionConnection::intPandelPDF_0_x(const double r, const double x) {
  const double rls = r/lambda_s_;
  const double rla = r/lambda_a_;
  double int_p = boost::math::gamma_p(rls, const_z_*x);
  int_p *= pow(tau_*const_z_,-rls);
  int_p *= exp(-rla);
  return int_p;
};
  
double PhotonDiffusionConnection::intPandelPDF_0_x_inv(const double r, const double prob_val) {
  const double rls = r/lambda_s_;
  const double rla = r/lambda_a_;
  double int_p = pow(tau_*const_z_,-rls) * exp(-rla);
  return boost::math::gamma_p_inv(rls, prob_val/int_p) /const_z_;
};

double PhotonDiffusionConnection::intPandelPDF_Quantile_inv(const double r, const double cont_quantile) {
  if (r==0.) 
    return 0.;
  
  if (cont_quantile <0 || cont_quantile >1 || r<0.) {
    log_error("computation error");
    return NAN;
  }
  //const double tot_Prob = intPandelPDF_0_inf (r);
  //return intPandelPDF_0_x_inv(r, tot_Prob*cont_quantile);
  
  //FASTER implementation: bypassing subfunction calls 
  const double rls = r/lambda_s_;
  return boost::math::gamma_p_inv(rls, cont_quantile)/const_z_;
};

bool PhotonDiffusionConnection::CorrectlyConfigured() const
{
  return (! (std::isnan(tresidual_early_) 
              || std::isnan(tresidual_late_)
              || std::isnan(lower_cont_quantile_)
              || std::isnan(upper_cont_quantile_)
              || std::isnan(min_pdfvalue_))
            && tresidual_early_>=0.
            && tresidual_late_>=0.
            && lower_cont_quantile_>=0.&& lower_cont_quantile_<1.
            && upper_cont_quantile_>0. && upper_cont_quantile_<=1 && upper_cont_quantile_>=lower_cont_quantile_
            && min_pdfvalue_>=0. && min_pdfvalue_<1.);
};

//make all these objects serializable
#if SERIALIZATION_ENABLED
  I3_SERIALIZABLE(Connection);
  I3_SERIALIZABLE(BoolConnection);
//   I3_SERIALIZABLE(DeltaTimeConnection);
//   I3_SERIALIZABLE(DynamicConnection);
//   I3_SERIALIZABLE(PhotonDiffusionConnection);
#endif //SERIALIZATION_ENABLED

