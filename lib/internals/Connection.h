/**
 * \file Connection.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: MapService.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#ifndef CONNECTION_H
#define CONNECTION_H

#include "IceHiveZ/__SERIALIZATION.h"

static const unsigned connection_version_ = 0;

#include <boost/make_shared.hpp>

#include "dataclasses/I3Constants.h"

#include "ToolZ/GCDinfo.h"
#include "ToolZ/OMKeyHash.h"
#include "ToolZ/Hitclasses.h"
#include "ToolZ/DistanceService.h"
#include "ToolZ/HashedGeometry.h"

//forward declarations for serialization
#if SERIALIZATION_ENABLED
class BoolConnection;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar, const BoolConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void load_construct_data(
    Archive & ar, BoolConnection * t, const unsigned int file_version);
}};
class DynamicConnection;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar, const DynamicConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void load_construct_data(
    Archive & ar, DynamicConnection * t, const unsigned int file_version);
}};

class DeltaTimeConnection;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar, const DeltaTimeConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void load_construct_data(
    Archive & ar, DeltaTimeConnection * t, const unsigned int file_version);
}};

class PhotonDiffusionConnection;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar, const PhotonDiffusionConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void load_construct_data(
    Archive & ar, PhotonDiffusionConnection * t, const unsigned int file_version);
}};

#endif //SERIALIZATION_ENABLED

//================== class defs =============================
//==========================================================

//================ base CLASS Connection ======================

/** baseclass from which all Connections derive
 * NOTE: Alle derived classes need to implement the AreConnected functions
 */
class Connection {
public: //exposed subinterface
  enum SpeedRating {
    NOTSET = 0,
    SLOW = 1,
    MEDIUM_SLOW = 2,
    MEDIUM = 3, 
    MEDIUM_FAST = 4,
    FAST = 4
  };  
private:  
#if SERIALIZATION_ENABLED
// NOTE: all derived classes need to implement the serialization of their objects
  friend class SERIALIZATION_NS::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED

protected: //property  
  /// the hashed Geometry
  HashedGeometryConstPtr hashedGeo_;
  
public: //constructor/destructor
  ///constructor: brings essential objects with it; hide it
  Connection(const HashedGeometryConstPtr& hashedGeo);  
protected: //hide these for class acting as an interface
  ///blank constructor; fields are not initialized; hide this constructor
  Connection();
public:  
  ///destructor: vrtual as it is rather trivial
  virtual ~Connection();
public: //methods
  ///get the Hasher
  CompactOMKeyHashServiceConstPtr GetHasher() const;  
  
public: //methods need to be overriden
  //the main interface to Connections (unfortunatelly cannot be templated)
  ///Are two hits causally connected
  ///\param h1 the one hit
  ///\param h2 the other hit
  virtual
  bool AreConnected (
    const AbsHit& h1,
    const AbsHit& h2) const =0;
  ///Are two hits causally connected (DAQ precision)
  ///\param h1 the one hit
  ///\param h2 the other hit
  virtual
  bool AreConnected (
    const AbsDAQHit& h1,
    const AbsDAQHit& h2) const =0;
  ///check if any derived class was provided with enough configuration to run correctly
  virtual 
  bool CorrectlyConfigured() const =0;
protected:
  ///configure this service with a hashedGeometry
  virtual
  void Configure(const HashedGeometryConstPtr& hashedGeo);
public:
  ///get the speedRating of this class
  virtual 
  SpeedRating GetSpeedRating() const =0;
};

//these are the base pointers
typedef boost::shared_ptr<Connection> ConnectionPtr;
typedef boost::shared_ptr<const Connection> ConnectionConstPtr;

#if SERIALIZATION_ENABLED
  SERIALIZATION_ASSUME_ABSTRACT(Connection);
#endif //SERIALIZATION_ENABLED

//=================== CLASS ConnectionBase =========  
  
template<class derived> //CRTP
class ConnectionBase : virtual public Connection {
#if SERIALIZATION_ENABLED
// NOTE: all derived classes need to implement the serialization of their objects
  friend class SERIALIZATION_NS::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
private:
  /// attribute an evaluation speed rating to the connection [0-not set, 1-slowest, 5-fastest]
  static const Connection::SpeedRating evalSpeedRating_;
protected: //constructor
  ///constructor: brings essential objects with it
  ConnectionBase();
protected:
  ConnectionBase(const HashedGeometryConstPtr& hashedGeo);
public:
  Connection::SpeedRating GetSpeedRating() const;
};  
template<class derived> const Connection::SpeedRating ConnectionBase<derived>::evalSpeedRating_(Connection::NOTSET);


//=================== CLASS BoolConnection =========

///everything is always connected ; dummy instrument
class BoolConnection : virtual public ConnectionBase<BoolConnection> {
#if SERIALIZATION_ENABLED
// NOTE: all derived classes need to implement the serialization of their objects
  friend class SERIALIZATION_NS::access;
  
  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar, const BoolConnection * t, const unsigned int file_version);
  
  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar, BoolConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
public:
  ///PARAM: should everything be connected?
  bool connect_everything_;
  
protected:
  /// constructor; hide it
  BoolConnection();
  /// fully initialized constructor; hide it
  BoolConnection(
    const bool connect_everything);
public:
  /// ctor
  BoolConnection (const HashedGeometryConstPtr& hashedGeo);
  /// fully initialized constructor
  BoolConnection(
    const HashedGeometryConstPtr& hashedGeo,
    const bool connect_everything);

public: // implements this function
  ///Are two hits causally connected
  ///\param h1 the one hit
  ///\param h2 the other hit
  bool AreConnected (
    const AbsHit& h1,
    const AbsHit& h2) const;
    
  ///Are two hits causally connected (DAQ precision)
  ///\param h1 the one hit
  ///\param h2 the other hit
  bool AreConnected (
    const AbsDAQHit& h1,
    const AbsDAQHit& h2) const;
public:    
  bool CorrectlyConfigured() const;  
public:  
  void Configure(const HashedGeometryConstPtr& hashedGeo);
};

typedef boost::shared_ptr<BoolConnection> BoolConnectionPtr;
typedef boost::shared_ptr<const BoolConnection> BoolConnectionConstPtr;

//=================== CLASS DTConnection =====================

/**
 * Class that implements AreConnected and transposes them to function 'Causal(dr, dt)'
 */
template <class derived>
class DTConnection : virtual public ConnectionBase<derived> {
protected: //properties
  /// hold this pointer closer
  DistanceServiceConstPtr distService_;
protected: //constructor
  ///constructor: brings essential objects with it
  DTConnection();
protected:
  DTConnection(const HashedGeometryConstPtr& hashedGeo);
public: //methods
  void Configure(const HashedGeometryConstPtr& hashedGeo);
public: // implements this function
  ///Are two hits causally connected
  ///\param h1 the one hit
  ///\param h2 the other hit
  bool AreConnected (
    const AbsHit& h1,
    const AbsHit& h2) const;
  ///Are two hits causally connected (DAQ precision)
  ///\param h1 the one hit
  ///\param h2 the other hit
  bool AreConnected (
    const AbsDAQHit& h1,
    const AbsDAQHit& h2) const;
protected:
  // defines this function
  virtual
  bool Causal(const double dr, const double dt) const =0;
};

#if SERIALIZATION_ENABLED
  SERIALIZATION_ASSUME_ABSTRACT(DTConnection);
#endif //SERIALIZATION_ENABLED
  
//=============== CLASS DeltaTimeConnection ===========

/**
 * A static connector probes if hits are with delta_t
 */
class DeltaTimeConnection : public DTConnection<DeltaTimeConnection>{
#if SERIALIZATION_ENABLED
// NOTE: all derived classes need to implement the serialization of their objects
  friend class SERIALIZATION_NS::access;
  
  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar, const DeltaTimeConnection * t, const unsigned int file_version);
  
  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar, DeltaTimeConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
protected:
  /// constructor
  DeltaTimeConnection();
  /// fully initialized constructor
  DeltaTimeConnection(
    const double tresidual_early,
    const double tresidual_late);
public: //constructor/destructor
  ///constructor: brings essential objects with it; hide it
  DeltaTimeConnection(const HashedGeometryConstPtr& hashedGeo);
  /// fully initialized constructor
  DeltaTimeConnection(
    const HashedGeometryConstPtr& hashedGeo,
    const double tresidual_early,
    const double tresidual_late);
public:
  /// the allowed time distance for h2 being earlier than h1
  double tresidual_early_;
  /// the allowed time distance for h2 being later than h1
  double tresidual_late_;

protected: //define these functions
  bool Causal(const double dr, const double dt) const;
  bool CorrectlyConfigured() const;
};

typedef boost::shared_ptr<DeltaTimeConnection> DeltaTimeConnectionPtr;
typedef boost::shared_ptr<const DeltaTimeConnection> DeltaTimeConnectionConstPtr;

//===================== CLASS DynamicConnection ===========

/*
 * A dynamic connector probes for the connection of hits by a time residual:
 * evaluated connected if tres in [tres_min, tres_plus]
 */
class DynamicConnection : public DTConnection<DynamicConnection> {
#if SERIALIZATION_ENABLED
// NOTE: all derived classes need to implement the serialization of their objects
  friend class SERIALIZATION_NS::access;
  
  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar, const DynamicConnection * t, const unsigned int file_version);
  
  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar, DynamicConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
public:
  ///the characteristic speed which governs the connection between DOMs
  double speed_;
  ///a permitted negative time residual, configured as a positive value 
  double tresidual_early_;
  ///a permitted positive time residual
  double tresidual_late_;
protected:
  /// constructor
  DynamicConnection();
public:
  DynamicConnection(
    const HashedGeometryConstPtr& hashedGeo);
protected: //define these functions 
  bool Causal(const double dr, const double dt) const;
  bool CorrectlyConfigured() const;  
};

typedef boost::shared_ptr<DynamicConnection> DynamicConnectionPtr;
typedef boost::shared_ptr<const DynamicConnection> DynamicConnectionConstPtr;


//=============== CLASS PhotonDiffusionConnection ===========

/// A diffuse connector probes if hits are with delta_t
class PhotonDiffusionConnection : public DTConnection<PhotonDiffusionConnection>{
#if SERIALIZATION_ENABLED
// NOTE: all derived classes need to implement the serialization of their objects
  friend class SERIALIZATION_NS::access;
  
  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar, const PhotonDiffusionConnection * t, const unsigned int file_version);
  
  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar, PhotonDiffusionConnection * t, const unsigned int file_version);
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED
private: //internal constants;
  /// the photon propagation speed in ice
  static const double c_ice_;  
  /// scattering time
  static const double tau_;
  /// scattering length
  static const double lambda_s_;
  /// absorption length
  static const double lambda_a_;
  /// helper==1./tau_ + I3Constants::c_ice/lambda_a_;
  static const double const_z_;
public: ///parameters
  ///param: permitted time residual early
  double tresidual_early_;
  ///param: permitted time residual late
  double tresidual_late_;
  /// lower containment quantile [0.05] (hits are on time)
  double lower_cont_quantile_;
  /// upper containment quantile [0.9] (hits are delayed)
  double upper_cont_quantile_;
  /// minimal required pdf-value ATM DISABLED
  double min_pdfvalue_;
  
protected: //constructors
  /// constructor; hidden
  PhotonDiffusionConnection();
public:
  PhotonDiffusionConnection(
    const HashedGeometryConstPtr& hashedGeo);  
protected: //define these functions
  bool Causal(const double dr, const double dt) const;
  bool CorrectlyConfigured() const;
private: //internal calls which are static
  //The pandel-function P is the ProbabilityDensityFunction of time-residual at a distant receiver
  //it's fixed parameter is the distance while the time shoudl be seen as teh running parameter
  //it also takes into account absorption, so that the dt-integral of this function is decreasing as r->inf
  /// evaluate the PandelPDF 
  /// @param distance distance to receiver
  /// @param tres runtime difference between emitter and receiver
  /// @return differential probability density to observe the photon
  static
  double PandelPDF(const double distance, const double tres);
  /// dt-integral of P on interval t=0..inf
  /// @return a arrival probability
  static
  double intPandelPDF_0_inf(const double r);
  /// time-integral of P on interval t=0..x
  static
  double intPandelPDF_0_x(const double r, const double x);
  /// returns x at where integral of P on interval t=0..x (the accumulative probability function left sided) is ==prop_val 
  static
  double intPandelPDF_0_x_inv(const double r, const double prob_val);
  /// returns x at where integral of P makes up this left sided quantile
  /// @return expected time residual at (left sided) quantile of time-residual distribution
  static
  double intPandelPDF_Quantile_inv(const double r, const double cont_quantile);
};

typedef boost::shared_ptr<PhotonDiffusionConnection> PhotonDiffusionConnectionPtr;
typedef boost::shared_ptr<const PhotonDiffusionConnection> PhotonDiffusionConnectionConstPtr;

#if SERIALIZATION_ENABLED
  SERIALIZATION_CLASS_VERSION(Connection, connection_version_);
  SERIALIZATION_CLASS_VERSION(BoolConnection, connection_version_);
  SERIALIZATION_CLASS_VERSION(DynamicConnection, connection_version_);  
  SERIALIZATION_CLASS_VERSION(DeltaTimeConnection, connection_version_);
  SERIALIZATION_CLASS_VERSION(PhotonDiffusionConnection, connection_version_);
#endif //SERIALIZATION_ENABLED

//==============================================================================
//========================== IMPLEMENTATIONS ===================================
//==============================================================================

  
//============== CLASS Connection ====================

inline
CompactOMKeyHashServiceConstPtr Connection::GetHasher() const
  { return hashedGeo_->GetHashService(); };
  
#if SERIALIZATION_ENABLED 
template <class Archive>
void Connection::serialize(Archive & ar, const unsigned version) {
  log_debug("(de)serialized Connection");
};
#endif //SERIALIZATION_ENABLED

//============== CLASS ConnectionBase =====================

#if SERIALIZATION_ENABLED
template <class derived> 
template <class Archive>
void ConnectionBase<derived>::serialize(Archive & ar, const unsigned version) {
  ar & SERIALIZATION_BASE_OBJECT_NVP( Connection );
  log_debug("(de)serialized BaseConnection");
};
#endif //SERIALIZATION_ENABLED

template <class derived>
ConnectionBase<derived>::ConnectionBase()
: Connection() 
{};

template <class derived>
ConnectionBase<derived>::ConnectionBase(const HashedGeometryConstPtr& hashedGeo) 
: Connection(hashedGeo)
{};

template <class derived>
Connection::SpeedRating ConnectionBase<derived>::GetSpeedRating() const
{ return evalSpeedRating_; };

//============== CLASS BoolConnection ====================
#if SERIALIZATION_ENABLED 
template <class Archive>
void BoolConnection::serialize(Archive & ar, const unsigned version) {
  ar & SERIALIZATION_BASE_OBJECT_NVP( ConnectionBase );
  ar & SERIALIZATION_NVP(connect_everything_);
  log_debug("(de)serialized BoolConnection");
};

template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar, const BoolConnection * t, const unsigned int file_version)
{};

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar, BoolConnection * t, const unsigned int file_version)
{
  ::new(t)BoolConnection();
};
#endif //SERIALIZATION_ENABLED

//============== CLASS DTConnection ============================== 

template<class derived>
DTConnection<derived>::DTConnection(
  const HashedGeometryConstPtr& hashedGeo)
: ConnectionBase<derived>(hashedGeo),
  distService_(hashedGeo->GetDistService())
{};

template<class derived>
DTConnection<derived>::DTConnection()
: ConnectionBase<derived>(),
  distService_()
{};

template<class derived>
bool DTConnection<derived>::AreConnected(
  const AbsHit& h1,
  const AbsHit& h2) const 
{
  const double dr= distService_->GetDistance(h1.GetDOMIndex(), h2.GetDOMIndex());
  const double dt=h1.TimeDiff(h2);
  return Causal(dr, dt);    
};

template<class derived>
bool DTConnection<derived>::AreConnected (
  const AbsDAQHit& h1,
  const AbsDAQHit& h2) const
{
  const double dr= distService_->GetDistance(h1.GetDOMIndex(), h2.GetDOMIndex());
  const double dt=h1.TimeDiff(h2);
  return Causal(dr, dt);    
};

template<class derived>
void DTConnection<derived>::Configure(
  const HashedGeometryConstPtr& hashedGeo)
{ 
  Connection::Configure(hashedGeo);
  distService_ = hashedGeo->GetDistService(); 
};

//============== CLASS DeltaTimeConnection ====================
#if SERIALIZATION_ENABLED 
template <class Archive>
void DeltaTimeConnection::serialize(Archive & ar, const unsigned version) {
  ar & SERIALIZATION_BASE_OBJECT_NVP( Connection );
  ar & SERIALIZATION_NVP(tresidual_early_);
  ar & SERIALIZATION_NVP(tresidual_late_);
  log_debug("(de)serialized DeltaTimeConnection");
};

template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar, const DeltaTimeConnection * t, const unsigned int file_version)
{ ar << SERIALIZATION_NS::make_nvp("hashedGeo", t->hashedGeo_); };

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar, DeltaTimeConnection * t, const unsigned int file_version)
{
  HashedGeometryConstPtr hashedGeo;
  ar >> SERIALIZATION_NS::make_nvp("hashedGeo", hashedGeo);
  ::new(t)DeltaTimeConnection(hashedGeo);
};
#endif //SERIALIZATION_ENABLED

//============== CLASS DynamicConnection ====================
#if SERIALIZATION_ENABLED 
template <class Archive>
void DynamicConnection::serialize(Archive & ar, const unsigned version) {
  ar & SERIALIZATION_BASE_OBJECT_NVP( Connection );
  ar & SERIALIZATION_NVP(speed_);
  ar & SERIALIZATION_NVP(tresidual_early_);
  ar & SERIALIZATION_NVP(tresidual_late_);
  log_debug("(de)serialized DynamicConnection");
};

template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar, const DynamicConnection * t, const unsigned int file_version)
{ 
  ar << SERIALIZATION_NS::make_nvp("hashedGeo", t->hashedGeo_);};

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar, DynamicConnection * t, const unsigned int file_version)
{
  HashedGeometryConstPtr hashedGeo;
  ar >> SERIALIZATION_NS::make_nvp("hashedGeo", hashedGeo);
  ::new(t)DynamicConnection(hashedGeo);
};
#endif //SERIALIZATION_ENABLED

//============== CLASS PhotonDiffusionConnection ====================
#if SERIALIZATION_ENABLED 
template <class Archive>
void PhotonDiffusionConnection::serialize(Archive & ar, const unsigned version) {
  ar & SERIALIZATION_BASE_OBJECT_NVP( Connection );
  ar & SERIALIZATION_NVP(tresidual_early_);
  ar & SERIALIZATION_NVP(tresidual_late_);
  ar & SERIALIZATION_NVP(lower_cont_quantile_);
  ar & SERIALIZATION_NVP(upper_cont_quantile_);
  ar & SERIALIZATION_NVP(min_pdfvalue_);
  log_debug("(de)serialized PhotonDiffusionConnection");
};

template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar, const PhotonDiffusionConnection * t, const unsigned int file_version)
{ ar << SERIALIZATION_NS::make_nvp("hashedGeo", t->hashedGeo_); };

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar, PhotonDiffusionConnection * t, const unsigned int file_version)
{
  HashedGeometryConstPtr hashedGeo;
  ar >> SERIALIZATION_NS::make_nvp("hashedGeo", hashedGeo);
  ::new(t)PhotonDiffusionConnection(hashedGeo);
};
#endif //SERIALIZATION_ENABLED

#endif //CONNECTION_H

