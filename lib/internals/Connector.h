/**
 * \file Connector.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: MapService.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#ifndef CONNECTOR_H
#define CONNECTOR_H

#include "IceHiveZ/__SERIALIZATION.h"

static const unsigned connector_version_ = 0;

#include <limits>
#include <list>
#include <map>
#include <iostream>
#include <boost/make_shared.hpp>

#include "dataclasses/I3Constants.h"

#include "ToolZ/OMKeyHash.h"
#include "ToolZ/Hitclasses.h"
#include "ToolZ/HashedGeometry.h"
#include "IceHiveZ/internals/Connection.h"
#include "IceHiveZ/internals/Relation.h"

//forward declarations for serialization
#if SERIALIZATION_ENABLED  
class Connector;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar,
    const Connector * t,
    const unsigned int version);

  template<class Archive>
  void load_construct_data(
    Archive & ar,
    Connector * t,
    const unsigned int version);
}};

class ConnectorBlock;

namespace SERIALIZATION_NS_BASE { namespace serialization {
  template<class Archive>
  void save_construct_data(
    Archive & ar,
    const ConnectorBlock * t,
    const unsigned int version);

  template<class Archive>
  void load_construct_data(
    Archive & ar,
    ConnectorBlock * t,
    const unsigned int version);
}};
#endif //SERIALIZATION_ENABLED  


//======================= Connector ============

/**
 * A service which tells you if hits are (causally and topologically) connected
 */
class Connector {
private:
  friend class ConnectorBlock;
#if SERIALIZATION_ENABLED
  friend class SERIALIZATION_NS::access;
  
  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar,
    const Connector * t,
    const unsigned int version);

  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar,
    Connector* t,
    const unsigned int version);
  
  template<class Archive>
  void serialize(Archive & ar, const unsigned int file_version);
#endif //SERIALIZATION_ENABLED  
private:
  ///a unique name for this service
  const std::string name_;
  ///pointer the OMKeyHasher
  const HashedGeometryConstPtr hashedGeo_;
  ///the connector for this service
  const ConnectionPtr connection_;
  ///the Relation for this Service
  const RelationPtr relation_;
  
public:
  ///constructor
  Connector(
    const std::string& name,
    const HashedGeometryConstPtr& hashedGeo,
    const ConnectionPtr& connection, 
    const RelationPtr& relation);
  
  std::string GetName() const;
  CompactOMKeyHashServiceConstPtr GetHashService() const;
  ConnectionPtr GetConnection() const;
  RelationPtr GetRelation() const;
    
  /// Are Hits h1 and h2 connected by being related and connected to each other?
  template <class Hitclass>
  bool Connected (const Hitclass& h1, const Hitclass& h2) const;
};

typedef boost::shared_ptr<Connector> ConnectorPtr;
typedef boost::shared_ptr<const Connector> ConnectorConstPtr;

#if SERIALIZATION_ENABLED
  SERIALIZATION_CLASS_VERSION(Connector, connector_version_);
#endif //SERIALIZATION_ENABLED

//============ CLASS ConnectorBlock ===========

///Holds a number of Connectors and and neccessary Services; is explicitly serializable
class ConnectorBlock {
#if SERIALIZATION_ENABLED  
  friend class SERIALIZATION_NS::access;
  
  template<class Archive>
  friend void SERIALIZATION_NS::save_construct_data(
    Archive & ar,
    const ConnectorBlock * t,
    const unsigned int version);

  template<class Archive>
  friend void SERIALIZATION_NS::load_construct_data(
    Archive & ar,
    ConnectorBlock * t,
    const unsigned int version);

  template<class Archive>
  void serialize(Archive & ar, const unsigned int version);
#endif //SERIALIZATION_ENABLED  
public:
  ///list of connectors
  typedef std::list<ConnectorPtr> ConnectorList;  
private: //property
  ///pointer the OMKeyHasher
  const HashedGeometryConstPtr hashedGeo_;
  ///list of connectors, which are 
  ConnectorList connectorlist_;
  ///the cumulative of all the Connectors of the connectorlist
  const RelationPtr cumulativeRel_;
  
public: //constructors
  /// blank constructor (need to fill this with AddConnector() calls)
  ConnectorBlock(
    const HashedGeometryConstPtr& hashedGeo);

public: //methods
  /// Add a Connector and add the Relation map to the cumRel Map
  void AddConnector (
    const ConnectorPtr& connector);
  
  ///check if to Hits are connected by the any of the connection services
  template <class Hitclass>
  bool Connected(
    const Hitclass& h1,
    const Hitclass& h2) const;
  
  ///diagnose the connections for these hits
  template <class Hitclass>
  void DiagnoseConnected(
    const Hitclass& h1,
    const Hitclass& h2) const;  
  
  ///pass Pointer to CompactOMKeyHasher to external
  CompactOMKeyHashServiceConstPtr GetHashService() const;
  ///retrieve a connector from the ConnectorList; 0 will pass the cumulative one
  ConnectorPtr GetConnector (const int index) const;
  ///Get the complete list of Relations  
  ConnectorList GetConnectorList() const;
};

typedef boost::shared_ptr<ConnectorBlock> ConnectorBlockPtr;
typedef boost::shared_ptr<const ConnectorBlock> ConnectorBlockConstPtr;

#if SERIALIZATION_ENABLED
  SERIALIZATION_CLASS_VERSION(ConnectorBlock, connector_version_);
#endif //SERIALIZATION_ENABLED

//==============================================================================
//========================== IMPLEMENTATIONS ===================================
//==============================================================================


//========================== CLASS Connector =========================
#if SERIALIZATION_ENABLED
template<class Archive>
void Connector::serialize(Archive & ar, const unsigned int file_version)
{
//   ar & SERIALIZATION_NS::make_nvp("name", const_cast<std::string &>(name_));  
//   ar & SERIALIZATION_NS::make_nvp("hasher", const_cast<CompactOMKeyHashServiceConstPtr &>(hasher_));
//   ar & SERIALIZATION_NS::make_nvp("connection", const_cast<ConnectionPtr &>(connection_));
//   ar & SERIALIZATION_NS::make_nvp("relation", const_cast<RelationPtr &>(relation_));
};

// //NOTE (de)serialization overrides
template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar,
  const Connector * t,
  const unsigned int version)
{
  ar << SERIALIZATION_NS::make_nvp("name", t->name_);  
  ar << SERIALIZATION_NS::make_nvp("hashedGeo", t->hashedGeo_);
  ar << SERIALIZATION_NS::make_nvp("connection", t->connection_);
  ar << SERIALIZATION_NS::make_nvp("relation", t->relation_);
};

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar,
  Connector * t,
  const unsigned int version)
{
  std::string name;
  ar >> SERIALIZATION_NS::make_nvp("name", name);
  HashedGeometryConstPtr hashedGeo;
  ar >> SERIALIZATION_NS::make_nvp("hashedGeo", hashedGeo);
  ConnectionPtr connection;
  ar >> SERIALIZATION_NS::make_nvp("connection", connection);
  RelationPtr relation;
  ar >> SERIALIZATION_NS::make_nvp("relation", relation);
  
  ::new(t) Connector(name, hashedGeo, connection, relation);
};
#endif //SERIALIZATION_ENABLED

inline
std::string Connector::GetName() const
{return name_;};

inline
CompactOMKeyHashServiceConstPtr Connector::GetHashService() const
{return hashedGeo_->GetHashService();};

inline 
ConnectionPtr Connector::GetConnection() const
{return connection_;};

inline 
RelationPtr Connector::GetRelation() const
{return relation_;};

template <class Hitclass>
bool Connector::Connected (
  const Hitclass& h1,
  const Hitclass& h2) const 
{
  bool related = relation_->AreRelated(h1.GetDOMIndex(), h2.GetDOMIndex());
  bool connected = connection_->AreConnected(h1, h2);
  if (h1.TimeDiff(h2)==0) { //NOTE this is probably unneccessary but required for generality
    //if hits occure at the exact same time, also check the reverse connection
    related |= relation_->AreRelated(h2.GetDOMIndex(), h1.GetDOMIndex());
    connected |= connection_->AreConnected(h2, h1);
  }
  const bool con = related && connected;
  log_debug_stream(name_<<": Hits are "<<(con ? "CONNECTED" : "NOT connected"));
  return con;
};

//========================== CLASS ConnectorBlock =========================
#if SERIALIZATION_ENABLED
template<class Archive>
void ConnectorBlock::serialize(Archive & ar, const unsigned int version)
{};

// //NOTE (de)serialization overrides
template<class Archive>
inline void SERIALIZATION_NS::save_construct_data(
  Archive & ar,
  const ConnectorBlock * t,
  const unsigned int version)
{
  ar << SERIALIZATION_NS::make_nvp("HashedGeo", t->hashedGeo_);
  ar << SERIALIZATION_NS::make_nvp("ConnectorList", t->connectorlist_);
};

template<class Archive>
inline void SERIALIZATION_NS::load_construct_data(
  Archive & ar,
  ConnectorBlock * t,
  const unsigned int version)
{
  HashedGeometryConstPtr hashedGeo;
  ar >> SERIALIZATION_NS::make_nvp("HashedGeo", hashedGeo);
  ConnectorBlock::ConnectorList connectorlist;
  ar >> SERIALIZATION_NS::make_nvp("ConnectorList", connectorlist);
  
  ::new(t) ConnectorBlock(hashedGeo);
  BOOST_FOREACH(ConnectorPtr& c, connectorlist) {
    log_trace_stream("Adding Connector "<<c->GetName(););
    t->AddConnector(c);
  };
};
#endif //SERIALIZATION_ENABLED

template <class Hitclass>
bool ConnectorBlock::Connected(
  const Hitclass& h1,
  const Hitclass& h2) const
{
  log_debug("Evaluating Connected()");

  bool any_relation = cumulativeRel_->AreRelated(h1.GetDOMIndex(), h2.GetDOMIndex());
  
  if (h1.TimeDiff(h2)==0.) { //NOTE this is probably unneccessary but required for generality
    //if hits occure at the exact same time, also check the reverse connection
    any_relation |= cumulativeRel_->AreRelated(h2.GetDOMIndex(), h1.GetDOMIndex());
    
    ConnectorList::const_iterator connectorlist_iter=connectorlist_.begin();
    const ConnectorList::const_iterator connectorlist_end=connectorlist_.end();
    while (connectorlist_iter != connectorlist_end){
      if ((*connectorlist_iter)->Connected(h1, h2)) {
        log_debug_stream("Hits are CONNECTED; evaluation of connector "<<(*connectorlist_iter)->GetName());
        return true;
      }
      if ((*connectorlist_iter)->Connected(h2, h1)) {
        log_debug_stream("Hits are CONNECTED; evaluation of connector "<<(*connectorlist_iter)->GetName());
        return true;
      }
      ++connectorlist_iter;
    }
  }
  
  if (! any_relation) {
    // none of the connectionServices hold a connection for this particular pair of DOMs
    log_debug("Hits are NOT connected; evaluation of cumulative connector");
    return false;
  }
  
  ConnectorList::const_iterator connectorlist_iter=connectorlist_.begin();
  const ConnectorList::const_iterator connectorlist_end=connectorlist_.end();
  while (connectorlist_iter != connectorlist_end){
    if ((*connectorlist_iter)->Connected(h1, h2)) {
      log_debug_stream("Hits are CONNECTED; evaluation of connector "<<(*connectorlist_iter)->GetName());
      return true;
    }
    ++connectorlist_iter;
  }

  log_debug_stream("Hits are NOT connected; evaluation of all connectors");
  return false;
};


template <class Hitclass>
void ConnectorBlock::DiagnoseConnected(
  const Hitclass& h1,
  const Hitclass& h2) const
{
  log_debug("Evaluating ConnectedDiagnose()");
  std::stringstream ss;
  
  bool connected = false;
 
  ConnectorList::const_iterator connectorlist_iter=connectorlist_.begin();
  const ConnectorList::const_iterator connectorlist_end=connectorlist_.end();
  while (connectorlist_iter != connectorlist_end){
    bool rel = (*connectorlist_iter)->GetRelation()->AreRelated(h1.GetDOMIndex(), h2.GetDOMIndex());
    bool con = (*connectorlist_iter)->GetConnection()->AreConnected(h1, h2);
    bool relcon = rel && con;

    ss<<" "<<(relcon?"+":"-")<<(*connectorlist_iter)->GetName()<<"("<<"R"<<(rel?"+":"-")<<"C"<<(con?"+":"-")<<")";
    connected |= relcon;
    ++connectorlist_iter;
  }
  ss<<(connected?" --Connected-- ":"");
  
  log_notice_stream(ss.str());
};


inline
CompactOMKeyHashServiceConstPtr
ConnectorBlock::GetHashService() const {
  return hashedGeo_->GetHashService();
};

inline
ConnectorBlock::ConnectorList
ConnectorBlock::GetConnectorList() const {
  return connectorlist_;
};

#endif //HIVECONNECTIONSERVICE_H
