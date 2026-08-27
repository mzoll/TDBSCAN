/**
 * \file Connector.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: DistanceMapService.cxx 129048 2015-02-13 14:42:33Z mzoll $
 * \version $Revision: 129048 $
 * \date $Date: 2015-02-13 15:42:33 +0100 (fre, 13 feb 2015) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 */

#include "IceHiveZ/internals/Connector.h"

#include <boost/foreach.hpp>

using namespace indexmatrix;

//============ CLASS Connector ===========

Connector::Connector(
  const std::string& name,
  const HashedGeometryConstPtr& hashedGeo,
  const ConnectionPtr& connection, 
  const RelationPtr& relation)
: name_(name),
  hashedGeo_(hashedGeo),
  connection_(connection),
  relation_(relation)
{
  //check hasher against the hashers of subclasses
//   if (!connection_->GetHasher())
//     log_fatal("connection_->GetHasher() not set");
//   if (!(relation_->GetHasher()))
//     log_fatal("(relation_->GetHasher()) not set");
//   
//   if (hasher_!=connection_->GetHasher()) //FIXME here need go a verify method
//     log_fatal_stream("trying to add potentially different hashers within the ConnectorBlock :: "
//       <<&(*hasher_)<<" vs "<<&(*(connection_->GetHasher())));
//   if (hasher_!=relation_->GetHasher())
//     log_fatal_stream("trying to add potentially different hashers within the ConnectorBlock :: "
//       <<&(*hasher_)<<" vs "<<&(*(relation_->GetHasher())));
};

#if SERIALIZATION_ENABLED
  #ifdef SERIALIZATON_SUPPORT_ICECUBE
I3_SERIALIZABLE(Connector);
  #endif // SERIALIZATON_SUPPORT_ICECUBE
  #ifdef SERIALIZATON_SUPPORT_BOOST
BOOST_CLASS_EXPORT_IMPLEMENT(Connector);
  #endif // SERIALIZATON_SUPPORT_BOOST
#endif //SERIALIZATON_ENABLED

//====================== CLASS ConnectorBlock ============

ConnectorBlock::ConnectorBlock(
  const HashedGeometryConstPtr& hashedGeo)
: hashedGeo_(hashedGeo),
  cumulativeRel_(boost::make_shared<Relation>(hashedGeo->GetHashService()))
{};

void ConnectorBlock::AddConnector (
  const ConnectorPtr& c) 
{
  log_info_stream("Adding Connector '"<<c->name_<<"' to ConnectorBlock";);
  //probe if the Hasher-object is the same
  //FIXME make a consistency test
  connectorlist_.push_back(c);
  cumulativeRel_->Join(*(c->relation_));
};

ConnectorPtr 
ConnectorBlock::GetConnector(const int index) const 
{
  if (!(index >= -1 && index < int(connectorlist_.size()))) {
    log_error("There is no pointer at that index in the ConnectorList; passing empty pointer");
    ConnectorPtr not_init;
    return not_init;
  }
  if (index == -1) {
    return boost::make_shared<Connector>("cumulative",
                                         hashedGeo_,
                                         boost::make_shared<BoolConnection>(hashedGeo_, true),
                                         cumulativeRel_);
  }
  else {
    ConnectorList::const_iterator connectorlist_iter = connectorlist_.begin(); //no indexing in lists :/
    for (size_t i = 0; i<(size_t)index && i<connectorlist_.size(); i++)
      ++connectorlist_iter;
    return *connectorlist_iter;
  }
};

#if SERIALIZATION_ENABLED
  I3_SERIALIZABLE(ConnectorBlock);
#endif //SERIALIZATON_ENABLED  