/**
 * \file ConnectionConfig.h
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveSplitter.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 * A confiuration of ringsettings that can be translated into Connectors
 */

#ifndef CONNECTIONCONFIG_H
#define CONNECTIONCONFIG_H

#include "IceHiveZ/__SERIALIZATION.h"

static const unsigned connectionconfig_version_ = 0;

#include <limits>
#include <list>
#include <map>

#include "ToolZ/HashedGeometry.h"

#include "dataclasses/I3Constants.h"
#include "IceHiveZ/internals/Connection.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>


//===================== CLASS ConnectionConfig =================

class ConnectionConfig : virtual public Connection {
public:
  virtual  
  ConnectionPtr BuildConnection(
    const HashedGeometryConstPtr& hashedGeo) =0;
};

typedef boost::shared_ptr<ConnectionConfig> ConnectionConfigPtr;
typedef boost::shared_ptr<const ConnectionConfig> ConnectionConfigConstPtr;

//SERIALIZATION_ASSUME_ABSTRACT(ConnectionConf);

//=================== helper template ==============

template<class Con>
class ConnectionConfig_Helper : public Con, virtual public ConnectionConfig {
  ConnectionPtr BuildConnection(
    const HashedGeometryConstPtr& hashedGeo) 
  {
    boost::shared_ptr<Con> connectionPtr = boost::make_shared<Con>(*this);
    connectionPtr->Configure(hashedGeo);
    return (ConnectionPtr)connectionPtr; //cast to basetype
  };
};

//=================== CLASS  BoolConnectionConfig ==============

typedef ConnectionConfig_Helper<BoolConnection> BoolConnectionConfig;

typedef boost::shared_ptr<BoolConnectionConfig> BoolConnectionConfigPtr;
typedef boost::shared_ptr<const BoolConnectionConfig> BoolConnectionConfigConstPtr;


//=================== CLASS  DeltaTimeConnectionConfig ==============

typedef ConnectionConfig_Helper<DeltaTimeConnection> DeltaTimeConnectionConfig;

typedef boost::shared_ptr<DeltaTimeConnectionConfig> DeltaTimeConnectionConfigPtr;
typedef boost::shared_ptr<const DeltaTimeConnectionConfig> DeltaTimeConnectionConfigConstPtr;


//=================== CLASS  DynamicConnectionConfig ==============

typedef ConnectionConfig_Helper<DynamicConnection> DynamicConnectionConfig;

typedef boost::shared_ptr<DynamicConnectionConfig> DynamicConnectionConfigPtr;
typedef boost::shared_ptr<const DynamicConnectionConfig> DynamicConnectionConfigConstPtr;


//=================== CLASS  PhotonDiffusionConnectionConfig ==============

typedef ConnectionConfig_Helper<PhotonDiffusionConnection> PhotonDiffusionConnectionConfig;

typedef boost::shared_ptr<PhotonDiffusionConnectionConfig> PhotonDiffusionConnectionConfigPtr;
typedef boost::shared_ptr<const PhotonDiffusionConnectionConfig> PhotonDiffusionConnectionConfigConstPtr;

#endif //CONNECTIONCONFIG_H
