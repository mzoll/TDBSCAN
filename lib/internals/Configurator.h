/**
 * \file Configuator.h
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

#ifndef CONFIGURATOR_H
#define CONFIGURATOR_H

#include "IceHiveZ/__SERIALIZATION.h"

static const unsigned icehiveconfig_version_ = 0;

#include <limits>
#include <list>
#include <map>

#include "dataclasses/geometry/I3Geometry.h"
#include "dataclasses/calibration/I3Calibration.h"
#include "dataclasses/status/I3DetectorStatus.h"
#include "ToolZ/OMKeyHash.h"
#include "ToolZ/Hitclasses.h"
#include "IceHiveZ/internals/Connection.h"
#include "IceHiveZ/internals/ConnectionConfig.h"
#include "IceHiveZ/internals/Relation.h"
#include "IceHiveZ/internals/RelationConfig.h"
#include "IceHiveZ/internals/Connector.h"

#include <boost/make_shared.hpp>
#include <boost/foreach.hpp>


///namespace in order to define parameters for IceHive

//===================== CLASS Configurator =========================

/** Stores the configuration so that a Relation and a Connection can be build
 */
class Configurator
{
  SET_LOGGER("Configurator");
private:
  /// a unique name
  std::string name_;
  
public:
  /// PARAM: The connector config to use
  ConnectionConfigPtr connectionConfig_; //NOTE polymorphic //TODO make into list
  /// PARAM: The Relation config to use; needs to have a caller <bool (const OMKey&, const OMKey&)>
  RelationConfigPtr relationConfig_;

public: //interface
  /// Constructor 
  Configurator(const std::string& name);
  
  /// Add (actually specify) a conneector configuration to build from
  void AddConnectionConfig(const ConnectionConfigPtr con);
  /// Add (actually specify) a relation configuration to build from
  void AddRelationConfig(const RelationConfigPtr rel);
  /// Get The name of this configurator
  std::string GetName() const;

  /// Build the Connector from the configurations contained in this object
  ///\param hashedgeo hashed Detector geometry
  Connector BuildConnector (
    const HashedGeometryConstPtr& hashedGeo) const;
};

typedef boost::shared_ptr<Configurator> ConfiguratorPtr;
typedef boost::shared_ptr<const Configurator> ConfiguratorConstPtr;

typedef std::list<Configurator> ConfiguratorList;


//====================== CLASS ConfiguratorBlock ==========

/**
 * Holds a block of COnfigurators and so can build a ConnectorBlock
 */
class ConfiguratorBlock {
public:
  ///holds all sub configurators
  ConfiguratorList config_list_;
  /// specify the OMKeys that should be hashed; a function of signature: bool (const OMKey&)
  boost::function<bool (const OMKey&)> hashOMKeys_;
public: //Setters/Manipulators
  /// add a sub-configurator
  void AddConfigurator(const Configurator& hc);
  /// set a set of OMKeys that are to be soelmny considered 
  void SetOMKeys(const boost::function<bool (const OMKey&)>& hashOMKeys);
public: //ctors
  ///constructor
  ConfiguratorBlock();
public: //methods
  ///Produce a ConnectorBlock by building all Connectors
  ConnectorBlock BuildConnectorBlock (
    const I3OMGeoMap& omgeo) const;
};

typedef boost::shared_ptr<ConfiguratorBlock> ConfiguratorBlockPtr;
typedef boost::shared_ptr<const ConfiguratorBlock> ConfiguratorBlockConstPtr;

#endif //CONFIGURATOR_H
