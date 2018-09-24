/**
 * \file Configurator.cxx
 *
 * (c) 2012 the IceCube Collaboration
 *
 * $Id: HiveSplitter.h 99900 2013-02-26 10:10:43Z mzoll $
 * \version $Revision: 99900 $
 * \date $Date: 2013-02-26 11:10:43 +0100 (Tue, 26 Feb 2013) $
 * \author Marcel Zoll <marcel.zoll@fysik.su.se>
 *
 */

#include "IceHiveZ/internals/Configurator.h"

#define MULTITHREADING_ENABLED 0 //cannot be activated if a python object is hooked up at ::RelationConfig::ConnectFrom/To

#if MULTITHREADING_ENABLED
  #include <thread>         // std::thread
  #include <mutex>          // std::mutex 
  static std::mutex mtx; 
#endif //MULTITHREADING_ENABLED
  
using namespace std;

//=============CLASS Configurator ========================

Configurator::Configurator(
  const std::string& name)
: name_(name)
{};

void Configurator::AddConnectionConfig (
  const ConnectionConfigPtr conconfig) 
{
  if (conconfig->CorrectlyConfigured()) {
    connectionConfig_ = conconfig;
    return;
  }
  log_error("ConnectionConfig is not sufficiently configured");
};

void Configurator::AddRelationConfig (
  const RelationConfigPtr relconfig) 
{
  relationConfig_ = relconfig;
};

std::string Configurator::GetName() const
{return name_;};

Connector Configurator::BuildConnector (
  const HashedGeometryConstPtr& hashedGeo) const 
{
  Connector con(name_, 
                hashedGeo,
                connectionConfig_->BuildConnection(hashedGeo),
                relationConfig_->BuildRelation(hashedGeo));
  return con;
}

// ================= CLASS class ConfiguratorBlock ============

bool AllTrue(const OMKey&) {
  return true;
};

ConfiguratorBlock::ConfiguratorBlock()
: config_list_(),
  hashOMKeys_(AllTrue)
{};

void ConfiguratorBlock::AddConfigurator(const Configurator& hc) {
  config_list_.push_back(hc);
};

void ConfiguratorBlock::SetOMKeys(const boost::function<bool (const OMKey&)>& hashOMKeys) {
  hashOMKeys_ = hashOMKeys;
};

ConnectorBlock ConfiguratorBlock::BuildConnectorBlock (
  const I3OMGeoMap& omgeo) const
{
  //evaluate which OMKeys need to be incorporated and build HasedGeo object
  std::set<OMKey> omkey_set;
  BOOST_FOREACH(const I3OMGeoMap::value_type& omkey_val, omgeo) {
    const OMKey& omkey = omkey_val.first;
    if (hashOMKeys_(omkey))
      omkey_set.insert(omkey_set.end(),omkey);
  }
  const HashedGeometryConstPtr hashedGeo = boost::make_shared<const HashedGeometry>(omgeo, omkey_set);
  
  ConnectorBlock cb(hashedGeo);
  
  //reorder the Configurator list so that the fastest evaluating Connectors are on top
  ConfiguratorList config_list_sorted(config_list_.begin(), config_list_.end());
  struct desc_speed {
    bool operator() (const Configurator& lhs, const Configurator& rhs) const
    { return lhs.connectionConfig_->GetSpeedRating() > rhs.connectionConfig_->GetSpeedRating(); };
  };
  config_list_sorted.sort(desc_speed());
  
#if MULTITHREADING_ENABLED
  struct Configurator_BuildConnector_Wrapper {
    ConnectorPtr con_; //return object of subcall
    const Configurator& config_; // wrapped object
    const HashedGeometryConstPtr hashedGeo_; //call parameter
    
    ///construct
    Configurator_BuildConnector_Wrapper (
      ConnectorPtr& con,
      const Configurator& config,
      const HashedGeometryConstPtr hashedGeo)
    : con_(con),
      config_(config),
      hashedGeo_(hashedGeo),
    {};
    
    /// parameterless call
    void operator() () {
      log_debug("BUILD");
      ConnectorPtr con = boost::make_shared<Connector>(config_.BuildConnector(hashedGeo_));
      log_debug("DONE");
      con_ = con;
    };
  }; //Configurator_BuildConnector_Wrapper

  //do the tedious hashing now, so there will be no I/O conflicts later between threads
  hashedGeo.distService_->HashAllDistances();
  
  //prepare some lists to store builders and their results
  std::vector<ConnectorPtr> connector_list;
  connector_list.reserve(config_list_sorted.size());
  std::vector<std::thread*> connector_builders;
  connector_builders.reserve(config_list_sorted.size());
  
  //spawn the threads
  BOOST_FOREACH(const Configurator& config, config_list_sorted) {
  //const Configurator& config = *config_list_sorted.begin(); {
    ConnectorPtr conptr;
    connector_list.push_back(conptr);
    Configurator_BuildConnector_Wrapper wrap(hashedGeo);
    log_info_stream("Spawn thread for "<<config.GetName());
    connector_builders.push_back(new std::thread(wrap));
  }

  //join all treads and add up the results
  for (size_t i=0; i< connector_builders.size(); i++) {
    log_warn_stream("waiting"<<i);
    connector_builders.at(i)->join();
    log_warn_stream("waited"<<i);
    //preserving order in list //NOTE also possible to do this asyncronous, I guess
    cb.AddConnector(connector_list[i]);
  }
#else
  BOOST_FOREACH(const Configurator& c, config_list_sorted) {
    const ConnectorPtr con= boost::make_shared<Connector>(c.BuildConnector(hashedGeo));
    cb.AddConnector(con);
  };
#endif
  return cb;
}